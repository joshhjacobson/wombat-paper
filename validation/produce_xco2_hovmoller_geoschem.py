import pandas as pd
import xarray as xr
import xesmf as xe


## setup
target_grid = xe.util.grid_global(1, 1, cf=True)
END_MONTH = "2017-04-01"


geoschem_spec_conc_glob = (
    "../1_transport/intermediates/GEOS_Chem/runs/run.v12.3.2.base/output/"
    "GEOSChem.SpeciesConc.*_0000z.nc4"
)
geoschem_level_edge_glob = (
    "../1_transport/intermediates/GEOS_Chem/runs/run.v12.3.2.base/output/"
    "GEOSChem.LevelEdgeDiags.*_0000z.nc4"
)


def compute_pressure_weights(ds, level_name: str):
    surface_pressure = ds[level_name].max(dim="ilev")
    weights = abs(ds[level_name].diff("ilev").rename(ilev="lev") / surface_pressure)
    # levels will take the value of the minuend’s ilev coordinate, instead we need to update
    # to the dataset's level coordinate
    weights["lev"] = ds["lev"]
    return weights


def compute_xco2(ds, mole_frac_name: str, weights_name: str):
    # compute pressure-weighted average over vertical column in units ppm
    return (ds[mole_frac_name] * ds[weights_name]).sum(dim="lev") * 1e6


print("Setup complete")

## produce xco2 datasets of monthly longitudinal averages


with xr.open_mfdataset(
    [
        "../1_transport/intermediates/GEOS_Chem/runs/run.v12.3.2.base/output/"
        "GEOSChem.LevelEdgeDiags.20140901_0000z.nc4",
        "../1_transport/intermediates/GEOS_Chem/runs/run.v12.3.2.base/output/"
        "GEOSChem.SpeciesConc.20140901_0000z.nc4",
    ]
) as ds:
    # precompute geos chem grid weights
    ds = ds.isel(time=0)
    ds["pressure_weights"] = compute_pressure_weights(ds, "Met_PEDGE")
    regridder_geoschem = xe.Regridder(
        ds[["SpeciesConc_CO2", "pressure_weights"]], target_grid, "conservative"
    )


def prep_geoschem_lev(ds):
    return ds["Met_PEDGE"]


da_geoschem_lev = xr.open_mfdataset(
    geoschem_level_edge_glob, preprocess=prep_geoschem_lev, chunks={"time": 1000}, parallel=True
)


def prep_geoschem(ds):
    return ds[["SpeciesConc_CO2"]]


with xr.open_mfdataset(
    geoschem_spec_conc_glob, preprocess=prep_geoschem, chunks={"time": 1000}, parallel=True
) as ds:
    ds["Met_PEDGE"] = da_geoschem_lev["Met_PEDGE"]
    ds = ds.where(ds["time"] < pd.to_datetime(END_MONTH), drop=True)
    ds["pressure_weights"] = compute_pressure_weights(ds, "Met_PEDGE")
    ds_geoschem = regridder_geoschem(ds[["SpeciesConc_CO2", "pressure_weights"]])
    da_geoschem_xco2 = compute_xco2(ds_geoschem, "SpeciesConc_CO2", "pressure_weights")
    da_geoschem_xco2_hov = da_geoschem_xco2.mean(dim="lon").resample(time="1M").mean()


print("Dataset configured")

da_geoschem_xco2_hov.to_netcdf("../hovmoller_array_geoschem.nc")

print("DONE")
