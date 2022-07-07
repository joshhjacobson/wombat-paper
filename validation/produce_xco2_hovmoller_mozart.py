import pandas as pd
import xarray as xr
import xesmf as xe

# TODO: setup dask distributed?


## setup
target_grid = xe.util.grid_global(1, 1, cf=True)
date_range = pd.date_range(start="2014-09", end="2017-04", freq="1M")
END_MONTH = "2017-04-01"

mozart_paths = []
for month in date_range:
    yyyy, mm, _ = str(month).split("-")
    mozart_paths.append(
        f"../1_transport/intermediates/MOZART/output/BasisFnsUpdated/{yyyy}{mm}/"
        f"BasisFnsUpdated.mz4.h0.{yyyy}-{mm}-01-03600.nc"
    )


def get_mozart_pressure_edges(ds):
    # compute pressure edges in units hPA; equation: https://www2.acom.ucar.edu/gcm/mozart-4
    return (ds["P0"] * ds["hyai"] + ds["PS"] * ds["hybi"]) * 0.01


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

with xr.open_dataset(mozart_paths[0], decode_times=False) as ds:
    # precompute mozart grid weights
    ds = ds.isel(time=0)
    ds["pressure_edge"] = get_mozart_pressure_edges(ds)
    ds["pressure_weights"] = compute_pressure_weights(ds, "pressure_edge")
    ds_mozart = ds[["CO2_VMR_avrg", "pressure_weights"]]
    regridder_mozart = xe.Regridder(ds_mozart, target_grid, "conservative")


def prep_mozart(ds):
    ds["time"] = pd.to_datetime(ds.date.values, format="%Y%m%d") + pd.to_timedelta(
        ds.datesec.values, unit="seconds"
    )
    print(ds["time"].values[0])
    ds["pressure_edge"] = get_mozart_pressure_edges(ds)
    ds["pressure_weights"] = compute_pressure_weights(ds, "pressure_edge")
    ds_mozart = regridder_mozart(ds[["CO2_VMR_avrg", "pressure_weights"]])
    da_mozart_xco2 = compute_xco2(ds_mozart, "CO2_VMR_avrg", "pressure_weights")
    return da_mozart_xco2.mean(dim="lon").resample(time="1M").mean()


with xr.open_mfdataset(
    mozart_paths[4:8],
    preprocess=prep_mozart,
    chunks={"time": 1000},
    parallel=True,
    decode_times=False,
) as da:
    da_mozart_xco2_hov = da.where(da["time"] < pd.to_datetime(END_MONTH), drop=True)


print("Dataset configured")

da_mozart_xco2_hov.to_netcdf("../data/hovmoller_array_mozart.nc")

print("DONE")
