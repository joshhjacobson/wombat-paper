import pandas as pd
import xarray as xr
import xesmf as xe
from matplotlib import gridspec
import matplotlib.pyplot as plt
import seaborn as sns

sns.set()


## setup
target_grid = xe.util.grid_global(1, 1, cf=True)
date_range = pd.date_range(start="2014-09", end="2017-04", freq="1M")
mozart_paths = []
for month in date_range:
    yyyy, mm, _ = str(month).split("-")
    mozart_paths.append(
        f"../1_transport/intermediates/MOZART/output/BasisFnsUpdated/{yyyy}{mm}/"
        f"BasisFnsUpdated.mz4.h0.{yyyy}-{mm}-01-03600.nc"
    )

geoschem_spec_conc_glob = (
    "../1_transport/intermediates/GEOS_Chem/runs/run.v12.3.2.base/output/"
    "GEOSChem.SpeciesConc.*_0000z.nc4"
)
geoschem_level_edge_glob = (
    "../1_transport/intermediates/GEOS_Chem/runs/run.v12.3.2.base/output/"
    "GEOSChem.LevelEdgeDiags.*_0000z.nc4"
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


## produce xco2 datasets of monthly longitudinal averages

with xr.open_dataset(mozart_paths[0], decode_times=False) as ds:
    # precompute mozart grid weights
    ds["pressure_edge"] = get_mozart_pressure_edges(ds)
    ds["pressure_weights"] = compute_pressure_weights(ds, "pressure_edge")
    ds_mozart = ds[["CO2_VMR_avrg", "pressure_weights"]]
    regridder_mozart = xe.Regridder(ds_mozart, target_grid, "conservative")


def prep_mozart(ds):
    ds["time"] = pd.to_datetime(ds.date.values, format="%Y%m%d") + pd.to_timedelta(
        ds.datesec.values, unit="seconds"
    )
    # ds = ds.resample(time="1M").mean()
    ds["pressure_edge"] = get_mozart_pressure_edges(ds)
    ds["pressure_weights"] = compute_pressure_weights(ds, "pressure_edge")
    ds_mozart = regridder_mozart(ds[["CO2_VMR_avrg", "pressure_weights"]])
    da_mozart_xco2 = compute_xco2(ds_mozart, "CO2_VMR_avrg", "pressure_weights")
    return da_mozart_xco2.mean(dim="lon")


with xr.open_mfdataset(
    mozart_paths,
    preprocess=prep_mozart,
    chunks={"time": 1000},
    parallel=True,
    decode_times=False,
) as da:
    da = da.where(da["time"] < pd.to_datetime("2017-04-01"), drop=True)
    da_mozart_xco2_hov = da.resample(time="1M").mean()
    # da_mozart_xco2_hov = da.where(da["time"] < pd.to_datetime("2017-04-01"), drop=True)


with xr.open_mfdataset(
    [
        "../1_transport/intermediates/GEOS_Chem/runs/run.v12.3.2.base/output/"
        "GEOSChem.LevelEdgeDiags.20140901_0000z.nc4",
        "../1_transport/intermediates/GEOS_Chem/runs/run.v12.3.2.base/output/"
        "GEOSChem.SpeciesConc.20140901_0000z.nc4",
    ]
) as ds:
    # precompute geos chem grid weights
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
    # NOTE: different order of operations...
    ds = (
        ds.where(ds["time"] < pd.to_datetime("2017-04-01"), drop=True)
        .resample(time="1M")
        .mean()
    )
    ds["pressure_weights"] = compute_pressure_weights(ds, "Met_PEDGE")
    ds_geoschem = regridder_geoschem(ds[["SpeciesConc_CO2", "pressure_weights"]])
    da_geoschem_xco2 = compute_xco2(ds_geoschem, "SpeciesConc_CO2", "pressure_weights")
    da_geoschem_xco2_hov = da_geoschem_xco2.mean(dim="lon")

da_xco2_hov_diff = da_geoschem_xco2_hov - da_mozart_xco2_hov


## create and save the plot

fig = plt.figure(figsize=(18, 8), constrained_layout=True)

gs = gridspec.GridSpec(2, 3, figure=fig, height_ratios=[1, 0.05])
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1], sharey=ax1)
ax3 = fig.add_subplot(gs[0, 2], sharey=ax1)
main_cbar_ax = fig.add_subplot(gs[1, :2])
diff_cbar_ax = fig.add_subplot(gs[1, 2])

# vmin, vmax = 395, 420
# sub_geoschem = da_geoschem_xco2_hov.plot(
#     ax=ax1, vmin=vmin, vmax=vmax, cmap="jet", add_colorbar=False
# )
# sub_mozart = da_mozart_xco2_hov.plot(
#     ax=ax2, vmin=vmin, vmax=vmax, cmap="jet", add_colorbar=False
# )
sub_geoschem = da_geoschem_xco2_hov.plot(ax=ax1, robust=True, cmap="jet", add_colorbar=False)
sub_mozart = da_mozart_xco2_hov.plot(ax=ax2, robust=True, cmap="jet", add_colorbar=False)
fig.colorbar(sub_geoschem, cax=main_cbar_ax, orientation="horizontal", label="CO2 [ppm]")

sub_diff = da_xco2_hov_diff.plot(ax=ax3, robust=True, cmap="RdYlBu_r", add_colorbar=False)
fig.colorbar(sub_diff, cax=diff_cbar_ax, orientation="horizontal", label="CO2 [ppm]")

ax1.set_title("GEOS Chem", fontsize=12)
ax1.set_ylabel("Month", fontsize=12)
ax2.set_title("MOZART", fontsize=12)
ax2.set_ylabel(None)
plt.setp(ax2.get_yticklabels(), visible=False)
ax3.set_title("Difference: GEOS Chem minus MOZART", fontsize=12)
ax3.set_ylabel(None)
plt.setp(ax3.get_yticklabels(), visible=False)
for ax in [ax1, ax2, ax3]:
    ax.set_xlabel("Latitude", fontsize=12)

fig.suptitle(
    "Hovmoller plots for latitude: monthly pressure-weighted vertical average (XCO2)",
    fontsize=14,
)
fig.savefig(f"../figures/hovmoller_lat_xco2.png", dpi=200)
