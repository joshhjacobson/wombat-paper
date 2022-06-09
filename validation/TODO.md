
## Data preparation

- [ ] create single input dataset for each model
    - [ ] MOZART: concatenate the "CO2_SRF_EMIS_avrg" variable across all files with format outputs/BasisFnsUpdated/yyyymm/*.h0.yyyy-mm-01-03600.nc
    - [ ] GEOS Chem: concatenate the "EmisCO2_Total" variable across all files with format runs/run.v12.3.2.base/output/HEMCO_diagnostics.yyyymm.nc, and compute the sum over the pressure levels

## Inputs

- [ ] annual aggregate flux: inputs could be mis-specified so check each model's diagnostics 
    - [ ] The input datasets both have units kg/m^2/s. For each dataset, compute total kg by multiplying through the correct units.
    - [ ] Produce a time series for of total global emitted mass of CO2 


## Outputs

merge first h1 (daily) file into one large dataset that covers the study period

focus on monthly analysis

plot maps of vertical average for the very first timestep (hourly version)
 - create difference plot

1. surface level
    - zonal averaging (hovmoller, two plots or difference)
        - don't need to do area weighting (fixed lat)
    - global average (time series, two lines on the same plot)

2. vertical averages
    - zonal averages (hovmoller)
    - global (time series)
        - should be same between mozart and geos chem


do this before moving on to a couple other sensitivity runs

## Notes

For vertical averaging, see profile.pdf
equations 1 & 2