
## Data preparation

- [x] create single input dataset for each model
    - [x] MOZART: concatenate the "CO2_SRF_EMIS_avrg" variable across all files with format outputs/BasisFnsUpdated/yyyymm/*.h0.yyyy-mm-01-03600.nc
    - [x] GEOS Chem: concatenate the "EmisCO2_Total" variable across all files with format runs/run.v12.3.2.base/output/HEMCO_diagnostics.yyyymm.nc, and compute the sum over the pressure levels

- [x] create single output dataset for each model
    - [x] MOZART: merge first h0 (hourly) file into one large dataset that covers the study period
        - [x] construct pressure edge variable; formula that involves surface pressure (PA), reference pressure (PO), and hybrid coefficients (P0*hyai + PS*hybi) / 100
        - [x] regrid to 1 x 1
    - [x] GEOS Chem: collect hourly dataset of mole fraction across study period
        - [x] regrid to 1 x 1

## Inputs

- [x] annual aggregate flux: inputs could be mis-specified so check each model's diagnostics 
    - [x] The input datasets both have units kg/m^2/s. For each dataset, compute total kg by multiplying through the correct units.
    - [x] Produce a time series for of total global emitted mass of CO2 

## Outputs

### Base run

- [x] compute vertical average (pressure levels) for the very first timestep (hourly version) for both models
    - [x] map the difference

- [x] surface level (average surface co2 to monthly)
    - [x] zonal averaging (hovmoller, two plots or difference)
        - don't need to do area weighting (fixed lat)
    - [ ] global average (time series, two lines on the same plot)

- [ ] vertical averages (monthly xco2 files)
    - [ ] longitudinal averages (hovmoller)
    - [ ] zonal averages (hovmoller)
    - [ ] global (time series)
        - should be same between mozart and geos chem


### Sensitivity runs

Repeat the above analysis for a couple sensitivity runs

## Notes

Focus on monthly analysis

For vertical averaging, see profile.pdf
equations 1 & 2

MOZART has surface as last element
GEOS Chem has surface as first element