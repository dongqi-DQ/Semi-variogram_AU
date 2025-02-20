# Semi-variogram_AU
Code and data for semi-variogram calculation using BoM daily rain gauges in Australia

Semi-variogram is based on Touma et al (2018) DOI: https://doi.org/10.1175/JCLI-D-18-0019.1


## Pre-processing
In `pre_process_BoM`:
- `Scraper_BoM.ipynb`: BoM daily rain gauge data scraper (csv files in zip files)
- `BoM_daily_gauge_to_nc.ipynb`: Convert BoM csv file to NetCDF file
- `Assign_percentile.ipynb`: Calculate the percentile of daily precipitation for each station (rainny days only)
- `Find_p90_all.ipynb`: Calculate the p90 value for each station. This creates an input file for semi-variogram

## Main script

`Calc_semi_variogram.ipynb`: calculates the semi-variogram for each station, and save the data for each station. The station CSVs have information for all neighbour stations, and the bins CSV are the binned semi-variogram <Equation 2 in Touma et al. (2018)>.


`Calc_scale_alpha.ipynb`: calculates the alpha value as shown in Equation 3 in Touma et al. (2018). This will output CSV files with scale of each extreme day for each station. Each station will end up with three CSV files marked by the station ID. Each file has all the extreme day info for the station.  
