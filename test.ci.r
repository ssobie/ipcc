library(ncdf4)
library(PCICt)

source('netcdf.calendar.R')

source('bisect.R')
source('io.R')
source('CI.R')
source('config.R')

tmp.dir <- '/local_temp/ssobie/citest/'

gcm.dir <- '/storage/data/climate/CMIP6/assembled/CanESM5/north_america/'
gcm.file <- 'tasmax_day_CanESM5_North_America_historical+ssp245_r1i1p1f1_gn_18500101-21801231.nc'

file.copy(from=paste0(gcm.dir,gcm.file),to=tmp.dir)

obs.dir <- '/storage/data/climate/observations/gridded/ANUSPLIN/ANUSPLIN_300ARCSEC/'
obs.file <- 'anusplin_tasmax_final.nc'
file.copy(from=paste0(obs.dir,obs.file),to=tmp.dir)

ci.file <- 'tasmax_day_BCCI_CanESM5_North_America_historical+ssp245_r1i1p1f1_gn_18500101-21801231.nc'

varname <- 'tasmax'

gcm.tmp <- paste0(tmp.dir,gcm.file)
obs.tmp <- paste0(tmp.dir,obs.file)
ci.tmp <- paste0(tmp.dir,ci.file)

test <- ci.netcdf.wrapper(gcm.tmp, obs.tmp, ci.tmp, varname)
