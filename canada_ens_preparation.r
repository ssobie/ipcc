##Create an RData object containing an ensemble average 
##of GCM data across Canada for plotting


library(rgdal)
library(raster)
library(ncdf4)
library(PCICt)

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

##--------------------------------------------------------

##********************************************************

var.name <- 'pr'
type <- 'annual'
ix <- 1

scenario <- 'ssp126'
interval <- '1971-2000'

experiment <- 'CMIP6'

##tmp.dir <- paste0(tmpdir,'/canada_map_prep/')
##if (!file.exists(tmp.dir)) {
##   dir.create(tmp.dir,recursive=TRUE)
##}

base.dir <- "/storage/data/projects/rci/data/cas/"
cmip.dir <- paste0(base.dir,tolower(experiment),'/one_degree/climatologies/')
save.dir <- paste0(base.dir,tolower(experiment),'/map_climatologies/')

##CMIP5
if (experiment == 'CMIP5') {

   omit.dirs <- c()

}
##CMIP6
if (experiment == 'CMIP6') {

}

all.dirs <- list.files(cmip.dir)
scen.dirs <- all.dirs[grep(scenario,all.dirs)]

var.stack <- stack()

for (gcm.dir in scen.dirs) {
   print(gcm.dir)
   read.dir <- paste0(cmip.dir,gcm.dir,'/',type,'/')
   var.files <- list.files(read.dir,var.name)
   int.file <- var.files[grep(interval,var.files)]
   var.brick <- brick(paste0(read.dir,int.file))
   var.subset <- subset(var.brick,ix)

   var.stack <- stack(var.stack,var.subset)
}

stack.file <- paste0(save.dir,var.name,'_',type,'_',scenario,'_',experiment,'_ensemble_stack_',interval,'.RData')
save(var.stack,file=stack.file)

var.ens.avg <- calc(var.stack,mean,na.rm=T)
save.file <- paste0(save.dir,var.name,'_',type,'_',scenario,'_',experiment,'_ensemble_average_',interval,'.RData')
save(var.ens.avg,file=save.file)
