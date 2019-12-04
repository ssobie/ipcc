##Script to calculate the standard derived variables from the 
##split CMIP6 downscaled files

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=TRUE)
source('/home/ssobie/code/repos/downscale_CMIP6/bccaqv2.derived.variable.support.r',chdir=TRUE)
source('/home/ssobie/code/repos/downscale_CMIP6/standard.variables.functions.r',chdir=TRUE)

##----

source('/home/ssobie/code/repos/ipcc/cmip6.canada.derived.files.r')

library(udunits2)
library(ncdf4)
library(PCICt)
library(foreach)
library(doParallel)
registerDoParallel(cores=2) # or some other number that you're comfortable with.

##--------------------------------------------------------------

gcm_unit_conversion <- function(varname,gcm.subset,gcm.nc) {

   units <- switch(varname,
                   pr='kg m-2 d-1',
                   tasmax='degC',
                   tasmin='degC')

   var.units <- ncatt_get(gcm.nc,varname,'units')$value
   if (var.units != units) {
      rv <- ud.convert(gcm.subset,var.units,units)
   } else {
      rv <- gcm.subset
   }
   return(rv)
}

##--------------------------------------------------------------

##****************************************************************
testing <- FALSE

res <- 'gr2'
if (testing) {
   tmpdir <- '/local_temp/ssobie'
   gcm <- 'BCC-CSM2-MR'
   scenario <- 'ssp585'
   run <- 'r1i1p1f1' 
   res <- NULL
   type <- 'annual'
   varname <- 'tasmax'
   pctl <- '000'
} else {
   args <- commandArgs(trailingOnly=TRUE)
   for(i in 1:length(args)){
       eval(parse(text=args[[i]]))
   }
}



tmp.dir <- paste0(tmpdir,'/',gcm,'_',scenario,'_',run,'_',type,'_',pctl,'_',varname,'/')
if (!file.exists(tmp.dir)) {
  dir.create(tmp.dir,recursive=TRUE)
}

base.dir <- '/storage/data/climate/CMIP6/'
if (is.null(res)) {
   write.dir <- paste0(base.dir,'/Derived/',gcm,'_',scenario,'_',run,'/')
} else {
   write.dir <- paste0(base.dir,'/Derived/',gcm,'_',scenario,'_',run,'_',res,'/')
}
if (!file.exists(write.dir))
  dir.create(write.dir,recursive=TRUE)

##Transfer template files for derived file creation
gcm.dir <- paste0(base.dir,'assembled/',gcm,'/')
var.files <- list.files(path=gcm.dir,pattern=varname)
rcp.files <- var.files[grep(scenario,var.files)]
if (is.null(res)) {
   gcm.file <- rcp.files[grep(run,rcp.files)]
} else {
   run.files <- rcp.files[grep(run,rcp.files)]
   gcm.file <- run.files[grep(res,run.files)] ##Include for GCMs with multiple resolutions
}

file.copy(from=paste0(gcm.dir,gcm.file),to=tmp.dir)
print('Done copying gcm file')

##---------------------------------------------------------------------------
##Annual Averages

if (type=='annual') {

  print('Ann Avg opening')
  out.dir <- paste0(tmp.dir,'annual/')
  dir.create(paste0(write.dir,'annual/'),recursive=TRUE,showWarnings=FALSE)
  ann.file <-  make_annual_file(varname,gcm,scenario,run,
                                gcm.file,out.dir,tmp.dir)
  ann.ncs <- nc_open(paste0(out.dir,ann.file),write=TRUE)  

}

##---------------------------------------------------------------------------
##Annual Block Maxima Files for writing

if (type=='annual_extremes') {

  print('Ann extremes opening')
  out.dir <- paste0(tmp.dir,'annual_extremes/')
  dir.create(paste0(write.dir,'annual_extremes/'),recursive=TRUE,showWarnings=FALSE)
  ext.file <- make_return_period_file(varname,gcm,scenario,run,
                                      gcm.file,out.dir,tmp.dir)
  ext.ncs <- nc_open(paste0(out.dir,ext.file),write=TRUE)  
}

##---------------------------------------------------------------------------
##Annual Quantile Files for writing

if (type == 'annual_quantiles') {

  print('Ann quantiles opening')
  out.dir <- paste0(tmp.dir,'annual_quantiles/')
  dir.create(paste0(write.dir,'annual_quantiles/'),recursive=TRUE,showWarnings=FALSE)
  quant.file <- make_quantile_file(varname,gcm,scenario,run,pctl,
                                   gcm.file,out.dir,tmp.dir)
  quant.ncs <- nc_open(paste0(out.dir,quant.file),write=TRUE)
}

##---------------------------------------------------------------------------
##Seasonal Average Files for writing

if (type=='seasonal') {
  print('Seasonal averages opening')
  out.dir <- paste0(tmp.dir,'seasonal/')
  dir.create(paste0(write.dir,'seasonal/'),recursive=TRUE,showWarnings=FALSE)
  seas.file <- make_seasonal_file(varname,gcm,scenario,run,
                                  gcm.file,out.dir,tmp.dir)
  seas.ncs <- nc_open(paste0(out.dir,seas.file),write=TRUE)
}

##---------------------------------------------------------------------------
##Monthly Average Files for writing

if (type=='monthly') {

  print('monthly avg opening')
  out.dir <- paste0(tmp.dir,'monthly/')
  dir.create(paste0(write.dir,'monthly/'),recursive=TRUE,showWarnings=FALSE)
  mon.file <- make_monthly_file(varname,gcm,scenario,run,
                                gcm.file,out.dir,tmp.dir)
  mon.ncs <- nc_open(paste0(out.dir,mon.file),write=TRUE)
}

##---------------------------------------------------------------------------
##---------------------------------------------------------------------------

##Iterate over the split files
gcm.nc <- nc_open(paste0(tmp.dir,gcm.file),write=FALSE)
gcm.dates <- netcdf.calendar(gcm.nc)
yearly.fac <- as.factor(format(gcm.dates,'%Y'))
seasonal.fac <- get_seasonal_fac(gcm.dates)
monthly.fac <- as.factor(format(gcm.dates,'%Y-%m'))  

lon <- ncvar_get(gcm.nc,'lon')
lat <- ncvar_get(gcm.nc,'lat')
n.lon <- length(lon)
n.lat <- length(lat)


for (j in 1:n.lat) {
   lat.ix <- j
   ltm <- proc.time()

   print(paste0('Latitude: ',j,' of ',n.lat))
   gcm.subset <- ncvar_get(gcm.nc,varname,start=c(1,j,1),count=c(-1,1,-1))
   flag <- is.na(gcm.subset[,1])
   gcm.list <- vector(mode='list',length=n.lon)
   gcm.list <- lapply(seq_len(nrow(gcm.subset)), function(k) gcm.subset[k,])
   rm(gcm.subset)

   ##----------------------------------------------------------
   ##Annual Averages 
   if (type=='annual') {
     annual_averages_for_model(varname,ann.ncs,lat.ix,n.lon,yearly.fac,flag,
                               gcm.list)
   }
   ##----------------------------------------------------------
   ##Seasonal Averages 
   if (type=='seasonal') {
     seasonal_averages_for_model(varname,seas.ncs,lat.ix,n.lon,seasonal.fac,flag,
                                 gcm.list)
   }
   ##----------------------------------------------------------
   ##Monthly Averages 
   if (type=='monthly') {
     monthly_averages_for_model(varname,mon.ncs,lat.ix,n.lon,monthly.fac,flag,
                                gcm.list)
   }
   ##----------------------------------------------------------
   ##Annual Extremes
   if (type=='annual_extremes') {
     annual_extremes_for_model(varname,ext.ncs,lat.ix,n.lon,yearly.fac,flag,
                               gcm.list)
   }
   ##----------------------------------------------------------
   if (type=='annual_quantiles') {
   ##Annual Quantiles
     annual_quantiles_for_model(varname,pctl,quant.ncs,lat.ix,n.lon,yearly.fac,flag,
                                gcm.list)
   }
   ##----------------------------------------------------------
   rm(gcm.list)
}##Longitude Loop
nc_close(gcm.nc)

file.remove(paste0(tmp.dir,"/",gcm.file))

if (type=='annual') {
   nc_close(ann.ncs)
   file.copy(from=paste0(out.dir,ann.file),to=paste0(write.dir,'annual/'),overwrite=TRUE)
}

if (type=='seasonal') {
   nc_close(seas.ncs)
   file.copy(from=paste0(out.dir,seas.file),to=paste0(write.dir,'seasonal/'),overwrite=TRUE)
}

if (type=='monthly') {
   nc_close(mon.ncs)
   file.copy(from=paste0(out.dir,mon.file),to=paste0(write.dir,'monthly/'),overwrite=TRUE)
}

if (type=='annual_extremes') {
   nc_close(ext.ncs)
   file.copy(from=paste0(out.dir,ext.file),to=paste0(write.dir,'annual_extremes/'),overwrite=TRUE)
}

if (type=='annual_quantiles') {
   nc_close(quant.ncs)
   file.copy(from=paste0(out.dir,quant.file),to=paste0(write.dir,'annual_quantiles/'),overwrite=TRUE)
}


