##Script to calculate and write the standard set of derived variables
##for the 800m data

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)
source('/home/ssobie/code/repos/downscale_CMIP6/bccaqv2.derived.variable.support.r',chdir=T)
source('/home/ssobie/code/repos/downscale_CMIP6/degree.day.variables.functions.r',chdir=T)

##----

source('/home/ssobie/code/repos/ipcc/cmip6.canada.derived.files.r')

library(udunits2)
library(ncdf4)
library(PCICt)
library(foreach)
library(climdex.pcic)

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

res <- NULL
if (testing) {
   tmpdir <- '/local_temp/ssobie'
   gcm <- 'GFDL-CM4'
   scenario <- 'ssp245'
   run <- 'r1i1p1f1'
   res <- 'gr1'
   type <- 'degree_days'
} else {
   args <- commandArgs(trailingOnly=TRUE)
   for(i in 1:length(args)){
       eval(parse(text=args[[i]]))
   }
}

tmp.dir <- paste0(tmpdir,'/',gcm,'_',scenario,'_',run,'_',type,'/')
if (!file.exists(tmp.dir)) {
  dir.create(tmp.dir,recursive=TRUE)
}

degree.names <- c('cdd','fdd','gdd','hdd')

base.dir <- '/storage/data/climate/CMIP6/'
###base.dir <-  '/storage/data/climate/CMIP5/daily/'

if (is.null(res)) {
   write.dir <- paste0(base.dir,'/Derived/',gcm,'_',scenario,'_',run,'/')
} else {
   write.dir <- paste0(base.dir,'/Derived/',gcm,'_',scenario,'_',run,'_',res,'/')
}
if (!file.exists(write.dir))
  dir.create(write.dir,recursive=TRUE)

##Transfer template files for derived file creation

gcm.dir <- paste0(base.dir,'assembled/',gcm,'/') ##CMIP6
###gcm.dir <- paste0(base.dir,gcm,'/') ##CMIP5
rcp.files <- list.files(path=gcm.dir,pattern=paste0('historical\\+',scenario))

if (is.null(res)) {
   run.files <- rcp.files[grep(run,rcp.files)]
} else {
   res.files <- rcp.files[grep(res,rcp.files)] ##Include for GCMs with multiple resolutions
   run.files <- res.files[grep(run,res.files)]
}

tasmax.file <- run.files[grep(paste0('tasmax_day_',gcm),run.files)]
tasmin.file <- run.files[grep(paste0('tasmin_day_',gcm),run.files)]

if (length(tasmax.file)!=1 | length(tasmin.file) !=1) {
   print(tasmax.file)
   print(tasmin.file)
   stop('More than one tasmax or tasmin file selected')
}

file.copy(from=paste0(gcm.dir,tasmax.file),to=tmp.dir)
print('Done copying gcm file')


if (type=='degree_days') {
  print('Degree days opening')
  freq <- 'Ann'
  out.dir <- paste0(tmp.dir,'degree_days/')
  dir.create(paste0(write.dir,'degree_days/'),recursive=TRUE,showWarnings=FALSE)
  dd.files <-  make_degree_day_files(degree.names,gcm,scenario,run,freq,
                                     tasmax.file,out.dir,tmp.dir)
  dd.ncs <- vector(mode='list',length=length(dd.files))
  for (d in seq_along(dd.files)) {
    dd.ncs[[d]] <- nc_open(paste0(out.dir,dd.files[d]),write=TRUE)
  }
  common.lat <- ncvar_get(dd.ncs[[1]],'lat')
  tas <- TRUE
  tasdiff <- FALSE

}

if (type=='gsl') {

  out.dir <- paste0(tmp.dir,'climdex/')
  dir.create(paste0(write.dir,'climdex/'),recursive=TRUE,showWarnings=FALSE)
  gsl.file <-  make_climdex_file('gsl',gcm,scenario,run,
                                 tasmax.file,out.dir,tmp.dir)
  gsl.ncs <- nc_open(paste0(out.dir,gsl.file),write=TRUE)
  common.lat <- ncvar_get(gsl.ncs,'lat')
  tas <- TRUE
  tasdiff <- FALSE

}

if (type=='dtr') {
  out.dir <- paste0(tmp.dir,'climdex/')
  dir.create(paste0(write.dir,'climdex/'),recursive=TRUE,showWarnings=FALSE)
  dtr.file <-  make_climdex_file('dtr',gcm,scenario,run,
                                 tasmax.file,out.dir,tmp.dir)
  dtr.ncs <- nc_open(paste0(out.dir,dtr.file),write=TRUE)
  common.lat <- ncvar_get(dtr.ncs,'lat')
  tas <- FALSE
  tasdiff <- TRUE
}



##---------------------------------------------------------------------------

##Iterate over the split files
tasmax.nc <- nc_open(paste0(tmp.dir,tasmax.file),write=FALSE)
tasmax.dates <- netcdf.calendar(tasmax.nc)
yearly.fac <- as.factor(format(tasmax.dates,'%Y'))
monthly.fac <- as.factor(format(tasmax.dates,'%Y-%m'))

print('TN Opening')
file.copy(paste0(gcm.dir,"/",tasmin.file),tmp.dir,overwrite=TRUE)
tasmin.nc <- nc_open(paste0(tmp.dir,tasmin.file),write=FALSE)

lon <- ncvar_get(tasmax.nc,'lon')
lat <- ncvar_get(tasmax.nc,'lat')
n.lon <- length(lon)
n.lat <- length(lat)

for (j in 1:n.lat) {
   lat.ix <- j
   print(paste0('Latitude: ',j,' of ',n.lat))
   tasmax.subset <- ncvar_get(tasmax.nc,'tasmax',start=c(1,j,1),count=c(-1,1,-1))
   tasmax.subset <- gcm_unit_conversion('tasmax',tasmax.subset,tasmax.nc) 

   tasmin.subset <- ncvar_get(tasmin.nc,'tasmin',start=c(1,j,1),count=c(-1,1,-1))    
   tasmin.subset <- gcm_unit_conversion('tasmin',tasmin.subset,tasmin.nc) 

   if (tas) {
     input.subset <- (tasmax.subset + tasmin.subset)/2
     input.list <- vector(mode='list',length=n.lon)
     input.list <- lapply(seq_len(nrow(input.subset)), function(k) input.subset[k,])
   }
   if (tasdiff) {
      input.subset <- tasmax.subset - tasmin.subset
      input.list <- vector(mode='list',length=n.lon)
      input.list <- lapply(seq_len(nrow(input.subset)), function(k) input.subset[k,])
   }
   flag <- is.na(input.subset[,1])
   rm(input.subset)
   rm(tasmax.subset)
   rm(tasmin.subset)

   ##----------------------------------------------------------
   if (type=='degree_days') {
     ##Degree Day 
     degree_days_for_model(degree.names,dd.ncs,lat.ix,n.lon,flag,
                           input.list,yearly.fac)
   }
   ##----------------------------------------------------------
   ##GSL
   if (type=='gsl') {
     gsl_for_model('gsl',gsl.ncs,lat.ix,n.lon,flag,
                   input.list,yearly.fac,tasmax.dates)
   }

   ##----------------------------------------------------------
   ##DTR 
   if (type=='dtr') {
     dtr_for_model('dtr',dtr.ncs,lat.ix,n.lon,flag,
                   input.list,monthly.fac)
   }
}##Latitude Loop
nc_close(tasmax.nc)
nc_close(tasmin.nc)

print('Removing lat band files')
file.remove(paste0(tmp.dir,"/",tasmax.file))
file.remove(paste0(tmp.dir,"/",tasmin.file))


##Move back

if (type=='degree_days') {
  for (d in seq_along(degree.names)) {
    nc_close(dd.ncs[[d]]) 
    file.copy(from=paste0(out.dir,dd.files[d]),to=paste0(write.dir,'degree_days/'),overwrite=TRUE)
  }
}

if (type=='gsl') {
  nc_close(gsl.ncs) 
  file.copy(from=paste0(out.dir,gsl.file),to=paste0(write.dir,'climdex/'),overwrite=TRUE)
}

if (type=='dtr') {
  nc_close(dtr.ncs) 
  file.copy(from=paste0(out.dir,dtr.file),to=paste0(write.dir,'climdex/'),overwrite=TRUE)
}




