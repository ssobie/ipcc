##Script to calculate and write the standard set of derived variables
##for the 800m data

##These three can remain the same, no dependence on input data
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)
source('/home/ssobie/code/repos/downscale_CMIP6/bccaqv2.simplified.climdex.support.r',chdir=T)
source('/home/ssobie/code/repos/downscale_CMIP6/climdex.variables.functions.r',chdir=T)

##----

source('/home/ssobie/code/repos/ipcc/cmip6.canada.derived.files.r')

library(udunits2)
library(ncdf4)
library(PCICt)
library(climdex.pcic)
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

res <- 'gr1'
if (testing) {
   tmpdir <- '/local_temp/ssobie'
   gcm <- 'BCC-CSM2-MR'
   scenario <- 'ssp585'
   run <- 'r1i1p1f1'
   res <- NULL
   type <- 'annual'
   climname <- 'prcptot'
} else {
   args <- commandArgs(trailingOnly=TRUE)
   for(i in 1:length(args)){
       eval(parse(text=args[[i]]))
   }
}

varname <- input.varname[[climname]]
climdex.name <- paste0('climdex.',climname)
climdex.info <- get.climdex.info(climdex.name)

tmp.dir <- paste0(tmpdir,'/',gcm,'_',scenario,'_',run,'_climdex_',type,'_',varname,'/')
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
if (type=='annual') {
  ##Annual Average Files for writing
  print('Ann Avg opening')
  print('Ann Avg opening')
  out.dir <- paste0(tmp.dir,'climdex/')
  dir.create(paste0(write.dir,'climdex/'),recursive=TRUE,showWarnings=FALSE)
  ann.file <-  make_climdex_file(climname,gcm,scenario,run,
                                 gcm.file,out.dir,tmp.dir)
  ann.ncs <- nc_open(paste0(out.dir,ann.file),write=TRUE)
  common.lat <- ncvar_get(ann.ncs,'lat')
}

##---------------------------------------------------------------------------
if (type=='monthly') {
  ##Monthly Average Files for writing
  print('Monthly avg opening')

  out.dir <- paste0(tmp.dir,'climdex/')
  dir.create(paste0(write.dir,'climdex/'),recursive=TRUE,showWarnings=FALSE)
  mon.file <- make_climdex_file(climname,gcm,scenario,run,
                                gcm.file,out.dir,tmp.dir)
  mon.ncs <- nc_open(paste0(out.dir,mon.file),write=TRUE)
  common.lat <- ncvar_get(mon.ncs,'lat')
}

if (type=='r9') {
  ##R95 Files for writing
  print('R9 opening')
  ann.names <- c(paste0(climname,'p'),
                 paste0(climname,'days'),
                 paste0(climname,'store'))
  ann.ncs <- vector(mode='list',length=length(ann.names))
  ann.files <- rep('A',length(ann.names))
  out.dir <- paste0(tmp.dir,'climdex/')
  dir.create(paste0(write.dir,'climdex/'),recursive=TRUE,showWarnings=FALSE)

  for (d in seq_along(ann.names)) {
     ann.files[d] <- make_climdex_file(ann.names[d],gcm,scenario,run,
                                   gcm.file,out.dir,tmp.dir) 
     ann.ncs[[d]] <- nc_open(paste0(out.dir,ann.files[d]),write=TRUE)  
  }
  names(ann.ncs) <- c('total','days','store')
  common.lat <- ncvar_get(ann.ncs[[1]],'lat')
}



##---------------------------------------------------------------------------
##---------------------------------------------------------------------------

print('Data opening')
gcm.nc <- nc_open(paste0(tmp.dir,gcm.file),write=FALSE)
gcm.dates <- netcdf.calendar(gcm.nc)
yearly.fac <- as.factor(format(gcm.dates,'%Y'))
monthly.fac <- as.factor(format(gcm.dates,'%Y-%m'))  

lon <- ncvar_get(gcm.nc,'lon')
lat <- ncvar_get(gcm.nc,'lat')
n.lon <- length(lon)
n.lat <- length(lat)

for (j in 1:n.lat) { 
   print(paste0('Latitude: ',j,' of ',n.lat))
   lat.ix <- j
   gcm.subset <- ncvar_get(gcm.nc,varname,start=c(1,j,1),count=c(-1,1,-1))
   gcm.subset <- gcm_unit_conversion(varname,gcm.subset,gcm.nc)

   flag <- is.na(gcm.subset[,1])
   gcm.list <- vector(mode='list',length=n.lon)
   gcm.list <- lapply(seq_len(nrow(gcm.subset)), function(k) gcm.subset[k,])
   rm(gcm.subset)

    ##----------------------------------------------------------
    ##Annual Averages 
    if (type=='annual') {
       annual_climdex_for_model(climname,ann.ncs,lat.ix,n.lon,yearly.fac,flag,
                                gcm.list)
    }
    ##----------------------------------------------------------
    ##Monthly Averages 
    if (type=='monthly') {
       monthly_climdex_for_model(climname,mon.ncs,lat.ix,n.lon,monthly.fac,flag,
                                 gcm.list)
    }
    ##----------------------------------------------------------
    ##Extreme Precip
    if (type=='r9') {
       r9_precip_for_model(climname,ann.names,ann.ncs,lat.ix,n.lon,yearly.fac,gcm.dates,flag,
                           gcm.list)
    }

  rm(gcm.list)
}##Latitude Loop
nc_close(gcm.nc)

if (type=='annual') {
   nc_close(ann.ncs)
   file.copy(from=paste0(out.dir,ann.file),to=paste0(write.dir,'climdex/'),overwrite=TRUE)
}

if (type=='monthly') {
  nc_close(mon.ncs)
  file.copy(from=paste0(out.dir,mon.file),to=paste0(write.dir,'climdex/'),overwrite=TRUE)
}

if (type=='r9') {
  for (d in seq_along(ann.names)) {
    nc_close(ann.ncs[[d]])  
    file.copy(from=paste0(out.dir,ann.files[d]),to=paste0(write.dir,'climdex/'),overwrite=TRUE)
  }
}


file.remove(paste0(tmp.dir,gcm.file))
