##Script to calculate CSDI and WSDI from the CMIP6 GCMs
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
registerDoParallel(cores=1) # or some other number that you're comfortable with.

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
   gcm <- 'CanESM5'
   scenario <- 'ssp245'
   run <- 'r1i1p2f1'
   res <- NULL
} else {
   args <- commandArgs(trailingOnly=TRUE)
   for(i in 1:length(args)){
       eval(parse(text=args[[i]]))
   }
}

base <- c(1971,2000)

tmp.dir <- paste0(tmpdir,'/',gcm,'_',scenario,'_',run,'_climdex_spell_quantiles/')
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
rcp.files <- list.files(path=gcm.dir,pattern=paste0('historical\\+',scenario))

if (is.null(res)) {
   gcm.files <- rcp.files[grep(run,rcp.files)]
} else {
   run.files <- rcp.files[grep(run,rcp.files)]
   gcm.files <- run.files[grep(res,run.files)] ##Include for GCMs with multiple resolutions
}

if (length(gcm.files) !=3) {
   print(gcm.files)
   stop('There need to be 3 input files, one for each variable')
}

print('Started copying these files')
print(gcm.files)
tasmax.file <- gcm.files[grep('tasmax',gcm.files)]
tasmin.file <- gcm.files[grep('tasmin',gcm.files)]

file.copy(from=paste0(gcm.dir,tasmax.file),to=tmp.dir)
file.copy(from=paste0(gcm.dir,tasmin.file),to=tmp.dir)
print('Done copying gcm files')

##---------------------------------------------------------------------------
  ##Monthly Average Files for writing
  print('Monthly avg opening')

  out.dir <- paste0(tmp.dir,'climdex/')
  dir.create(paste0(write.dir,'climdex/'),recursive=TRUE,showWarnings=FALSE)
  wsdi.file <- make_climdex_file('wsdi',gcm,scenario,run,
                                   tasmax.file,out.dir,tmp.dir)
  wsdi.ncs <- nc_open(paste0(out.dir,wsdi.file),write=TRUE)

  csdi.file <- make_climdex_file('csdi',gcm,scenario,run,
                                   tasmin.file,out.dir,tmp.dir)
  csdi.ncs <- nc_open(paste0(out.dir,csdi.file),write=TRUE)

  print('Done with making Climdex files')
##---------------------------------------------------------------------------
##---------------------------------------------------------------------------

print('Data opening')
tasmax.nc <- nc_open(paste0(tmp.dir,tasmax.file),write=FALSE)
tasmax.dates <- netcdf.calendar(tasmax.nc)

tasmin.nc <- nc_open(paste0(tmp.dir,tasmin.file),write=FALSE)
tasmin.dates <- netcdf.calendar(tasmin.nc)

yearly.fac <- as.factor(format(tasmax.dates,'%Y'))
monthly.fac <- as.factor(format(tasmax.dates,'%Y-%m'))  

lon <- ncvar_get(tasmax.nc,'lon')
lat <- ncvar_get(tasmax.nc,'lat')
n.lon <- length(lon)
n.lat <- length(lat)

for (j in 1:n.lat) { 
   print(paste0('Latitude: ',j,' of ',n.lat))
   print('Extracting latitude band and converting to list')
   tasmax.subset <- ncvar_get(tasmax.nc,'tasmax',start=c(1,j,1),count=c(-1,1,-1))
   tasmax.subset <- gcm_unit_conversion('tasmax',tasmax.subset,tasmax.nc)
   tasmax.flag <- is.na(tasmax.subset[,1])
   flen <- length(tasmax.subset[,1])
   if (any(tasmax.flag)) {
      stop('NA values present in the TX GCM')
   }
   tasmax.list <- vector(mode='list',length=n.lon)
   tasmax.list <- lapply(seq_len(nrow(tasmax.subset)), function(k) tasmax.subset[k,])
   rm(tasmax.subset)

   tasmin.subset <- ncvar_get(tasmin.nc,'tasmin',start=c(1,j,1),count=c(-1,1,-1))
   tasmin.subset <- gcm_unit_conversion('tasmin',tasmin.subset,tasmin.nc)
   tasmin.flag <- is.na(tasmin.subset[,1])
   if (any(tasmin.flag)) {
      stop('NA values present in the TN GCM')
   }
   tasmin.list <- vector(mode='list',length=n.lon)
   tasmin.list <- lapply(seq_len(nrow(tasmin.subset)), function(k) tasmin.subset[k,])
   rm(tasmin.subset)

   n.col <- length(levels(yearly.fac))
   wsdi.matrix <- csdi.matrix <- matrix(NA,nrow=n.lon,ncol=n.col)

   ##----------------------------------------------------------
   print('Calculating Climdex Object')
   climdex.objects <- foreach(tasmax=tasmax.list,tasmin=tasmin.list,
                        .export=c('climdexInput.raw','tasmax.dates','tasmin.dates','base') 
                        ) %do% {
                           objects <- climdexInput.raw(tmax=tasmax,tmin=tasmin,
                                      tmax.dates=tasmax.dates,tmin.dates=tasmin.dates,base=base)
                        }
   rm(tasmax.list)   
   rm(tasmin.list)

   ##----------------------------------------------------------
   print('Calculating and writing the Wsdi index')
   wsdi.values <- foreach(obj=climdex.objects,.export='climdex.wsdi') %do% {
                           climdex.values <- climdex.wsdi(obj)}
   wsdi.sub.matrix <- matrix(unlist(wsdi.values),nrow=flen,ncol=n.col,byrow=TRUE)

   rm(wsdi.values)
   wsdi.matrix[!tasmax.flag,] <- wsdi.sub.matrix
   rm(wsdi.sub.matrix)      
   
   ncvar_put(wsdi.ncs,varid='wsdiETCCDI',vals=wsdi.matrix,
             start=c(1,j,1),count=c(-1,1,-1))
   rm(wsdi.matrix)

   ##---------------------------------------------------------
   print('Calculating and writing the Csdi index')
   csdi.values <- foreach(obj=climdex.objects,.export='climdex.csdi') %do% {
                           climdex.values <- climdex.csdi(obj)}
   csdi.sub.matrix <- matrix(unlist(csdi.values),nrow=flen,ncol=n.col,byrow=TRUE)
   rm(csdi.values)
   csdi.matrix[!tasmax.flag,] <- csdi.sub.matrix
   rm(csdi.sub.matrix)
   ncvar_put(csdi.ncs,varid='csdiETCCDI',vals=csdi.matrix,
             start=c(1,j,1),count=c(-1,1,-1))
   rm(csdi.matrix)


   rm(climdex.objects)
   gc()

}##Latitude Loop

nc_close(tasmax.nc)
nc_close(tasmin.nc)

nc_close(wsdi.ncs)
file.copy(from=paste0(out.dir,wsdi.file),to=paste0(write.dir,'climdex/'),overwrite=TRUE)

nc_close(csdi.ncs)
file.copy(from=paste0(out.dir,csdi.file),to=paste0(write.dir,'climdex/'),overwrite=TRUE)

file.remove(paste0(tmp.dir,tasmax.file))
file.remove(paste0(tmp.dir,tasmin.file))

