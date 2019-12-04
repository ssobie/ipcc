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

res <- NULL
if (testing) {
   tmpdir <- '/local_temp/ssobie'
   gcm <- 'BCC-CSM2-MR'
   scenario <- 'ssp585'
   run <- 'r1i1p1f1'
   res <- NULL
} else {
   args <- commandArgs(trailingOnly=TRUE)
   for(i in 1:length(args)){
       eval(parse(text=args[[i]]))
   }
}

base <- c(1971,2000)

tmp.dir <- paste0(tmpdir,'/',gcm,'_',scenario,'_',run,'_climdex_temp_quantiles/')
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
rcp.files <- list.files(path=gcm.dir,pattern=scenario)

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
  tx90p.file <- make_climdex_file('tx90p',gcm,scenario,run,
                                   tasmax.file,out.dir,tmp.dir)
  tx90p.ncs <- nc_open(paste0(out.dir,tx90p.file),write=TRUE)

  tx10p.file <- make_climdex_file('tx10p',gcm,scenario,run,
                                   tasmax.file,out.dir,tmp.dir)
  tx10p.ncs <- nc_open(paste0(out.dir,tx10p.file),write=TRUE)

  tn90p.file <- make_climdex_file('tn90p',gcm,scenario,run,
                                   tasmin.file,out.dir,tmp.dir)
  tn90p.ncs <- nc_open(paste0(out.dir,tn90p.file),write=TRUE)

  tn10p.file <- make_climdex_file('tn10p',gcm,scenario,run,
                                   tasmin.file,out.dir,tmp.dir)
  tn10p.ncs <- nc_open(paste0(out.dir,tn10p.file),write=TRUE)

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


   n.col <- length(levels(monthly.fac))
   tx90p.matrix <- tx10p.matrix <- tn90p.matrix <- tn10p.matrix <-matrix(NA,nrow=n.lon,ncol=n.col)

   ##----------------------------------------------------------
   print('Calculating Climdex Object')
   climdex.objects <- foreach(tasmax=tasmax.list,tasmin=tasmin.list,
                        .export=c('climdexInput.raw','tasmax.dates','tasmin.dates','base') 
                        ) %dopar% {
                           objects <- climdexInput.raw(tmax=tasmax,tmin=tasmin,
                                      tmax.dates=tasmax.dates,tmin.dates=tasmin.dates,base=base)
                        }
   rm(tasmax.list)   
   rm(tasmin.list)

   ##----------------------------------------------------------
   print('Calculating and writing the TX90p index')
   tx90p.values <- foreach(obj=climdex.objects,.export='climdex.tx90p') %dopar% {
                           climdex.values <- climdex.tx90p(obj)}
   tx90p.sub.matrix <- matrix(unlist(tx90p.values),nrow=flen,ncol=n.col,byrow=TRUE)
   rm(tx90p.values)
   tx90p.matrix[!tasmax.flag,] <- tx90p.sub.matrix
   rm(tx90p.sub.matrix)      
   ncvar_put(tx90p.ncs,varid='tx90pETCCDI',vals=tx90p.matrix,
             start=c(1,j,1),count=c(-1,1,-1))
   rm(tx90p.matrix)

   ##---------------------------------------------------------
   print('Calculating and writing the TX10p index')
   tx10p.values <- foreach(obj=climdex.objects,.export='climdex.tx10p') %dopar% {
                           climdex.values <- climdex.tx10p(obj)}
   tx10p.sub.matrix <- matrix(unlist(tx10p.values),nrow=flen,ncol=n.col,byrow=TRUE)
   rm(tx10p.values)
   tx10p.matrix[!tasmax.flag,] <- tx10p.sub.matrix
   rm(tx10p.sub.matrix)
   ncvar_put(tx10p.ncs,varid='tx10pETCCDI',vals=tx10p.matrix,
             start=c(1,j,1),count=c(-1,1,-1))
   rm(tx10p.matrix)

   ##---------------------------------------------------------
   print('Calculating and writing the TN90p index')
   tn90p.values <- foreach(obj=climdex.objects,.export='climdex.tn90p') %dopar% {
                           climdex.values <- climdex.tn90p(obj)}
   tn90p.sub.matrix <- matrix(unlist(tn90p.values),nrow=flen,ncol=n.col,byrow=TRUE)
   rm(tn90p.values)
   tn90p.matrix[!tasmax.flag,] <- tn90p.sub.matrix
   rm(tn90p.sub.matrix)
   ncvar_put(tn90p.ncs,varid='tn90pETCCDI',vals=tn90p.matrix,
             start=c(1,j,1),count=c(-1,1,-1))
   rm(tn90p.matrix)

   ##---------------------------------------------------------
   print('Calculating and writing the TN10p index')
   tn10p.values <- foreach(obj=climdex.objects,.export='climdex.tn10p') %dopar% {
                           climdex.values <- climdex.tn10p(obj)}
   tn10p.sub.matrix <- matrix(unlist(tn10p.values),nrow=flen,ncol=n.col,byrow=TRUE)
   rm(tn10p.values)
   tn10p.matrix[!tasmax.flag,] <- tn10p.sub.matrix
   rm(tn10p.sub.matrix)
   ncvar_put(tn10p.ncs,varid='tn10pETCCDI',vals=tn10p.matrix,
             start=c(1,j,1),count=c(-1,1,-1))
   rm(tn10p.matrix)

   rm(climdex.objects)
   gc()

}##Latitude Loop

nc_close(tasmax.nc)
nc_close(tasmin.nc)


nc_close(tx90p.ncs)
file.copy(from=paste0(out.dir,tx90p.file),to=paste0(write.dir,'climdex/'),overwrite=TRUE)

nc_close(tx10p.ncs)
file.copy(from=paste0(out.dir,tx10p.file),to=paste0(write.dir,'climdex/'),overwrite=TRUE)

nc_close(tn90p.ncs)
file.copy(from=paste0(out.dir,tn90p.file),to=paste0(write.dir,'climdex/'),overwrite=TRUE)

nc_close(tn10p.ncs)
file.copy(from=paste0(out.dir,tn10p.file),to=paste0(write.dir,'climdex/'),overwrite=TRUE)


file.remove(paste0(tmp.dir,tasmax.file))
file.remove(paste0(tmp.dir,tasmin.file))

