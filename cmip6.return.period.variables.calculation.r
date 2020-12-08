##Script to calculate and write the standard set of derived variables
##for the 800m data

##These three can remain the same, no dependence on input data
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

##----

library(ncdf4)
library(PCICt)
library(extRemes)
library(ismev)

##--------------------------------------------------------------

make_return_period_file <- function(varname,gcm,scenario,run,rp,interval,deg=NULL,
                                    gcm.file,tmp.dir) {

   gcm.nc <- nc_open(paste0(tmp.dir,gcm.file),write=TRUE)
   gcm.series <- netcdf.calendar(gcm.nc)
   gcm.time <- ncvar_get(gcm.nc,'time')
   gcm.cal <- attr(gcm.series,'cal')
   nc_close(gcm.nc)
   
   bnds <- strsplit(interval,'-')[[1]]
   yix <- gcm.series >= as.PCICt(paste0(bnds[1],'-01-01'),cal=gcm.cal) &
          gcm.series <= as.PCICt(paste0(bnds[2],'-12-30'),cal=gcm.cal)

   rp.one <- gsub(paste0(varname,'_'),paste0(varname,'_RP',rp,'_'),gcm.file)
   if (grepl('one|two|three',deg)) {
      rp.deg <- gsub(pattern=run,replacement=paste0(run,'_',deg,'_deg'),rp.one)
      rp.file <- gsub(pattern='[0-9]{4}-[0-9]{4}',replacement=interval,rp.deg)
   } else {
      rp.file <- gsub(pattern='[0-9]{4}-[0-9]{4}',replacement=interval,rp.one)
   }
   work <- paste0('cdo -O timmean ',tmp.dir,gcm.file,' ',tmp.dir,rp.file)
   system(work)

   Sys.sleep(1)
   print('Created empty RP file')
   nc <- nc_open(paste0(tmp.dir,rp.file),write=TRUE)

   ##Fix the time step
   gcm.date <- gcm.time[yix][1]  
   ncvar_put(nc,'time',gcm.date)   
   stnd_name <- ncatt_get(nc,varname,'standard_name')$value

   ncatt_put(nc,varname,'standard_name',paste0('RP',rp,' for ',stnd_name))
   long_name <- ncatt_get(nc,varname,'long_name')$value
   ncatt_put(nc,varname,'long_name',paste0('RP',rp,' for ',long_name))
   ncatt_put(nc,0,'history','')

   nc_close(nc)

   return(list(file=rp.file,index=yix))

}

##--------------------------------------------------------------

calc_return_periods <- function(ts.yearly,var.name,rperiod) {

  if (sum(is.na(ts.yearly)) == length(ts.yearly)) {
    return(NA)
  } else {

    inf.flag <- is.infinite(ts.yearly)
    ts.yearly[inf.flag] <- NA

    na.flag <- is.na(ts.yearly)
    if(sum(na.flag)>0) {
      ts.to.fit <- as.vector(ts.yearly[-which(na.flag)])
    } else {
      ts.to.fit <- as.vector(ts.yearly)
    }
    if (var.name=='tasmin') {
      ts.to.fit <- -ts.to.fit
      u.len <- length(unique(ts.to.fit))
      f.len <- length(ts.to.fit)
      if (u.len < (f.len*2/3))
        ts.to.fit <- jitter(ts.to.fit,amount=3)
    }

###    ts.fit <- fevd(ts.to.fit,type='GEV')
###    ts.old <- gev.fit(ts.to.fit,show=FALSE)
###    ts.fit$results$par <- ts.old$mle
###    names(ts.fit$results$par) <- c('location','scale','shape')
###    ts.rps <- return.level(ts.fit,return.period=as.numeric(rperiod),make.plot=F)

    ##Include Try-Catch for errors
if (1==1) {
    ts.rps <- tryCatch({
       ts.fit <- fevd(ts.to.fit,type='GEV')
       ts.old <- gev.fit(ts.to.fit,show=FALSE)
       ts.fit$results$par <- ts.old$mle
       names(ts.fit$results$par) <- c('location','scale','shape')
       ts.rps <- return.level(ts.fit,return.period=as.numeric(rperiod),make.plot=F)
       return(ts.rps)
    }, warning=function(war){
       ts.rps <- NA
       return(ts.rps)
       ##message('Warning from fevd fit')
       ##browser()
       ##ts.fit <- gev.fit(ts.to.fit,show=FALSE)
       ts.rps <- NA ##return.level(ts.fit,return.period=as.numeric(rperiod),make.plot=F)
    }, error=function(err){
       print('Error in fevd')
       ###browser()
       ###ts.fit <- gev.fit(ts.to.fit,show=FALSE)
       ts.rps <- NA ###return.level(ts.fit,return.period=as.numeric(rperiod),make.plot=F)
       return(ts.rps)
    },finally={
    })

}

    rv <- as.numeric(ts.rps)

    if (var.name=='tasmin')
      rv <- -as.numeric(ts.rps) ##$return.level
    return(rv)
  }
}

##--------------------------------------------------------------

check_rp_outliers <- function(data,var.name) {
  rv <- data
  if (var.name=='tasmax')
    flags <- which(data > 75)
  if (var.name=='tasmin')
    flags <- which(data < -75)
  if (var.name=='pr')
    flags <- which(data > 1500)
  if (length(flags)!=0) {
    rv[flags] <- NA
  }
  return(rv)
}



##--------------------------------------------------------------

return_periods_for_model <- function(varname,rp.ncs,lat.ix,n.lon,flag,rperiod,
                                     ann.list) {
   ##Variables
       flen <- sum(!flag)
       rp.vector <- rep(NA,length=n.lon)

       if (flen!=0) { ##Some Real Values
          sub.list <- ann.list[!flag]
          rp.vals <- rep(0,flen)
          for (i in 1:flen) {
             ts.yearly <- sub.list[[i]]
             rp.vals[i] <- calc_return_periods(ts.yearly,varname,rperiod)            
          }
          rp.checked <- check_rp_outliers(rp.vals,varname)
          ##print(rp.checked)
          rp.vector[!flag] <- rp.checked
          rm(rp.vals)
          rm(sub.list)
      } else {
         print('All NA values')
      }

      ncvar_put(rp.ncs,varid=varname,vals=rp.vector,
                start=c(1,lat.ix,1),count=c(-1,1,1))
      rm(rp.vector)
      gc()
} 



##--------------------------------------------------------------

##****************************************************************
testing <- FALSE

res <- NULL
if (testing) {
   tmpdir <- '/local_temp/ssobie'
   gcm <- 'IPSL-CM6A-LR'
   scenario <- 'ssp245'
   run <- 'r1i1p1f1'
   res <- NULL
   rp  <- 5
   varname <- 'pr'
   interval <- '1971-2000'
   deg <- "NULL"
   base.dir <- '/storage/data/climate/CMIP6/Derived/'

   gcm.dir <- paste0(base.dir,gcm,'_',scenario,'_',run,'/annual_extremes/')
   gcm.files <- list.files(path=gcm.dir,pattern=paste0('_',gcm,'_'))
   run.files <-  gcm.files[grep(run,gcm.files)]
   scen.files <- run.files[grep(scenario,run.files)]
   gcm.file <- scen.files[grep(varname,scen.files)]
   write.dir <- paste0(base.dir,gcm,'_',scenario,'_',run,'/return_periods/') ###degree_anomaly_rps/')

} else {
   args <- commandArgs(trailingOnly=TRUE)
   for(i in 1:length(args)){
       eval(parse(text=args[[i]]))
   }

   gcm.file <- gcmfile
   gcm.dir <- gcmdir
   write.dir <- writedir
}

##-----------------------------------------------

tmp.dir <- paste0(tmpdir,'/',gcm,'_',scenario,'_',run,'_return_period_',rp,'_',varname,'/')
if (!file.exists(tmp.dir)) {
  dir.create(tmp.dir,recursive=TRUE)
}

if (!file.exists(write.dir)) {
  dir.create(write.dir,recursive=TRUE)
}

file.copy(from=paste0(gcm.dir,gcm.file),to=tmp.dir)
print('Done copying annual extremes file')


##---------------------------------------------------------------------------
##Return Period Files for writing
print('Return Period opening')
print('Return Period opening')

rp.file <-  make_return_period_file(varname,gcm,scenario,run,rp,interval,deg,
                                    gcm.file,tmp.dir)

###file.copy(from=paste0(tmp.dir,rp.file),to=write.dir,overwrite=TRUE)

rp.ncs <- nc_open(paste0(tmp.dir,rp.file$file),write=TRUE)


##---------------------------------------------------------------------------
##---------------------------------------------------------------------------


print('Data opening')
gcm.nc <- nc_open(paste0(tmp.dir,gcm.file),write=FALSE)
gcm.dates <- netcdf.calendar(gcm.nc)

lon <- ncvar_get(gcm.nc,'lon')
lat <- ncvar_get(gcm.nc,'lat')
n.lon <- length(lon)
n.lat <- length(lat)

for (j in 1:n.lat) { 
   print(paste0('Latitude: ',j,' of ',n.lat))
   lat.ix <- j
   gcm.subset <- ncvar_get(gcm.nc,varname,start=c(1,j,1),count=c(-1,1,-1))[,rp.file$index]
   flag <- is.na(gcm.subset[,1])
   ann.list <- lapply(seq_len(nrow(gcm.subset)), function(k) gcm.subset[k,])
   rm(gcm.subset)

   return_periods_for_model(varname,rp.ncs,lat.ix,n.lon,flag,rp,
                            ann.list)      
   rm(ann.list)
}##Latitude Loop

nc_close(gcm.nc)
nc_close(rp.ncs)

file.copy(from=paste0(tmp.dir,rp.file$file),to=write.dir,overwrite=TRUE)

file.remove(paste0(tmp.dir,gcm.file))
file.remove(paste0(tmp.dir,rp.file$file))
