#Calculates the regional average climatology value at the intervals
#when the global average temperature has increased by one degree
#increments.



library(rgdal)
library(raster)
library(zoo)
library(ncdf4)

##----------------------------------------------------------
##Read a shapefile

get_shape_file <- function(region,shape.dir) {
   ##shape.dir <- '/storage/data/projects/rci/data/assessments/shapefiles/bc_common/'
   reg.shp <- readOGR(shape.dir,region,stringsAsFactors=F)
   return(reg.shp)
}

##----------------------------------------------------------
## Find the overlapping cells
##This gets the masked values including any overlapping cells
##Normal mask omits the partially covered cells
## **This step takes a while so it should be done once, then
## saved for future use **
##Use a spatial copy of the file with one time step for speed

cell_mask <- function(file.one,shape) {

   ras <- rasterize(shape,file.one,getCover=T)
   ras[ras==0] <- NA
   return(ras)
}

##----------------------------------------------------------
##Find the area weighting for the average

area_weighting <- function(file.brick,area.shape) {

   brick.one <- subset(file.brick,10)
   area.overlay <- cell_mask(brick.one,area.shape)
   area.masked <- mask(brick.one,area.overlay)
   cc <- area(area.masked) ##Returns area including blank regions
   cc.masked <- mask(cc,area.overlay)
   area.weights <- cc.masked / cellStats(cc.masked,sum)

   rv <- list(weights=area.weights,mask=area.overlay)
   return(rv)
}

##----------------------------------------------------------
area_average_time_series <- function(file.brick,area.cover) {

   area.series <- mask(file.brick,area.cover$mask)
   area.weighted.series <- area.series * area.cover$weights
   area.time.series <- cellStats(area.weighted.series,sum)
   return(area.time.series)
}

##---------------------------------------------------------
##Create time series
create_time_series <- function(tasmax.file,tasmin.file,tmp.dir,map.shape=NULL) {
   tasmax.fix <- gsub("tasmax_","tasmax_180_",tasmax.file)
   work <- paste0("cdo sellonlatbox,-180,180,-90,90 ",tmp.dir,tasmax.file," ",tmp.dir,tasmax.fix)
   system(work)
   Sys.sleep(1)
   tasmin.fix <- gsub("tasmin_","tasmin_180_",tasmin.file)
   work <- paste0("cdo sellonlatbox,-180,180,-90,90 ",tmp.dir,tasmin.file," ",tmp.dir,tasmin.fix)
   system(work)
   Sys.sleep(1)

   tasmax <- brick(paste0(tmp.dir,tasmax.fix))
   tasmin <- brick(paste0(tmp.dir,tasmin.fix))
   tas <- (tasmax + tasmin) / 2

   area.cover <- area_weighting(tas,map.shape)
   tas.series <- area_average_time_series(tas,area.cover)

   dates <- tasmax@z$Date
   tas.zoo <- zoo(tas.series,dates)
   tas.roll <- rollmean(tas.zoo,31)

   file.remove(paste0(tmp.dir,tasmax.fix))
   file.remove(paste0(tmp.dir,tasmin.fix))        
   return(list(roll=tas.roll,raw=tas.zoo))
}

##----------------------------------------------------------
##Find the area weighting for the average

global_weighting <- function(file.brick) {

   brick.one <- subset(file.brick,10)
   cc <- area(brick.one)
   area.weights <- cc / cellStats(cc,sum)
   rv <- area.weights
   return(rv)
}

##----------------------------------------------------------
global_average_time_series <- function(file.brick,area.weights) {

   weighted.series <- file.brick * area.weights
   time.series <- cellStats(weighted.series,sum)
   return(time.series)
}

##---------------------------------------------------------
##Create time series
global_time_series <- function(tasmax.file,tasmin.file,tmp.dir) {
   tasmax <- brick(paste0(tmp.dir,tasmax.file))
   tasmin <- brick(paste0(tmp.dir,tasmin.file))
   tas <- (tasmax + tasmin) / 2

   area.cover <- global_weighting(tas)
   tas.series <- global_average_time_series(tas,area.cover)

   dates <- tasmax@z$Date
   tas.zoo <- zoo(tas.series,dates)
   tas.roll <- rollmean(tas.zoo,31)
   return(list(roll=tas.roll,raw=tas.zoo))
}


##---------------------------------------------------------

find_anomaly_timing <- function(tas.series,base,anoms) {

   tas.raw <- tas.series$raw
   raw.dates <- index(tas.raw)
   raw.years <- format(raw.dates,'%Y')
   raw.series <- as.numeric(tas.raw)

   rst <- head(grep(base[1],raw.years),1)
   ren <- tail(grep(base[2],raw.years),2)
   raw.baseline <- mean(raw.series[rst:ren])
   raw.anomalies <- raw.series - raw.baseline

   tas.roll <- tas.series$roll
   dates <- index(tas.roll)
   years <- format(dates,'%Y')
   series <- as.numeric(tas.roll)

   bst <- head(grep(base[1],years),1)
   ben <- tail(grep(base[2],years),2)
   if (length(bst)==0) {browser()}
   if (length(ben)==0) {browser()}
   baseline <- mean(series[bst:ben])
   anomalies <- series - baseline

   anom.years <- rep(0,length(anoms))
   anom.values <- rep(0,length(anoms))

   for (i in seq_along(anoms)) {
      ix <- which.min(abs(anoms[i] - anomalies))
      anom.years[i] <- years[ix]

      rx <- which(years[ix] == raw.years)      
      anom.values[i] <- mean(raw.anomalies[(rx-15):(rx+15)])
      if (is.na(anom.values[i])) { browser()}
   }   

   flag <- anoms - max(anomalies,na.rm=T) > 0.5
   anom.years[flag] <- NA
   anom.values[flag] <- NA
   flag <- as.numeric(anom.years) > 2100
   anom.years[flag] <- NA
   anom.values[flag] <- NA

##   lines(as.numeric(years),anomalies)
##   points(anom.years,anom.values,pch=18)

   return(list(years=anom.years,values=anom.values,anoms=raw.anomalies,rolled=anomalies))
}


anomaly_timing_climatology <- function(tas.roll,anoms,base) {

   series <- tas.roll$raw
   years <- format(index(tas.roll$raw),'%Y')
   bst <- head(grep(base[1],years),1)
   ben <- tail(grep(base[2],years),2)
   if (length(bst)==0) {browser()}
   if (length(ben)==0) {browser()}
   baseline <- mean(series[bst:ben])
   anomalies <- series - baseline

   clims <- rep(NA,length(anoms))
   for (i in seq_along(anoms)) {
      if (is.na(anoms[i])) { 
         clims[i] <- NA
      } else {
         ix <- which(anoms[i] == years)
         if (length(ix) == 0) {browser()}
         if (ix+15 > length(anomalies)) {browser()}
         clims[i] <- mean(anomalies[(ix-15):(ix+15)])
      }
   }   

   return(clims)
}

##---------------------------------------------------------

experiment <- 'cmip5'
scenario <- 'rcp85' ##'ssp585' ##

###data.dir <- paste0("/storage/data/projects/rci/data/cas/",experiment,"/one_degree/")
if (experiment=='cmip5') {
   data.dir <- "/storage/data/climate/CMIP5/daily/Derived/"
} 
if (experiment=='cmip6') {
   data.dir <- "/storage/data/climate/CMIP6/Derived/"
}

write.dir <- paste0("/storage/data/projects/rci/data/cas/",experiment,"/degree_climatologies/")

tmp.dir <- "/local_temp/ssobie/cmip_temp_timing/"

if (!file.exists(tmp.dir)) {
   dir.create(tmp.dir,recursive=T)
}



cmip5.list <- list(c('ACCESS1-0','r1i1p1'),
                   c('bcc-csm1-1','r1i1p1'),
                   c('bcc-csm1-1-m','r1i1p1'),
                   c('BNU-ESM','r1i1p1'),
                   c('CanESM2','r1i1p1'),
                   c('CCSM4','r2i1p1'),
                   c('CESM1-CAM5','r1i1p1'),                  
                   c('CNRM-CM5','r1i1p1'),
                   c('CSIRO-Mk3-6-0','r1i1p1'),
                   c('FGOALS-g2','r1i1p1'),
                   c('GFDL-CM3','r1i1p1'),
                   c('GFDL-ESM2G','r1i1p1'),
                   c('GFDL-ESM2M','r1i1p1'),
                   c('HadGEM2-AO','r1i1p1'),
                   c('HadGEM2-CC','r1i1p1'),
                   c('HadGEM2-ES','r1i1p1'),
                   c('inmcm4','r1i1p1'),
                   c('IPSL-CM5A-LR','r1i1p1'),
                   c('IPSL-CM5A-MR','r1i1p1'),
                   c('MIROC-ESM','r1i1p1'),
                   c('MIROC-ESM-CHEM','r1i1p1'),
                   c('MIROC5','r3i1p1'),
                   c('MPI-ESM-LR','r3i1p1'),
                   c('MPI-ESM-MR','r1i1p1'),
                   c('MRI-CGCM3','r1i1p1'),
                   c('NorESM1-M','r1i1p1'))

cmip6.list <- list(c('ACCESS-CM2','r1i1p1f1'),
                   c('ACCESS-ESM1-5','r1i1p1f1'),
                   c('BCC-CSM2-MR','r1i1p1f1'),
                   c('CanESM5','r1i1p2f1'),
                   c('CNRM-CM6-1','r1i1p1f2'),
                   c('CNRM-ESM2-1','r1i1p1f2'),
                   c('EC-Earth3','r4i1p1f1'),
                   c('EC-Earth3-Veg','r1i1p1f1'),
                   c('FGOALS-g3','r1i1p1f1'),
                   c('GFDL-ESM4','r1i1p1f1'),
                   c('HadGEM3-GC31-LL','r1i1p1f3'),
                   c('INM-CM4-8','r1i1p1f1'),
                   c('INM-CM5-0','r1i1p1f1'),
                   c('IPSL-CM6A-LR','r1i1p1f1'),
                   c('KACE-1-0-G','r2i1p1f1'),
                   c('KIOST-ESM','r1i1p1f1'),
                   c('MIROC6','r1i1p1f1'),
                   c('MIROC-ES2L','r1i1p1f2'),
                   c('MPI-ESM1-2-HR','r1i1p1f1'),
                   c('MPI-ESM1-2-LR','r1i1p1f1'),
                   c('MRI-ESM2-0','r1i1p1f1'),
                   c('NESM3','r1i1p1f1'),
                   c('NorESM2-LM','r1i1p1f1'),
                   c('NorESM2-MM','r1i1p1f1'),
                   c('UKESM1-0-LL','r1i1p1f2'))


pcic12.list <- list(c('ACCESS1-0','r1i1p1'),
                   c('CanESM2','r1i1p1'),
                   c('CCSM4','r2i1p1'),
                   c('CNRM-CM5','r1i1p1'),
                   c('CSIRO-Mk3-6-0','r1i1p1'),
                   c('GFDL-ESM2G','r1i1p1'),
                   c('HadGEM2-CC','r1i1p1'),
                   c('HadGEM2-ES','r1i1p1'),
                   c('inmcm4','r1i1p1'),
                   c('MIROC5','r3i1p1'),
                   c('MPI-ESM-LR','r3i1p1'),
                   c('MRI-CGCM3','r1i1p1'))


gcm.list <- switch(experiment,
                   pcic12=pcic12.list,
                   cmip5=cmip5.list,
                   cmip6=cmip6.list)

shape.dir <- '/storage/data/projects/rci/data/assessments/shapefiles/bc_common/'
region <- 'canada_boundary'
map.shape <- get_shape_file(region,shape.dir)

anoms <- c(1,2,3,4)
base <- c(1971,2000)

raw.anomalies <- vector(mode='list',length=length(gcm.list))
rolled.anomalies <- vector(mode='list',length=length(gcm.list))

anomaly.years <- matrix(0,nrow=length(gcm.list),ncol=length(anoms))
anomaly.values <- matrix(0,nrow=length(gcm.list),ncol=length(anoms))
tas.clims <- shp.clims <- matrix(0,nrow=length(gcm.list),ncol=length(anoms))

##--------------------------------------------------------------------
##Global Anomalies at GCM native resolution

for (j in seq_along(gcm.list)) {

   gcm.info <- gcm.list[[j]]
   gcm <- gcm.info[1]
   run <- gcm.info[2]

   print(gcm)
   gcm.dir <- paste0(data.dir,gcm,'_',scenario,'_',run,'/annual/')
   
   tasmax.file <- list.files(path=gcm.dir,pattern="tasmax_annual_average_")
   tasmin.file <- list.files(path=gcm.dir,pattern="tasmin_annual_average_")

   if (length(tasmax.file) > 1 | length(tasmax.file) == 0) {browser()}
   if (length(tasmin.file) > 1 | length(tasmin.file) == 0) {browser()}

   print(tasmax.file)
   print(tasmin.file)

   file.copy(from=paste0(gcm.dir,tasmax.file),to=tmp.dir,overwrite=TRUE)
   file.copy(from=paste0(gcm.dir,tasmin.file),to=tmp.dir,overwrite=TRUE)
   print('Done copying')
   
   ##Deal with irregularly spaced latitude values near the poles
   if (grepl('(FGOALS-g2|FGOALS-g3|GFDL-ESM2G|GFDL-ESM2M)',gcm)) { 
      nc <- nc_open(paste0(tmp.dir,tasmax.file))
      nlon <- nc$dim$lon$len
      nlat <- nc$dim$lat$len
      nc_close(nc)
      tasmax.reg <- gsub("tasmax_","tasmax_reg_",tasmax.file)
      work <- paste0("cdo remapcon,r",nlon,"x",nlat," ",tmp.dir,tasmax.file," ",tmp.dir,tasmax.reg)
      system(work)
      Sys.sleep(1)

      tasmin.reg <- gsub("tasmin_","tasmin_reg_",tasmin.file)
      work <- paste0("cdo remapcon,r",nlon,"x",nlat," ",tmp.dir,tasmin.file," ",tmp.dir,tasmin.reg)
      system(work)
      Sys.sleep(1)
      tas.roll <- global_time_series(tasmax.reg,tasmin.reg,tmp.dir)
      shp.roll <- create_time_series(tasmax.reg,tasmin.reg,tmp.dir,map.shape=map.shape) 
      file.remove(paste0(tmp.dir,tasmax.reg))
      file.remove(paste0(tmp.dir,tasmin.reg))
   } else {
       tas.roll <- global_time_series(tasmax.file,tasmin.file,tmp.dir)
       shp.roll <- create_time_series(tasmax.file,tasmin.file,tmp.dir,map.shape=map.shape) 
   }
   tas.timing <- find_anomaly_timing(tas.roll,base,anoms)
   anomaly.years[j,] <- tas.timing$years
   anomaly.values[j,] <- tas.timing$values
   raw.anomalies[[j]] <- tas.timing$anoms
   rolled.anomalies[[j]] <- tas.timing$rolled

   tas.clims[j,] <- anomaly_timing_climatology(tas.roll,tas.timing$years,base)
   shp.clims[j,] <- anomaly_timing_climatology(shp.roll,tas.timing$years,base)

   file.remove(paste0(tmp.dir,tasmax.file))
   file.remove(paste0(tmp.dir,tasmin.file))

}

gcms <- unlist(lapply(gcm.list,function(x){x[1]}))
runs <- unlist(lapply(gcm.list,function(x){x[2]}))
rv <- rbind(c('Model','Run','1Degrees','2Degrees','3Degrees','4Degrees'),cbind(gcms,runs,round(shp.clims,2)))

write.file <- paste0(write.dir,"canada_",experiment,"_",scenario,"_tas_degree_anomaly_values.csv")
write.table(rv,file=write.file,sep=',',quote=F,row.name=F,col.name=F)





##------------------------------------------------------------
##Temperature anomaly timing using a shape file to clip the
##model data firrst
if (1==0) {

##plot(c(),xlim=c(1850,2100),ylim=c(-2,8),xlab='Year',ylab='Temp Anomaly',main='Baseline 1971-2000')

for (j in seq_along(gcm.list)) {

   gcm.info <- gcm.list[[j]]
   gcm <- gcm.info[1]
   run <- gcm.info[2]

   print(gcm)
   gcm.dir <- paste0(data.dir,gcm,'_',scenario,'_',run,'/annual/')
   tasmax.file <- list.files(path=gcm.dir,pattern=paste0("tasmax_annual_average_",gcm))
   tasmin.file <- list.files(path=gcm.dir,pattern=paste0("tasmin_annual_average_",gcm))

   file.copy(from=paste0(gcm.dir,tasmax.file),to=tmp.dir,overwrite=TRUE)
   file.copy(from=paste0(gcm.dir,tasmin.file),to=tmp.dir,overwrite=TRUE)

   tas.roll <- create_time_series(tasmax.file,tasmin.file,map.shape,tmp.dir)
   tas.timing <- find_anomaly_timing(tas.roll,base,anoms)
   anomaly.years[j,] <- tas.timing$years
   anomaly.values[j,] <- tas.timing$values
   raw.anomalies[[j]] <- tas.timing$anoms
   rolled.anomalies[[j]] <- tas.timing$rolled

   file.remove(paste0(tmp.dir,tasmax.file))
   file.remove(paste0(tmp.dir,tasmin.file))
}

rv <- rbind(c('Model','1Degrees','2Degrees','3Degrees','4Degrees'),cbind(gcm.list,anomaly.years))


##write.table(rv,file='/storage/data/projects/rci/data/winter_sports/temperature.anomaly.years.csv',
##            sep=',',quote=FALSE,row.name=FALSE,col.name=F)
 
##plot(1865:2085,rolled.anomalies[[1]],type='l',ylim=c(-2,8),xlab='Year',ylab='Temp Anomaly')

##for (j in seq_along(gcm.list)) {
##   lines(1865:2085,rolled.anomalies[[j]])
##}


}