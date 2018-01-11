##Script to compare Arctic Anomalies with Global Anomalies for NASA GISTEMP

library(ncdf4)
library(PCICt)

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

get.areas <- function(lon,lat,dm) {
  lond <- diff(lon)[1]/2
  latd <- diff(lat)[1]/2      
  
  lon.up <- lon + lond
  lon.dn <- lon - lond
  lat.up <- lat + latd
  lat.dn <- lat - latd
  
  areas <- matrix(NA,nrow=dm[2],ncol=dm[1])

  for (i in 1:length(lon)) {
    for (j in 1:length(lat)) {
     areas[j,i] <- 2*(pi/180)*6371^2 * abs(sin(lat.up[j]/180*pi)-sin(lat.dn[j]/180*pi)) *abs(lon.up[i]-lon.dn[i])
    }
  }

  return(areas)
}

##---------------------------------------------------------------------
##GISTEMP

get.obs.means <- function(file.name,var.name,yst,yen) {

  data.ann <- file.name
 
  nc <- nc_open(data.ann)
  lon <- ncvar_get(nc,'lon')
  lat <- ncvar_get(nc,'lat')
  time <- netcdf.calendar(nc)
  time.st <- grep(yst,time)
  time.en <- grep(yen,time)

  base.st <- grep('1861',time)
  base.en <- grep('1880',time)  

  offset.st <- grep('1951',time)
  offset.en <- grep('1980',time)  

  tas.anom <- ncvar_get(nc,var.name)
  obs.mask <- is.na(tas.anom[,,time.st:time.en])
  
  ar <- get.areas(lon,lat,dim(tas.anom))

  wts <- ar[,1]/sum(ar[,1])

  ac.ix <- lat > 60

  arctic.wts <- ar[ac.ix,1]/sum(ar[ac.ix,1])
  global.wts <- ar[,1]/sum(ar[,1])

  arctic.subset <- tas.anom[,ac.ix,time.st:time.en]
  arctic.mid <- apply(arctic.subset,c(1,3),weighted.mean,arctic.wts,na.rm=T)
  arctic.mean <- apply(apply(arctic.subset,c(1,3),weighted.mean,arctic.wts,na.rm=T),2,mean,na.rm=T)

  global.subset <- tas.anom[,,time.st:time.en]
  global.mean <- apply(apply(global.subset,c(1,3),weighted.mean,global.wts,na.rm=T),2,mean,na.rm=T)


  base.subset <- tas.anom[,,base.st:base.en]
  base.mean <- apply(apply(base.subset,c(1,3),weighted.mean,global.wts,na.rm=T),2,mean,na.rm=T)
  offset.subset <- tas.anom[,,offset.st:offset.en]
  offset.mean <- apply(apply(offset.subset,c(1,3),weighted.mean,global.wts,na.rm=T),2,mean,na.rm=T)
  
  global.offset <- mean(offset.mean)-mean(base.mean)
  arctic.offset <- global.offset * 1.89
  
  nc_close(nc)

  rv <- list(arctic=arctic.mean+arctic.offset,global=global.mean+global.offset,mask=obs.mask)
  return(rv)
}




