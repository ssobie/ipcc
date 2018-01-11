##Script to compare Arctic Anomalies with Global Anomalies for NASA GISTEMP

library(ncdf4)
library(raster)
library(PCICt)

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

data.dir <- '/storage/data/projects/rci/data/nrcan/nasa_gistemp/'
data.ann <- paste0(data.dir,'gistemp1200_ERSSTv4_annual.nc')

nc <- nc_open(data.ann)
lon <- ncvar_get(nc,'lon')
lat <- ncvar_get(nc,'lat')
time <- netcdf.calendar(nc)
time.st <- grep('1900',time)
time.en <- length(time)

var.name <- 'tempanomaly'
tas.anom <- ncvar_get(nc,var.name)

rs <- raster(data.ann)
ar <- as.matrix(area(rs))
wts <- ar[,1]/sum(ar[,1])

ac.ix <- lat > 65

arctic.wts <- ar[ac.ix,1]/sum(ar[ac.ix,1])
global.wts <- ar[,1]/sum(ar[,1])

arctic.subset <- tas.anom[,ac.ix,time.st:time.en]
arctic.mid <- apply(arctic.subset,c(1,3),weighted.mean,arctic.wts,na.rm=T)
arctic.mean <- apply(apply(arctic.subset,c(1,3),weighted.mean,arctic.wts,na.rm=T),2,mean)

global.subset <- tas.anom[,,time.st:time.en]
global.mean <- apply(apply(global.subset,c(1,3),weighted.mean,global.wts,na.rm=T),2,mean)
nc_close(nc)

if (1==0) {
  plot.file <- '/storage/data/projects/rci/data/nrcan/nasa_gistemp/gistemp_arctic_global_anomaly_comparison.png'
  png(plot.file,width=800,height=800)
  par(mar=c(5,5,4,2))
  plot(global.mean,arctic.mean,pch=18,cex=3,main='NASA GISTEMP Anomalies',
                   xlab='Global Avg TAS Anomalies (degC)',ylab='Arctic Avg TAS Anomalies (degC)',
                   cex.lab=2,cex.axis=2,cex.main=2.5)
  abline(h=0,v=0,col='gray',lty=2,lwd=2)
  box(which='plot')
  dev.off()
}


