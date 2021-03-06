##Script to compare Arctic Anomalies with Global Anomalies for NASA GISTEMP

library(ncdf4)
library(scales)
library(PCICt)

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)


calcSSE <- function(par,mydata){
 
  loessMod <- try(loess(arctic ~ global, data=mydata, span=x), silent=T)
  res <- try(loessMod$residuals, silent=T)
  if(class(res)!="try-error"){
    if((sum(res, na.rm=T) > 15)){
      sse <- sum(res^2)  
    }
  }else{
    sse <- 99999
  }
  return(sse)
}

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

get.arctic.average <- function(file.name) {

   nc <- nc_open(file.name)
   lon <- ncvar_get(nc,'lon')
   lat <- ncvar_get(nc,'lat')
   time <- netcdf.calendar(nc)
   obs.st <- grep('1900',time)
   obs.en <- grep('2016',time)

   clim.st <- grep('1951',time[obs.st:obs.en])
   clim.en <- grep('1980',time[obs.st:obs.en])

   time.en <- length(time)

   var.name <- 'tas'
   tas.series <- ncvar_get(nc,var.name) - 273
   tas.subset <- tas.series[,,obs.st:obs.en]

   dim(tas.subset)
   tas.clim <- apply(tas.subset[,,clim.st:clim.en],c(1,2),mean,na.rm=T)
   tas.avg <- aperm(apply(tas.clim,c(1,2),function(x){rep(x,each=dim(tas.subset)[3])}),c(2,3,1))

   tas.anom <- tas.subset - tas.avg

   ##rs <- raster(file.name)
   ##ar2 <- as.matrix(area(rs))
   ar <- get.areas(lon,lat,dim(tas.anom))
   wts <- ar[,1]/sum(ar[,1])

   ac.ix <- lat > 60

   arctic.wts <- ar[ac.ix,1]/sum(ar[ac.ix,1])
   global.wts <- ar[,1]/sum(ar[,1])

   arctic.subset <- tas.anom[,ac.ix,]
   arctic.mid <- apply(arctic.subset,c(1,3),weighted.mean,arctic.wts,na.rm=T)
   arctic.mean <- apply(apply(arctic.subset,c(1,3),weighted.mean,arctic.wts,na.rm=T),2,mean)

   global.subset <- tas.anom[,,]
   global.mean <- apply(apply(global.subset,c(1,3),weighted.mean,global.wts,na.rm=T),2,mean)
   nc_close(nc)
   
   if (length(arctic.mean) == 151) {
     arctic.mean <- arctic.mean[-1]
     global.mean <- global.mean[-1]
   }
   rv <- list(arctic=arctic.mean,       
              global=global.mean)
   return(rv)
}


##GISTEMP
data.dir <- '/storage/data/projects/rci/data/nrcan/nasa_gistemp/'
data.ann <- paste0(data.dir,'gistemp1200_ERSSTv4_annual.nc')
 
nc <- nc_open(data.ann)
lon <- ncvar_get(nc,'lon')
lat <- ncvar_get(nc,'lat')
time <- netcdf.calendar(nc)
time.st <- grep('1900',time)
time.en <- grep('2016',time)

var.name <- 'tempanomaly'
tas.anom <- ncvar_get(nc,var.name)

##rs <- raster(data.ann)
##ar2 <- as.matrix(area(rs))

ar <- get.areas(lon,lat,dim(tas.anom))

wts <- ar[,1]/sum(ar[,1])

ac.ix <- lat > 60

arctic.wts <- ar[ac.ix,1]/sum(ar[ac.ix,1])
global.wts <- ar[,1]/sum(ar[,1])

arctic.subset <- tas.anom[,ac.ix,time.st:time.en]
arctic.mid <- apply(arctic.subset,c(1,3),weighted.mean,arctic.wts,na.rm=T)
arctic.mean <- apply(apply(arctic.subset,c(1,3),weighted.mean,arctic.wts,na.rm=T),2,mean)

global.subset <- tas.anom[,,time.st:time.en]
global.mean <- apply(apply(global.subset,c(1,3),weighted.mean,global.wts,na.rm=T),2,mean)
nc_close(nc)


##HadCRU4
data.dir <- '/storage/data/projects/rci/data/nrcan/hadcru4/'
data.ann <- paste0(data.dir,'HadCRUT.4.5.0.0.median.annual.nc')
 
nc <- nc_open(data.ann)
lon <- ncvar_get(nc,'lon')
lat <- ncvar_get(nc,'lat')
time <- netcdf.calendar(nc)
time.st <- grep('1900',time)
time.en <- grep('2016',time)

var.name <- 'tas'
tas.anom <- ncvar_get(nc,var.name)

##rs <- raster(data.ann)
##ar2 <- as.matrix(area(rs))
ar <- get.areas(lon,lat,dim(tas.anom))
wts <- ar[,1]/sum(ar[,1])

ac.ix <- lat > 60

arctic.wts <- ar[ac.ix,1]/sum(ar[ac.ix,1])
global.wts <- ar[,1]/sum(ar[,1])

arctic.subset <- tas.anom[,ac.ix,time.st:time.en]
arctic.mid <- apply(arctic.subset,c(1,3),weighted.mean,arctic.wts,na.rm=T)
had.arctic.mean <- apply(apply(arctic.subset,c(1,3),weighted.mean,arctic.wts,na.rm=T),2,mean,na.rm=T)

global.subset <- tas.anom[,,time.st:time.en]
had.global.mean <- apply(apply(global.subset,c(1,3),weighted.mean,global.wts,na.rm=T),2,mean,na.rm=T)
nc_close(nc)

  gissdata <- list(arctic=arctic.mean,global=global.mean)
  gissfit <- lm(arctic~global,data=gissdata)

  haddata <- list(arctic=had.arctic.mean,global=had.global.mean)
  hadfit <- lm(arctic~global,data=haddata)

  giss.arctic.sub <- arctic.mean[104:117]
  giss.global.sub <- global.mean[104:117]

  had.arctic.sub <- had.arctic.mean[104:117]
  had.global.sub <- had.global.mean[104:117]


if (1==0) {


gcm.list <- c('ACCESS1-0',
              'ACCESS1-3',
              'bcc-csm1-1',
              'bcc-csm1-1-m',
              'BNU-ESM',
              'CanESM2',
              'CCSM4',
              'CNRM-CM5',
              'CSIRO-Mk3-6-0',
              'FGOALS-g2',
              'GFDL-CM3',
              'GFDL-ESM2G',
              'GFDL-ESM2M',
              'HadGEM2-AO',
              'HadGEM2-CC',
              'HadGEM2-ES',
              'inmcm4',
              'IPSL-CM5A-LR',
              'IPSL-CM5A-MR',
              'IPSL-CM5B-LR',
              'MIROC5',
              'MIROC-ESM',
              'MIROC-ESM-CHEM',
              'MPI-ESM-MR',
              'MPI-ESM-LR',
              'MRI-CGCM3',
              'NorESM1-M')

rcp26.list <- c('bcc-csm1-1',
              'bcc-csm1-1-m',
              'BNU-ESM',
              'CanESM2',
              'CCSM4',
              'CNRM-CM5',
              'CSIRO-Mk3-6-0',
              'FGOALS-g2',
              'GFDL-CM3',
              'GFDL-ESM2G',
              'GFDL-ESM2M',
              'HadGEM2-AO',
              'HadGEM2-ES',
              'IPSL-CM5A-LR',
              'IPSL-CM5A-MR',
              'MIROC5',
              'MIROC-ESM',
              'MIROC-ESM-CHEM',
              'MPI-ESM-MR',
              'MPI-ESM-LR',
              'MRI-CGCM3',
              'NorESM1-M')


arctic.26 <- vector(length=length(rcp26.list),mode='list')
global.26 <- vector(length=length(rcp26.list),mode='list')

arctic.45 <- vector(length=length(gcm.list),mode='list')
global.45 <- vector(length=length(gcm.list),mode='list')

arctic.85 <- vector(length=length(gcm.list),mode='list')
global.85 <- vector(length=length(gcm.list),mode='list')

for (g in seq_along(gcm.list)) {
  gcm <- gcm.list[g]
  print(gcm)
  gcm.dir <- paste0('/storage/data/climate/downscale/BCCAQ2/CMIP5/global/annual/',gcm,'/')
  rcp45.file <- paste0(gcm.dir,'tas_ann_',gcm,'_historical+rcp45_r1i1p1_19500101-21001231.nc')
  rcp85.file <- paste0(gcm.dir,'tas_ann_',gcm,'_historical+rcp85_r1i1p1_19500101-21001231.nc')
  ag.45 <- get.arctic.average(rcp45.file)
  ag.85 <- get.arctic.average(rcp85.file)
  arctic.45[[g]] <- ag.45$arctic
  global.45[[g]] <- ag.45$global
  arctic.85[[g]] <- ag.85$arctic
  global.85[[g]] <- ag.85$global

}

for (g in seq_along(rcp26.list)) {
  gcm <- rcp26.list[g]
  print(gcm)
  gcm.dir <- paste0('/storage/data/climate/downscale/BCCAQ2/CMIP5/global/annual/',gcm,'/')
  rcp26.file <- paste0(gcm.dir,'tas_ann_',gcm,'_historical+rcp26_r1i1p1_19500101-21001231.nc')
  ag.26 <- get.arctic.average(rcp26.file)
  arctic.26[[g]] <- ag.26$arctic
  global.26[[g]] <- ag.26$global

}



  my.rcp85 <- list(arctic=unlist(arctic.85),global=unlist(global.85))      
  rcp85.fit <- lm(arctic~global,data=my.rcp85)
  my.rcp45 <- list(arctic=unlist(arctic.45),global=unlist(global.45))      
  rcp45.fit <- lm(arctic~global,data=my.rcp45)
  my.rcp26 <- list(arctic=unlist(arctic.26),global=unlist(global.26))      
  rcp26.fit <- lm(arctic~global,data=my.rcp26)
##browser()
  all.rcps <- list(arctic=c(unlist(arctic.26),unlist(arctic.45),unlist(arctic.85)),
                   global=c(unlist(global.26),unlist(global.45),unlist(global.85)))

  rcp26.loess <- loess(arctic~global,data=my.rcp26,span=0.25)  
  rcp45.loess <- loess(arctic~global,data=my.rcp45,span=0.25)  
  rcp85.loess <- loess(arctic~global,data=my.rcp85,span=0.25)  
  rcps.loess <- loess(arctic~global,data=all.rcps,span=0.25)

  global.avg.85 <- apply(matrix(unlist(global.85),nrow=27,ncol=150,byrow=T),2,mean)
  global.avg.45 <- apply(matrix(unlist(global.45),nrow=27,ncol=150,byrow=T),2,mean)
  global.avg.26 <- apply(matrix(unlist(global.26),nrow=22,ncol=150,byrow=T),2,mean)     

  arctic.avg.85 <- apply(matrix(unlist(arctic.85),nrow=27,ncol=150,byrow=T),2,mean)
  arctic.avg.45 <- apply(matrix(unlist(arctic.45),nrow=27,ncol=150,byrow=T),2,mean)
  arctic.avg.26 <- apply(matrix(unlist(arctic.26),nrow=22,ncol=150,byrow=T),2,mean)     
  

##browser()
  smoothed.rcp26 <- predict(rcp26.loess)
  smoothed.rcp45 <- predict(rcp45.loess)
  smoothed.rcp85 <- predict(rcp85.loess)
  smoothed.rcps <- predict(rcps.loess)



  ##plot.file <- paste0('/storage/data/projects/rci/data/nrcan/nasa_gistemp/gcm_smoothed_gistemp_hadcru4_arctic_global_anomaly_comparison.png')
  plot.file <- paste0('/storage/data/projects/rci/data/nrcan/nasa_gistemp/gcm_smoothed_hadcru4_arctic_global_anomaly_comparison.png')
  png(plot.file,width=800,height=800)
  par(mar=c(5,3,4,5))
  plot(had.global.mean,had.arctic.mean,xlim=c(-1,6),ylim=c(-2,16),pch=0,col='black',cex=2,main='GCM,  HadCRU4 Anomalies', ##NASA GISTEMP,
                   xlab='Global Avg TAS Anomalies (degC)',ylab='',
                   cex.lab=2,cex.axis=2,cex.main=2.5,axes=F)
  
  axis(1,seq(-1,6,1),seq(-1,6,1),cex.axis=2)
  axis(4,seq(0,15,5),seq(0,15,5),cex.axis=2)
  mtext("Arctic Avg TAS Anomalies (degC)",side=4,line=3,cex=2)

  test <- chull(unlist(global.85),unlist(arctic.85))
  hpts <- c(test,test[1])
  polygon(unlist(global.85)[hpts],unlist(arctic.85)[hpts],col=alpha('red',0.3),border=alpha('red',0.3))
  test <- chull(unlist(global.45),unlist(arctic.45))
  hpts <- c(test,test[1])
  polygon(unlist(global.45)[hpts],unlist(arctic.45)[hpts],col=alpha('green',0.3),border=alpha('green',0.3))
  test <- chull(unlist(global.26),unlist(arctic.26))
  hpts <- c(test,test[1])
  polygon(unlist(global.26)[hpts],unlist(arctic.26)[hpts],col=alpha('blue',0.3),border=alpha('blue',0.3))

  ##points(global.mean,arctic.mean,pch=0,cex=3)
  points(had.global.mean,had.arctic.mean,pch=0,col='black',cex=2)
  abline(coef=hadfit$coefficients,col='black',lwd=4)
  g <- order(all.rcps$global)
  i <- order(my.rcp26$global)
  j <- order(my.rcp45$global)
  k <- order(my.rcp85$global)
  lines(all.rcps$global[g],smoothed.rcps[g],col='darkgray',lwd=5) 
  abline(h=0,v=0,col='gray',lty=2,lwd=2)

  yrs <- 1951:2100
  ix <- seq(0,150,50)
  ix[1] <- 1
  for (j in c(4)) {
    points(x=global.avg.26[ix[j]],y=arctic.avg.26[ix][j],pch=19,cex=3,col='lightblue')
    text(x=global.avg.26[ix[j]],y=arctic.avg.26[ix][j]+0.5,yrs[ix][j],col='lightblue',cex=1.5)
  }
  for (j in 3:4) {
    points(x=global.avg.45[ix[j]],y=arctic.avg.45[ix][j],pch=19,cex=3,col='green')
    text(x=global.avg.45[ix[j]],y=arctic.avg.45[ix][j]+0.5,yrs[ix][j],col='green',cex=1.5)
  }
  for(j in 3:4) {
    points(x=global.avg.85[ix[j]],y=arctic.avg.85[ix][j],pch=19,cex=3,col='red')
    text(x=global.avg.85[ix[j]],y=arctic.avg.85[ix][j]+0.5,yrs[ix][j],col='red',cex=1.5)
  }
  legend('topleft',leg=c('HadCRU4','Obs. Fit','RCP2.6','RCP4.5','RCP8.5','Model Fit'),
          col=c('black','black',alpha('blue',0.3),alpha('green',0.3),alpha('red',0.3),'darkgray'),
          pch=c(0,15,15,15,15,15),cex=2)
  box(which='plot')
  dev.off()
}


if (1==0) {
##  for (g in seq_along(gcm.list)) {
##     points(global.85[[g]],arctic.85[[g]],col=alpha('red',0.25),pch=2,cex=2.5) 
##     points(global.45[[g]],arctic.45[[g]],col=alpha('green',0.25),pch=2,cex=2.5) 
##  }
##  for (g in seq_along(rcp26.list)) {
##     points(global.26[[g]],arctic.26[[g]],col=alpha('lightblue',0.25),pch=2,cex=2.5) 
##  }

  ##abline(coef=rcp26.fit$coefficients,col='blue',lwd=3)
  ##abline(coef=rcp45.fit$coefficients,col='green',lwd=3)
  ##abline(coef=rcp85.fit$coefficients,col='red',lwd=3)
  ##abline(coef=gissfit$coefficients,col='black',lwd=3)

#  lines(my.rcp26$global[i],smoothed.rcp26[i],col='blue',lwd=3)
#  lines(my.rcp45$global[j],smoothed.rcp45[j],col='darkgreen',lwd=3)
#  lines(my.rcp85$global[k],smoothed.rcp85[k],col='darkred',lwd=3)  

}



if (1==0) {
  plot.file <- '/storage/data/projects/rci/data/nrcan/nasa_gistemp/gistemp_hadcru4_arctic_global_anomaly_comparison.png'
  png(plot.file,width=800,height=800)
  par(mar=c(5,5,4,2))
  plot(global.mean,arctic.mean,xlim=c(-1,2),ylim=c(-2,4),pch=18,cex=3,main='NASA GISTEMP and HadCRU4 Anomalies',
                   xlab='Global Avg TAS Anomalies (degC)',ylab='Arctic Avg TAS Anomalies (degC)',
                   cex.lab=2,cex.axis=2,cex.main=2.5)
  points(global.mean,arctic.mean,pch=18,cex=3)
  points(had.global.mean,had.arctic.mean,pch=16,col='red',cex=2)
  
  points(giss.global.sub,giss.arctic.sub,pch=18,cex=3,col='green')
  points(had.global.sub,had.arctic.sub,pch=16,col='blue',cex=2)

  abline(coef=gissfit$coefficients,col='black',lwd=3)
  abline(coef=hadfit$coefficients,col='red',lwd=3)
  
  abline(h=0,v=0,col='gray',lty=2,lwd=2)
  legend('topleft',leg=c('GISTEMP (1900-2002)','GISTEMP (2003-2016)','HadCRU4 (1900-2002)','HadCRU4 (2003-2016)'),col=c('black','green','red','blue'),pch=c(18,18,16,16),cex=2)
  box(which='plot')
  dev.off()
}

if (1==0) {
##HadCRU4 Individual Simulations
data.dir <- '/storage/data/projects/rci/data/nrcan/HadCRU4/annual/'
slen <- 100 
hadsim.arctic.mean <- vector(mode='list',length=slen)
hadsim.global.mean <- vector(mode='list',length=slen)

  plot.file <- '/storage/data/projects/rci/data/nrcan/nasa_gistemp/gistemp_hadcru4_all_simulations_arctic_global_anomaly.png'
  png(plot.file,width=800,height=800)
  par(mar=c(5,5,4,2))
  plot(global.mean,arctic.mean,xlim=c(-1,2),ylim=c(-2,4),pch=18,cex=3,main='NASA GISTEMP and HadCRU4 Anomalies',
                   xlab='Global Avg TAS Anomalies (degC)',ylab='Arctic Avg TAS Anomalies (degC)',
                   cex.lab=2,cex.axis=2,cex.main=2.5)

for (i in 1:slen) {
  print(i)
  data.ann <- paste0(data.dir,'HadCRUT.4.5.0.0.anomalies.annual.',i,'.nc')
  nc <- nc_open(data.ann)
  var.name <- 'temperature_anomaly'
  tas.anom <- ncvar_get(nc,var.name)

  arctic.subset <- tas.anom[,ac.ix,time.st:time.en] 
  hadsim.arctic.mean[[i]] <- apply(apply(arctic.subset,c(1,3),weighted.mean,arctic.wts,na.rm=T),2,mean,na.rm=T)

  global.subset <- tas.anom[,,time.st:time.en]
  hadsim.global.mean[[i]] <- apply(apply(global.subset,c(1,3),weighted.mean,global.wts,na.rm=T),2,mean,na.rm=T)
  nc_close(nc)

  points(hadsim.global.mean[[i]],hadsim.arctic.mean[[i]],pch=1,col='red',cex=2)

}


  points(global.mean,arctic.mean,pch=18,cex=3)
##  points(had.global.mean,had.arctic.mean,pch=16,col='red',cex=2)
##  abline(coef=gissfit$coefficients,col='black',lwd=3)
##  abline(coef=hadfit$coefficients,col='red',lwd=3)

  abline(h=0,v=0,col='gray',lty=2,lwd=2)
  legend('topleft',leg=c('GISTEMP','HadCRU4'),col=c('black','red'),pch=c(18,16),cex=2)
  box(which='plot')
  dev.off()

  plot.file <- '/storage/data/projects/rci/data/nrcan/nasa_gistemp/gistemp_hadcru4_all_simulations_arctic_global_time_series.png'
  png(plot.file,width=800,height=800)
  par(mar=c(5,5,4,2))
  plot(c(),xlim=c(1900,2020),ylim=c(-2,3),xlab='Year',ylab='Temp. Anomaly (degC)',main='HadCRU4 Ensemble Anomalies',
                cex.lab=2,cex.axis=2,cex.main=2.5)
  lapply(hadsim.arctic.mean,function(x){lines(1900:2017,x)})
  lapply(hadsim.global.mean,function(x){lines(1900:2017,x,col='red')})
  lines(1900:2017,had.global.mean,col='orange',lwd=3)
  lines(1900:2017,had.arctic.mean,col='blue',lwd=3)
  abline(h=0,col='gray')
  legend('topleft',leg=c('Arctic Anomalies','Global Anomalies'),col=c('black','red'),pch=c(18,16),cex=2)
  box(which='plot')
  dev.off()  
}