##Script to compare Arctic Anomalies with Global Anomalies for NASA GISTEMP

library(ncdf4)
library(scales)
library(PCICt)
library(nlme)

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

get.arctic.average <- function(file.name,obs.mask) {

   nc <- nc_open(file.name)
   lon <- ncvar_get(nc,'lon')
   lat <- ncvar_get(nc,'lat')
   time <- netcdf.calendar(nc)
   obs.st <- grep('1900',time)
   obs.en <- grep('2099',time)

   clim.st <- grep('1951',time[obs.st:obs.en])
   clim.en <- grep('1980',time[obs.st:obs.en])
  
   time.en <- length(time)

   var.name <- 'tas'
   tas.series <- ncvar_get(nc,var.name) - 273
   tas.subset <- tas.series[,,obs.st:obs.en]

   tas.subset[obs.mask] <- NA

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
   arctic.mean <- apply(apply(arctic.subset,c(1,3),weighted.mean,arctic.wts,na.rm=T),2,mean,na.rm=T)

   global.subset <- tas.anom[,,]
   global.mean <- apply(apply(global.subset,c(1,3),weighted.mean,global.wts,na.rm=T),2,mean,na.rm=T)

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
time.en <- time.en <- grep('2002',time)

var.name <- 'tempanomaly'
tas.anom <- ncvar_get(nc,var.name)
gistemp.mask <- is.na(tas.anom[,,time.st:time.en])

##rs <- raster(data.ann)
##ar2 <- as.matrix(area(rs))

ar <- get.areas(lon,lat,dim(tas.anom))

wts <- ar[,1]/sum(ar[,1])

ac.ix <- lat > 60

arctic.wts <- ar[ac.ix,1]/sum(ar[ac.ix,1])
global.wts <- ar[,1]/sum(ar[,1])

arctic.subset <- tas.anom[,ac.ix,time.st:time.en]
arctic.mid <- apply(arctic.subset,c(1,3),weighted.mean,arctic.wts,na.rm=T)
gis.arctic.mean <- apply(apply(arctic.subset,c(1,3),weighted.mean,arctic.wts,na.rm=T),2,mean)

global.subset <- tas.anom[,,time.st:time.en]
gis.global.mean <- apply(apply(global.subset,c(1,3),weighted.mean,global.wts,na.rm=T),2,mean)
nc_close(nc)


##HadCRU4
data.dir <- '/storage/data/projects/rci/data/nrcan/hadcru4/'
data.ann <- paste0(data.dir,'HadCRUT.4.5.0.0.median.annual.nc')
 
nc <- nc_open(data.ann)
lon <- ncvar_get(nc,'lon')
lat <- ncvar_get(nc,'lat')
time <- netcdf.calendar(nc)
time.st <- grep('1900',time)
time.en <- grep('2002',time)

var.name <- 'tas'
tas.anom <- ncvar_get(nc,var.name)
hadcru4.mask <- is.na(tas.anom[,,time.st:time.en])

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


gistemp.data <- list(arctic=gis.arctic.mean,global=gis.global.mean)
gistemp.fit <- lm(arctic~global,data=gistemp.data)

hadcru4.data <- list(arctic=had.arctic.mean,global=had.global.mean)
hadcru4.fit <- lm(arctic~global,data=hadcru4.data)

##************************************************************************
obs.grid <- 'gistemp'
chose.mask <- gistemp.mask


obs.mask <- array(FALSE,c(dim(chose.mask)[c(1,2)],200))
obs.mask[,,1:dim(chose.mask)[3]] <- gistemp.mask


gcm.list <- c('ACCESS1-0',
              'ACCESS1-3',
              'bcc-csm1-1',
              'bcc-csm1-1-m',
              'CanESM2',
              'CCSM4',
              'CNRM-CM5',
              'CSIRO-Mk3-6-0',
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

##              'FGOALS-g2',

rcp26.list <- c('bcc-csm1-1',
              'bcc-csm1-1-m',
              'CanESM2',
              'CCSM4',
              'CNRM-CM5',
              'CSIRO-Mk3-6-0',
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

rg.45.slope <- rep(0,length(gcm.list))
rg.85.slope <- rep(0,length(gcm.list))
rg.26.slope <- rep(0,length(rcp26.list))


for (g in seq_along(gcm.list)) {
  gcm <- gcm.list[g]
  print(gcm)
  gcm.dir <- paste0('/storage/data/climate/downscale/BCCAQ2/CMIP5/global/annual/',gcm,'/')
  rcp45.file <- paste0(gcm.dir,'tas_ann_',obs.grid,'_grid_',gcm,'_historical+rcp45_r1i1p1_19000101-21001231.nc')
  rcp85.file <- paste0(gcm.dir,'tas_ann_',obs.grid,'_grid_',gcm,'_historical+rcp85_r1i1p1_19000101-21001231.nc')
  ag.45 <- get.arctic.average(rcp45.file,obs.mask)
  ag.85 <- get.arctic.average(rcp85.file,obs.mask)
  arctic.45[[g]] <- ag.45$arctic
  global.45[[g]] <- ag.45$global
  arctic.85[[g]] <- ag.85$arctic
  global.85[[g]] <- ag.85$global

  my.rcp45 <- list(arctic=ag.45$arctic,global=ag.45$global)      

  rcp45.fit <- lm(arctic~global,data=my.rcp45)
  rg.45.slope[g] <- rcp45.fit$coefficients[2]            

  my.rcp85 <- list(arctic=ag.85$arctic,global=ag.85$global)      
  rcp85.fit <- lm(arctic~global,data=my.rcp85)
  rg.85.slope[g] <- rcp85.fit$coefficients[2]            

}

for (g in seq_along(rcp26.list)) {
  gcm <- rcp26.list[g]
  print(gcm)
  gcm.dir <- paste0('/storage/data/climate/downscale/BCCAQ2/CMIP5/global/annual/',gcm,'/')
  rcp26.file <- paste0(gcm.dir,'tas_ann_',obs.grid,'_grid_',gcm,'_historical+rcp26_r1i1p1_19000101-21001231.nc')
  ag.26 <- get.arctic.average(rcp26.file,obs.mask)
  arctic.26[[g]] <- ag.26$arctic
  global.26[[g]] <- ag.26$global
  my.rcp26 <- list(arctic=ag.26$arctic,global=ag.26$global)      
  rcp26.fit <- lm(arctic~global,data=my.rcp26)
  rg.26.slope[g] <- rcp26.fit$coefficients[2]            

}

ratios.26 <- mapply('/',arctic.26,global.26)
ratios.45 <- mapply('/',arctic.45,global.45)
ratios.85 <- mapply('/',arctic.85,global.85)

plot.file <- '/storage/data/projects/rci/data/nrcan/plots/gcm_ratios_arctic_global_time_series2.png'
png(plot.file,width=800,height=800)
  par(mar=c(5,5,4,2))
  par(mfrow=c(3,1))
  plot(c(),xlim=c(2001,2099),ylim=c(-2,6),xlab='Year',ylab='Arctic/Global Anomaly Ratio',main='RCP2.6 GCM Arctic Amplification',
                cex.lab=2,cex.axis=2,cex.main=2.5)
  apply(ratios.26,2,function(y,x){lines(x,y)},1901:2100)
  lines(1901:2100,apply(ratios.26,1,mean),col='red',lwd=3)
  abline(h=2,col='gray')
  box(which='plot')
  plot(c(),xlim=c(2001,2099),ylim=c(-2,6),xlab='Year',ylab='Arctic/Global Anomaly Ratio',main='RCP4.5 GCM Arctic Amplification',
                cex.lab=2,cex.axis=2,cex.main=2.5)
  apply(ratios.45,2,function(y,x){lines(x,y)},1901:2100)
  lines(1901:2100,apply(ratios.45,1,mean),col='red',lwd=3)
  abline(h=2,col='gray')
  box(which='plot')
  plot(c(),xlim=c(2001,2099),ylim=c(-2,6),xlab='Year',ylab='Arctic/Global Anomaly Ratio',main='RCP8.5 GCM Arctic Amplification',
                cex.lab=2,cex.axis=2,cex.main=2.5)
  apply(ratios.85,2,function(y,x){lines(x,y)},1901:2100)
  lines(1901:2100,apply(ratios.85,1,mean),col='red',lwd=3)
  abline(h=2,col='gray')
  box(which='plot')
  dev.off()



browser()

##  my.rcp85 <- list(arctic=unlist(arctic.85),global=unlist(global.85))      
##  rcp85.fit <- lm(arctic~global,data=my.rcp85)
##  my.rcp45 <- list(arctic=unlist(arctic.45),global=unlist(global.45))      
##  rcp45.fit <- lm(arctic~global,data=my.rcp45)
##  my.rcp26 <- list(arctic=unlist(arctic.26),global=unlist(global.26))      
##  rcp26.fit <- lm(arctic~global,data=my.rcp26)

model.slopes <- matrix(NA,nrow=length(gcm.list),ncol=3)
model.slopes[gcm.list %in% rcp26.list,1] <- round(rg.26.slope,2)
model.slopes[,2] <- round(rg.45.slope,2)
model.slopes[,3] <- round(rg.85.slope,2)

model.slopes <- cbind(c(gcm.list,'5%','Mean','95%'),
                      rbind(model.slopes,
                            round(apply(model.slopes,2,quantile,0.05,na.rm=T),2),
                            round(apply(model.slopes,2,mean,na.rm=T),2),
                            round(apply(model.slopes,2,quantile,0.95,na.rm=T),2)))

##write.file <- paste0('/storage/data/projects/rci/data/nrcan/nasa_gistemp/slopes/',obs.grid,'_slopes_2002.csv')

##write.table(model.slopes,file=write.file,quote=F,row.name=F,col.name=F,sep=',')


if (1==0) {
  ccsm4.rcp26 <- list(arctic=arctic.26[[5]],global=global.26[[5]])      
  rcp26.fit <- lm(arctic~global,data=ccsm4.rcp26)
  ccsm4.rcp45 <- list(arctic=arctic.45[[7]],global=global.45[[7]])      
  rcp45.fit <- lm(arctic~global,data=ccsm4.rcp45)
  ccsm4.rcp85 <- list(arctic=arctic.85[[7]],global=global.85[[7]])      
  rcp85.fit <- lm(arctic~global,data=ccsm4.rcp85)

##abline(rcp26.fit$coefficients, col='green',lwd=4)
##abline(rcp45.fit$coefficients, col='orange',lwd=4)
##abline(rcp85.fit$coefficients, col='red',lwd=4)
}


##Autocorrelation for obs linear fit
##Lag-1 autocorrelation

gistemp.ac <- gls(arctic~global,data=gistemp.data,correlation = corARMA(p=1,q=1))
print('GISTEMP Fit')
print('Standard Fit')
print(coef(gistemp.fit))
print(confint(gistemp.fit,level=0.9))
print('Lag-1 Fit')
print(coef(gistemp.ac))
print(confint(gistemp.ac,level=0.9))

print('HadCRU4 Fit')
hadcru4.ac <- gls(arctic~global,data=hadcru4.data,correlation = corARMA(p=1,q=1))
print('Standard fit')
print(coef(hadcru4.fit))
print(confint(hadcru4.fit,level=0.9))
print('Lag-1 Fit')
print(coef(hadcru4.ac))
print(confint(hadcru4.ac,level=0.9))

