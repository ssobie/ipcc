##Script to compare Arctic Anomalies with Global Anomalies for NASA GISTEMP

library(ncdf4)
library(scales)
library(PCICt)
library(nlme)
library(zyp)


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
   obs.en <- length(time)
   obs.st <- obs.en-105

   print(obs.st)
   print(obs.en)
   clim.st <- 51 ##grep('1951',time[obs.st:obs.en])
   clim.en <- 80 ##grep('1980',time[obs.st:obs.en])
  
   time.en <- length(time)

   var.name <- 'tas'
   tas.series <- ncvar_get(nc,var.name) - 273
   
   tas.subset <- tas.series[,,obs.st:obs.en]
   rm(tas.series)              
   tas.subset[obs.mask] <- NA
   
   print(dim(tas.subset))

   tas.clim <- apply(tas.subset[,,clim.st:clim.en],c(1,2),mean,na.rm=T)
   tas.avg <- aperm(apply(tas.clim,c(1,2),function(x){rep(x,each=dim(tas.subset)[3])}),c(2,3,1))

   tas.anom <- tas.subset - tas.avg

   rm(tas.subset)
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
time.en <- time.en <- grep('2005',time)

var.name <- 'tempanomaly'
tas.anom <- ncvar_get(nc,var.name)
gistemp.mask <- is.na(tas.anom[,,time.st:time.en])

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
time.en <- grep('2005',time)

var.name <- 'tas'
tas.anom <- ncvar_get(nc,var.name)
hadcru4.mask <- is.na(tas.anom[,,time.st:time.en])
obs.mask <- hadcru4.mask
obs.grid <- 'hadcru4'

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

##*************************************************************************
##GCMs

gcm.list <- c('ACCESS1-0','ACCESS1-3','BNU-ESM','CanESM2','CSIRO-Mk3-6-0','GFDL-ESM2M','IPSL-CM5A-MR','IPSL-CM5B-LR',
              'MPI-ESM-LR','MPI-ESM-MR','NorESM1-M','bcc-csm1-1','bcc-csm1-1-m','CNRM-CM5','inmcm4',
              'HadGEM2-CC','HadGEM2-ES','IPSL-CM5A-LR','MIROC5','MIROC-ESM','MIROC-ESM-CHEM')

arctic.past <- vector(length=length(gcm.list),mode='list')
global.past <- vector(length=length(gcm.list),mode='list')

rg.past.slope <- rep(0,length(gcm.list))
rg.past.upper <- rep(0,length(gcm.list))
rg.past.lower <- rep(0,length(gcm.list))

global.trends <- matrix(0,nrow=length(gcm.list),ncol=3)
arctic.trends <- matrix(0,nrow=length(gcm.list),ncol=3)

for (g in seq_along(gcm.list)) {
  gcm <- gcm.list[g]
  print(gcm)
  gcm.dir <- paste0('/storage/data/climate/downscale/BCCAQ2/CMIP5/global/control/annual/',gcm,'/')
  
  past.file <- list.files(path=gcm.dir,pattern=paste0('tas_ann_',obs.grid,'_grid_',gcm,'_piControl*'),full.name=T)
  print(past.file)

  ag.past <- get.arctic.average(past.file,obs.mask)

  arctic.past[[g]] <- ag.past$arctic
  global.past[[g]] <- ag.past$global

  my.past <- list(arctic=ag.past$arctic,global=ag.past$global)      
  my.global <- list(time=1:length(ag.past$global),data=ag.past$global)
  my.arctic <- list(time=1:length(ag.past$arctic),data=ag.past$arctic)
  past.fit <- lm(arctic~global,data=my.past)
  
  global.zyp <- zyp.sen(data~time,data=my.global)
  global.trends[g,1] <- global.zyp$coefficients[2]
  global.trends[g,c(2,3)] <- confint(global.zyp,level=0.9)[1,]

  arctic.zyp <- zyp.sen(data~time,data=my.arctic)
  arctic.trends[g,1] <- arctic.zyp$coefficients[2]
  arctic.trends[g,c(2,3)] <- confint(arctic.zyp,level=0.9)[1,]

  rg.past.slope[g] <- past.fit$coefficients[2]            
  conf.ints <- confint(past.fit,level=0.9)
  rg.past.upper[g] <- conf.ints[2,2]
  rg.past.lower[g] <- conf.ints[2,1]

}

model.slopes <- round(rg.past.slope,2)
model.lower <- round(rg.past.lower,2)
model.upper <- round(rg.past.upper,2)
model.matrix <- cbind(model.slopes,model.lower,model.upper)

model.slopes <- cbind(c(gcm.list,'5%','Mean','95%'),
                      rbind(model.matrix,
                            round(apply(model.matrix,2,quantile,0.05,na.rm=T),2),
                            round(apply(model.matrix,2,mean,na.rm=T),2),
                            round(apply(model.matrix,2,quantile,0.95,na.rm=T),2)))

write.file <- paste0('/storage/data/projects/rci/data/nrcan/nasa_gistemp/slopes/',obs.grid,'_slopes_Control.csv')
write.table(model.slopes,file=write.file,quote=F,row.name=F,col.name=F,sep=',')

global.tm <- cbind(c(gcm.list,'5%','Mean','95%'),
                      rbind(global.trends,
                            apply(global.trends,2,quantile,0.05,na.rm=T),
                            apply(global.trends,2,mean,na.rm=T),
                            apply(global.trends,2,quantile,0.95,na.rm=T)))

write.file <- paste0('/storage/data/projects/rci/data/nrcan/nasa_gistemp/slopes/',obs.grid,'_global_trends_Control.csv')
write.table(global.tm,file=write.file,quote=F,row.name=F,col.name=F,sep=',')

arctic.tm <- cbind(c(gcm.list,'5%','Mean','95%'),
                      rbind(arctic.trends,
                            apply(arctic.trends,2,quantile,0.05,na.rm=T),
                            apply(arctic.trends,2,mean,na.rm=T),
                            apply(arctic.trends,2,quantile,0.95,na.rm=T)))

write.file <- paste0('/storage/data/projects/rci/data/nrcan/nasa_gistemp/slopes/',obs.grid,'_arctic_trends_Control.csv')
write.table(arctic.tm,file=write.file,quote=F,row.name=F,col.name=F,sep=',')





if (1==0) {
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

}


if (1==0) {
plot.file <- '/storage/data/projects/rci/data/nrcan/nasa_gistemp/gistemp_CONTROL_arctic_global_anomaly.png'
  png(plot.file,width=800,height=800)
  par(mar=c(5,5,4,2))
  plot(gis.global.mean,gis.arctic.mean,xlim=c(-1,1),ylim=c(-2,2),cex=2,main='NASA GISTEMP and piControl GCM Anomalies',
                   xlab='Global Avg TAS Anomalies (degC)',ylab='Arctic Avg TAS Anomalies (degC)',
                   cex.lab=2,cex.axis=2,cex.main=2.5,pch=16,col='red')
  mapply(FUN=points,global.past,arctic.past,pch=18,col='black',cex=2)
  points(gis.global.mean,gis.arctic.mean,cex=2,pch=16,col='red')
   
  abline(h=0,v=0,col='gray',lty=2,lwd=2)
  legend('topleft',leg=c('GISTEMP','GCM'),col=c('red','black'),pch=c(16,18),cex=2)
  box(which='plot')
  dev.off()
}
if(1==1) {
  plot.file <- '/storage/data/projects/rci/data/nrcan/nasa_gistemp/hadcru4_CONTROL_arctic_global_anomaly.png'
  png(plot.file,width=800,height=800)
  par(mar=c(5,5,4,2))
  plot(had.global.mean,had.arctic.mean,xlim=c(-1,1),ylim=c(-2,2),cex=2,main='HadCRU4 and piControl GCM Anomalies',
                   xlab='Global Avg TAS Anomalies (degC)',ylab='Arctic Avg TAS Anomalies (degC)',
                   cex.lab=2,cex.axis=2,cex.main=2.5,pch=16,col='red')
  mapply(FUN=points,global.past,arctic.past,pch=18,col='black',cex=2)
  points(had.global.mean,had.arctic.mean,cex=2,pch=16,col='red')
   
  abline(h=0,v=0,col='gray',lty=2,lwd=2)
  legend('topleft',leg=c('HadCRU4','GCM'),col=c('red','black'),pch=c(16,18),cex=2)
  box(which='plot')
  dev.off()
}

