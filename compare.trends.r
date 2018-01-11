##Script to compare Arctic Anomalies with Global Anomalies for NASA GISTEMP

library(ncdf4)
library(scales)
library(PCICt)
library(nlme)
library(zyp)


source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)
source('/storage/home/ssobie/code/repos/ipcc/obs.read.r',chdir=T)
source('/storage/home/ssobie/code/repos/ipcc/arctic.average.r',chdir=T)


##---------------------------------------------------------------------
##Gridded Obs
gis.file <- '/storage/data/projects/rci/data/nrcan/nasa_gistemp/gistemp1200_ERSSTv4_annual.nc'
gistemp.data <- get.obs.means(file.name=gis.file,var.name='tempanomaly','1900','2005')
##gistemp.fit <- lm(arctic~global,data=gistemp.data)
gistemp.mask <- gistemp.data$mask

had.file <- '/storage/data/projects/rci/data/nrcan/hadcru4/HadCRUT.4.5.0.0.median.annual.nc'
hadcru4.data <- get.obs.means(file.name=had.file,var.name='tas','1900','2005')
##hadcru4.fit <- lm(arctic~global,data=hadcru4.data)
hadcru4.mask <- hadcru4.data$mask

##---------------------------------------------------------------------
##*************************************************************************
##GCMs
gcm.list <- c('ACCESS1-3','CanESM2','CSIRO-Mk3-6-0','GFDL-ESM2M','IPSL-CM5A-MR',
              'MIROC-ESM','NorESM1-M','bcc-csm1-1','CNRM-CM5',
              'HadGEM2-ES','IPSL-CM5A-LR','MIROC-ESM-CHEM')


obs.data <- gistemp.data
chose.mask <- gistemp.mask
obs.grid <- 'gistemp'

obs.mask <- array(FALSE,c(dim(chose.mask)[c(1,2)],106))
obs.mask[,,1:dim(chose.mask)[3]] <- gistemp.mask



##---------------------------------------------------------------------

control.past <- vector(length=length(gcm.list),mode='list')
natural.past <- vector(length=length(gcm.list),mode='list')
rcp85.past <- vector(length=length(gcm.list),mode='list')

global.trends <- matrix(0,nrow=length(gcm.list),ncol=3)
arctic.trends <- matrix(0,nrow=length(gcm.list),ncol=3)

for (g in seq_along(gcm.list)) {
  gcm <- gcm.list[g]
  print(gcm)

  ##Control Runs
  base.dir <- '/storage/data/climate/downscale/BCCAQ2/CMIP5/global/control/annual/'
  scenario <- 'piControl' 
  gcm.dir <- paste0(base.dir,gcm)  
  file.control <- list.files(path=gcm.dir,pattern=paste0('tas_ann_',obs.grid,'_grid_',gcm,'_',scenario,'*'),full.name=T)
  print(file.control)
  nc.control <- nc_open(file.control)
  time <- netcdf.calendar(nc.control)
    obs.en <- length(time)
    obs.st <- obs.en-105
    clim.st <- 51       
    clim.en <- 80
  control.past[[g]] <- get.arctic.average(nc.control,obs.mask,obs.st,obs.en,clim.st,clim.en)
  
  ##Natural Runs
  base.dir <- '/storage/data/climate/downscale/BCCAQ2/CMIP5/global/historicalNAT/annual/'
  scenario <- 'historicalNat' 
  gcm.dir <- paste0(base.dir,gcm)  
  file.natural <- list.files(path=gcm.dir,pattern=paste0('tas_ann_',obs.grid,'_grid_',gcm,'_',scenario,'*'),full.name=T)
  print(file.natural)
  nc.natural <- nc_open(file.natural)
  time <- netcdf.calendar(nc.natural)
    obs.st <- grep('1900',time)
    obs.en <- grep('2005',time)
    clim.st <- grep('1951',time[obs.st:obs.en])
    clim.en <- grep('1980',time[obs.st:obs.en])
  natural.past[[g]] <- get.arctic.average(nc.natural,obs.mask,obs.st,obs.en,clim.st,clim.en)

  ##RCP85 Runs
  base.dir <- '/storage/data/climate/downscale/BCCAQ2/CMIP5/global/annual/'
  scenario <- 'rcp85'
  gcm.dir <- paste0(base.dir,gcm)  
  file.rcp85 <- list.files(path=gcm.dir,pattern=paste0('tas_ann_',obs.grid,'_grid_',gcm,'_historical\\+',scenario,'*'),full.name=T)
  print(file.rcp85)
  nc.rcp85 <- nc_open(file.rcp85)
  time <- netcdf.calendar(nc.rcp85)
    obs.st <- grep('1900',time)
    obs.en <- grep('2005',time)
    clim.st <- grep('1951',time[obs.st:obs.en])
    clim.en <- grep('1980',time[obs.st:obs.en])
  rcp85.past[[g]] <- get.arctic.average(nc.rcp85,obs.mask,obs.st,obs.en,clim.st,clim.en)

  ##my.past <- list(arctic=ag.past$arctic,global=ag.past$global)      
  ##my.global <- list(time=1:length(ag.past$global),data=ag.past$global)
  ##my.arctic <- list(time=1:length(ag.past$arctic),data=ag.past$arctic)
  ##past.fit <- lm(arctic~global,data=my.past)

}

  control.arctic <- matrix(unlist(lapply(control.past,function(x){return(x$arctic)})),nrow=12,ncol=106,byrow=T) 
  control.arctic.df <- list(time=1:106,data=apply(control.arctic,2,mean))
  control.arctic.zyp <- zyp.sen(data~time,data=control.arctic.df)
  control.arctic.trend <- c(round(control.arctic.zyp$coefficients[2]*100,2),round(confint(control.arctic.zyp,level=0.9)[2,]*100,2))

  control.global <- matrix(unlist(lapply(control.past,function(x){return(x$global)})),nrow=12,ncol=106,byrow=T)
  control.global.df <- list(time=1:106,data=apply(control.global,2,mean))
  control.global.zyp <- zyp.sen(data~time,data=control.global.df)
  control.global.trend <- c(round(control.global.zyp$coefficients[2]*100,2),round(confint(control.global.zyp,level=0.9)[2,]*100,2))

  natural.arctic <- matrix(unlist(lapply(natural.past,function(x){return(x$arctic)})),nrow=12,ncol=106,byrow=T) 
  natural.arctic.df <- list(time=1:106,data=apply(natural.arctic,2,mean))
  natural.arctic.zyp <- zyp.sen(data~time,data=natural.arctic.df)
  natural.arctic.trend <- c(round(natural.arctic.zyp$coefficients[2]*100,2),round(confint(natural.arctic.zyp,level=0.9)[2,]*100,2))

  natural.global <- matrix(unlist(lapply(natural.past,function(x){return(x$global)})),nrow=12,ncol=106,byrow=T)
  natural.global.df <- list(time=1:106,data=apply(natural.global,2,mean))
  natural.global.zyp <- zyp.sen(data~time,data=natural.global.df)
  natural.global.trend <- c(round(natural.global.zyp$coefficients[2]*100,2),round(confint(natural.global.zyp,level=0.9)[2,]*100,2))

  rcp85.arctic <- matrix(unlist(lapply(rcp85.past,function(x){return(x$arctic)})),nrow=12,ncol=106,byrow=T) 
  rcp85.arctic.df <- list(time=1:106,data=apply(rcp85.arctic,2,mean))
  rcp85.arctic.zyp <- zyp.sen(data~time,data=rcp85.arctic.df)
  rcp85.arctic.trend <- c(round(rcp85.arctic.zyp$coefficients[2]*100,2),round(confint(rcp85.arctic.zyp,level=0.9)[2,]*100,2))

  rcp85.global <- matrix(unlist(lapply(rcp85.past,function(x){return(x$global)})),nrow=12,ncol=106,byrow=T)
  rcp85.global.df <- list(time=1:106,data=apply(rcp85.global,2,mean))
  rcp85.global.zyp <- zyp.sen(data~time,data=rcp85.global.df)
  rcp85.global.trend <- c(round(rcp85.global.zyp$coefficients[2]*100,2),round(confint(rcp85.global.zyp,level=0.9)[2,]*100,2))

  obs.global <- list(time=1:length(obs.data$global),data=obs.data$global)
  obs.global.zyp <- zyp.sen(data~time,data=obs.global)
  obs.global.trends <- c(round(obs.global.zyp$coefficients[2]*100,2),
                  round(confint(obs.global.zyp,level=0.9)[2,]*100,2))


  obs.arctic <- list(time=1:length(obs.data$arctic),data=obs.data$arctic)
  obs.arctic.zyp <- zyp.sen(data~time,data=obs.arctic)
  obs.arctic.trends <- c(round(obs.arctic.zyp$coefficients[2]*100,2),
                         round(confint(obs.arctic.zyp,level=0.9)[2,]*100,2))




plot.file <- paste0('/storage/data/projects/rci/data/nrcan/plots/',obs.grid,'_gcm_control_natural_trends_arctic_global_time_series.png')
png(plot.file,width=800,height=800)
  par(mar=c(5,5,4,2))
  par(mfrow=c(4,1))
  plot(c(),xlim=c(1900,2005),ylim=c(-2,3),xlab='Year',ylab='Temperature Anomaly (degC)',main='Control GCM Series',
                                                                          cex.lab=2,cex.axis=2,cex.main=2.5)
  apply(control.global,1,function(y,x){lines(x,y,col='gray')},1900:2005)
  apply(control.arctic,1,function(y,x){lines(x,y,col='gray')},1900:2005)
  lines(1900:2005,apply(control.global,2,mean),col='green',lwd=5)
  lines(1900:2005,apply(control.arctic,2,mean),col='blue',lwd=5)
  abline(h=0,lty=2)
  text(1920,2.5,paste0('Global Trend: ',control.global.trend[1],' (',control.global.trend[2],' to ',control.global.trend[3],')'),cex=2)
  points(1902,2.5,cex=3,pch=15,col='green')
  text(1980,2.5,paste0('Arctic Trend: ',control.arctic.trend[1],' (',control.arctic.trend[2],' to ',control.arctic.trend[3],')'),cex=2)
  points(1960,2.5,cex=3,pch=15,col='blue')
  box(which='plot')

  plot(c(),xlim=c(1900,2005),ylim=c(-2,3),xlab='Year',ylab='Temperature Anomaly (degC)',main='HistoricalNAT GCM Series',
                                                                          cex.lab=2,cex.axis=2,cex.main=2.5)
  apply(natural.global,1,function(y,x){lines(x,y,col='gray')},1900:2005)
  apply(natural.arctic,1,function(y,x){lines(x,y,col='gray')},1900:2005)
  lines(1900:2005,apply(natural.global,2,mean),col='green',lwd=5)
  lines(1900:2005,apply(natural.arctic,2,mean),col='blue',lwd=5)
  abline(h=0,lty=2)
  text(1920,2.5,paste0('Global Trend: ',natural.global.trend[1],' (',natural.global.trend[2],' to ',natural.global.trend[3],')'),cex=2)
  points(1902,2.5,cex=3,pch=15,col='green')
  text(1980,2.5,paste0('Arctic Trend: ',natural.arctic.trend[1],' (',natural.arctic.trend[2],' to ',natural.arctic.trend[3],')'),cex=2)
  points(1960,2.5,cex=3,pch=15,col='blue')
  box(which='plot')

  plot(c(),xlim=c(1900,2005),ylim=c(-2,3),xlab='Year',ylab='Temperature Anomaly (degC)',main='HistoricalALL GCM Series',
                                                                          cex.lab=2,cex.axis=2,cex.main=2.5)
  apply(rcp85.global,1,function(y,x){lines(x,y,col='gray')},1900:2005)
  apply(rcp85.arctic,1,function(y,x){lines(x,y,col='gray')},1900:2005)
  lines(1900:2005,apply(rcp85.global,2,mean),col='green',lwd=5)
  lines(1900:2005,apply(rcp85.arctic,2,mean),col='blue',lwd=5)
  abline(h=0,lty=2)
  text(1920,2.5,paste0('Global Trend: ',rcp85.global.trend[1],' (',rcp85.global.trend[2],' to ',rcp85.global.trend[3],')'),cex=2)
  points(1902,2.5,cex=3,pch=15,col='green')
  text(1980,2.5,paste0('Arctic Trend: ',rcp85.arctic.trend[1],' (',rcp85.arctic.trend[2],' to ',rcp85.arctic.trend[3],')'),cex=2)
  points(1960,2.5,cex=3,pch=15,col='blue')
  box(which='plot')

##Observations
  plot(c(),xlim=c(1900,2005),ylim=c(-2,3),xlab='Year',ylab='Temperature Anomaly (degC)',main=paste0(toupper(obs.grid),' Series'),
                                                                          cex.lab=2,cex.axis=2,cex.main=2.5)
  lines(1900:2005,obs.data$global,col='green',lwd=5)
  lines(1900:2005,obs.data$arctic,col='blue',lwd=5)

  abline(h=0,lty=2)
  text(1920,2.5,paste0('Global Trend: ',obs.global.trends[1],' (',obs.global.trends[2],' to ',obs.global.trends[3],')'),cex=2)
  points(1902,2.5,cex=3,pch=15,col='green')
  text(1980,2.5,paste0('Arctic Trend: ',obs.arctic.trends[1],' (',obs.arctic.trends[2],' to ',obs.arctic.trends[3],')'),cex=2)
  points(1960,2.5,cex=3,pch=15,col='blue')
  box(which='plot')

  dev.off()
