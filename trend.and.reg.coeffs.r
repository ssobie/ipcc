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
gistemp.to.fit <- list(arctic=gistemp.data$arctic,global=gistemp.data$global)
gistemp.fit <- lm(arctic~global,data=gistemp.to.fit)
gistemp.slopes <- c(gistemp.fit$coefficients[2],confint(gistemp.fit,level=0.9)[2,])
gistemp.mask <- gistemp.data$mask

had.file <- '/storage/data/projects/rci/data/nrcan/hadcru4/HadCRUT.4.5.0.0.median.annual.nc'
hadcru4.data <- get.obs.means(file.name=had.file,var.name='tas','1900','2005')
hadcru4.to.fit <- list(arctic=hadcru4.data$arctic,global=hadcru4.data$global)
hadcru4.fit <- lm(arctic~global,data=hadcru4.to.fit)
hadcru4.slopes <- c(hadcru4.fit$coefficients[2],confint(hadcru4.fit,level=0.9)[2,])
hadcru4.mask <- hadcru4.data$mask

##---------------------------------------------------------------------
##*************************************************************************
##GCMs

##Control List
##gcm.list <- c('ACCESS1-0','ACCESS1-3','BNU-ESM','CanESM2','CSIRO-Mk3-6-0','GFDL-ESM2M','IPSL-CM5A-MR','IPSL-CM5B-LR',
##              'MPI-ESM-LR','MPI-ESM-MR','NorESM1-M','bcc-csm1-1','bcc-csm1-1-m','CNRM-CM5','inmcm4',
##              'HadGEM2-CC','HadGEM2-ES','IPSL-CM5A-LR','MIROC5','MIROC-ESM','MIROC-ESM-CHEM')

##Natural List
##gcm.list <- c('ACCESS1-3','CanESM2','CESM1-CAM5','CSIRO-Mk3-6-0','GFDL-ESM2M','IPSL-CM5A-MR',
##              'MIROC-ESM','NorESM1-M','bcc-csm1-1','CCSM4','CNRM-CM5','GFDL-CM3',
##              'HadGEM2-ES','IPSL-CM5A-LR','MIROC-ESM-CHEM','MRI-CGCM3')


##RCP26 List
##gcm.list <- c('bcc-csm1-1','bcc-csm1-1-m','CanESM2','CCSM4','CNRM-CM5','CSIRO-Mk3-6-0','GFDL-CM3','GFDL-ESM2G',
##              'GFDL-ESM2M','HadGEM2-AO','HadGEM2-ES','IPSL-CM5A-LR','IPSL-CM5A-MR','MIROC5','MIROC-ESM','MIROC-ESM-CHEM',
##              'MPI-ESM-MR','MPI-ESM-LR','MRI-CGCM3','NorESM1-M')

##gcm.list <- c('ACCESS1-0','ACCESS1-3','bcc-csm1-1','bcc-csm1-1-m','CanESM2','CCSM4','CNRM-CM5','CSIRO-Mk3-6-0',
##              'GFDL-CM3','GFDL-ESM2G','GFDL-ESM2M','HadGEM2-AO','HadGEM2-CC','HadGEM2-ES','inmcm4','IPSL-CM5A-LR',
##              'IPSL-CM5A-MR','IPSL-CM5B-LR','MIROC5','MIROC-ESM','MIROC-ESM-CHEM','MPI-ESM-MR', 'MPI-ESM-LR',
##              'MRI-CGCM3','NorESM1-M')


##Common List
gcm.list <- c('ACCESS1-3','CanESM2','CSIRO-Mk3-6-0','GFDL-ESM2M','IPSL-CM5A-MR',
              'MIROC-ESM','NorESM1-M','bcc-csm1-1','CNRM-CM5',
              'HadGEM2-ES','IPSL-CM5A-LR','MIROC-ESM-CHEM')

##Control/Natural Runs
##base.dir <- '/storage/data/climate/downscale/BCCAQ2/CMIP5/global/historicalNAT/annual/'
##scenario <- 'historicalNat' ##'piControl'

##RCP Runs
base.dir <- '/storage/data/climate/downscale/BCCAQ2/CMIP5/global/annual/'
scenario <- 'rcp85'

obs.data <- gistemp.data
obs.grid <- 'gistemp'
obs.slopes <- gistemp.slopes
chose.mask <- gistemp.mask

obs.mask <- array(FALSE,c(dim(chose.mask)[c(1,2)],106))
obs.mask[,,1:dim(chose.mask)[3]] <- gistemp.mask

##---------------------------------------------------------------------

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
  gcm.dir <- paste0(base.dir,gcm,'/')  
  file.name <- list.files(path=gcm.dir,pattern=paste0('tas_ann_',obs.grid,'_grid_',gcm,'_historical\\+',scenario,'*'),full.name=T)
  ##file.name <- list.files(path=gcm.dir,pattern=paste0('tas_ann_',obs.grid,'_grid_',gcm,'_',scenario,'*'),full.name=T)
  
  nc <- nc_open(file.name)
  time <- netcdf.calendar(nc)

  if (scenario =='piControl') {
    obs.en <- length(time)
    obs.st <- obs.en-105
    clim.st <- 51       
    clim.en <- 80
  } else {     
    obs.st <- grep('1900',time)
    obs.en <- grep('2005',time)
    clim.st <- grep('1951',time[obs.st:obs.en])
    clim.en <- grep('1980',time[obs.st:obs.en])
  }

  ag.past <- get.arctic.average(nc,obs.mask,obs.st,obs.en,clim.st,clim.en)

  arctic.past[[g]] <- ag.past$arctic
  global.past[[g]] <- ag.past$global

  my.past <- list(arctic=ag.past$arctic,global=ag.past$global)      
  my.global <- list(time=1:length(ag.past$global),data=ag.past$global)
  my.arctic <- list(time=1:length(ag.past$arctic),data=ag.past$arctic)
  past.fit <- lm(arctic~global,data=my.past)
  
  global.zyp <- zyp.sen(data~time,data=my.global)
  global.trends[g,1] <- global.zyp$coefficients[2]
  global.trends[g,c(2,3)] <- confint(global.zyp,level=0.9)[2,]

  arctic.zyp <- zyp.sen(data~time,data=my.arctic)
  arctic.trends[g,1] <- arctic.zyp$coefficients[2]
  arctic.trends[g,c(2,3)] <- confint(arctic.zyp,level=0.9)[2,]

  rg.past.slope[g] <- past.fit$coefficients[2]            
  conf.ints <- confint(past.fit,level=0.9)
  rg.past.upper[g] <- conf.ints[2,2]
  rg.past.lower[g] <- conf.ints[2,1]

}

all.global <- unlist(global.past)
all.arctic <- unlist(arctic.past)
all.data <- list(arctic=all.arctic,global=all.global)
all.fit <- lm(arctic~global,data=all.data)
all.slopes <- c(all.fit$coefficients[2],confint(all.fit,level=0.9)[2,])
print(all.slopes)

model.slopes <- round(rg.past.slope,2)
model.lower <- round(rg.past.lower,2)
model.upper <- round(rg.past.upper,2)
model.matrix <- cbind(model.slopes,model.lower,model.upper)

model.slopes <- rbind(c('GCM','Slope','5%','95%'),
                      cbind(c(gcm.list,'5%','Mean','95%',' ',obs.grid),
                      rbind(model.matrix,
                            round(apply(model.matrix,2,quantile,0.05,na.rm=T),2),
                            round(apply(model.matrix,2,mean,na.rm=T),2),
                            round(apply(model.matrix,2,quantile,0.95,na.rm=T),2),
                            c(' ',' ',' '),
                            round(obs.slopes,2))))

##write.file <- paste0('/storage/data/projects/rci/data/nrcan/nasa_gistemp/slopes/',obs.grid,'_slopes_',scenario,'.2005.matching.csv')
##write.table(model.slopes,file=write.file,quote=F,row.name=F,col.name=F,sep=',')

 
  obs.global <- list(time=1:length(obs.data$global),data=obs.data$global)
  obs.global.zyp <- zyp.sen(data~time,data=obs.global)
  obs.global.trends <- c(round(obs.global.zyp$coefficients[2]*100,2),
                  round(confint(obs.global.zyp,level=0.9)[2,]*100,2))

global.tm <- rbind(c('GCM','Trend','5%','95%'),
                      cbind(c(gcm.list,'5%','Mean','95%',' ',obs.grid),
                      rbind(round(global.trends*100,2),
                            round(apply(global.trends*100,2,quantile,0.05,na.rm=T),2),
                            round(apply(global.trends*100,2,mean,na.rm=T),2),
                            round(apply(global.trends*100,2,quantile,0.95,na.rm=T),2),  
                            c(' ',' ',' '),
                            obs.global.trends)))

##write.file <- paste0('/storage/data/projects/rci/data/nrcan/nasa_gistemp/slopes/',obs.grid,'_global_trends_',scenario,'.2005.csv')
##write.table(global.tm,file=write.file,quote=F,row.name=F,col.name=F,sep=',')

  obs.arctic <- list(time=1:length(obs.data$arctic),data=obs.data$arctic)
  obs.arctic.zyp <- zyp.sen(data~time,data=obs.arctic)
  obs.arctic.trends <- c(round(obs.arctic.zyp$coefficients[2]*100,2),
                         round(confint(obs.arctic.zyp,level=0.9)[2,]*100,2))

arctic.tm <- rbind(c('GCM','Trend','5%','95%'),
                      cbind(c(gcm.list,'5%','Mean','95%',' ',obs.grid),
                      rbind(round(arctic.trends*100,2),
                            round(apply(arctic.trends*100,2,quantile,0.05,na.rm=T),2),
                            round(apply(arctic.trends*100,2,mean,na.rm=T),2),
                            round(apply(arctic.trends*100,2,quantile,0.95,na.rm=T),2),  
                            c(' ',' ',' '),
                            obs.arctic.trends)))

##write.file <- paste0('/storage/data/projects/rci/data/nrcan/nasa_gistemp/slopes/',obs.grid,'_arctic_trends_',scenario,'.2005.csv')
##write.table(arctic.tm,file=write.file,quote=F,row.name=F,col.name=F,sep=',')



if(1==0) {
plot.file <- '/storage/data/projects/rci/data/nrcan/plots/gistemp_control_arctic_global_anomaly_matching.png'
  png(plot.file,width=800,height=800)
  par(mar=c(5,5,4,2))
  plot(obs.data$global,obs.data$arctic,xlim=c(-1,2),ylim=c(-2,4),cex=2,main='NASA GISTEMP and piControl GCM Anomalies',
                   xlab='Global Avg TAS Anomalies (degC)',ylab='Arctic Avg TAS Anomalies (degC)',
                   cex.lab=2,cex.axis=2,cex.main=2.5,pch=16,col='red')
  mapply(FUN=points,global.past,arctic.past,pch=18,col='black',cex=2)
  points(obs.data$global,obs.data$arctic,cex=2,pch=16,col='red')

  abline(h=0,v=0,col='gray',lty=2,lwd=2)
  legend('topleft',leg=c('GISTEMP','GCM'),col=c('red','black'),pch=c(16,18),cex=2)
  box(which='plot')
  dev.off()
}



