##Script to compare Arctic Anomalies with Global Anomalies for NASA GISTEMP

library(ncdf4)
library(scales)
library(PCICt)
library(nlme)
library(zyp)


source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)
source('/storage/home/ssobie/code/repos/ipcc/obs.read.r',chdir=T)
source('/storage/home/ssobie/code/repos/ipcc/arctic.average.r',chdir=T)


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

tas.data <- function(gcm.list,base.dir,obs.grid,scenario,chose.mask) {

  arctic.past <- vector(length=length(gcm.list),mode='list')
  global.past <- vector(length=length(gcm.list),mode='list')
  ag.slopes <- rep(0,length(gcm.list))
  ag.lower <- rep(0,length(gcm.list))
  ag.upper <- rep(0,length(gcm.list))
  ag.cors <- rep(0,length(gcm.list))
  for (g in seq_along(gcm.list)) {
    gcm <- gcm.list[g]
    print(gcm)
    gcm.dir <- paste0(base.dir,gcm,'/')
    if (scenario=='rcp85') {
      file.name <- list.files(path=gcm.dir,pattern=paste0('tas_ann_',obs.grid,'_grid_',gcm,'_historical\\+',scenario,'*'),full.name=T)
      print(file.name)
    } else {
        file.name <- list.files(path=gcm.dir,pattern=paste0('tas_ann_',obs.grid,'_grid_',gcm,'_',scenario,'*'),full.name=T)
    }
    nc <- nc_open(file.name)
    time <- netcdf.calendar(nc)

    if (scenario=='piControl') {
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

    obs.mask <- array(FALSE,c(dim(chose.mask)))
    obs.mask[,,1:dim(chose.mask)[3]] <- chose.mask
    obs.mask <- obs.mask[,,1:106]

    ag.gcm <- get.arctic.average(nc,obs.mask,obs.st,obs.en,clim.st,clim.en)

    arctic.past[[g]] <- ag.gcm$arctic
    global.past[[g]] <- ag.gcm$global

    my.past <- list(arctic=ag.gcm$arctic,global=ag.gcm$global)      
    past.fit <- lm(arctic~global,data=my.past)
    ag.slopes[g] <- past.fit$coefficients[2]  
    ag.cors[g] <- cor(ag.gcm$global,ag.gcm$arctic)
    ag.lower[g] <- confint(past.fit,level=0.9)[2,1]
    ag.upper[g] <- confint(past.fit,level=0.9)[2,2]

  }
  rv <- list(arctic=arctic.past,
             global=global.past,
             slopes=ag.slopes,
             lower=ag.lower,
             upper=ag.upper,
             cors=ag.cors)
  return(rv)
}


get.random.arctic.average <- function(nc,obs.mask,rlen) {

   lon <- ncvar_get(nc,'lon')
   lat <- ncvar_get(nc,'lat')
   time <- netcdf.calendar(nc)
   time.en <- length(time)

   var.name <- 'tas'
   tas.series <- ncvar_get(nc,var.name) - 273

   ar <- get.areas(lon,lat,c(dim(tas.series)[c(1,2)],rlen))
   wts <- ar[,1]/sum(ar[,1])

   ac.ix <- lat > 60

   arctic.wts <- ar[ac.ix,1]/sum(ar[ac.ix,1])
   global.wts <- ar[,1]/sum(ar[,1])
   N <- 1000
   ag.slopes <- rep(NA,N)
   arctic.vals <- c()
   global.vals <- c()

   arctic.means <- rep(NA,N)
   arctic.sds <- rep(NA,N)
   global.means <- rep(NA,N)
   global.sds <- rep(NA,N)
   ag.cor <- rep(NA,N)

   for (i in 1:N) {
     print(i)
     rix <- sample.int(time.en,rlen,replace=F)
     tas.subset <- tas.series[,,rix]
     tas.subset[obs.mask] <- NA
     tas.clim <- apply(tas.subset[,,51:80],c(1,2),mean,na.rm=T)
     tas.avg <- aperm(apply(tas.clim,c(1,2),function(x){rep(x,each=dim(tas.subset)[3])}),c(2,3,1))
     tas.anom <- tas.subset - tas.avg

     arctic.subset <- tas.anom[,ac.ix,]
     arctic.mean <- apply(apply(arctic.subset,c(1,3),weighted.mean,arctic.wts,na.rm=T),2,mean,na.rm=T)
     arctic.vals <- c(arctic.vals,arctic.mean)
     arctic.means[i] <- mean(arctic.mean)
     arctic.sds[i] <- sd(arctic.mean)

     global.subset <- tas.anom[,,]
     global.mean <- apply(apply(global.subset,c(1,3),weighted.mean,global.wts,na.rm=T),2,mean,na.rm=T)
     global.vals <- c(global.vals,global.mean)
     global.means[i] <- mean(global.mean)
     global.sds[i] <- sd(global.mean)

     my.past <- list(arctic=arctic.mean,global=global.mean)      
     past.fit <- lm(arctic~global,data=my.past)
     ag.slopes[i] <- past.fit$coefficients[2]  
     ag.cor[i] <- cor(global.mean,arctic.mean)
          
   }

   rv <- list(slopes=ag.slopes,
              cors = ag.cor,
              gm=global.means,
              gs=global.sds,
              am=arctic.means,
              as=arctic.sds)   
   return(rv)
}




random.tas.data <- function(gcm.list,base.dir,obs.grid,scenario,chose.mask) {

  control.slopes <- vector(length=length(gcm.list),mode='list')
  control.cors <- vector(length=length(gcm.list),mode='list') 
  control.arctic.means <- vector(length=length(gcm.list),mode='list') 
  control.arctic.sds <- vector(length=length(gcm.list),mode='list') 
  control.global.means <- vector(length=length(gcm.list),mode='list') 
  control.global.sds <- vector(length=length(gcm.list),mode='list') 

  rlen <- 106
  for (g in seq_along(gcm.list)) {
    gcm <- gcm.list[g]
    print(gcm)
    gcm.dir <- paste0(base.dir,gcm,'/')
    file.name <- list.files(path=gcm.dir,pattern=paste0('tas_ann_',obs.grid,'_grid_',gcm,'_',scenario,'*'),full.name=T)

    nc <- nc_open(file.name)
    time <- netcdf.calendar(nc)

    obs.mask <- array(FALSE,c(dim(chose.mask)))
    obs.mask[,,1:dim(chose.mask)[3]] <- chose.mask
    obs.mask <- obs.mask[,,1:rlen]
    ag.gcm <- get.random.arctic.average(nc,obs.mask,rlen)

    control.slopes[[g]] <- ag.gcm$slopes
    control.cors[[g]] <- ag.gcm$cors
    control.arctic.means[[g]] <- ag.gcm$am
    control.arctic.sds[[g]] <- ag.gcm$as
    control.global.means[[g]] <- ag.gcm$gm
    control.global.sds[[g]] <- ag.gcm$gs

    nc_close(nc)
  }
  rv <- list(slopes=control.slopes,
             cors=control.cors,
             ameans=control.arctic.means,
             asds=control.arctic.sds,
             gmeans=control.global.means,
             gsds=control.global.sds)
                 
  return(rv)
}






##*************************************************************************
##Obs
had.file <- '/storage/data/projects/rci/data/nrcan/hadcru4/HadCRUT.4.5.0.0.median.annual.nc'
hadcru4.data <- get.obs.means(file.name=had.file,var.name='tas','1900','2016')
hadcru4.to.fit <- list(arctic=hadcru4.data$arctic,global=hadcru4.data$global)
hadcru4.fit <- lm(arctic~global,data=hadcru4.to.fit)
hadcru4.slopes <- c(hadcru4.fit$coefficients[2],confint(hadcru4.fit,level=0.9)[2,])
hadcru4.mask <- hadcru4.data$mask

chose.mask <- hadcru4.mask
obs.grid <- 'hadcru4'

##GCMs

gcm.list <- c('ACCESS1-0','ACCESS1-3','CanESM2','CSIRO-Mk3-6-0','GFDL-ESM2M','IPSL-CM5A-MR','IPSL-CM5B-LR',
              'MPI-ESM-LR','MPI-ESM-MR','NorESM1-M','bcc-csm1-1','bcc-csm1-1-m','CNRM-CM5','inmcm4',
              'HadGEM2-CC','HadGEM2-ES','IPSL-CM5A-LR','MIROC5','MIROC-ESM','MIROC-ESM-CHEM')
nat.list <- c('ACCESS1-3','CanESM2','CSIRO-Mk3-6-0','GFDL-ESM2M','IPSL-CM5A-MR',
              'NorESM1-M','bcc-csm1-1','CNRM-CM5',
              'HadGEM2-ES','IPSL-CM5A-LR','MIROC-ESM','MIROC-ESM-CHEM')
lab.list <- c('ACC0','ACC3','Can','CSI','GFDLM','IPSLM','IPSLB',
              'MPIL','MPIM','Nor','bcc','bccm','CNRM','inm',
              'HCC','HES','IPSLL','MIR5','MIRE','MIRC')


base.dir <-'/storage/data/climate/downscale/BCCAQ2+PRISM/CMIP5/global/control/annual/'
control.data <- random.tas.data(gcm.list,base.dir,obs.grid,scenario='piControl',chose.mask)

base.dir <-'/storage/data/climate/downscale/BCCAQ2+PRISM/CMIP5/global/annual/'
rcp85.data <- tas.data(gcm.list,base.dir,obs.grid,scenario='rcp85',chose.mask)

base.dir <-'/storage/data/climate/downscale/BCCAQ2+PRISM/CMIP5/global/control/annual/'
control.run <- tas.data(gcm.list,base.dir,obs.grid,scenario='piControl',chose.mask)

base.dir <-'/storage/data/climate/downscale/BCCAQ2+PRISM/CMIP5/global/historicalNAT/annual/'
natural.run <- tas.data(nat.list,base.dir,obs.grid,scenario='historicalNat',chose.mask)
nat.sub <- gcm.list %in% nat.list

x <- seq(3,60,3)

png('/storage/data/projects/rci/data/nrcan/plots/slope.comparison.png',height=700,width=1000)
plot(c(),xlim=c(0,60),ylim=c(0,7),axes=F,cex=2,ylab='Slope',xlab='',cex.lab=1.5)
par(oma=c(10,5,2,2))
axis(1,at=seq(3,60,3),lab.list,las=2,cex=1.5,cex.axis=1.5)
axis(2,at=0:9,0:9,cex=1.5,cex.axis=1.5)
for (i in 1:20) {
 lx <- x[i]   
 boxplot(at=lx,control.data$slopes[[i]],add=T,axes=F)

 points(lx-0.5,rcp85.data$upper[i],pch=2,col='red',cex=1)
 points(lx-0.5,rcp85.data$slopes[i],pch='-',col='red',cex=3)
 points(lx-0.5,rcp85.data$lower[i],pch=6,col='red',cex=1)
 lines(c(lx-0.5,lx-0.5),c(rcp85.data$lower[i],rcp85.data$upper[i]),lwd=3,col='red')

 points(lx-1.5,control.run$upper[i],pch=2,col='blue',cex=1)
 points(lx-1.5,control.run$slopes[i],pch='-',col='blue',cex=3)
 points(lx-1.5,control.run$lower[i],pch=6,col='blue',cex=1)
 lines(c(lx-1.5,lx-1.5),c(control.run$lower[i],control.run$upper[i]),lwd=3,col='blue')
}
for (j in 1:length(nat.list)) {
 lc <- x[nat.sub]
 points(lc[j]-1,natural.run$upper[j],pch=2,col='green',cex=1)
 points(lc[j]-1,natural.run$slopes[j],pch='-',col='green',cex=3)
 points(lc[j]-1,natural.run$lower[j],pch=6,col='green',cex=1)
 lines(c(lc[j]-1,lc[j]-1),c(natural.run$lower[j],natural.run$upper[j]),lwd=3,col='green')
}
abline(v=seq(1,60,3))
legend('topleft',leg=c('Control','Natural','Historical'),col=c('blue','green','red'),pch=15,cex=1.5)
box(which='plot')
dev.off()
browser()

arctic.tm <- cbind(c(gcm.list,'5%','Mean','95%'),
                      rbind(arctic.trends,
                            apply(arctic.trends,2,quantile,0.05,na.rm=T),
                            apply(arctic.trends,2,mean,na.rm=T),
                            apply(arctic.trends,2,quantile,0.95,na.rm=T)))

write.file <- paste0('/storage/data/projects/rci/data/nrcan/nasa_gistemp/slopes/',obs.grid,'_arctic_trends_Control.csv')
write.table(arctic.tm,file=write.file,quote=F,row.name=F,col.name=F,sep=',')




if(1==0) {
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

