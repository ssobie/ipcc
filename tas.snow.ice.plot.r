##Script to make a multiplot figure of Snow and Ice data in the Canadian Arctic

library(ncdf4)
library(scales)
library(PCICt)

source('/storage/home/ssobie/code/repos/ipcc/obs.read.r',chdir=T)
source('/storage/home/ssobie/code/repos/ipcc/arctic.average.r',chdir=T)

ice.plot <- function(tas,ice.data,colour,bounds=FALSE) {
  glen <- ncol(tas)-1
  newx <- seq(-100, 400, by=0.5)
  ice.fit <- matrix(nrow=length(newx),ncol=glen)

  tas.anoms <- vector(mode='list',length=glen)
  ice.anoms <- vector(mode='list',length=glen)

  for (i in 2:glen) {
     ##tas.mean <- mean(tas[c(11,12,13),i])
     tas.anom <- tas[,i] ##-tas.mean         
     tas.anoms[[i]] <- tas.anom
     ix <- order(tas.anom)

     ##ice.mean <- mean(ice.data[c(11,12,13),i])
     ice.anom <- ice.data[,i] ## - ice.mean) ##/ice.mean
     ice.anoms[[i]] <- ice.anom
     ##points(ice.data[ix,i],tas.anom[ix],pch=16,col='gray')
     ##points(ice.anom[ix],tas.anom[ix],pch=16,col='gray')
     ##lines(ice.data[ix,i],tas.anom[ix])

     ##x <- as.vector(ice.data[ix,i])
     x <- as.vector(ice.anom[ix])
     y <- as.vector(tas.anom[ix])
     df <- data.frame(x = x,
                 y = y)
     mod <- lm(y ~ x, data = df)
     ##abline(mod,col=colour)
     preds <- predict(mod, newdata = data.frame(x=newx))                     
     ice.fit[,i] <- preds
  }

  ice.unlist <- unlist(ice.anoms)
  ix <- which(is.na(ice.unlist))
  ice.unlist[ix] <- mean(ice.unlist,na.rm=T)
  tas.unlist <- unlist(tas.anoms)
  tx <- which(is.na(tas.unlist))
  tas.unlist[tx] <- mean(tas.unlist,na.rm=T)
  
  test <- chull(ice.unlist,tas.unlist)
  hpts <- c(test,test[1])
  polygon(ice.unlist[hpts],tas.unlist[hpts],col=alpha(colour,0.3),border=alpha(colour,0.3))

  ice.lower <- apply(ice.fit,1,quantile,0.05,na.rm=T)
  ice.upper <- apply(ice.fit,1,quantile,0.95,na.rm=T)
  ice.mean <- apply(ice.fit,1,mean,na.rm=T)

if (bounds) {
##  lines(newx,ice.lower,col='blue',lwd=3)
##  lines(newx,ice.mean,col=colour,lwd=3)
##  lines(newx,ice.upper,col='blue',lwd=3)
  }
}##Ice plot


##Arctic/Global Averages

tas.data <- function(gcm.list,base.dir,obs.grid,scenario,chose.mask) {

  arctic.past <- vector(length=length(gcm.list),mode='list')
  global.past <- vector(length=length(gcm.list),mode='list')

  for (g in seq_along(gcm.list)) {
    gcm <- gcm.list[g]
    print(gcm)
    gcm.dir <- paste0(base.dir,gcm,'/')
    all.files <- list.files(path=gcm.dir,pattern=paste0('tas_ann_',obs.grid,'_grid_',gcm,'_historical\\+',scenario),full.name=T)
    file.name <- all.files[grep('18500101',all.files)]
    print(file.name)
    nc <- nc_open(file.name)
    time <- netcdf.calendar(nc)

    obs.st <- grep('1861',time)   ##grep('1900',time) ##
    if (length(obs.st) ==0) {
     obs.st <- 1
    }
    
    obs.en <- grep('2099',time)
    clim.st <- grep('1951',time[obs.st:obs.en])## grep('1861',time[obs.st:obs.en]) ##
    if (length(clim.st) ==0) {
      obs.st <- 1
    }
    clim.en <- grep('1980',time[obs.st:obs.en]) 
    nst <- grep('1900',time[obs.st:obs.en]) 
    nen <- length(time[obs.st:obs.en])
    obs.mask <- array(FALSE,c(dim(chose.mask)[c(1,2)],obs.en-obs.st+1))
    obs.mask[,,1:dim(chose.mask)[3]] <- chose.mask

    ag.gcm <- get.arctic.average(nc,obs.mask,obs.st,obs.en,clim.st,clim.en)
    
    arctic.offset <- mean(ag.gcm$arctic[grep('1951',time[obs.st:obs.en]):grep('1980',time[obs.st:obs.en])]) - 
          mean(ag.gcm$arctic[grep('1862',time[obs.st:obs.en]):grep('1880',time[obs.st:obs.en])])
    print(arctic.offset)      
    global.offset <- mean(ag.gcm$global[grep('1951',time[obs.st:obs.en]):grep('1980',time[obs.st:obs.en])]) - 
          mean(ag.gcm$global[grep('1862',time[obs.st:obs.en]):grep('1880',time[obs.st:obs.en])])
    print(global.offset)

    arctic.past[[g]] <- ag.gcm$arctic[nst:nen] + arctic.offset
    global.past[[g]] <- ag.gcm$global[nst:nen] + global.offset


  } 
  rv <- list(arctic=arctic.past,
               global=global.past)

  return(rv)                  
}


##Read the snow/ice data
snow.ice.csv <- function(rcp) {

  sic.temps <- read.csv(paste0(read.dir,'tas_sic_rcp',rcp,'_arctic_1860.csv'),as.is=T,header=F)
  ice.free.nep <- read.csv(paste0(read.dir,'ice_free_nep_rcp',rcp,'_1860.csv'),as.is=T,header=F)
  ice.free.caa <- read.csv(paste0(read.dir,'ice_free_caa_rcp',rcp,'_1860.csv'),as.is=T,header=F)

  snc.temps <- read.csv(paste0(read.dir,'tas_snc_rcp',rcp,'_arctic_1860.csv'),as.is=T,header=F)
  snow.cover.ea <- read.csv(paste0(read.dir,'snow_cover_eurasia_rcp',rcp,'_1860.csv'),as.is=T,header=F)
  snow.cover.na <- read.csv(paste0(read.dir,'snow_cover_north_america_rcp',rcp,'_1860.csv'),as.is=T,header=F)

  models <- arctic.temps[1,-1]
  runs <- arctic.temps[2,-1]

  gcms <- paste0(models,'-r',runs)
  ex <- 242
  ice.tas.data <- sic.temps[43:ex,]
  ice.nep.data <- ice.free.nep[43:ex,]
  ice.caa.data <- ice.free.caa[43:ex,]

  snow.tas.data <- snc.temps[43:ex,]
  snow.ea.data <- snow.cover.ea[43:ex,]
  snow.na.data <- snow.cover.na[43:ex,]

  ice.tas <- round(apply(ice.tas.data,2,as.numeric),1)
  ice.nep <- round(apply(ice.nep.data,2,as.numeric),1)
  ice.caa <- round(apply(ice.caa.data,2,as.numeric),1)

  snow.tas <- round(apply(snow.tas.data,2,as.numeric),1)
  snow.ea <- round(apply(snow.ea.data,2,as.numeric),1)
  snow.na <- round(apply(snow.na.data,2,as.numeric),1)

  rv <- list(icetas=ice.tas,nep=ice.nep,caa=ice.caa,
             snowtas=snow.tas,sea=snow.ea,sna=snow.na)
  return(rv)           

}



##----------------------------------------------------------------
##Observations

had.file <- '/storage/data/projects/rci/data/nrcan/hadcru4/HadCRUT.4.5.0.0.median.annual.nc'
hadcru4.data <- get.obs.means(file.name=had.file,var.name='tas','1901','2016')
hadcru4.to.fit <- list(arctic=hadcru4.data$arctic[1:116],global=hadcru4.data$global[1:116])
hadcru4.fit <- lm(arctic~global,data=hadcru4.to.fit)
hadcru4.slopes <- c(hadcru4.fit$coefficients[2],confint(hadcru4.fit,level=0.9)[2,])
hadcru4.mask <- hadcru4.data$mask

chose.mask <- hadcru4.mask
obs.grid <- 'hadcru4'

##----------------------------------------------------------------
##GCM data
base.dir <-'/storage/data/climate/downscale/BCCAQ2+PRISM/CMIP5/global/annual/'

##ACCESS1-0
gcm.list <- c('bcc-csm1-1','ACCESS1-0','ACCESS1-3','bcc-csm1-1-m','CanESM2','CCSM4','CNRM-CM5','CSIRO-Mk3-6-0',
              'GFDL-CM3','GFDL-ESM2G','GFDL-ESM2M','HadGEM2-AO','HadGEM2-ES','inmcm4','IPSL-CM5A-LR',
              'IPSL-CM5A-MR','IPSL-CM5B-LR','MIROC5','MIROC-ESM','MIROC-ESM-CHEM','MPI-ESM-MR', 'MPI-ESM-LR',
              'MRI-CGCM3','NorESM1-M','GISS-E2-H','GISS-E2-R','NorESM1-ME')

gcm.len <- length(gcm.list)

rcp85.data <- tas.data(gcm.list,base.dir,obs.grid,scenario='rcp85',chose.mask)
rcp45.data <- tas.data(gcm.list,base.dir,obs.grid,scenario='rcp45',chose.mask)


##RCP26 List
rcp26.list <- c('bcc-csm1-1','bcc-csm1-1-m','CanESM2','CCSM4','CNRM-CM5','CSIRO-Mk3-6-0','GFDL-CM3','GFDL-ESM2G',
              'GFDL-ESM2M','HadGEM2-AO','HadGEM2-ES','IPSL-CM5A-LR','IPSL-CM5A-MR','MIROC5','MIROC-ESM','MIROC-ESM-CHEM',
              'MPI-ESM-MR','MPI-ESM-LR','MRI-CGCM3','NorESM1-M','NorESM1-ME','GISS-E2-H','GISS-E2-R')
gcm26.len <- length(rcp26.list)
rcp26.data <- tas.data(rcp26.list,base.dir,obs.grid,scenario='rcp26',chose.mask)

all.rcps <- list(arctic=c(unlist(rcp26.data$arctic),unlist(rcp45.data$arctic),unlist(rcp85.data$arctic)),
                 global=c(unlist(rcp26.data$global),unlist(rcp45.data$global),unlist(rcp85.data$global)))
all.rcp85 <- list(arctic=unlist(rcp85.data$arctic),global=unlist(rcp85.data$global))
rcps.loess <- loess(arctic~global,data=all.rcps,span=0.75)
smoothed.rcps <- predict(rcps.loess)
ey <- 200
  global.avg.85 <- apply(matrix(unlist(rcp85.data$global),nrow=gcm.len,ncol=ey,byrow=T),2,mean)
  global.avg.45 <- apply(matrix(unlist(rcp45.data$global),nrow=gcm.len,ncol=ey,byrow=T),2,mean)
  global.avg.26 <- apply(matrix(unlist(rcp26.data$global),nrow=gcm26.len,ncol=ey,byrow=T),2,mean)

  arctic.avg.85 <- apply(matrix(unlist(rcp85.data$arctic),nrow=gcm.len,ncol=ey,byrow=T),2,mean)
  arctic.avg.45 <- apply(matrix(unlist(rcp45.data$arctic),nrow=gcm.len,ncol=ey,byrow=T),2,mean)
  arctic.avg.26 <- apply(matrix(unlist(rcp26.data$arctic),nrow=gcm26.len,ncol=ey,byrow=T),2,mean)

##----------------------------------------------------------------
##Ice data - RCP85

proj.dir <- '/storage/data/projects/rci/data/nrcan/'

read.dir <- paste0(proj.dir,'ice_snow_data/')
plot.dir <- paste0(proj.dir,'plots/')

arctic.temps <- read.csv(paste0(read.dir,'arctic_temps.csv'),as.is=T,header=F)
ice.free.nep <- read.csv(paste0(read.dir,'ice_free_nep.csv'),as.is=T,header=F)
ice.free.caa <- read.csv(paste0(read.dir,'ice_free_caa.csv'),as.is=T,header=F)

snow.cover.ea <- read.csv(paste0(read.dir,'snow_cover_eurasia.csv'),as.is=T,header=F)
snow.cover.na <- read.csv(paste0(read.dir,'snow_cover_north_america.csv'),as.is=T,header=F)

models <- arctic.temps[1,-1]
runs <- arctic.temps[2,-1]

ice.tas.data <- arctic.temps[3:27,]
ice.nep.data <- ice.free.nep[3:27,]
ice.caa.data <- ice.free.caa[3:27,]

snow.ea.data <- snow.cover.ea[3:27,]
snow.na.data <- snow.cover.na[3:27,]

ice.tas <- round(apply(ice.tas.data,2,as.numeric),1)-273
ice.nep <- round(apply(ice.nep.data,2,as.numeric),1)
ice.caa <- round(apply(ice.caa.data,2,as.numeric),1)

snow.ea <- round(apply(snow.ea.data,2,as.numeric),1)
snow.na <- round(apply(snow.na.data,2,as.numeric),1)

##----------------------------------------------------------------
##Ice data - rcp26,45

si.26 <- snow.ice.csv('26')
si.45 <- snow.ice.csv('45')
si.85 <- snow.ice.csv('85')

##****************************************************************************
##Create the two figure plot

  plot.file <- paste0(plot.dir,'arctic_tas_and_ice_free_days_all_rcps_4figure_new2.png')
  png(plot.file,width=1400,height=800)
  ##par(mfrow=c(1,2))
  layout(mat=cbind(matrix(1,nrow=2,ncol=2),matrix(2:5,nrow=2,ncol=2,byrow=T)))
  par(mar=c(5,5,4,1))
  plot(hadcru4.data$global,hadcru4.data$arctic,xlim=c(-1,6.2),ylim=c(-2,16),pch=0,col='black',cex=2,
                   main='GCM, HadCRU4 Temperature Anomalies',
                   xlab='Global Temperature Anomalies (\u00B0C)',ylab='Arctic Temperature Anomalies (\u00B0C)',
                   cex.lab=2,cex.axis=2,cex.main=2,axes=F)
  axis(1,seq(-1,6,1),seq(-1,6,1),cex.axis=2)
  axis(2,seq(0,15,5),seq(0,15,5),cex.axis=2)

  test <- chull(unlist(rcp85.data$global),unlist(rcp85.data$arctic))
  hpts <- c(test,test[1])
  polygon(unlist(rcp85.data$global)[hpts],unlist(rcp85.data$arctic)[hpts],col=alpha('red',0.3),border=alpha('red',0.3))
  test <- chull(unlist(rcp45.data$global),unlist(rcp45.data$arctic))
  hpts <- c(test,test[1])
  polygon(unlist(rcp45.data$global)[hpts],unlist(rcp45.data$arctic)[hpts],col=alpha('green',0.3),border=alpha('green',0.3))
  test <- chull(unlist(rcp26.data$global),unlist(rcp26.data$arctic))
  hpts <- c(test,test[1])
  polygon(unlist(rcp26.data$global)[hpts],unlist(rcp26.data$arctic)[hpts],col=alpha('blue',0.3),border=alpha('blue',0.3))
  points(hadcru4.data$global,hadcru4.data$arctic,pch=0,col='black',cex=2)
  abline(coef=hadcru4.fit$coefficients,col='black',lwd=4)
  g <- order(all.rcps$global)
  lines(all.rcps$global[g],smoothed.rcps[g],col='darkgray',lwd=5)
  abline(h=seq(-5,20,5),col='gray',lty=2,lwd=2)

  yrs <- 1901:2100
  ix <- seq(50,2100,50)
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

  ##-----------------------------------------------------------------------
  ##Snow/Ice Plots

  par(mar=c(5,1,4,5))
  plot(c(),xlim=c(-50,300),ylim=c(-2,16),cex=2,main='Northeast Passage',
                   xlab='Ice Free Day Anomaly',ylab='',
                   cex.lab=2,cex.axis=2,cex.main=2,axes=F)
  axis(1,seq(-50,400,50),seq(-50,400,50),cex.axis=2)
  axis(4,seq(0,15,5),seq(0,15,5),cex.axis=2)
  mtext("Arctic Temperature Anomalies (\u00B0C)",side=4,line=3,cex=1.5)

  ice.plot(si.85$icetas,si.85$nep,'red',bounds=TRUE)
  ice.plot(si.45$icetas,si.45$nep,'green') 
  ice.plot(si.26$icetas,si.26$nep,'blue')
  
  abline(h=seq(-5,20,5),col='gray',lty=2,lwd=2)
  legend('topleft',leg=c('RCP2.6','RCP4.5','RCP8.5'),
          col=c(alpha('blue',0.3),alpha('green',0.3),alpha('red',0.3)),
          pch=c(15,15,15),cex=2)

  ## text(x=50,y=11,'North East\n Passage',cex=1.5) 
  ##text(x=-1,y=11,'North East\n Passage',cex=2)
  box(which='plot')

  plot(c(),xlim=c(-50,300),ylim=c(-2,16),cex=2,main='Canadian Arctic Archipelago',
                   xlab='Ice Free Day Anomaly',ylab='',
                   cex.lab=2,cex.axis=2,cex.main=2,axes=F)
  axis(1,seq(-50,400,50),seq(-50,400,50),cex.axis=2)
  axis(4,seq(0,15,5),seq(0,15,5),cex.axis=2)
  mtext("Arctic Temperature Anomalies (\u00B0C)",side=4,line=3,cex=1.5)

  ice.plot(si.85$icetas,si.85$caa,'red')
  ice.plot(si.45$icetas,si.45$caa,'green')
  ice.plot(si.26$icetas,si.26$caa,'blue')

  abline(h=seq(-5,20,5),col='gray',lty=2,lwd=2)
  legend('topleft',leg=c('RCP2.6','RCP4.5','RCP8.5'),
          col=c(alpha('blue',0.3),alpha('green',0.3),alpha('red',0.3)),
          pch=c(15,15,15),cex=2)
  box(which='plot')

  plot(c(),xlim=c(-100,25),ylim=c(-2,16),cex=2,main='North America',
                   xlab='Snow Cover Day Anomaly',ylab='Arctic Temperature Anomalies (\u00B0C)',
                   cex.lab=2,cex.axis=2,cex.main=2,axes=F)
  axis(1,seq(-100,40,20),seq(-100,40,20),cex.axis=2)
  axis(4,seq(0,15,5),seq(0,15,5),cex.axis=2)
  mtext("Arctic Temperature Anomalies (\u00B0C)",side=4,line=3,cex=1.5)

  ice.plot(si.85$snowtas,si.85$sna,'red')
  ice.plot(si.45$snowtas,si.45$sna,'green')
  ice.plot(si.26$snowtas,si.26$sna,'blue')

  abline(h=seq(-5,20,5),col='gray',lty=2,lwd=2)
  legend('topright',leg=c('RCP2.6','RCP4.5','RCP8.5'),
          col=c(alpha('blue',0.3),alpha('green',0.3),alpha('red',0.3)),
          pch=c(15,15,15),cex=2)


  ##text(x=300,y=11,'North America',cex=1.5)
  ##text(x=0.12,y=11,'North America',cex=2)
  box(which='plot')

  plot(c(),xlim=c(-100,25),ylim=c(-2,16),cex=2,main='Eurasia',
                   xlab='Snow Cover Day Anomaly',ylab='Arctic Temperature Anomalies (\u00B0C)',
                   cex.lab=2,cex.axis=2,cex.main=2,axes=F)
  axis(1,seq(-100,40,20),seq(-100,40,20),cex.axis=2)
  axis(4,seq(0,15,5),seq(0,15,5),cex.axis=2)
  mtext("Arctic Temperature Anomalies (\u00B0C)",side=4,line=3,cex=1.5)

  ice.plot(si.85$snowtas,si.85$sea,'red')
  ice.plot(si.45$snowtas,si.45$sea,'green')
  ice.plot(si.26$snowtas,si.26$sea,'blue')

  abline(h=seq(-5,20,5),col='gray',lty=2,lwd=2)
  legend('topright',leg=c('RCP2.6','RCP4.5','RCP8.5'),
          col=c(alpha('blue',0.3),alpha('green',0.3),alpha('red',0.3)),
          pch=c(15,15,15),cex=2)

  ##text(x=300,y=11,'Eurasia',cex=1.5)
  ##text(x=0.16,y=11,'Eurasia',cex=2)
  box(which='plot')


  dev.off()



