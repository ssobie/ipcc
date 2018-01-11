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

get.area.average <- function(file.name,var.name,bnds) {

   nc <- nc_open(file.name)
   lon <- ncvar_get(nc,'lon')
   lon <- ((lon + 180) %% 360) - 180
   lat <- ncvar_get(nc,'lat')
   time <- netcdf.calendar(nc)

   obs.st <- head(grep('1861',time),1)
   if (length(obs.st) ==0) {
     obs.st <- 1
   }
   obs.en <- length(time) ##tail(grep('2100',time),1)
   print(obs.st)
   print(obs.en)
   clim.st <- grep('1861',time[obs.st:obs.en])
   clim.en <- grep('1880',time[obs.st:obs.en])

   lon.st <- which.min(abs(bnds$lon[1]-lon))
   lon.en <- which.min(abs(bnds$lon[2]-lon))
   lat.st <- which.min(abs(bnds$lat[1]-lat))
   lat.en <- which.min(abs(bnds$lat[2]-lat))

   data.series <- ncvar_get(nc,var.name)
   data.subset <- data.series[,,obs.st:obs.en]

   print(dim(data.subset))

   data.clim <- apply(data.subset[,,clim.st:clim.en],c(1,2),mean,na.rm=T)
   data.avg <- aperm(apply(data.clim,c(1,2),function(x){rep(x,each=dim(data.subset)[3])}),c(2,3,1))      
   data.anom <- data.subset - data.avg
   anom.sub <- data.anom[,lat.st:lat.en,]

   ar <- get.areas(lon,lat[lat.st:lat.en],dim(anom.sub))
   arctic.wts <- ar[,1]/sum(ar[,1])
   
   arctic.mid <- apply(anom.sub,c(1,3),weighted.mean,arctic.wts,na.rm=T)
   arctic.mean <- apply(apply(anom.sub,c(1,3),weighted.mean,arctic.wts,na.rm=T),2,mean,na.rm=T)

   nc_close(nc)
   
   rv <- arctic.mean 
   return(rv)
}

##*************************************************************************

var.name <- 'tas'
scenario <- 'rcp85'

##Arctic Area
bnds <- list(lon=c(-180,180),
             lat=c(60,90))   


##Snow cover models
##Snow Cover Models
gcm.list <- c('bcc-csm1-1-m+r1i1p1',
              'CanESM2+r1i1p1','CanESM2+r2i1p1','CanESM2+r3i1p1','CanESM2+r4i1p1','CanESM2+r5i1p1',
              'CCSM4+r1i1p1','CCSM4+r2i1p1','CCSM4+r3i1p1','CCSM4+r4i1p1','CCSM4+r5i1p1',
              'CNRM-CM5+r1i1p1',
              'CSIRO-Mk3-6-0+r1i1p1','CSIRO-Mk3-6-0+r2i1p1','CSIRO-Mk3-6-0+r3i1p1',
              'CSIRO-Mk3-6-0+r4i1p1','CSIRO-Mk3-6-0+r5i1p1',
              'GISS-E2-H+r1i1p1','GISS-E2-H+r2i1p1','GISS-E2-H+r3i1p1',
              'GISS-E2-R+r1i1p1',
              'MIROC5+r1i1p1','MIROC5+r2i1p1','MIROC5+r3i1p1',
              'MIROC-ESM+r1i1p1','MIROC-ESM-CHEM+r1i1p1',
              'MPI-ESM-LR+r1i1p1','MPI-ESM-LR+r2i1p1','MPI-ESM-LR+r3i1p1',
              'MPI-ESM-MR+r1i1p1',
              'MRI-CGCM3+r1i1p1',
              'NorESM1-M+r1i1p1')

##GCMs
##Sea Ice Models

gcm.list <- c('bcc-csm1-1','bcc-csm1-1-m','CanESM2','CCSM4','CNRM-CM5',
              'GISS-E2-H','GISS-E2-R','HadGEM2-AO','HadGEM2-ES','IPSL-CM5A-LR','IPSL-CM5A-MR',
              'MIROC5','MIROC-ESM','MIROC-ESM-CHEM','MPI-ESM-LR','MRI-CGCM3','NorESM1-M')

gcm.list <- c('ACCESS1-0','bcc-csm1-1','bcc-csm1-1-m','CanESM2','CCSM4','CNRM-CM5',
              'GISS-E2-H','GISS-E2-R','HadGEM2-AO','HadGEM2-CC','HadGEM2-ES','inmcm4','IPSL-CM5A-LR','IPSL-CM5A-MR',
              'MIROC5','MIROC-ESM','MIROC-ESM-CHEM','MPI-ESM-LR','MRI-CGCM3','NorESM1-M')

##gcm.list <- c('bcc-csm1-1','bcc-csm1-1-m','CanESM2','CCSM4','CNRM-CM5',
##              'CSIRO-Mk3-6-0','MIROC5','MIROC-ESM','MIROC-ESM-CHEM','MPI-ESM-LR','MRI-CGCM3','NorESM1-M')



##series.data <- vector(length=length(gcm.list),mode='list')
series.data <- matrix(NA,nrow=240,ncol=length(gcm.list))
gcms <- c()
runs <- c()
for (g in seq_along(gcm.list)) {
  gcm.name <- strsplit(gcm.list[g],'\\+')[[1]]
  gcm <- gcm.name[1]
  gcms[g] <- gcm
  run <- 'r1i1p1' ##gcm.name[2]
  runs[g] <- substr(run,2,2)

##  gcm <- gcm.list[g]    
##  gcms[g] <- gcm
##  runs[g] <- 1 


  gcm.dir <- paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/CMIP5/global/annual/',gcm,'/')
  
  all.files <- list.files(path=gcm.dir,pattern=paste0('tas_ann_',gcm,'_*'),full.name=T)
  past.files <- all.files[grep('18500101',all.files)]
  run.files <- past.files[grep('r1i1p1',past.files)]
  scen.file <- run.files[grep(scenario,run.files)]

  print(scen.file)

  ag.past <- get.area.average(scen.file,var.name,bnds)
  if (length(ag.past)==241) {
     ag.past <- ag.past[2:241]
  }
  if (length(ag.past)==231) {
     ag.past <- c(rep(NA,10),ag.past[12:241])
  }
  if (length(ag.past)==232) {
     ag.past <- c(rep(NA,10),ag.past[12:241])
  }
  if (length(ag.past)==239) {
     ag.past <- c(ag.past,NA)
  }

  series.data[,g] <- ag.past 
}

header1 <- c('Model',gcms)
header2 <- c('Realization',runs)
years <- 1861:2100
full.data <- rbind(header1,header2,cbind(years,round(series.data,1)))

write.file <- paste0('/storage/data/projects/rci/data/nrcan/ice_snow_data/tas_sic_',scenario,'_arctic_1860.csv')
write.table(full.data,file=write.file,quote=F,row.name=F,col.name=F,sep=',')


