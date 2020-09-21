##Methods to prepare data for KKZ scenario selection
##
## 1) Read in the climatologies for past and future
## 2) Regrid (conservative interpolation) to common grid (144x73 for now to match past work)
## 3) Compute absolute/relative anomalies between past and future
## 4) Clip anomalies to the subsetting region
## 5) Save files as RData 

##-------------------------------------------------------------------

library(raster)
library(ncdf4)


get_region_boundaries <- function(region) {

   regions <- list(
      AUS=list(lat=c(-45, -11), lon=c(110, 155), name='Australia'),
      AMZ=list(lat=c(-20, 12), lon=c(278, 326), name='Amazon Basin'),
      SSA=list(lat=c(-56, -20), lon=c(284, 320), name='Southern South America'),
      CAM=list(lat=c(10, 30), lon=c(244, 277), name='Central America'),
      WNA=list(lat=c(30, 60), lon=c(230, 257), name='Western North America'),
      CNA=list(lat=c(30, 50), lon=c(257, 275), name='Central North America'),
      ENA=list(lat=c(25, 50), lon=c(275, 300), name='Eastern North America'),
      ALA=list(lat=c(60, 72), lon=c(190, 257), name='Alaska'),
      GRL=list(lat=c(50, 85), lon=c(257, 350), name='Greenland'),
      MED=list(lat=c(30, 48), lon=c(350, 40), name='Mediterranean Basin'),
      NEU=list(lat=c(48, 75), lon=c(350, 40), name='Northern Europe'),
      WAF=list(lat=c(-12, 18), lon=c(340, 22), name='Western Africa'),
      EAF=list(lat=c(-12, 18), lon=c(22, 52), name='Eastern Africa'),
      SAF=list(lat=c(-35, -12), lon=c(350, 52), name='Southern Africa'),
      SAH=list(lat=c(18, 30), lon=c(340, 65), name='Sahara'),
      SEA=list(lat=c(-11, 20), lon=c(95, 155), name='Southeast Asia'),
      EAS=list(lat=c(20, 50), lon=c(100, 145), name='East Asia'),
      SAS=list(lat=c(5, 30), lon=c(65, 100), name='South Asia'),
      CAS=list(lat=c(30, 50), lon=c(40, 75), name='Central Asia'),
      TIB=list(lat=c(30, 50), lon=c(75, 100), name='Tibet'),
      NAS=list(lat=c(50, 70), lon=c(40, 180), name='North Asia'))

   region.info <- regions[[region]]
   return(region.info)

}

is_ratio <- function(var.name) {
   ratios <- list(cddETCCDI=FALSE,csdiETCCDI=FALSE,cwdETCCDI=FALSE,dtrETCCDI=FALSE,
                  fdETCCDI=FALSE,gslETCCDI=FALSE,idETCCDI=FALSE,
                  prcptotETCCDI=TRUE,
                  r10mmETCCDI=FALSE,r1mmETCCDI=FALSE,r20mmETCCDI=FALSE,
                  r95pETCCDI=TRUE,r99pETCCDI=TRUE,r95daysETCCDI=TRUE,r99daysETCCDI=TRUE,
                  r95storeETCCDI=TRUE,r99storeETCCDI=TRUE,
                  rx1dayETCCDI=TRUE,rx2dayETCCDI=TRUE,rx5dayETCCDI=TRUE,
                  sdiiETCCDI=TRUE,
                  suETCCDI=FALSE,su30ETCCDI=FALSE,
                  tn10pETCCDI=FALSE,tn90pETCCDI=FALSE,
                  tnnETCCDI=FALSE,tnxETCCDI=FALSE,
                  trETCCDI=FALSE,
                  tx10pETCCDI=FALSE,tx90pETCCDI=FALSE,
                  txnETCCDI=FALSE,txxETCCDI=FALSE,
                  wsdiETCCDI=FALSE,     
                  pr=TRUE,tasmax=FALSE,tasmin=FALSE)
   rv <- ratios[[var.name]]
   return(rv)
}

##----------------------------------------------------------------

gcm.list <- list(list(gcm='ACCESS1-0',runs='r1i1p1'),
                 list(gcm='bcc-csm1-1',runs='r1i1p1'),
                 list(gcm='bcc-csm1-1-m',runs='r1i1p1'),
                 list(gcm='BNU-ESM',runs='r1i1p1'),
                 list(gcm='CanESM2',runs=c('r1i1p1','r2i1p1','r3i1p1','r4i1p1','r5i1p1')),
                 list(gcm='CCSM4',runs=c('r1i1p1','r2i1p1','r6i1p1')),
                 list(gcm='CMCC-CM',runs='r1i1p1'),
                 list(gcm='CMCC-CMS',runs='r1i1p1'),
                 list(gcm='CNRM-CM5',runs='r1i1p1'),
                 list(gcm='CSIRO-Mk3-6-0',runs=c('r1i1p1','r2i1p1','r3i1p1','r4i1p1','r5i1p1',
                                                 'r6i1p1','r7i1p1','r8i1p1','r9i1p1','r10i1p1')),
                 list(gcm='FGOALS-g2',runs=c('r1i1p1')),
                 list(gcm='FGOALS-s2',runs=c('r1i1p1','r2i1p1','r3i1p1')),
                 list(gcm='GFDL-ESM2G',runs='r1i1p1'),
                 list(gcm='GFDL-ESM2M',runs='r1i1p1'),
                 list(gcm='GISS-E2-R',runs=c('r6i1p1','r6i1p3')),
                 list(gcm='HadGEM2-CC',runs='r1i1p1'),
                 list(gcm='HadGEM2-ES',runs=c('r2i1p1','r3i1p1','r4i1p1')),
                 list(gcm='inmcm4',runs='r1i1p1'),
                 list(gcm='IPSL-CM5A-LR',runs=c('r1i1p1','r2i1p1','r3i1p1','r4i1p1')),
                 list(gcm='IPSL-CM5B-LR',runs='r1i1p1'),
                 list(gcm='MIROC-ESM-CHEM',runs='r1i1p1'),
                 list(gcm='MIROC-ESM',runs='r1i1p1'),
                 list(gcm='MIROC5',runs=c('r1i1p1','r2i1p1','r3i1p1')),
                 list(gcm='MPI-ESM-LR',runs=c('r1i1p1','r2i1p1','r3i1p1')),
                 list(gcm='MPI-ESM-MR',runs=c('r1i1p1','r2i1p1','r3i1p1')),
                 list(gcm='MRI-CGCM3',runs='r1i1p1'),
                 list(gcm='NorESM1-M',runs='r1i1p1'))

##----------------------------------------------------------------

tmp.dir <- '/local_temp/ssobie/cmip5_prep/'
if (!file.exists(tmp.dir)) {
    dir.create(tmp.dir,recursive=T)
}

base.dir <- '/storage/data/climate/CLIMDEX/CMIP5_new/climatologies/'
save.dir <- '/storage/data/climate/CLIMDEX/CMIP5_new/KKZ/'

##-------------------
scenario <- 'rcp45'
type <- 'climdex'
past.int <- '1971-2000'
proj.int <- '2071-2100'
region <- 'WNA'

##-------------------

for (g in seq_along(gcm.list)) {
   gcm.info <- gcm.list[[g]]
   gcm <- gcm.info$gcm
   cat(paste0(gcm,'..'))
   runs <- gcm.info$runs

   for (j in seq_along(runs)) {
      run <- runs[j]

      hist.dir <- paste0(base.dir,gcm,'_historical_',run,'/')
      hist.files <- list.files(path=hist.dir,pattern='ETCCDI_climatology')
      past.files <- hist.files[grep(past.int,hist.files)]

      rcp.dir <-  paste0(base.dir,gcm,'_',scenario,'_',run,'/')
      rcp.files <- list.files(path=rcp.dir,pattern='ETCCDI_climatology')
      proj.files <- rcp.files[grep(proj.int,rcp.files)]

      stopifnot(length(past.files)==length(proj.files))

      var.names <- as.vector(sapply(past.files,function(x){strsplit(x,'_')[[1]][1]}))

      if (type=='climdex') {
         clipped.anomaly.median <- old.anomalies <- rep(0,length(past.files))
      }
      if (type=='seasonal') {
         clipped.anomaly.median <- matrix(0,nrow=length(past.files),ncol=4)
      }

      for (i in seq_along(past.files)) {
         var.name <- var.names[i]
         print(var.name)
         past.file <- past.files[i]
         proj.file <- proj.files[i]

         ##Copy to temp
         file.copy(from=paste0(hist.dir,past.file),to=tmp.dir,overwrite=TRUE)
         file.copy(from=paste0(rcp.dir,proj.file),to=tmp.dir,overwrite=TRUE)

         ##--------------------------------------------

         past.remap.file <- gsub(pattern=var.name,replacement=paste0(var.name,'_regrid_'),past.file)
         proj.remap.file <- gsub(pattern=var.name,replacement=paste0(var.name,'_regrid_'),proj.file)
         anomaly.remap.file <- gsub(pattern=proj.int,replacement=paste0(proj.int,'-',past.int),proj.remap.file)
         clipped.anomaly.file <- paste0(region,'_subset_',anomaly.remap.file)

         work <- paste0("cdo remapcon,r144x73 ",tmp.dir,past.file," ",tmp.dir,past.remap.file)
         system(work)
         Sys.sleep(1)
         work <- paste0("cdo remapcon,r144x73 ",tmp.dir,proj.file," ",tmp.dir,proj.remap.file)
         system(work)
         Sys.sleep(1)

         ##Average regionally first, then calculate the change

         region.bounds <- get_region_boundaries(region)
         map.extent <- extent(c(region.bounds$lon,region.bounds$lat))
if (1==0) {
         past.brick <- brick(paste0(tmp.dir,past.remap.file))
         clipped.past <- crop(past.brick,map.extent)
         clipped.past.mean <- cellStats(clipped.past,mean,na.rm=T)

         proj.brick <- brick(paste0(tmp.dir,proj.remap.file))
         clipped.proj <- crop(proj.brick,map.extent)
         clipped.proj.mean <- cellStats(clipped.proj,mean,na.rm=T)

         if (!is_ratio(var.name)) {
            ###clipped.anomaly.median[i] <-  clipped.proj.mean - clipped.past.mean
            old.anomalies[i] <-  clipped.proj.mean - clipped.past.mean
         } else {
            ###clipped.anomaly.median[i] <-  (clipped.proj.mean - clipped.past.mean)/clipped.past.mean*100
            old.anomalies[i] <-  (clipped.proj.mean - clipped.past.mean)/clipped.past.mean*100
         }
}
        ##Alex's method
         nc <- nc_open(paste0(tmp.dir,past.remap.file))
         proj.nc <- nc_open(paste0(tmp.dir,proj.remap.file))
         lon <- ncvar_get(nc,'lon')
         lat <- ncvar_get(nc,'lat')
         lon.x <- which((lon >= region.bounds$lon[1]) & (lon <= region.bounds$lon[2]))
         lat.y <- which((lat >= region.bounds$lat[1]) & (lat <= region.bounds$lat[2]))

         clipped.past.median <- mean(ncvar_get(nc,var.name)[lon.x,lat.y],na.rm=T)                  
         clipped.proj.median <- mean(ncvar_get(proj.nc,var.name)[lon.x,lat.y],na.rm=T)                  

         if (!is_ratio(var.name)) {
            clipped.anomaly.median[i] <-  clipped.proj.median - clipped.past.median
         } else {
            clipped.anomaly.median[i] <-  (clipped.proj.median - clipped.past.median)/clipped.past.median*100
         }
         nc_close(nc)
         nc_close(proj.nc)

if (1==0) {
         anomaly.brick <- brick(paste0(tmp.dir,anomaly.remap.file))
         clipped.anomaly <- crop(anomaly.brick,map.extent)


         ##Compute Anomalies
         if (is_ratio(var.name)) {
            offset <- 100
            work <- paste0("cdo sub ",tmp.dir,proj.remap.file," ",tmp.dir,past.remap.file," ",tmp.dir,"tmp_diff.nc")
            system(work)
            Sys.sleep(1)
            work <- paste0("cdo div ",tmp.dir,"tmp_diff.nc ",tmp.dir,past.remap.file," ",tmp.dir,anomaly.remap.file)
            system(work)
            Sys.sleep(1)
            file.remove(paste0(tmp.dir,"tmp_diff.nc"))
         } else{
            offset <- 1
            work <- paste0("cdo sub ",tmp.dir,proj.remap.file," ",tmp.dir,past.remap.file," ",tmp.dir,anomaly.remap.file)
            system(work)
            Sys.sleep(1)
         }

    
         if (type=='climdex') {
            clipped.anomaly.median[i] <- cellStats(clipped.anomaly,median,na.rm=T)*offset
         }
         if (type=='seasonal') {
            clipped.anomaly.median[i,] <- cellStats(clipped.anomaly,median,na.rm=T)*offset
         }
     }
         ##file.remove(paste0(tmp.dir,anomaly.remap.file))
         file.remove(paste0(tmp.dir,proj.remap.file))
         file.remove(paste0(tmp.dir,past.remap.file))
         file.remove(paste0(tmp.dir,proj.file))    
         file.remove(paste0(tmp.dir,past.file))

     }

     if (type=='climdex') {
        names(clipped.anomaly.median) <- var.names
     }

     if (type=='seasonal') {
        rownames(clipped.anomaly.median) <- var.names
     }

     save.file <- paste0(region,'_nc_mean_first_',gcm,'_',scenario,'_',run,'_',type,'_',proj.int,'_',past.int,'.RData')
     save(clipped.anomaly.median,file=paste0(save.dir,save.file))
 
   }
}