##Methods to prepare data for KKZ scenario selection
##
## 1) Read in the climatologies for past and future
## 2) Regrid (conservative interpolation) to common grid (144x73 for now to match past work)
## 3) Compute absolute/relative anomalies between past and future
## 4) Clip anomalies to the subsetting region
## 5) Save files as RData 

##Revised to use the IPCC-WGI Version 4 regions

##-------------------------------------------------------------------

library(raster)
library(rgdal)

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

##                 list(gcm='GFDL-CM4',runs='r1i1p1f1'),

gcm.list <- list(list(gcm='ACCESS-CM2',runs='r1i1p1f1'),
                 list(gcm='ACCESS-ESM1-5',runs=c('r1i1p1f1','r2i1p1f1','r3i1p1f1')),
                 list(gcm='BCC-CSM2-MR',runs='r1i1p1f1'),
                 list(gcm='CanESM5',runs=c('r1i1p2f1','r2i1p2f1','r3i1p2f1',
                             'r4i1p2f1','r5i1p2f1','r6i1p2f1',
                             'r7i1p2f1','r8i1p2f1','r9i1p2f1',
                             'r10i1p2f1')),
                 list(gcm='CNRM-CM6-1',runs='r1i1p1f2'),
                 list(gcm='CNRM-ESM2-1',runs='r1i1p1f2'),
                 list(gcm='EC-Earth3',runs='r4i1p1f1'),
                 list(gcm='EC-Earth3-Veg',runs='r1i1p1f1'),
                 list(gcm='FGOALS-g3',runs='r1i1p1f1'),
                 list(gcm='GFDL-ESM4',runs='r1i1p1f1'),
                 list(gcm='HadGEM3-GC31-LL',runs='r1i1p1f3'),
                 list(gcm='INM-CM4-8',runs='r1i1p1f1'),
                 list(gcm='INM-CM5-0',runs='r1i1p1f1'),
                 list(gcm='IPSL-CM6A-LR',runs='r1i1p1f1'),
                 list(gcm='KACE-1-0-G',runs='r2i1p1f1'),
                 list(gcm='KIOST-ESM',runs='r1i1p1f1'),
                 list(gcm='MIROC6',runs='r1i1p1f1'),
                 list(gcm='MIROC-ES2L',runs='r1i1p1f2'),
                 list(gcm='MPI-ESM1-2-HR',runs='r1i1p1f1'),
                 list(gcm='MPI-ESM1-2-LR',runs='r1i1p1f1'),
                 list(gcm='MRI-ESM2-0',runs='r1i1p1f1'),
                 list(gcm='NESM3',runs='r1i1p1f1'),
                 list(gcm='NorESM2-LM',runs='r1i1p1f1'),
                 list(gcm='NorESM2-MM',runs='r1i1p1f1'),
                 list(gcm='UKESM1-0-LL',runs='r1i1p1f2'))

###Add GCFL-CM4 separately to include the resolution difference
gcm.list <- list(list(gcm='CNRM-CM6-1',runs='r1i1p1f2'),
                 list(gcm='CNRM-ESM2-1',runs='r1i1p1f2'),
                 list(gcm='EC-Earth3',runs='r4i1p1f1'),
                 list(gcm='EC-Earth3-Veg',runs='r1i1p1f1'),
                 list(gcm='FGOALS-g3',runs='r1i1p1f1'),
                 list(gcm='GFDL-ESM4',runs='r1i1p1f1'),
                 list(gcm='HadGEM3-GC31-LL',runs='r1i1p1f3'),
                 list(gcm='INM-CM4-8',runs='r1i1p1f1'),
                 list(gcm='INM-CM5-0',runs='r1i1p1f1'),
                 list(gcm='IPSL-CM6A-LR',runs='r1i1p1f1'),
                 list(gcm='KACE-1-0-G',runs='r2i1p1f1'),
                 list(gcm='KIOST-ESM',runs='r1i1p1f1'),
                 list(gcm='MIROC6',runs='r1i1p1f1'),
                 list(gcm='MIROC-ES2L',runs='r1i1p1f2'),
                 list(gcm='MPI-ESM1-2-HR',runs='r1i1p1f1'),
                 list(gcm='MPI-ESM1-2-LR',runs='r1i1p1f1'),
                 list(gcm='MRI-ESM2-0',runs='r1i1p1f1'),
                 list(gcm='NESM3',runs='r1i1p1f1'),
                 list(gcm='NorESM2-LM',runs='r1i1p1f1'),
                 list(gcm='NorESM2-MM',runs='r1i1p1f1'),
                 list(gcm='UKESM1-0-LL',runs='r1i1p1f2'))


gcm.list <- list(list(gcm='CanESM5',runs=c('r7i1p2f1','r8i1p2f1','r9i1p2f1',
                          'r10i1p2f1')))


##--------------------------------------------------------------------

testing <- FALSE

if (testing) {
   tmpdir <- '/local_temp/ssobie/cmip6_prep'
} else {
   args <- commandArgs(trailingOnly=TRUE)
   for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
   }
}

tmp.dir <- paste0(tmpdir,'/kkz_prep/') ##
if (!file.exists(tmp.dir)) {
    dir.create(tmp.dir,recursive=T)
}

base.dir <- '/storage/data/climate/CMIP6/Derived/'

##-------------------
##gcm <- 'MRI-ESM2-0'
scenario <- 'ssp245'
##run <- 'r1i1p1f1'
res <- NULL
type <- 'climdex'
past.int <- '1986-2005'
proj.int <- '2081-2100'
region <- 'CAN'
res <- NULL

##-------------------
##Load the IPCC-WGI shapfiles
##shape.dir <- '/storage/data/projects/rci/data/assessments/shapefiles/ipcc_regions_v4' 
##ipcc.v4 <- readOGR(shape.dir,'IPCC-WGI-reference-regions-v4')

shape.dir <- '/storage/data/projects/rci/data/assessments/shapefiles/bc_common/' 
ipcc.v4 <- readOGR(shape.dir,"canada_boundary")
##-------------------

for (g in seq_along(gcm.list)) {
   gcm.info <- gcm.list[[g]]
   gcm <- gcm.info$gcm
   cat(paste0(gcm,'..'))
   runs <- gcm.info$runs

   for (j in seq_along(runs)) {
      run <- runs[j]
      cat(paste0(run,'..'))

      if (!is.null(res)) {
         read.dir <- paste0(base.dir,gcm,'_',scenario,'_',run,'_',res,'/',type,'/climatologies/')
      } else {
         read.dir <- paste0(base.dir,gcm,'_',scenario,'_',run,'/',type,'/climatologies/')
      }

      save.dir <- '/storage/data/climate/CMIP6/KKZ/'

      files <- list.files(path=read.dir,pattern='_annual_')
      past.files <- files[grep(past.int,files)]
      if (length(past.files) != 33) {
         print(paste0('Length of files ',length(files)))
         print(past.files)
         stop('Not the right number of files')
      }
      proj.files <- files[grep(proj.int,files)]

      stopifnot(length(past.files)==length(proj.files))

      var.names <- as.vector(sapply(past.files,function(x){strsplit(x,'_')[[1]][1]}))

      if (type=='climdex') {
         clipped.anomaly.median <- rep(0,length(past.files))
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
         file.copy(from=paste0(read.dir,past.file),to=tmp.dir,overwrite=TRUE)
         file.copy(from=paste0(read.dir,proj.file),to=tmp.dir,overwrite=TRUE)

         ##--------------------------------------------

         past.remap.file <- gsub(pattern=var.name,replacement=paste0(var.name,'_regrid_'),past.file)
         proj.remap.file <- gsub(pattern=var.name,replacement=paste0(var.name,'_regrid_'),proj.file)
         anomaly.remap.file <- gsub(pattern=proj.int,replacement=paste0(proj.int,'-',past.int),proj.remap.file)
         clipped.anomaly.file <- paste0(region,'_subset_',anomaly.remap.file)

         ###work <- paste0("cdo remapcon,r144x73 ",tmp.dir,past.file," ",tmp.dir,past.remap.file)
         ##Revise to 1 degree grid following the latest CMIP6 comparison of Climdex
         work <- paste0("cdo remapcon,r360x181 ",tmp.dir,past.file," ",tmp.dir,past.remap.file)
         system(work)
         Sys.sleep(1)
         ###work <- paste0("cdo remapcon,r144x73 ",tmp.dir,proj.file," ",tmp.dir,proj.remap.file)
         work <- paste0("cdo remapcon,r360x181 ",tmp.dir,proj.file," ",tmp.dir,proj.remap.file)   
         system(work)
         Sys.sleep(1)

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

         ###reg.ix <- which(ipcc.v4$Acronym==region)
         region.shape <- ipcc.v4 ###[reg.ix,]        
         anomaly.brick <- brick(paste0(tmp.dir,anomaly.remap.file))
         print(extent(anomaly.brick))
         if (xmax(anomaly.brick) > 200) {
            x.offset <- -180 - xmin(anomaly.brick)
            xmin(anomaly.brick) <- xmin(anomaly.brick) + x.offset
            xmax(anomaly.brick) <- xmax(anomaly.brick) + x.offset
         }
         clipped.anomaly <- mask(anomaly.brick,region.shape)
    
         if (type=='climdex') {
            clipped.anomaly.median[i] <- cellStats(clipped.anomaly,median,na.rm=T)*offset
         }
         if (type=='seasonal') {
            clipped.anomaly.median[i,] <- cellStats(clipped.anomaly,median,na.rm=T)*offset
         }

         file.remove(paste0(tmp.dir,anomaly.remap.file))
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

      save.file <- paste0(region,'_',gcm,'_',scenario,'_',run,'_',type,'_',proj.int,'_',past.int,'.RData')
      save(clipped.anomaly.median,file=paste0(save.dir,save.file))
   }
}