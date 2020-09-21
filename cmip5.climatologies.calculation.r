##Script to calculate climatologies from the derived variable files
##for CMIP6

##---------------------------------------------------------
##GCMs to match Alex's scenario selection paper

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


##---------------------------------------------------------
##Averaging method for non-annual parameters

climdex.names <- c("cddETCCDI","cwdETCCDI","dtrETCCDI","fdETCCDI",
                   "gslETCCDI","idETCCDI","prcptotETCCDI","r10mmETCCDI",   
                   "r20mmETCCDI","r95daysETCCDI","r95pETCCDI","r95storeETCCDI",
                   "r99daysETCCDI","r99pETCCDI","r99storeETCCDI","rx1dayETCCDI",
                   "rx2dayETCCDI","rx5dayETCCDI","sdiiETCCDI","su30ETCCDI",
                   "suETCCDI","tn10pETCCDI","tn90pETCCDI","tnnETCCDI",
                   "tnxETCCDI","trETCCDI","tx10pETCCDI","tx90pETCCDI",  
                   "txnETCCDI","txxETCCDI")

climdex.freqs <- c('ann','mon')

get_agg_fxn <- function(var.name) {
   fx <- switch(var.name,
                pr='mean',
                tasmax='mean',
                tasmin='mean',
                rx1dayETCCDI='max',
                rx2dayETCCDI='max',
                rx5dayETCCDI='max',
                txxETCCDI='max',
                tnxETCCDI='max',
                txnETCCDI='min',
                tnnETCCDI='min',
                tn10pETCCDI='mean',
                tx10pETCCDI='mean',
                tn90pETCCDI='mean',
                tx90pETCCDI='mean',
                dtrETCCDI='mean')
   if (is.null(fx)) {
      fx <- mean
   }
   return(fx)
}

##---------------------------------------------------------
##Aggregate non-annual parameters to annual

aggregate_to_annual <- function(input.file,var.name,freq,
                                read.dir,write.dir) {

   agg.fxn <- get_agg_fxn(var.name)
   
   agg.file <- gsub(pattern=paste0('_',freq,'_'),replacement='_ann_',x=input.file)
   work <- paste0('cdo year',agg.fxn,' ',read.dir,input.file,' ',
                  write.dir,agg.file)
   system(work)
   Sys.sleep(1)                     
   return(agg.file)
}  

##---------------------------------------------------------
##Subset the derived file by the given time interval

subset_by_time <- function(input.file,interval,freq,read.dir,write.dir) {
   yrs <- strsplit(interval,'-')[[1]]
   subset.file <- gsub(pattern='[0-9]{4}-[0-9]{4}',replacement=interval,input.file)
   work <- paste0('cdo seldate,',yrs[1],'-01-01T00:00:00,',yrs[2],'-12-31T23:59:59 ',
                  read.dir,input.file,' ',write.dir,subset.file)
   system(work)
   Sys.sleep(1)                     
   return(subset.file)
}

##---------------------------------------------------------
##Subset the derived file by the given time interval

make_climatology <- function(input.file,clim.file,clim.fx,read.dir,write.dir) {

   work <- paste0('cdo ',clim.fx,' ',
                  read.dir,input.file,' ',write.dir,clim.file)
   system(work)
   Sys.sleep(1)                     
   ##return(clim.file)
}

##---------------------------------------------------------

##*********************************************************

##args <- commandArgs(trailingOnly=TRUE)
##for(i in 1:length(args)){
##    eval(parse(text=args[[i]]))
##}


tmp.dir <- '/local_temp/ssobie/cmip5_clims/' ##paste0(tmpdir,'/') ##
if (!file.exists(tmp.dir)) {
    dir.create(tmp.dir,recursive=T)
}

base.dir <- '/storage/data/climate/CLIMDEX/CMIP5_new/'

##-------------------

scenario <- 'rcp45'
type <- 'climdex'
intervals <- '2071-2100' ##c('1986-2005','2081-2100')


##-------------------

for (g in seq_along(gcm.list)) {
   gcm.info <- gcm.list[[g]]
   gcm <- gcm.info$gcm
   cat(paste0(gcm,'..'))
   runs <- gcm.info$runs

   for (j in seq_along(runs)) {
      run <- runs[j]
      cat(paste0(run,'..'))
      version <- list.files(path=paste0(base.dir,scenario,'/',gcm,'/',run,'/'))
      last.dir <- list.files(path=paste0(base.dir,scenario,'/',gcm,'/',run,'/',version,'/'))
      read.dir <- paste0(base.dir,scenario,'/',gcm,'/',run,'/',version,'/',last.dir,'/')
      write.dir <- paste0(base.dir,'climatologies/',gcm,'_',scenario,'_',run,'/')

      if (!file.exists(write.dir)) {
         dir.create(write.dir,recursive=T)
      }

      files <- list.files(path=read.dir,pattern='ETCCDI_yr')
      files <- files[-grep('alt',files)]

      var.names <- as.vector(sapply(files,function(x){strsplit(x,'_')[[1]][1]}))

      for (i in seq_along(files)) {
         var.name <- var.names[i]
         cat(paste0(var.name,'..'))
         print(var.name)
         var.file <- files[i]
         
         ##Copy to temp
         file.copy(from=paste0(read.dir,var.file),to=tmp.dir,overwrite=TRUE)
 
         ##--------------------------------------------
         ##Climdex Variables
         if (type=='climdex') {
            ##Subset by time
            for (interval in intervals) {
               time.file <- subset_by_time(var.file,interval,var.freq,
                                           tmp.dir,tmp.dir)
               clim.file <- gsub(pattern='_yr_',replacement='_climatology_',time.file)
               make_climatology(time.file,clim.file,'timmean',tmp.dir,tmp.dir)
               file.copy(from=paste0(tmp.dir,clim.file),to=write.dir,overwrite=TRUE)
               Sys.sleep(1)
               file.remove(paste0(tmp.dir,clim.file))
               file.remove(paste0(tmp.dir,time.file))
            } 
            file.remove(paste0(tmp.dir,var.file))
         }
         ##--------------------------------------------
         ##Annual
         if (type=='annual') {
            ##Subset by time
            for (interval in intervals) {
               time.file <- subset_by_time(var.file,interval,var.freq,
                                           tmp.dir,tmp.dir)
               clim.file <- gsub(pattern=paste0(var.freq,'_',var.avg),replacement=paste0(var.freq,'_',var.avg,'_climatology'),time.file)
               make_climatology(time.file,clim.file,'timmean',tmp.dir,tmp.dir)
               file.copy(from=paste0(tmp.dir,clim.file),to=write.dir,overwrite=TRUE)
               Sys.sleep(1)
               file.remove(paste0(tmp.dir,clim.file))
               file.remove(paste0(tmp.dir,time.file))
            }
         }
         ##--------------------------------------------
         ##Seasonal
         if (type=='seasonal') {
            ##Subset by time
            for (interval in intervals) {
               time.file <- subset_by_time(var.file,interval,var.freq,
                                           tmp.dir,tmp.dir)
               clim.file <- gsub(pattern=paste0(var.freq,'_',var.avg),replacement=paste0(var.freq,'_',var.avg,'_climatology'),time.file)
               make_climatology(time.file,clim.file,'yseasmean',tmp.dir,tmp.dir)
               file.copy(from=paste0(tmp.dir,clim.file),to=write.dir,overwrite=TRUE)
               Sys.sleep(1)
               file.remove(paste0(tmp.dir,clim.file))
               file.remove(paste0(tmp.dir,time.file))
            }
         }
         ##--------------------------------------------
         ##Monthly
         if (type=='monthly') {
            ##Subset by time
            for (interval in intervals) {
               time.file <- subset_by_time(var.file,interval,var.freq,
                                          tmp.dir,tmp.dir)
               clim.file <- gsub(pattern=paste0(var.freq,'_',var.avg),replacement=paste0(var.freq,'_',var.avg,'_climatology'),time.file)
               make_climatology(time.file,clim.file,'ymonmean',tmp.dir,tmp.dir)
               file.copy(from=paste0(tmp.dir,clim.file),to=write.dir,overwrite=TRUE)
               Sys.sleep(1)
               file.remove(paste0(tmp.dir,clim.file))
               file.remove(paste0(tmp.dir,time.file))
            }
         }

      }
   cat('\n')
   }
}
