##Script to calculate climatologies from the derived variable files
##for CMIP6

##---------------------------------------------------------
##GCMs based on the available CMIP6/Derived set of GCMs

##                 list(gcm='GFDL-CM4',runs='r1i1p1f1'),

##list(gcm='ACCESS-ESM1-5',runs=c('r1i1p1f1','r2i1p1f1','r3i1p1f1')),
##                 list(gcm='BCC-CSM2-MR',runs='r1i1p1f1'),
##                 list(gcm='CanESM5',runs=c('r1i1p2f1','r2i1p2f1','r3i1p2f1','r4i1p2f1',
##                                           'r5i1p2f1','r6i1p2f1','r7i1p2f1','r8i1p2f1'),                 
##                                           'r9i1p2f1','r10i1p2f1'),

gcm.list <- list(list(gcm='ACCESS-CM2',runs='r1i1p1f1'),
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
                 list(gcm='MIROC6',runs='r1i1p1f1'),
                 list(gcm='MPI-ESM1-2-HR',runs='r1i1p1f1'),
                 list(gcm='MPI-ESM1-2-LR',runs='r1i1p1f1'),
                 list(gcm='MRI-ESM2-0',runs='r1i1p1f1'),
                 list(gcm='NESM3',runs='r1i1p1f1'),
                 list(gcm='NorESM2-LM',runs='r1i1p1f1'),
                 list(gcm='UKESM1-0-LL',runs='r1i1p1f2'))

##gcm.list <- list(list(gcm='CanESM5',runs=c('r1i1p2f1','r2i1p2f1','r3i1p2f1','r4i1p2f1',
##                                           'r5i1p2f1','r6i1p2f1','r7i1p2f1','r8i1p2f1'),                 
##                                           'r9i1p2f1','r10i1p2f1'))          

gcm.list <- list(list(gcm='KACE-1-0-G',runs=c('r2i1p1f1')))

##---------------------------------------------------------
##Averaging method for non-annual parameters

climdex.names <- c("cddETCCDI","csdiETCCDI","cwdETCCDI","dtrETCCDI","fdETCCDI",
                   "gslETCCDI","idETCCDI","prcptotETCCDI","r1mmETCCDI","r10mmETCCDI",   
                   "r20mmETCCDI","r95daysETCCDI","r95pETCCDI","r95storeETCCDI",
                   "r99daysETCCDI","r99pETCCDI","r99storeETCCDI","rx1dayETCCDI",
                   "rx2dayETCCDI","rx5dayETCCDI","sdiiETCCDI","su30ETCCDI",
                   "suETCCDI","tn10pETCCDI","tn90pETCCDI","tnnETCCDI",
                   "tnxETCCDI","trETCCDI","tx10pETCCDI","tx90pETCCDI",  
                   "txnETCCDI","txxETCCDI","wsdiETCCDI")

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

testing <- TRUE

if (!testing) {
   args <- commandArgs(trailingOnly=TRUE)
   for(i in 1:length(args)){
       eval(parse(text=args[[i]])) 
   }
} else {
   tmpdir <- '/local_temp/ssobie/cmip6_clims/' ##
}

tmp.dir <-  paste0(tmpdir,'/') ##
if (!file.exists(tmp.dir)) {
    dir.create(tmp.dir,recursive=T)
}

base.dir <- '/storage/data/climate/CMIP6/Derived/'

##-------------------

scenario <- 'ssp245'
type <- 'climdex'
intervals <- c('1971-2000','1981-2010','1986-2005','2041-2070','2071-2100','2081-2100')
res <- 'NULL'

##-------------------

for (g in seq_along(gcm.list)) {
   gcm.info <- gcm.list[[g]]
   gcm <- gcm.info$gcm
   cat(paste0(gcm,'..'))
   runs <- gcm.info$runs
   
   if (gcm=='GFDL-CM4') {
      ##stop('Deal with the different resolutions')
   }

   for (j in seq_along(runs)) {
      run <- runs[j]
      cat(paste0(run,'..'))
      read.dir <- paste0(base.dir,gcm,'_',scenario,'_',run,'/',type,'/')
      if (gcm=='GFDL-CM4') {
         read.dir <- paste0(base.dir,gcm,'_',scenario,'_',run,'_',res,'/',type,'/')
      }
      write.dir <- paste0(read.dir,'climatologies/')

      if (!file.exists(write.dir)) {
         dir.create(write.dir,recursive=T)
      }

      files <- list.files(path=read.dir,pattern=gcm)
###      files <- files[grep('r1mmETCCDI',files)]
      print(files)
      if (gcm=='GFDL-CM4') {
         files <- files[grep(res,files)]
         ##stop('Deal with the different resolutions')
      }

      var.names <- as.vector(sapply(files,function(x){strsplit(x,'_')[[1]][1]})) ##'dtrETCCDI' ##
      var.freqs <- as.vector(sapply(files,function(x){strsplit(x,'_')[[1]][2]})) ##'ann' ##

      for (i in seq_along(files)) {
         var.name <- var.names[i]
         var.freq <- var.freqs[i]
         cat(paste0(var.name,'..'))
         print(var.name)
         var.file <- files[i]

         ##Copy to temp
         file.copy(from=paste0(read.dir,var.file),to=tmp.dir,overwrite=TRUE)
 
         ##--------------------------------------------
         ##Climdex Variables
         if (type=='climdex') {
            if (var.freq=='mon') {
               agg.file <- aggregate_to_annual(var.file,var.name,var.freq,
                             tmp.dir,tmp.dir)
               var.file <- agg.file
               var.freq <- 'ann'                              

            }

            ##Subset by time
            for (interval in intervals) {

               time.file <- subset_by_time(var.file,interval,var.freq,
                                           tmp.dir,tmp.dir)
               clim.file <- gsub(pattern=paste0('_',var.freq,'_'),replacement='_climatology_',time.file)
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
