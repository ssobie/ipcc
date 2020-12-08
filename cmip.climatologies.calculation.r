##Script to calculate climatologies from the derived variable files
##for CMIP6

library(ncdf4)


##---------------------------------------------------------
##GCMs based on the available CMIP6/Derived set of GCMs

cmip6.list <- list(list(gcm='ACCESS-CM2',runs='r1i1p1f1'),
                 list(gcm='ACCESS-ESM1-5',runs=c('r1i1p1f1','r2i1p1f1','r3i1p1f1')),     
                 list(gcm='BCC-CSM2-MR',runs='r1i1p1f1'),
                 list(gcm='CNRM-CM6-1',runs='r1i1p1f2'),
                 list(gcm='CNRM-ESM2-1',runs='r1i1p1f2'),
                 list(gcm='EC-Earth3',runs='r4i1p1f1'),
                 list(gcm='EC-Earth3-Veg',runs='r1i1p1f1'),
                 list(gcm='FGOALS-g3',runs='r1i1p1f1'),
                 list(gcm='GFDL-ESM4',runs='r1i1p1f1'),
                 list(gcm='HadGEM3-GC31-LL',runs='r1i1p1f3'),
                 list(gcm='KACE-1-0-G',runs=c('r2i1p1f1')),
                 list(gcm='KIOST-ESM',runs=c('r1i1p1f1')),
                 list(gcm='INM-CM4-8',runs='r1i1p1f1'),
                 list(gcm='INM-CM5-0',runs='r1i1p1f1'),
                 list(gcm='IPSL-CM6A-LR',runs='r1i1p1f1'),
                 list(gcm='MIROC6',runs='r1i1p1f1'),
                 list(gcm='MIROC-ES2L',runs='r1i1p1f2'),
                 list(gcm='MPI-ESM1-2-HR',runs='r1i1p1f1'),
                 list(gcm='MPI-ESM1-2-LR',runs='r1i1p1f1'),
                 list(gcm='MRI-ESM2-0',runs='r1i1p1f1'),
                 list(gcm='NESM3',runs='r1i1p1f1'),
                 list(gcm='NorESM2-LM',runs='r1i1p1f1'),
                 list(gcm='NorESM2-MM',runs='r1i1p1f1'),                 
                 list(gcm='UKESM1-0-LL',runs='r1i1p1f2'),
                 list(gcm='CanESM5',runs=c('r1i1p2f1','r2i1p2f1','r3i1p2f1','r4i1p2f1',
                                           'r5i1p2f1','r6i1p2f1','r7i1p2f1','r8i1p2f1'),                 
                                           'r9i1p2f1','r10i1p2f1'))          


cmip5.list <- list(list(gcm='bcc-csm1-1',runs='r1i1p1'),
                   list(gcm='bcc-csm1-1-m',runs='r1i1p1'),
                   list(gcm='BNU-ESM',runs='r1i1p1'),
                   list(gcm='CanESM2',runs='r1i1p1'),
                   list(gcm='CCSM4',runs='r2i1p1'),
                   list(gcm='CESM1-CAM5',runs='r1i1p1'),
                   list(gcm='CNRM-CM5',runs='r1i1p1'),
                   list(gcm='CSIRO-Mk3-6-0',runs='r1i1p1'),
                   list(gcm='FGOALS-g2',runs='r1i1p1'),
                   list(gcm='GFDL-CM3',runs='r1i1p1'),
                   list(gcm='GFDL-ESM2G',runs='r1i1p1'),
                   list(gcm='GFDL-ESM2M',runs='r1i1p1'),
                   list(gcm='HadGEM2-AO',runs='r1i1p1'),
                   list(gcm='HadGEM2-ES',runs='r1i1p1'),
                   list(gcm='IPSL-CM5A-LR',runs='r1i1p1'),
                   list(gcm='IPSL-CM5A-MR',runs='r1i1p1'),
                   list(gcm='MIROC-ESM',runs='r1i1p1'),
                   list(gcm='MIROC-ESM-CHEM',runs='r1i1p1'),
                   list(gcm='MIROC5',runs='r3i1p1'),
                   list(gcm='MPI-ESM-LR',runs='r3i1p1'),
                   list(gcm='MPI-ESM-MR',runs='r1i1p1'),
                   list(gcm='MRI-CGCM3',runs='r1i1p1'),
                   list(gcm='NorESM1-M',runs='r1i1p1'),
                   list(gcm='NorESM1-ME',runs='r1i1p1'),
                   list(gcm='ACCESS1-0',runs='r1i1p1'), 
                   list(gcm='HadGEM2-CC',runs='r1i1p1'),
                   list(gcm='inmcm4',runs='r1i1p1'))



##---------------------------------------------------------
##Averaging method for non-annual parameters

##Not used, files are taken from list.files
climdex.names <- c("cddETCCDI","csdiETCCDI","cwdETCCDI","dtrETCCDI","fdETCCDI",
                   "gslETCCDI","idETCCDI","prcptotETCCDI","r1mmETCCDI","r10mmETCCDI",   
                   "r20mmETCCDI","r95daysETCCDI","r95pETCCDI",
                   "r99daysETCCDI","r99pETCCDI","rx1dayETCCDI",
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
   work <- paste0('cdo -s year',agg.fxn,' ',read.dir,input.file,' ',
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
   work <- paste0('cdo -s seldate,',yrs[1],'-01-01T00:00:00,',yrs[2],'-12-31T23:59:59 ',
                  read.dir,input.file,' ',write.dir,subset.file)
   system(work)
   Sys.sleep(1)                     
   return(subset.file)
}

##---------------------------------------------------------
##Subset the derived file by the given time interval

make_climatology <- function(input.file,clim.file,clim.fx,read.dir,write.dir) {

   work <- paste0('cdo -s ',clim.fx,' ',
                  read.dir,input.file,' ',write.dir,clim.file)
   system(work)
   Sys.sleep(1)                     
##   return(clim.file)
}

make_average_file <- function(file,fxn,pattern,replacement,tmp.dir) {
   avg.file <- gsub(pattern=pattern,
                    replacement=replacement,
                    file)
   work <- paste0('cdo -O -s ',fxn,' ',tmp.dir,file,' ',tmp.dir,avg.file)
   system(work)
   return(avg.file)
}

##---------------------------------------------------------

fix_gcm_coords <- function(coord.file,clim.file,tmp.dir) {


   if (grepl('(FGOALS-g2|GFDL-ESM2G|GFDL-ESM2M|FGOALS-g3)',coord.file)) {
      nc <- nc_open(paste0(tmp.dir,coord.file))
      nlon <- nc$dim$lon$len
      nlat <- nc$dim$lat$len
      nc_close(nc)
      file.reg <- paste0("reg_fix_",coord.file)
      work <- paste0("cdo -s remapcon,r",nlon,"x",nlat," ",tmp.dir,coord.file," ",tmp.dir,file.reg)
      system(work)
      Sys.sleep(1)

      work <- paste0("cdo -s sellonlatbox,-180,180,-90,90 ",tmp.dir,file.reg," ",tmp.dir,clim.file)
      system(work)
      Sys.sleep(1)
      file.remove(paste0(tmp.dir,coord.file))
      file.remove(paste0(tmp.dir,file.reg))
   } else {

      work <- paste0("cdo -s sellonlatbox,-180,180,-90,90 ",tmp.dir,coord.file," ",tmp.dir,clim.file)
      system(work)
      Sys.sleep(1)
      file.remove(paste0(tmp.dir,coord.file))
   }
   ###return(file.fix)
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
   cmip <- 'CMIP6'
   scenario <- 'ssp126'
   type='monthly'
   cmip6.list <- list(list(gcm='CNRM-CM6-1',runs='r1i1p1f2'),
                      list(gcm='CNRM-ESM2-1',runs='r1i1p1f2'),
                      list(gcm='EC-Earth3-Veg',runs='r1i1p1f1'))

}

tmp.dir <-  paste0(tmpdir,'/',cmip,'_',scenario,'_',type,'/') ##
if (!file.exists(tmp.dir)) {
    dir.create(tmp.dir,recursive=T)
}

if (cmip=='CMIP6') {
   base.dir <- '/storage/data/climate/CMIP6/Derived/'
   gcm.list <- cmip6.list
} 

if (cmip=='CMIP5') {
   base.dir <- '/storage/data/climate/CMIP5/daily/Derived/'
   gcm.list <- cmip5.list
} 


##-------------------


intervals <- c('1971-2000','1981-2010','1986-2005','2011-2040','2041-2070','2071-2100','2081-2100')
##intervals <- '2081-2100'
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

      clean.files <- list.files(write.dir,pattern='ETCCDI_climatology',full.name=T)
      file.remove(clean.files)

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
            ##Subset by time
            for (interval in intervals) {

               time.file <- subset_by_time(var.file,interval,var.freq,
                                           tmp.dir,tmp.dir)

               if (var.freq=='mon') {
                  ##Monthly climatologies
                  
                  mon.file <- gsub(pattern=paste0(var.name,'_mon'),
                                   replacement=paste(var.name,'_monthly_average_climatology',sep=''),time.file)
                  mon.coord.file <- gsub(pattern=paste0(var.name,'_mon'),replacement=paste0(var.name,'_monthly_coord_fix'),time.file)                   
                  make_climatology(time.file,mon.coord.file,clim.fx='ymonmean',
                                                tmp.dir,tmp.dir)

                  fix_gcm_coords(mon.coord.file,mon.file,tmp.dir)
                  file.copy(from=paste0(tmp.dir,mon.file),to=write.dir,overwrite=TRUE)
                  file.remove(paste0(tmp.dir,mon.file))

                  ##Seasonal Climatologies
                  seas.fx <- get_agg_fxn(var.name)
                  seas.file <- gsub(pattern=paste0(var.name,'_mon'),
                                   replacement=paste(var.name,'_seasonal_average_climatology',sep=''),time.file)
                  seas.coord.file <- gsub(pattern=paste0(var.name,'_mon'),replacement=paste0(var.name,'_seasonal_coord_fix'),time.file)                   
                  avg.file <- make_average_file(time.file,fxn=paste0('seas',seas.fx),
                                   pattern=paste0(var.name,'_mon'),replacement=paste0(var.name,'_seasonal'),tmp.dir)
                  make_climatology(avg.file,seas.coord.file,clim.fx='yseasmean',
                                                tmp.dir,tmp.dir)
                  fix_gcm_coords(seas.coord.file,seas.file,tmp.dir)

                  file.copy(from=paste0(tmp.dir,seas.file),to=write.dir,overwrite=TRUE)
                  file.remove(paste0(tmp.dir,seas.file))

                  ##Annual Climatologies
                  ann.fx <- get_agg_fxn(var.name)
                  ann.file <- gsub(pattern=paste0(var.name,'_mon'),
                                   replacement=paste(var.name,'_annual_average_climatology',sep=''),time.file)
                  ann.coord.file <- gsub(pattern=paste0(var.name,'_mon'),replacement=paste0(var.name,'_annual_coord_fix'),time.file)                   
                  avg.file <- make_average_file(time.file,fxn=paste0('year',ann.fx),
                                   pattern=paste0(var.name,'_mon'),replacement=paste0(var.name,'_annual'),tmp.dir)
                  make_climatology(avg.file,ann.coord.file,clim.fx='timmean',
                                                tmp.dir,tmp.dir)
                  fix_gcm_coords(ann.coord.file,ann.file,tmp.dir)
                  file.copy(from=paste0(tmp.dir,ann.file),to=write.dir,overwrite=TRUE)
                  file.remove(paste0(tmp.dir,ann.file))

   
               } else { ##Annual climdex variables

                  clim.file <- gsub(pattern=paste0('_',var.freq,'_'),replacement='_annual_average_climatology_',time.file)
                  coord.file <- gsub(pattern=paste0('_',var.freq,'_'),replacement='_annual_average_coord_fix_',time.file)
                  make_climatology(time.file,coord.file,'timmean',tmp.dir,tmp.dir)
                  fix_gcm_coords(coord.file,clim.file,tmp.dir)
                  file.copy(from=paste0(tmp.dir,clim.file),to=write.dir,overwrite=TRUE)
                  Sys.sleep(1)
                  file.remove(paste0(tmp.dir,clim.file))
               }
               
               file.remove(paste0(tmp.dir,time.file))

            } 
            file.remove(paste0(tmp.dir,var.file))
         }

         ##--------------------------------------------
         ##Annual
         if (type=='annual') {
            var.avg <- switch(var.name,
                              pr='total',tasmax='average',tasmin='average')
         
            ##Subset by time
            for (interval in intervals) {
               time.file <- subset_by_time(var.file,interval,var.freq,
                                           tmp.dir,tmp.dir)
               coord.file <- gsub(pattern=paste0(var.freq,'_',var.avg),replacement=paste0(var.freq,'_',var.avg,'_coord_fix'),time.file)
               clim.file <- gsub(pattern=paste0(var.freq,'_',var.avg),replacement=paste0(var.freq,'_',var.avg,'_climatology'),time.file)
               make_climatology(time.file,coord.file,'timmean',tmp.dir,tmp.dir)

               fix_gcm_coords(coord.file,clim.file,tmp.dir)

               file.copy(from=paste0(tmp.dir,clim.file),to=write.dir,overwrite=TRUE)
               Sys.sleep(1)
               file.remove(paste0(tmp.dir,clim.file))
               file.remove(paste0(tmp.dir,time.file))
            }
         }
         ##--------------------------------------------
         ##Seasonal
         if (type=='seasonal') {
            var.avg <- switch(var.name,
                              pr='total',tasmax='average',tasmin='average')

            ##Subset by time
            for (interval in intervals) {
               time.file <- subset_by_time(var.file,interval,var.freq,
                                           tmp.dir,tmp.dir)
               coord.file <- gsub(pattern=paste0(var.freq,'_',var.avg),replacement=paste0(var.freq,'_',var.avg,'_coord_fix'),time.file)
               clim.file <- gsub(pattern=paste0(var.freq,'_',var.avg),replacement=paste0(var.freq,'_',var.avg,'_climatology'),time.file)
               make_climatology(time.file,coord.file,'yseasmean',tmp.dir,tmp.dir)
               fix_gcm_coords(coord.file,clim.file,tmp.dir)
               file.copy(from=paste0(tmp.dir,clim.file),to=write.dir,overwrite=TRUE)
               Sys.sleep(1)
               file.remove(paste0(tmp.dir,clim.file))
               file.remove(paste0(tmp.dir,time.file))
            }
         }
         ##--------------------------------------------
         ##Monthly
         if (type=='monthly') {
            var.avg <- switch(var.name,
                              pr='total',tasmax='average',tasmin='average')

            ##Subset by time
            for (interval in intervals) {
               time.file <- subset_by_time(var.file,interval,var.freq,
                                          tmp.dir,tmp.dir)
               coord.file <- gsub(pattern=paste0(var.freq,'_',var.avg),replacement=paste0(var.freq,'_',var.avg,'_coord_fix'),time.file)
               clim.file <- gsub(pattern=paste0(var.freq,'_',var.avg),replacement=paste0(var.freq,'_',var.avg,'_climatology'),time.file)
               make_climatology(time.file,coord.file,'ymonmean',tmp.dir,tmp.dir)
               fix_gcm_coords(coord.file,clim.file,tmp.dir)
               file.copy(from=paste0(tmp.dir,clim.file),to=write.dir,overwrite=TRUE)
               Sys.sleep(1)
               file.remove(paste0(tmp.dir,clim.file))
               file.remove(paste0(tmp.dir,time.file))
            }
         }
         ##--------------------------------------------
         ##Degree Days
         if (type=='degree_days') {
            ##Subset by time
            for (interval in intervals) {
               time.file <- subset_by_time(var.file,interval,var.freq,
                                          tmp.dir,tmp.dir)
               coord.file <- gsub(pattern=paste0(var.freq,'_'),replacement=paste0(var.freq,'_dd_coord_fix'),time.file)
               clim.file <- gsub(pattern=paste0(var.freq,'_'),replacement=paste0(var.freq,'_climatology_'),time.file)
               make_climatology(time.file,coord.file,'timmean',tmp.dir,tmp.dir)
               fix_gcm_coords(coord.file,clim.file,tmp.dir)
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
