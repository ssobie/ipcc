##Set up and submit jobs to run BCCS and BCCA for the specified models
##

submit_rps_job <- function(var.name,gcm,run,scenario,
                           interval,deg,rp,
                           gcm.file,
                           gcmdir,writedir,
                           rps.template) {     

      rps.array <- readLines(rps.template)

      ##RPS
      ##Output located at line 10
      rps.array[6] <- paste0("#PBS -o ./out_dir/",gcm,"-",var.name,"-",run,"-rcp",substr(scenario,15,17),"-",interval,".out")
      rps.array[7] <- paste0("#PBS -e ./out_dir/",gcm,"-",var.name,"-",run,"-rcp",substr(scenario,15,17),"-",interval,".err")
      rps.array[8] <- paste0("#PBS -N ",tolower(substr(gcm,1,3)),".",substr(run,1,3),".",substr(scenario,15,17))

      rps.array[14] <- paste0('gcm="',gcm,'"')
      rps.array[15] <- paste0('run="',run,'"')
      rps.array[16] <- paste0('scenario="',scenario,'"')
      rps.array[17] <- paste0('rp="',rp,'"')
      rps.array[18] <- paste0('interval="',interval,'"')
      rps.array[19] <- paste0('deg="',deg,'"')
      rps.array[20] <- paste0('varname="',var.name,'"')
      rps.array[21] <- paste0('gcmfile="',gcm.file,'"')
      rps.array[22] <- paste0('gcmdir="',gcmdir,'"')
      rps.array[23] <- paste0('writedir="',writedir,'"')

      rps.file <- paste0("./sub-RP.",var.name,".pbs")
      writeLines(rps.array,rps.file)
      work <- paste0("qsub ",rps.file)
      print(work)
      system(work)
      Sys.sleep(2)

}

##---------------------------------------------------------
deg_anomaly_intervals <- function(gcm,run,deg.anoms) {

   dgx <- deg.anoms$Model == gcm & deg.anoms$Run == run

   deg.years <- deg.anoms[dgx,3:5] ##1,2,3 degree warming

   deg.intervals <- sapply(deg.years,function(x){paste0( (x-14),'-',(x+15))})

   return(deg.intervals)

}
##---------------------------------------------------------


##*****************************************************************************
##---------------------------------------------------------------------------

cmip.dir <- '/storage/data/climate/CMIP6/Derived/'
deg.file <- paste0("/storage/data/projects/rci/data/cas/cmip6/degree_climatologies/",
                    "global_cmip6_ssp585_tas_degree_anomaly_years.csv")

###cmip.dir <- '/storage/data/climate/CMIP5/daily/Derived/'

deg.anoms <- read.csv(deg.file,header=T,as.is=T)

rps.template <- "/storage/home/ssobie/code/repos/ipcc/template.run.return.periods.pbs"

##---------------------------------------------------------------------------

all.gcms <- list.files(path=cmip.dir)

##scen.gcms <- all.gcms[grep('ssp585',all.gcms)]  ##For Degree Anomaly Return Periods
##gcm <- 'ACCESS-CM2'
##gcm.list <- scen.gcms[grep(gcm,scen.gcms)]

cmip6.list <- list(list(gcm='ACCESS-CM2',runs='r1i1p1f1'),
                 list(gcm='ACCESS-ESM1-5',runs='r1i1p1f1'),
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
                 list(gcm='CanESM5',runs='r1i1p2f1'))


gcm.list <- list(list(gcm='bcc-csm1-1',runs='r1i1p1'),
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

gcm.list <- list(list(gcm='CNRM-CM6-1',runs='r1i1p1f2'),
                 list(gcm='CNRM-ESM2-1',runs='r1i1p1f2'))

scenario <- "ssp126"

varnames <- c('pr','tasmax','tasmin') ##

intervals <- c('1971-2000','1981-2010','2011-2040','2041-2070','2071-2100') ## ###

##---------------------------------------------------------------------------

for (gcm.info in gcm.list) {
   print(gcm.info)
   file.parts <- gcm.info ##strsplit(gcm.info,'_')[[1]]

   gcm <- file.parts$gcm ##[1]
   ###scenario <-  ##file.parts[2]
   run <- file.parts$runs ##[3]

   ##Use intervals from global degree anomaly timing
   ###intervals <- deg_anomaly_intervals(gcm,run,deg.anoms)  

   gcmdir <- paste0(cmip.dir,gcm,'_',scenario,'_',run,'/annual_extremes/')
   writedir <- paste0(cmip.dir,gcm,'_',scenario,'_',run,'/return_periods/') ###degree_anomaly_rps/')

   for (varname in varnames) {
      print(varname)
      var.files <- list.files(gcmdir,patter=paste0(varname,'_annual'))

      for (i in seq_along(var.files)) {
         var.file <- var.files[i]

         if (varname=='pr') {
            rps <- c('5','20','50')
         } else {
            rps <- c('5','20')
         }
         for (rp in rps) {
            print(rp)
            for (interval in intervals) {
               print(interval)
               deg <- NULL ###c('one','two','three')[which(intervals %in% interval)]

               ##Check for existing file
               all.files <- list.files(path=writedir,pattern=paste0(varname,'_RP',rp,'_'))
               check.file <- NULL ###all.files[grep(interval,all.files)]

               if (length(check.file) == 1) {
                  print('Existing file:')
                  print(check.file)
               } else if (length(check.file) > 1) {     
                  browser()
               } else {
                  submit_rps_job(varname,gcm,run,scenario,
                              interval,deg,rp,
                              var.file,
                              gcmdir,writedir,   
                              rps.template)
               }
            }##Intervals
         }##RPs
     }##Var.files
  }##Var.names

}##GCMs



