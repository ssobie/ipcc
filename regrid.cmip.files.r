##Code to calculate properly weighted averages (by latitude)
##for gridded data in Canada and the provinces

library(rgdal)
library(raster)
library(ncdf4)
library(PCICt)

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

##----------------------------------------------------------
##Create an average temperature file
make_tas_file <- function(tasmax.file,tasmin.file,tas.file,tmp.dir) {

   work1 <- paste0("cdo add ",
                   tmp.dir,tasmax.file," ",
                   tmp.dir,tasmin.file," ",tmp.dir,"tmp.nc")
   system(work1)
   Sys.sleep(1)
   work2 <- paste0("cdo divc,2 ",tmp.dir,"tmp.nc ",
                   tmp.dir,tas.file)
   system(work2)                   
   Sys.sleep(1)
   file.remove(paste0(tmp.dir,"tmp.nc"))
}                 

##----------------------------------------------------------
##Regrid GCM file to one degree grid
regrid_remapcon <- function(input.dir,input.file,output.dir,output.file) {
###/storage/home/ssobie/grid_files/canada.one.deg.grid.txt "
   work <- paste0("cdo -s remapcon,/storage/home/ssobie/grid_files/north.hemis.one.deg.grid.txt ",
                   input.dir,input.file," ",
                   output.dir,output.file)
   system(work)
   Sys.sleep(1)
}                 

##----------------------------------------------------------
##**********************************************************

testing <- TRUE

if (testing) {
   tmpdir <- '/local_temp/ssobie'

   experiment <- 'CMIP6'
   type <- 'return_periods'
   scenario <- 'ssp126'
} else {

   args <- commandArgs(trailingOnly=TRUE)
   for(i in 1:length(args)){
       eval(parse(text=args[[i]]))
   }
   
}


   ##CMIP5
if (experiment == 'CMIP5') {
   one.dir <- "/storage/data/projects/rci/data/cas/cmip5/one_degree/climatologies/"
}

##CMIP6
if (experiment == 'CMIP6') {
   

   gcm.list <- list(c('ACCESS-CM2','r1i1p1f1'),
                   c('ACCESS-ESM1-5','r1i1p1f1'),
                   c('BCC-CSM2-MR','r1i1p1f1'),
                   c('CanESM5','r1i1p2f1'),
                   c('CNRM-CM6-1','r1i1p1f2'),
                   c('CNRM-ESM2-1','r1i1p1f2'),
                   c('EC-Earth3','r4i1p1f1'),
                   c('EC-Earth3-Veg','r1i1p1f1'),
                   c('FGOALS-g3','r1i1p1f1'),
                   c('GFDL-ESM4','r1i1p1f1'),
                   c('HadGEM3-GC31-LL','r1i1p1f3'),
                   c('INM-CM4-8','r1i1p1f1'),
                   c('INM-CM5-0','r1i1p1f1'),
                   c('IPSL-CM6A-LR','r1i1p1f1'),
                   c('KACE-1-0-G','r2i1p1f1'),
                   c('KIOST-ESM','r1i1p1f1'),
                   c('MIROC6','r1i1p1f1'),
                   c('MIROC-ES2L','r1i1p1f2'),
                   c('MPI-ESM1-2-HR','r1i1p1f1'),
                   c('MPI-ESM1-2-LR','r1i1p1f1'),
                   c('MRI-ESM2-0','r1i1p1f1'),
                   c('NESM3','r1i1p1f1'),
                   c('NorESM2-LM','r1i1p1f1'),
                   c('NorESM2-MM','r1i1p1f1'),
                   c('UKESM1-0-LL','r1i1p1f2'))

  gcm.list <- list(c('HadGEM3-GC31-LL','r1i1p1f3'))

   ##Omit the P1 CanESM5
   one.dir <- "/storage/data/projects/rci/data/cas/cmip6/one_degree/climatologies/"

}


   tmp.dir <- paste0(tmpdir,'/',experiment,'_',type,'_regridding/')
   if (!file.exists(tmp.dir)) {
      dir.create(tmp.dir,recursive=TRUE)
   }


   for (gcm.info in gcm.list) {

      gcm <- gcm.info[[1]]
      run <- gcm.info[[2]]
      print(paste0(gcm,'_',scenario,'_',run))

      if (experiment=='CMIP6') {
           if (type == 'return_periods') {
              gcm.dir <- paste0("/storage/data/climate/CMIP6/Derived/",gcm,"_",scenario,"_",run,"/return_periods/")
           } else {
              gcm.dir <- paste0("/storage/data/climate/CMIP6/Derived/",gcm,"_",scenario,"_",run,"/",type,"/climatologies/")
           }

           write.dir <- paste0(one.dir,gcm,"_",scenario,"_",run,"/",type,"/")
      } 
      if (experiment=='CMIP5') {
           gcm.dir <- paste0("/storage/data/climate/CMIP5/daily/Derived/",gcm,"_",scenario,"_",run,"/",type,"/")
           write.dir <- paste0(one.dir,gcm,"_",scenario,"_",run,"/",type,"/")
      }


      if (!file.exists(write.dir)) {
         dir.create(write.dir,recursive=T)
      }

      clim.files <- list.files(path=gcm.dir,pattern='.nc')

      if (length(clim.files) == 0) {
         print('Missing one of the files selected')
         browser()
      } 
 
      for (file in clim.files) {
         ##Remap to common grid
         one.deg.file <- gsub(paste0("_",gcm,"_"),paste0("_",gcm,"_regridded_"),file)
         regrid_remapcon(gcm.dir,file,write.dir,one.deg.file)


      }
   }

 


