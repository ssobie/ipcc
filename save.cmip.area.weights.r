##Script to produce tables of projected temperature and precipitation changes
##The output format is the same as the MOTI tables 

library(raster)
library(ncdf4)
library(rgeos)
library(PCICt)
library(rgdal)

##-----------------------------------------------------------------------------------------------
##Area average functions

get_shape_file <- function(region,shape.dir) {
   reg.shp <- readOGR(shape.dir,region,stringsAsFactors=F)
   return(reg.shp)
}

cell_mask <- function(file.brick,shape) {
   ras <- rasterize(shape,file.brick,getCover=T)
   ras[ras==0] <- NA
   return(ras)
}

##--------------------------------------------------------------------------------------

fix_gcm_coords <- function(gcm.file,tmp.dir) {

   file.name <- basename(gcm.file)
   file.fix <- paste0("coord_fix_",basename(gcm.file))
   if (length(file.name) > 1 | length(file.name)==0) { browser()}
   
   if (grepl('(FGOALS-g2|GFDL-ESM2G|GFDL-ESM2M|FGOALS-g3)',file.name)) {
      nc <- nc_open(paste0(tmp.dir,file.name))
      nlon <- nc$dim$lon$len
      nlat <- nc$dim$lat$len
      nc_close(nc)
      file.reg <- paste0("reg_fix_",file.name)
      work <- paste0("cdo -s remapcon,r",nlon,"x",nlat," ",tmp.dir,file.name," ",tmp.dir,file.reg)
      system(work)
      Sys.sleep(1)

      work <- paste0("cdo -s sellonlatbox,-180,180,-90,90 ",tmp.dir,file.reg," ",tmp.dir,file.fix)
      system(work)
      Sys.sleep(1)
      file.remove(paste0(tmp.dir,file.name))
      file.remove(paste0(tmp.dir,file.reg))  
   } else {
  
      work <- paste0("cdo -s sellonlatbox,-180,180,-90,90 ",tmp.dir,file.name," ",tmp.dir,file.fix)
      system(work)
      Sys.sleep(1)
      file.remove(paste0(tmp.dir,file.name))
   }
   return(file.fix)
}
 

##--------------------------------------------------------------------------------------
##Split apart area masked to remove area weight calculation to happen once

area_weights <- function(clim.file,area.shape,tmp.dir) {

   file.copy(from=clim.file,to=tmp.dir,overwrite=TRUE)
   fix.file <- basename(clim.file) ###fix_gcm_coords(clim.file,tmp.dir)

   file.brick <- brick(paste0(tmp.dir,fix.file))
   file.one <- subset(file.brick,1)
   area.overlay <- cell_mask(file.one,area.shape)
   cc <- area(file.one) ##Returns area including blank regions
   cc.masked <- mask(cc,area.overlay)
   area.weights <- cc.masked / cellStats(cc.masked,sum)

   rv <- list(overlay=area.overlay,weights=area.weights)
   return(rv)
}

##--------------------------------------------------------------------------------------

##Read in the regional time series
calculate_area_weights <- function(read.dir,tmp.dir,
                                   clip.shp) {
                                   

  ##GCM and Regional averages grouping monthly and seasonal values - needs pre-computed seasonal and monthly files
  ann.files <- list.files(path=paste(read.dir,'annual/climatologies',sep=''),pattern=paste0('tasmax_annual'),full.name=TRUE)
  ###ann.files <- list.files(path=paste(read.dir,'annual/climatologies',sep=''),pattern=paste0('pr_annual'),full.name=TRUE)
  ann.file <- ann.files[grep('1971-2000',ann.files)]
  print(ann.file)

  rv <- area_weights(ann.file,clip.shp,tmp.dir)

  return(rv)
}


##---------------------------------------------------------------------
##*********************************************************************

cmip <- "CMIP5"
scenario <- "rcp26"

regions <- c("quebec","ontario","bc","prairies","northern_canada","atlantic_canada","canada_boundary")
tmpdir <- "/local_temp/ssobie/area_regions"

tmp.dir <- paste0(tmpdir,'/areas_',cmip,'_',scenario,'/')


if (!file.exists(tmp.dir)) {
   dir.create(tmp.dir,recursive=TRUE)
}


if (cmip=='CMIP6') {
   read.dir <- "/storage/data/climate/CMIP6/Derived/"
}
if (cmip=='CMIP5') {
   read.dir <- "/storage/data/climate/CMIP5/daily/Derived/"
}


shape.dir <- "/storage/data/projects/rci/data/assessments/shapefiles/canada_regions/"

save.base <- paste0("/storage/data/projects/rci/data/cas/",tolower(cmip),"/area_info/")

cmip5.list <- list(c('ACCESS1-0','r1i1p1'),
                   c('bcc-csm1-1','r1i1p1'),
                   c('bcc-csm1-1-m','r1i1p1'),
                   c('BNU-ESM','r1i1p1'),
                   c('CanESM2','r1i1p1'),
                   c('CCSM4','r2i1p1'),
                   c('CESM1-CAM5','r1i1p1'),
                   c('CNRM-CM5','r1i1p1'),
                   c('CSIRO-Mk3-6-0','r1i1p1'),
                   c('FGOALS-g2','r1i1p1'),
                   c('GFDL-CM3','r1i1p1'),
                   c('GFDL-ESM2G','r1i1p1'),
                   c('GFDL-ESM2M','r1i1p1'),
                   c('HadGEM2-AO','r1i1p1'),
                   c('HadGEM2-CC','r1i1p1'),
                   c('HadGEM2-ES','r1i1p1'),
                   c('inmcm4','r1i1p1'),
                   c('IPSL-CM5A-LR','r1i1p1'),
                   c('IPSL-CM5A-MR','r1i1p1'),
                   c('MIROC-ESM','r1i1p1'),
                   c('MIROC-ESM-CHEM','r1i1p1'),
                   c('MIROC5','r3i1p1'),
                   c('MPI-ESM-LR','r3i1p1'),
                   c('MPI-ESM-MR','r1i1p1'),
                   c('MRI-CGCM3','r1i1p1'),
                   c('NorESM1-M','r1i1p1'),
                   c('NorESM1-ME','r1i1p1'))
##cmip5.list <- list(c('ACCESS1-0','r1i1p1'),
##                   c('bcc-csm1-1','r1i1p1'),
##                   c('bcc-csm1-1-m','r1i1p1'))
cmip5.list <- list(c('CESM1-CAM5','r1i1p1'))


cmip6.list <- list(c('ACCESS-CM2','r1i1p1f1'),
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

##cmip6.list <- list(c('ACCESS-CM2','r1i1p1f1'),
##                   c('ACCESS-ESM1-5','r1i1p1f1'),
##                   c('BCC-CSM2-MR','r1i1p1f1'))

model.dirs <- switch(cmip,
                     CMIP6=cmip6.list,
                     CMIP5=cmip5.list) 

model.list <- unlist(lapply(model.dirs,function(x,y){paste(x[1],y,x[2],sep='_')},scenario))

past.int <- "1971-2000"


for (region in regions) {
   clip.shp <- get_shape_file(region,shape.dir)

   save.dir <- paste0(save.base,region,'/')
   if (!file.exists(save.dir)) {
      dir.create(save.dir,recursive=TRUE)
   }


   for (m in seq_along(model.list)) {
       model <- model.list[m]
       print(model)
       gcm <- strsplit(model,'_')[[1]][1]

       
       model.dir <- paste0(read.dir,model,'/')
       save.file <- paste0(save.dir,region,'_',gcm,'_rcp26_area_info.RData')

       area.info <- calculate_area_weights(model.dir,tmp.dir,clip.shp)
       save(area.info,file=save.file)

   }    

}