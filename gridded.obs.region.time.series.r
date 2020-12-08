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

   file.name <- gcm.file ##basename(gcm.file)
   file.fix <- paste0("coord_fix_",gcm.file)
   if (length(file.name) > 1 ) { browser()}
   
      work <- paste0("cdo sellonlatbox,-180,180,-90,90 ",tmp.dir,file.name," ",tmp.dir,file.fix)
      system(work)
      Sys.sleep(1)
      file.remove(paste0(tmp.dir,file.name))

   return(file.fix)
}
 

##--------------------------------------------------------------------------------------

area_weighted_average <- function(clim.file,area.shape,tmp.dir) {


   fix.file <- basename(clim.file) ##fix_gcm_coords(clim.file,tmp.dir)

   file.brick <- brick(paste0(tmp.dir,fix.file))
   file.one <- subset(file.brick,1)

   area.overlay <- cell_mask(file.one,area.shape)
   area.masked <- mask(file.brick,area.overlay)
   cc <- area(area.masked) ##Returns area including blank regions
   cc.masked <- mask(cc,area.overlay)
   area.weights <- cc.masked / cellStats(cc.masked,sum)

   area.weighted <- area.masked * area.weights
   area.average <- cellStats(area.weighted,sum)

   dates <- file.brick@z$Date
   rm(file.brick)
   file.remove(paste0(tmp.dir,fix.file))

   rv <- list(dates=dates,series=area.average)
   return(rv)
}

##--------------------------------------------------------------------------------------

##Read in the regional time series
tas_average_time_series <- function(var.name,input.file,period,
                                    read.dir,tmp.dir,
                                    clip.shp,
                                    past.int) {
 
  file.copy(from=input.file,to=tmp.dir,overwrite=TRUE)
  file <- basename(input.file)

  ##-------------------------------------------------
  ##Extract subset of data for region
  reg.avg <- area_weighted_average(file,clip.shp,tmp.dir)

  ##Calculate anomalies
  years <- format(reg.avg$dates,'%Y')
  bounds <- strsplit(past.int,'-')[[1]]
  yst <- head(grep(bounds[1],years),1)
  yen <- tail(grep(bounds[2],years),1)

  base <- mean(reg.avg$series[yst:yen],na.rm=T)
  anomalies <- reg.avg$series - base

  rv <- list(dates=reg.avg$dates,anoms=anomalies)

  return(rv)
}


##*********************************************************************
##*********************************************************************

make_obs_time_series <- function(var.name,model,period,input.file,
                             region,scenario,clip.shp,
                             read.dir,write.dir,tmp.dir,past.int) {
   
    ##Read in the regional time series
    region.avgs <- tas_average_time_series(var.name=var.name,input.file,period,
                               read.dir=model.dir,tmp.dir=tmp.dir,
                               clip.shp=clip.shp,past.int)
    
    output <- rbind(c('Dates',toupper(model)),
                    cbind(as.character(region.avgs$dates),round(region.avgs$anoms,2)))

    bounds <- paste0(format(range(region.avgs$dates),'%Y'),collapse='-')
    print(bounds)
    write.file <- paste0(var.name,'_',model,'_',region,'_time_series_anomalies_obs_',bounds,'.csv')

    write.table(output,
                file=paste0(write.dir,write.file),
                quote=F,row.name=F,col.name=F,sep=',')

}

##---------------------------------------------------------------------
##*********************************************************************

tmp.dir <- "/local_temp/ssobie/obs_time_series/"

if (!file.exists(tmp.dir)) {
   dir.create(tmp.dir,recursive=TRUE)
}

##------------------------------------

regions <- c("northern_canada",'prairies','bc','ontario','quebec','atlantic_canada')

regions <-  'canada_boundary'

for (region in regions) {

period <- 'annual'
var.name <- 'tas'
past.int <- '1971-2000'

##------------------------------------

shape.dir <- "/storage/data/projects/rci/data/assessments/shapefiles/canada_regions/"
clip.shp <- get_shape_file(region,shape.dir)

write.dir <- paste0("/storage/data/projects/rci/data/cas/gridded_observations/time_series/")

if (!file.exists(write.dir))
  dir.create(write.dir,recursive=TRUE)

##HadCRU4
if (1==0) {
read.dir <- "/storage/data/projects/rci/data/nrcan/hadcru4/"
model <- 'hadcru4'
input.file <- paste0(read.dir,"HadCRUT.4.5.0.0.median.annual.nc")

make_obs_time_series(var.name,model,period,input.file,
                     region,scenario,clip.shp,
                     read.dir,write.dir,tmp.dir,past.int)
}

##GISTEMP

read.dir <- "/storage/data/projects/rci/data/nrcan/nasa_gistemp/"
model <- 'gistemp'
input.file <- paste0(read.dir,"gistemp1200_ERSSTv4_annual.nc")

make_obs_time_series(var.name,model,period,input.file,
                     region,scenario,clip.shp,
                     read.dir,write.dir,tmp.dir,past.int)

}