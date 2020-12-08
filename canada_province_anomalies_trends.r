##Code to calculate properly weighted averages (by latitude)
##for gridded data in Canada and the provinces

library(rgdal)
library(raster)
library(ncdf4)
library(PCICt)

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

##----------------------------------------------------------
##Read a shapefile

get_shape_file <- function(region,shape.dir) {
   ##shape.dir <- '/storage/data/projects/rci/data/assessments/shapefiles/bc_common/'
   reg.shp <- readOGR(shape.dir,region,stringsAsFactors=F)
   return(reg.shp)
}

##----------------------------------------------------------
## Find the overlapping cells
##This gets the masked values including any overlapping cells
##Normal mask omits the partially covered cells
## **This step takes a while so it should be done once, then
## saved for future use **
##Use a spatial copy of the file with one time step for speed

cell_mask <- function(file.one,shape) {

   ras <- rasterize(shape,file.one,getCover=T)
   ras[ras==0] <- NA
   return(ras)
}

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
regrid_remapcon <- function(input.file,output.file,tmp.dir) {

   work <- paste0("cdo remapcon,/storage/home/ssobie/grid_files/canada.one.deg.grid.txt ",
                   tmp.dir,input.file," ",
                   tmp.dir,output.file)
   system(work)
   Sys.sleep(1)
}                 

##----------------------------------------------------------
##Find the area weighting for the average

area_weighting <- function(file.brick,area.shape) {

   brick.one <- subset(file.brick,10)
   area.overlay <- cell_mask(brick.one,area.shape)
   area.masked <- mask(brick.one,area.overlay)
   cc <- area(area.masked) ##Returns area including blank regions
   cc.masked <- mask(cc,area.overlay)
   area.weights <- cc.masked / cellStats(cc.masked,sum)

   rv <- list(weights=area.weights,mask=area.overlay)
   return(rv)       
}

##----------------------------------------------------------
area_average_time_series <- function(file.brick,area.cover) {

   area.series <- mask(file.brick,area.cover$mask)
   area.weighted.series <- area.series * area.cover$weights
   area.time.series <- cellStats(area.weighted.series,sum)
   return(area.time.series)
}

##----------------------------------------------------------

trend_and_projection <- function(one.deg.file,tmp.dir,area.shape,
                                 trend.file,anom.file,series.file,save.dir) {

      file.brick <- brick(paste0(tmp.dir,one.deg.file))
      time <- file.brick@z$Date
      yr.dates <- format(time,'%Y')
      yr.fac <- as.factor(yr.dates)

      area.cover <- area_weighting(file.brick,area.shape)

      area.series <- area_average_time_series(file.brick,area.cover)
 
      year.series <- tapply(area.series,yr.fac,mean)
      yrs <- as.numeric(levels(yr.fac))
      base.ix <- yrs >= 1970 & yrs <= 2010
      past.ix <- yrs >= 1981 & yrs <= 2010
      proj.ix <- yrs >= 2071 & yrs <= 2100
      base.series <- data.frame(time=as.numeric(levels(yr.fac))[base.ix],ts=year.series[base.ix])
      lm.fit <- lm(ts~time,base.series)

      sm.fit <- summary(lm.fit)

      lm.info <- list(coeffs=sm.fit$coefficients,
                 series=base.series)
   
      save(lm.info,file=paste0(save.dir,'trends/',trend.file))

      ##Projected change
      past.series <- data.frame(time=as.numeric(levels(yr.fac))[past.ix],ts=year.series[past.ix])
      proj.series <- data.frame(time=as.numeric(levels(yr.fac))[proj.ix],ts=year.series[proj.ix])
      anomaly <- round(mean(year.series[proj.ix],na.rm=T) - mean(year.series[past.ix],na.rm=T),2)
      percent <- round(anomaly / mean(year.series[past.ix],na.rm=T) * 100,2)

      proj.info <- list(anom=anomaly,
                        past=past.series,
                        proj=proj.series)

      save(proj.info,file=paste0(save.dir,'anomalies/',anom.file))

      series.info <- list(years=yrs,
                          series=year.series)
      save(series.info,file=paste0(save.dir,'series/',series.file))

}
##----------------------------------------------------------

##**********************************************************
testing <- FALSE

if (testing) {
   experiment <- 'CMIP5'   
   tmpdir <- '/local_temp/ssobie/can_at'
   region <- 'canada_boundary'

} else {

 
}

tmp.dir <- paste0(tmpdir,'trend_anom_series_',region,'_experiment/')
  
if (!file.exists(tmp.dir)) {
   dir.create(tmp.dir,recursive=TRUE)
}



shape.dir <- '/storage/data/projects/rci/data/cas/canada/shapefiles/'

area.shape <- get_shape_file(region,shape.dir)

##Model trends and anomalies

##CMIP5
if (experiment == 'CMIP5') {
   files <- list.files(path="/storage/data/climate/CMIP5/daily/Derived/")
   files <- files[-grep("HadCM3",files)]
   files <- files[-grep("CanCM4",files)]
   
   save.dir <- "/storage/data/projects/rci/data/cas/cmip5/"
}

##CMIP6
if (experiment == 'CMIP6') {
   ##gcm <- 'CanESM5'
   ##ssp <- 'ssp245'
   ##run <- 'r2i1p2f1'
   ##Omit the P1 CanESM5
   files <- list.files(path="/storage/data/climate/CMIP6/Derived/")
   files <- files[-grep('CanESM5_ssp245_r[0-9]i1p1f1',files)]
   save.dir <- "/storage/data/projects/rci/data/cas/cmip6/"
   ##SKip over CNRM-CM6-1-HR for now
   ##files <- files[58:60]
   ##files <- files[90:length(files)]
}


for (file in files) {
   print(file)
   file.split <- strsplit(file,'_')[[1]]
   gcm <- file.split[1]
   ssp <- file.split[2]
   run <- file.split[3]

   if (any(grepl('(gr1|gr2)',file.split))) {
      res <- file.split[4]
      gcm.dir <- paste0("/storage/data/climate/CMIP6/Derived/",gcm,"_",ssp,"_",run,"_",res,"/annual/")
   } else {
      gcm.dir <- switch(experiment,
                 CMIP6=paste0("/storage/data/climate/CMIP6/Derived/",gcm,"_",ssp,"_",run,"/annual/"),
                 CMIP5=paste0("/storage/data/climate/CMIP5/daily/Derived/",gcm,"_",ssp,"_",run,"/annual/"))
   }
   tasmax.file <- list.files(path=gcm.dir,pattern='tasmax_annual')
   tasmin.file <- list.files(path=gcm.dir,pattern='tasmin_annual')
   pr.file <- list.files(path=gcm.dir,pattern='pr_annual')

   if (length(tasmax.file) > 1 | length(tasmin.file) > 1 | length(pr.file) > 1) {
      print(tasmax.file)
      print('-')
      print(tasmin.file)
      print('-')
      print(pr.file)
      print('Too many GCM files selected')
   } else if (length(tasmax.file) == 0 | length(tasmin.file) == 0 | length(pr.file) == 0) {
      print(tasmax.file)
      print('-')
      print(tasmin.file)
      print('-')
      print(pr.file)
      print('Missing one of the files selected')
   } else {
   
      file.copy(from=paste0(gcm.dir,tasmax.file),to=tmp.dir,overwrite=TRUE)
      file.copy(from=paste0(gcm.dir,tasmin.file),to=tmp.dir,overwrite=TRUE)
      file.copy(from=paste0(gcm.dir,pr.file),to=tmp.dir,overwrite=TRUE)
 
      tas.file <- gsub('tasmax_annual','tas_annual',tasmax.file) 
      make_tas_file(tasmax.file,tasmin.file,tas.file,tmp.dir)

      ##Remap to common grid
      tas.one.deg.file <- paste0("tas_annual_",gcm,"_",ssp,"_",run,"_one_degree.nc")
      pr.one.deg.file <- paste0("pr_annual_",gcm,"_",ssp,"_",run,"_one_degree.nc")

      regrid_remapcon(tas.file,tas.one.deg.file,tmp.dir)
      regrid_remapcon(pr.file,pr.one.deg.file,tmp.dir)

      pr.trend.file <- paste0("pr_annual_trend_",gcm,"_",ssp,"_",run,"_1971-2010.RData")
      pr.anom.file <- paste0("anomaly_pr_annual_",region,"_",gcm,"_",ssp,"_",run,"_1971-2010.RData")
      pr.series.file <- paste0("pr_annual_series_",region,"_",gcm,"_",ssp,"_",run,"_1971-2010.RData")
      trend_and_projection(pr.one.deg.file,tmp.dir,area.shape,
                           pr.trend.file,pr.anom.file,pr.series.file,save.dir)

      tas.trend.file <- paste0("tas_annual_trend_",region,"_",gcm,"_",ssp,"_",run,"_1971-2010.RData")
      tas.anom.file <- paste0("anomaly_tas_annual_",region,"_",gcm,"_",ssp,"_",run,"_1971-2010.RData")
      tas.series.file <- paste0("tas_annual_series_",region,"_",gcm,"_",ssp,"_",run,"_1971-2010.RData")
      trend_and_projection(tas.one.deg.file,tmp.dir,area.shape,
                           tas.trend.file,tas.anom.file,tas.series.file,save.dir)


      ##Clean up
      file.remove(paste0(tmp.dir,tasmax.file))
      file.remove(paste0(tmp.dir,tasmin.file))
      file.remove(paste0(tmp.dir,tas.file))
      file.remove(paste0(tmp.dir,pr.file))      
      file.remove(paste0(tmp.dir,tas.one.deg.file)) 
      file.remove(paste0(tmp.dir,pr.one.deg.file)) 
   }
###   browser()

}



