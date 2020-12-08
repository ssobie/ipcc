##Code to calculate properly weighted averages (by latitude)
##for gridded data in Canada and the provinces

library(rgdal)
library(raster)
library(ncdf4)
library(PCICt)
library(zoo)

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
##Correct for 0-360 longitude

correct_longitude_dimension <- function(input.file,output.file) {
   print('Fixing longitude dimension')                            
   work <- paste0("cdo sellonlatbox,-180,180,-90,90 ",input.file," ",output.file)
   system(work)
   Sys.sleep(1)
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

region_time_series <- function(input.file,tmp.dir,area.shape,season) {
                               

      file.brick <- brick(paste0(tmp.dir,input.file))
      time <- file.brick@z$Date
      yr.dates <- format(time,'%Y')
      yr.fac <- as.factor(yr.dates)

      seas.ix <- grep(season,time)
      file.sub <- subset(file.brick,seas.ix)

      area.cover <- area_weighting(file.sub,area.shape)

      area.series <- area_average_time_series(file.sub,area.cover)

      seas.dates <- time[seas.ix]


      past.ix <- seas.dates >= as.Date('1961-01-01') & seas.dates <= as.Date('1990-12-31')
      proj.ix <- seas.dates >= as.Date('2011-01-01') & seas.dates <= as.Date('2040-12-31')

      anom <- round( mean(area.series[proj.ix]) - mean(area.series[past.ix]), 1)
      
      return(list(series=area.series,anom=anom))

}
##----------------------------------------------------------

##**********************************************************

tmpdir <- '/local_temp/ssobie/'

region <- 'taiga_plains'
type <- 'monthly'

tmp.dir <- paste0(tmpdir,'taiga_test_',region,'_experiment/')
  
if (!file.exists(tmp.dir)) {
   dir.create(tmp.dir,recursive=TRUE)
}

gcm <- 'ACCESS1-0'
ssp <- 'rcp45'
run <- 'r1i1p1'

shape.dir <- '/storage/data/projects/rci/data/cas/canada/shapefiles/'
eco.dir <- '/storage/data/projects/rci/data/assessments/shapefiles/eco_provinces'

area.shape <- get_shape_file(region,eco.dir)

##---------------------------------------------------------------------------------------------------
##CMIP5 Raw
cmip5.dir <- paste0("/storage/data/climate/CMIP5/daily/Derived/",gcm,"_",ssp,"_",run,"/",type,"/")
gcm.dir <- cmip5.dir 

tasmax.file <- list.files(path=gcm.dir,pattern=paste0('tasmax_',type))
tasmin.file <- list.files(path=gcm.dir,pattern=paste0('tasmin_',type))

file.copy(from=paste0(gcm.dir,tasmax.file),to=tmp.dir,overwrite=TRUE)
file.copy(from=paste0(gcm.dir,tasmin.file),to=tmp.dir,overwrite=TRUE)
 
tas.file <- gsub(paste0('tasmax_',type),paste0('tas_',type),tasmax.file) 
make_tas_file(tasmax.file,tasmin.file,tas.file,tmp.dir)

lon.file <- gsub(paste0('tas_',type),paste0('tas_corr_',type),tas.file) 
correct_longitude_dimension(paste0(tmp.dir,tas.file),paste0(tmp.dir,lon.file))
 
cmip5.series <- region_time_series(lon.file,tmp.dir,area.shape,season='*-12-01')
croll <- rollmean(cmip5.series$series,31)
dec.dates <- seq(from=as.Date('1850-12-01'),by='year',to=as.Date('2100-12-01'))

##---------------------------------------------------------------------------------------------------
##BCCAQv2 Downscaled
bccaqv2.dir <- "/storage/data/climate/downscale/BCCAQ2/bccaqv2_derived/monthly/"

gcm.files <- list.files(path=bccaqv2.dir,pattern=gcm)
scen.files <- gcm.files[grep(ssp,gcm.files)]
tasmax.file <- scen.files[grep('tasmax',scen.files)]
tasmin.file <- scen.files[grep('tasmin',scen.files)]

file.copy(from=paste0(bccaqv2.dir,tasmax.file),to=tmp.dir,overwrite=TRUE)
file.copy(from=paste0(bccaqv2.dir,tasmin.file),to=tmp.dir,overwrite=TRUE)
 
tas.file <- gsub(paste0('tasmax_',type),paste0('tas_',type),tasmax.file) 
make_tas_file(tasmax.file,tasmin.file,tas.file,tmp.dir)

bccaqv2.series <- region_time_series(tas.file,tmp.dir,area.shape,season='*-12-01')


bd <- seq(from=as.Date('1950-12-01'),by='year',to=as.Date('2100-12-01'))
broll <- rollmean(bccaqv2.series$series,31)


plot(dec.dates,cmip5.series$series,col='red',type='l',lwd=2,ylim=c(-35,-5),ylab='Dec Average Temperature (decC)')
lines(dec.dates[15:235],croll,col='red',lwd=3)
lines(bd,bccaqv2.series$series,col='blue',lwd=2)
lines(bd[15:135],broll,col='blue',lwd=3)
box(which='plot')

browser()


##Remap to common grid
##tas.one.deg.file <- paste0("tas_seasonal_",gcm,"_",ssp,"_",run,"_one_degree.nc")
##regrid_remapcon(tas.file,tas.one.deg.file,tmp.dir)

##tas.series.file <- paste0("tas_annual_series_",region,"_",gcm,"_",ssp,"_",run,"_1971-2010.RData")
##trend_and_projection(tas.one.deg.file,tmp.dir,area.shape,
##                     tas.trend.file,tas.anom.file,tas.series.file,save.dir)

##Clean up
file.remove(paste0(tmp.dir,tasmax.file))
file.remove(paste0(tmp.dir,tasmin.file))
file.remove(paste0(tmp.dir,tas.file))
file.remove(paste0(tmp.dir,tas.one.deg.file)) 




