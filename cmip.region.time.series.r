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

area_weighted_average <- function(clim.file,area.shape,tmp.dir) {

   if (length(clim.file) > 1 ) { browser()}
   
   fix.file <- fix_gcm_coords(clim.file,tmp.dir)
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

   rv <- list(dates=dates,series=area.average)
   return(rv)
}



new_area_weighted_average <- function(clim.file,area.shape,area.info,tmp.dir) {

   
   fix.file <- fix_gcm_coords(clim.file,tmp.dir)

   file.brick <- brick(paste0(tmp.dir,fix.file))
   dates <- file.brick@z$Date

   dates.ix <- dates >= as.Date('1850-01-01') & dates <= as.Date('2100-12-31') 

   area.masked <- mask(file.brick,area.info$overlay)

   area.weighted <- area.masked * area.info$weights
   area.average <- cellStats(area.weighted,sum)

   rm(file.brick)
   file.remove(paste0(tmp.dir,fix.file))

   rv <- list(dates=dates[dates.ix],series=area.average[dates.ix])
   return(rv)

}


##--------------------------------------------------------------------------------------

##Read in the regional time series
t_p_average_time_series <- function(var.name,period,
                                    read.dir,tmp.dir,
                                    clip.shp,area.info,
                                    past.int) {
 
  if (var.name == 'tas') {                                   
     print('Making tas file')    
     tasmax.file <- list.files(path=paste0(read.dir,period,'/'),pattern=paste0('tasmax_',period,'_'),full.name=TRUE)
       if (length(tasmax.file) > 1 | length(tasmax.file) == 0) { browser()}
     tas.file <- gsub('tasmax_','tas_',basename(tasmax.file))
     tasmin.file <- list.files(path=paste0(read.dir,period,'/'),pattern=paste0('tasmin_',period,'_'),full.name=TRUE)
       if (length(tasmin.file) > 1 | length(tasmin.file) == 0) { browser()}
     file.copy(from=tasmax.file,to=tmp.dir,overwrite=TRUE)
     file.copy(from=tasmin.file,to=tmp.dir,overwrite=TRUE)
     work <- paste0('cdo -s add ',tmp.dir,basename(tasmax.file),' ',tmp.dir,basename(tasmin.file),' ',tmp.dir,'tas_add.nc')
     system(work)
     Sys.sleep(1)
     work <- paste0('cdo -s divc,2 ',tmp.dir,'tas_add.nc ',tmp.dir,tas.file)
     system(work)
     Sys.sleep(1)
     file <- tas.file
     file.remove(paste0(tmp.dir,'tas_add.nc'))
  } else { 

     ##GCM and Regional averages grouping monthly and seasonal values - needs pre-computed seasonal and monthly files
     input.file <- list.files(path=paste0(read.dir,period,'/'),pattern=paste0(var.name,'_',period,'_'),full.name=TRUE)
       if (length(input.file) > 1 | length(input.file) == 0) { browser()}
     file.copy(from=input.file,to=tmp.dir,overwrite=TRUE)
     file <- basename(input.file)
  }

  ##-------------------------------------------------
  ##Extract subset of data for region
  reg.avg <- new_area_weighted_average(file,clip.shp,area.info,tmp.dir)

  ##Calculate anomalies
  years <- format(reg.avg$dates,'%Y')
  bounds <- strsplit(past.int,'-')[[1]]
  yst <- head(grep(bounds[1],years),1)
  yen <- tail(grep(bounds[2],years),1)

  base <- mean(reg.avg$series[yst:yen],na.rm=T)
  anomalies <- reg.avg$series - base

  if (grepl('HadGEM2-ES',read.dir)) {
     ydates <- as.Date(paste0(format(reg.avg$dates,'%Y'),'-01-01'))
     reg.avg$dates <- ydates
  }

  rv <- list(dates=reg.avg$dates,anoms=anomalies)

  return(rv)
}

##-------------------------------------------------

load_area_info <- function(region,model,var.name,scenario,load.dir) {

  pr.vars <- c('pr',
               'prcptotETCCDI','sdiiETCCDI',
               'r1mmETCCDI','r10mmETCCDI','r20mmETCCDI',
               'rx1dayETCCDI','rx2dayETCCDI','rx5dayETCCDI',
               'r95pETCCDI','r95daysETCCDI','r99pETCCDI','r99daysETCCDI',
               'pr_rp5','pr_rp20','pr_rp50',
               'cddETCCDI','cdd90ETCCDI','cddmaxETCCDI','cwdETCCDI')
   var.pr <- FALSE
   var.pr <- any(pr.vars %in% var.name)

   if (model=='KIOST-ESM' & var.pr) {
      load.file <- paste0(load.dir,region,'_KIOST-ESM_pr_area_info.RData')
   } else if (model == 'CESM1-CAM5' & scenario == 'rcp26') {
      load.file <- paste0(load.dir,region,'_CESM1-CAM5_rcp26_area_info.RData')
   } else {
      load.file <- paste0(load.dir,region,'_',model,'_area_info.RData')
   }

   load(load.file)
   return(area.info)
}


##*********************************************************************
##*********************************************************************

make_time_series <- function(model.list,var.name,period,ts_function,
                             cmip,region,scenario,clip.shp,
                             read.dir,write.dir,area.dir,tmp.dir,past.int) {
   
    region.avgs <- vector(mode='list',length=length(model.list))
    gcm.list <- run.list <- rep('',length(model.list))

    years <- 1850:2100
    seas <- as.vector(sapply(years,function(x){paste(x,c('-01-15','-04-15','-07-15','-10-15'),sep='')}))
    months <- paste0('-',sprintf('%02d',1:12),'-01')
    mons <- as.vector(sapply(years,function(x){paste(x,months,sep='')}))
    years <- paste0(years,'-01-01')

    dates <- as.Date(switch(period,annual=years,seasonal=seas,monthly=mons))
    nrow <- length(dates)

    series.matrix <- matrix(NA,ncol=length(model.list),nrow=nrow)

    for (m in seq_along(model.list)) {

       model <- model.list[m]
       print(model)
       model.info <- strsplit(model,'_')[[1]]
       gcm.list[m] <- model.info[1]
       run.list[m] <- model.info[3]

       model.dir <- paste0(read.dir,model,'/')

       area.info <- load_area_info(region,model.info[1],var.name,scenario,area.dir)

       region.avgs <- ts_function(var.name=var.name,period,
                                       read.dir=model.dir,tmp.dir=tmp.dir,
                                       clip.shp=clip.shp,area.info=area.info,
                                       past.int=past.int)

       date.ix <-  dates %in% as.Date(format(region.avgs$dates,'%Y-%m-%d'))

       series.matrix[date.ix,m] <- round(region.avgs$anoms,2)

    }

    write.file <- paste0(var.name,'_',region,'_time_series_anomalies_',cmip,'_',scenario,'_1850-2100.csv')
    header <- c('Dates',model.list)
    write.table(rbind(header,cbind(as.character(dates),series.matrix)),
                file=paste0(write.dir,write.file),
                quote=F,row.name=F,col.name=F,sep=',')

}

##---------------------------------------------------------------------
##*********************************************************************

testing <- TRUE

if (testing) {

  tmpdir <- "/local_temp/ssobie"
  cmip <- "CMIP5"
  scenario <- "rcp26"
  ##regions <-  'canada_boundary'
  regions <- c('bc','prairies','ontario','quebec',
               'northern_canada','atlantic_canada')
} else {
   args <- commandArgs(trailingOnly=TRUE)
   for(i in 1:length(args)){
       eval(parse(text=args[[i]]))
   }

}

for (region in regions) {

tmp.dir <- paste0(tmpdir,"/cmip_",region,"_",scenario,"_time_series/")

area.dir <- paste0("/storage/data/projects/rci/data/cas/",tolower(cmip),"/area_info/",region,"/")

if (!file.exists(tmp.dir)) {
   dir.create(tmp.dir,recursive=TRUE)
}

##------------------------------------
##Testing


period <- 'annual'
var.name <- 'tas'
past.int <- '1971-2000'

##------------------------------------

if (cmip=='CMIP6') {
   read.dir <- "/storage/data/climate/CMIP6/Derived/"
}
if (cmip=='CMIP5') {
   read.dir <- "/storage/data/climate/CMIP5/daily/Derived/"
}


shape.dir <- "/storage/data/projects/rci/data/assessments/shapefiles/canada_regions/"
clip.shp <- get_shape_file(region,shape.dir)

write.dir <- paste0("/storage/data/projects/rci/data/cas/",tolower(cmip),"/time_series/",region,"/")

if (!file.exists(write.dir))
  dir.create(write.dir,recursive=TRUE)



scenario.models <- list.files(read.dir,pattern=scenario)
###                   c('CESM1-CAM5','r1i1p1'),
cmip5.list <- list(c('ACCESS1-0','r1i1p1'),
                   c('bcc-csm1-1','r1i1p1'),
                   c('bcc-csm1-1-m','r1i1p1'),
                   c('BNU-ESM','r1i1p1'),
                   c('CanESM2','r1i1p1'),
                   c('CCSM4','r2i1p1'),
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

##cmip5.list <- list(c('HadGEM2-ES','r1i1p1'),
##                   c('bcc-csm1-1','r1i1p1'),
##                   c('bcc-csm1-1-m','r1i1p1'))


                   
if (scenario=='rcp26') {
   omx <- c(grep('ACCESS1-0',cmip5.list),
            grep('HadGEM2-CC',cmip5.list),
            grep('inmcm4',cmip5.list))
   cmip5.list <- cmip5.list[-omx]
}


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


cmip6.models <- c("ACCESS-CM2",paste0("ACCESS-ESM1-5_",scenario,"_r1i1p1f1"),
            "BCC-CSM2-MR",
            paste0("CanESM5_",scenario,"_r1i1p2f1"),
            "CNRM-CM6-1","CNRM-ESM2-1",
            "EC-Earth3_","EC-Earth3-Veg","FGOALS-g3","GFDL-ESM4","HadGEM3-GC31-LL",
            "INM-CM4-8","INM-CM5-0","IPSL-CM6A-LR","KACE-1-0-G","KIOST-ESM",
            "MIROC6","MIROC-ES2L","MPI-ESM1-2-HR","MPI-ESM1-2-LR","MRI-ESM2-0",
            "NESM3","NorESM2-LM","NorESM2-MM",
            paste0("UKESM1-0-LL_",scenario,"_r1i1p1f2"))


##model.list <- scenario.models[grepl(paste0(models,collapse='|'),scenario.models)]

model.dirs <- switch(cmip,
                     CMIP6=cmip6.list,
                     CMIP5=cmip5.list) 

model.list <- unlist(lapply(model.dirs,function(x,y){paste(x[1],y,x[2],sep='_')},scenario))

##Temperature and Precipitation
ts_function <- t_p_average_time_series
make_time_series(model.list,var.name,period,ts_function,
                 cmip,region,scenario,clip.shp,
                 read.dir,write.dir,area.dir,tmp.dir,past.int)


}