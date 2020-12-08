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
##**********************************************************
tmp.dir <- '/local_temp/ssobie/cmip5_rankings/'
if (!file.exists(tmp.dir)) {
   dir.create(tmp.dir,recursive=TRUE)
}

save.dir <- "/storage/data/projects/rci/data/cas/cmip5/trends/"

shape.dir <- '/storage/data/projects/rci/data/assessments/shapefiles/bc_common/'
region <- 'canada_boundary'
area.shape <- get_shape_file(region,shape.dir)



##gcm <- 'CanESM5'
##ssp <- 'ssp245'
##run <- 'r2i1p2f1'

##files <- list.files(path="/storage/data/climate/CMIP5/daily/Derived/")
files <- list.files(path="/storage/data/climate/CMIP6/Derived/")
##Omit the P1 CanESM5
files <- files[-grep('CanESM5_ssp245_r[0-9]i1p1f1',files)]

##Only SSP245
files <- files[grep('ssp245',files)]

##SKip over CNRM-CM6-1-HR for now
##files <- files[58:60]
##files <- files[90:length(files)]

pr.vars <- list(annual="pr_annual",annual_extremes="pr_annual",
                climdex=c("cddETCCDI_ann","cwdETCCDI_ann","prcptotETCCDI_ann",
                          "r10mmETCCDI_ann","r1mmETCCDI_ann",
                          "r20mmETCCDI_ann","r95pETCCDI_ann","r99pETCCDI_ann",
                          "sdiiETCCDI_ann"))
tas.vars <- list(annual=c("tasmax_annual","tasmin_annual"),annual_extremes=c("tasmax_annual","tasmin_annual"),
                 climdex=c("csdiETCCDI_ann","fdETCCDI_ann","gslETCCDI_ann","idETCCDI_ann",
                           "suETCCDI_ann","trETCCDI_ann","wsdiETCCDI_ann"))

pr.matrix <- matrix(NA,nrow=length(files),ncol=length(unlist(pr.vars))+3)
tas.matrix <- matrix(NA,nrow=length(files),ncol=length(unlist(tas.vars)))

f.i <- 1
v.i <- 1

for (file in files) {
   print(file)
   file.split <- strsplit(file,'_')[[1]]
   gcm <- file.split[1]
   ssp <- file.split[2]
   run <- file.split[3]
   
   if (any(grepl('(gr1|gr2)',file.split))) {
      res <- file.split[4]
      run <- paste0(run,'_',res)
   }

   pr.matrix[f.i,1] <- gcm
   pr.matrix[f.i,2] <- ssp  
   pr.matrix[f.i,3] <- run
   v.i <- 4
   for (type in names(pr.vars)) {

      if (any(grepl('(gr1|gr2)',file.split))) {
         res <- file.split[4]
         gcm.dir <- paste0("/storage/data/climate/CMIP6/Derived/",gcm,"_",ssp,"_",run,"_",res,"/",type,"/")
      } else {
         gcm.dir <- paste0("/storage/data/climate/CMIP6/Derived/",gcm,"_",ssp,"_",run,"/",type,"/")
         ###gcm.dir <- paste0("/storage/data/climate/CMIP5/daily/Derived/",gcm,"_",ssp,"_",run,"/annual/")
      }
      var.list <- pr.vars[[type]]
      for (var.name in var.list) {
        print(paste0(type," ",var.name))
        var.file <- list.files(path=gcm.dir,pattern=var.name)

        if (length(var.file) > 1) {
           print(var.file)
           print('Too many GCM files selected')
        } else if (length(var.file) == 0) {
           print('Missing one of the files selected')
        } else {
 
            file.copy(from=paste0(gcm.dir,var.file),to=tmp.dir,overwrite=TRUE)
            ##Remap to common grid
            one.deg.file <- paste0(var.name,"_",gcm,"_",ssp,"_",run,"_one_degree.nc")
            regrid_remapcon(var.file,one.deg.file,tmp.dir)
  
            file.brick <- brick(paste0(tmp.dir,one.deg.file))
            time <- file.brick@z$Date
            yr.dates <- format(time,'%Y')

            area.cover <- area_weighting(file.brick,area.shape)
            year.series <- area_average_time_series(file.brick,area.cover)
 
            yrs <- as.numeric(yr.dates)
            base.ix <- yrs >= 1970 & yrs <= 2010
            past.ix <- yrs >= 1981 & yrs <= 2010
            proj.ix <- yrs >= 2071 & yrs <= 2100

            ##Projected change
            anomaly <- round(mean(year.series[proj.ix],na.rm=T) - mean(year.series[past.ix],na.rm=T),2)
            percent <- round((mean(year.series[proj.ix],na.rm=T) - mean(year.series[past.ix],na.rm=T)) / 
                              mean(year.series[past.ix],na.rm=T) *100,2)

            proj.info <- list(anom=anomaly,
                              prct=percent,
                              years=yrs,
                              series=year.series)
            print(paste0("Elements: ",f.i,", ",v.i))
            pr.matrix[f.i,v.i] <- anomaly
            v.i <- v.i + 1
            proj.file <- paste0(save.dir,var.name,"_",gcm,"_",ssp,"_",run,"_pr_info_2071-2100_1981-2010.RData")

            ##save(proj.info,file=proj.file)
            file.remove(paste0(tmp.dir,var.file))
            file.remove(paste0(tmp.dir,one.deg.file))           
         }
      }
   }
   f.i <- f.i + 1
}
 


