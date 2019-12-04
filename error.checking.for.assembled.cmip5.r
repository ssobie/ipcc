##Script to verify the assembled GCM files match the downloaded files

library(ncdf4)
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

##-----------------------------------------------------------------------------
##Find the matching download files that correspond to the assembled file

find_downloaded_files <- function(gcm,centre,esgf.dir) {

   gcm.dir <- paste0(esgf.dir,centre,'/',gcm)
   gcm.files <- list.files(path=gcm.dir,
                           pattern=paste0(varname,'_day_',gcm),
                           recursive=TRUE,full.name=TRUE)
   return(gcm.files)
}

##-----------------------------------------------------------------------------
##Find the downloaded files that correspond to the concatenated file

find_matching_download_files <- function(file,download.files) {

   file.split <- strsplit(file,'_')[[1]]
   run <- file.split[grepl('r[0-9]{1,2}i1p[0-9]{1}',file.split)]
   scenario <- file.split[grepl('historical',file.split)]
   future <- grepl('rcp',scenario)
   
   if (future) {
      rcp <- strsplit(scenario,'\\+')[[1]][2]
      hist.files <- download.files[grepl(paste0('historical_',run),download.files)]
      proj.files <- download.files[grepl(paste0(rcp,'_',run),download.files)]
      match.files <- sort(c(hist.files,proj.files))
   } else {
      hist.files <- download.files[grepl(paste0(scenario,'_',run),download.files)]
      match.files <- sort(hist.files)
   }
   return(match.files)
}

##-----------------------------------------------------------------------------
##Check the time series
check_time_series_from_file <- function(files,dir) {

   time.matrix <- matrix(NA,nrow=length(files),ncol=3)
   colnames(time.matrix) <- c('File','Monotonically Increasing?','No NA Values?')
   for (i in seq_along(files)) {
      file <- files[i]
##      print(file)
      time.matrix[i,1] <- file
      nc <- nc_open(paste0(dir,file))
      time.series <- netcdf.calendar(nc)
      nc_close(nc)
      diff.range <- unique(range(diff(time.series)))
      diff.test <- diff.range == 1 
      if (length(diff.range) != 1) {diff.test <- FALSE}
      time.matrix[i,2] <- diff.test     
      time.na <- is.na(time.series)
      time.matrix[i,3] <- sum(time.na)==0      
   }
   return(time.matrix)
}

##-----------------------------------------------------------------------------
#Dimension Check - confirm the spatial dimensions and coordinates are the same

check_file_spatial_information <- function(files,download.files,dir) {

   space.matrix <- matrix(NA,nrow=length(files),ncol=3)
   colnames(space.matrix) <- c('File','Same Longitudes?','Same Latitudes?')

   for (i in seq_along(files)) {
      file <- files[i]
##      print(file)

      match.files <- find_matching_download_files(file,download.files)
##      print(match.files)
      space.matrix[i,1] <- file
      nc <- nc_open(paste0(dir,file))
      cat.lon <- ncvar_get(nc,'lon')
      cat.lat <- ncvar_get(nc,'lat')
      nc_close(nc)
      
      r.len <- length(match.files)
      ref.lons <- rep(NA,r.len)
      ref.lats <- rep(NA,r.len)
      for (j in 1:r.len) {
         ref.nc <- nc_open(match.files[j])
         ref.lons[j] <- sum(cat.lon - ncvar_get(ref.nc,'lon')) == 0
         ref.lats[j] <- sum(cat.lat - ncvar_get(ref.nc,'lat')) == 0
         nc_close(ref.nc)
      }
      space.matrix[i,2] <- sum(ref.lons == 0) == 0
      space.matrix[i,3] <- sum(ref.lats == 0) == 0
   }
   return(space.matrix)
}



##-----------------------------------------------------------------------------
##Random cell check - extract ~10 cells at random and confirm they are the 
##same in the downloaded and assembled files

check_random_cells <- function(var.name,files,download.files,dir) {
   t.len <- 5 
   cell.matrix <- matrix(NA,nrow=length(files),ncol=t.len)

   for (i in seq_along(files)) {
      file <- files[i]
      ##print(paste0('File number ',i))
      match.files <- find_matching_download_files(file,download.files)
      cell.matrix[i,1] <- file
      nc <- nc_open(paste0(dir,file))
      cat.lon <- ncvar_get(nc,'lon')
      cat.lat <- ncvar_get(nc,'lat')
      cat.time <- as.Date(as.character(netcdf.calendar(nc)))
      cat.one <- ncvar_get(nc,var.name,start=c(1,1,1000),count=c(-1,-1,1))
      
      rand.na <- TRUE
      while (rand.na) {
         lon.ix <- sample(1:length(cat.lon),t.len)
         lat.ix <- sample(1:length(cat.lat),t.len)  
         cat.test <- cat.one[lon.ix,lat.ix]
         rand.na <- any(is.na(cat.test))
      }
      ##print('Random cells selected')
      ##print(paste0('Lons ',lon.ix))
      ##print(paste0('Lats ',lat.ix))

      cat.series <- matrix(NA,nrow=t.len,ncol=length(cat.time))
      for (k in 1:t.len) {
         cat.series[k,] <- ncvar_get(nc,var.name,start=c(lon.ix[k],lat.ix[k],1),count=c(1,1,-1))
      }      
      nc_close(nc)
      ##print('Extracted cells from concatenated series')      
      r.len <- length(match.files)
      
      for (k in 1:t.len) {
        ref.values <- c() 
        ref.time <- as.Date('1000-01-01')
        for (j in 1:r.len) {
            ref.nc <- nc_open(match.files[j])
            ref.time <- c(ref.time,as.Date(as.character(netcdf.calendar(ref.nc))))
            ref.values <- c(ref.values,ncvar_get(ref.nc,var.name,start=c(lon.ix[k],lat.ix[k],1),count=c(1,1,-1)))
            nc_close(ref.nc)
         }
         ref.time <- ref.time[-1]
         flag <- duplicated(ref.time)
         ref.values <- ref.values[!flag]

         cell.matrix[i,k] <- sum(cat.series[k,] - ref.values) == 0
         if (sum(cat.series[k,] - ref.values) != 0)
            browser()
      }
      print(paste0('Done with comparison for file ',i))
   }
   return(cell.matrix)
}



##-----------------------------------------------------------------------------
##*****************************************************************************

##args <- commandArgs(trailingOnly=TRUE)
##for(i in 1:length(args)){
##    eval(parse(text=args[[i]]))
##}
varname <- 'tasmin'
gcm <- 'inmcm4'
centre <- 'INM'

base.dir <- '/storage/data/climate/CMIP5/daily/'
esgf.dir <- '/storage/data/climate/CMIP5/incoming/output1/'
gcm.dir <- paste0(base.dir,gcm,'/')
var.files <- list.files(path=gcm.dir,pattern=varname)

download.files <- find_downloaded_files(gcm,centre,esgf.dir)

time.check <- check_time_series_from_file(var.files,gcm.dir)
print('Time check')
print(time.check)

space.check <- check_file_spatial_information(var.files,download.files,gcm.dir)
print('Space check')
print(space.check)

cell.check <- check_random_cells(varname,var.files,download.files,gcm.dir)
print('Random cell check')
print(cell.check)
  
##Create a formatted xlsx file with the checking results