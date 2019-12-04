##Clip the model simulations to north america and remove the unlimited time dimension

library(ncdf4)
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

get_dates_from_file <- function(file) {
   nc <- nc_open(file)
   time.range <- as.Date(as.character(range(netcdf.calendar(nc))))
   nc_close(nc)
   return(time.range)
}

subset_dates_by_ncks <- function(input.file,tmp.file) {

   nc <- nc_open(input.file)
   times <- as.Date(as.character(netcdf.calendar(nc)))
   dst <- head(which(times >= as.Date('1900-01-01')),1)-1
   den <- tail(which(times <= as.Date('2100-12-31')),1)-1

   work <- paste0('ncks -d time,',dst,',',den,' ',input.file,' ',tmp.file)
   print(work)
   system(work)
   new.dates <- get_dates_from_file(tmp.file)
   output.file <- gsub(pattern='[0-9]{8}-[0-9]{8}',
                      replacement=paste0(format(new.dates[1],'%Y%m%d'),'-',format(new.dates[2],'%Y%m%d')),input.file)
   file.copy(from=tmp.file,to=output.file,overwrite=TRUE) 
   return(output.file)
}

subset_dates <- function(input.file,tmp.file) {
 
  work <- paste0('cdo -O seldate,1900-01-01T00:00:00,2100-12-31T23:59:59 ',input.file,' ',tmp.file)
  print(work)

  system(work)
  new.dates <- get_dates_from_file(tmp.file)
  output.file <- gsub(pattern='[0-9]{8}-[0-9]{8}',
                    replacement=paste0(format(new.dates[1],'%Y%m%d'),'-',format(new.dates[2],'%Y%m%d')),input.file)
  file.copy(from=tmp.file,to=output.file,overwrite=TRUE) 
  return(output.file)

}

subset_to_north_america <- function(input.file,output.file) {
  work <- paste('ncks -O -d lon,180.,310. -d lat,5.,90. ',input.file,' ',output.file,sep='')
  system(work)
}

remove_unlimited_time_dim <- function(input.file,output.file) {
   work <- paste0('nccopy -u ',input.file,' ',output.file)
   print(work)
   system(work)
}


##*************************************************************************

##args <- commandArgs(trailingOnly=TRUE)
##for(i in 1:length(args)){
##    eval(parse(text=args[[i]]))
##}


varname <- 'tasmin'
scenario <- 'ssp126'
gcm <- 'CanESM5'

tmpdir <- '/local_temp/ssobie/cmip6/'
tmp.dir <- paste0(tmpdir,varname,'/') ##'/local_temp/ssobie/cmip5/' ##tmpdir

if (!file.exists(tmp.dir)) {
  dir.create(tmp.dir,recursive=TRUE)
}

read.dir <- '/storage/data/climate/CMIP6/assembled/'
write.dir <- paste0(read.dir,gcm,'/north_america/')

file.pattern <- paste0(varname,'_day_',gcm,'_historical\\+',scenario)
gcm.files <- list.files(path=paste0(read.dir,gcm,'/'),pattern=file.pattern)

for (g in seq_along(gcm.files)) {
   file <- gcm.files[g]
   print(file)
   file.copy(from=paste0(read.dir,gcm,'/',file),to=tmp.dir) ##,overwrite=T)
   print('Copied file')
   sub.tmp <- paste0('SUBSET_',file)
   sub.file <- gsub(gcm,paste0(gcm,'_North_America'),file)      
   subset_to_north_america(paste0(tmp.dir,file),paste0(tmp.dir,sub.tmp))
   print('Subset to North America')
   remove_unlimited_time_dim(paste0(tmp.dir,sub.tmp),paste0(tmp.dir,sub.file))   
   print('Removed Unlimited Time Dimension')
   print('Subset by time')
   ###new.file <- subset_dates(paste0(tmp.dir,sub.file),paste0(tmp.dir,'time.tmp.nc'))   
   new.file <- subset_dates_by_ncks(paste0(tmp.dir,sub.file),paste0(tmp.dir,'time.tmp.nc'))   

   print('Done subsetting by time')
   print('Copying back')
   file.copy(from=new.file,to=write.dir,overwrite=TRUE)
   print('Done copying back')
   file.remove(paste0(tmp.dir,sub.tmp))
   file.remove(paste0(tmp.dir,file))
   file.remove(paste0(tmp.dir,'time.tmp.nc'))
   file.remove(new.file)
}
   
   