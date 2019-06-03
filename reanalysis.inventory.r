###Script to assemble a human readable inventory of Reanalysis
## data at PCIC

##------------------------------------------------------------
##Create table of ERA5 details

library(ncdf4)

create_era5_inventory <- function(files,read.dir,write.dir) {

   header <- c('Variable','Netcdf Variable','Standard Name','Frequency','Start','End')

   file.matrix <- matrix(NA,nrow=length(files),ncol=length(header))
   
   for (i in seq_along(files)) {
      file <- files[i]
      print(file)
      fsplit <- strsplit(file,'_')[[1]]
      vend <- grep('ERA5',fsplit)
      file.matrix[i,1] <- paste(fsplit[1:(vend-2)],collapse='_') ##Variable Name
      nc <- nc_open(paste0(read.dir,file))
      var.name <- names(nc$var)
      file.matrix[i,2] <- var.name
      long.name <- nc$var[[var.name]]$longname
      file.matrix[i,3] <- long.name
      file.matrix[i,4] <- fsplit[vend-1] ##Frequency
      years <- unlist(regmatches(file,gregexpr('[0-9]{8}-[0-9]{8}',file)))
      file.matrix[i,5] <- substr(years,1,4)
      file.matrix[i,6] <- substr(years,10,13)            
      nc_close(nc)
   }

   write.table(rbind(header,file.matrix),file=paste0(write.dir,'ERA5_reanalysis_inventory.csv'),quote=FALSE,row.name=F,col.name=F,sep=',')
}

##------------------------------------------------------------
##ERA5 downloaded data
if (1==0) {
era5.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/ERA5/concat/'
write.dir <- '/storage/data/projects/rci/stat.downscaling/inventories/'

hourly.files <- list.files(path=era5.dir,pattern='hour_ERA5')
daily.files <- list.files(path=era5.dir,pattern='day_ERA5')

files <- sort(c(hourly.files,daily.files))
create_era5_inventory(files,era5.dir,write.dir)
}



##------------------------------------------------------------
##ERA-Interim downloaded data
