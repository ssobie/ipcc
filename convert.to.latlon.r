##Script to convert the gcms to the same grid as observations.

##******************************************************

tmp.dir <- '/local_temp/ssobie/convert/'

if (!file.exists(tmp.dir)) {
  dir.create(tmp.dir,recursive=TRUE)
}


gcm <- 'NorESM1-M'

##------------------------------------

base.dir <- '/storage/data/projects/CanSISE/n2w/CMIP5/monthly/'
gcm.dir <- paste(base.dir,gcm,sep='')
all.files <- list.files(path=gcm.dir,pattern=paste0('sic_count_',gcm))
scen.files <- all.files[grep('r1i1p1',all.files)]
sic.files <- scen.files[grep('1850',scen.files)]  
files <- sic.files

for (file in files) {
    print(file)
  stereo.file <- gsub(pattern=paste0(varname,'_count'),replacement=paste0(varname,'_stereo'),file)
  move.over <- paste0('mv ',gcm.dir,'/',file,' ',gcm.dir,'/',stereo.file)
  system(move.over)
          
  move.to <- paste("rsync -av --progress ",gcm.dir,'/',stereo.file, " " ,tmp.dir,sep='')
  print(move.to)
  system(move.to)

  data.file <- paste0(tmp.dir,stereo.file)
  tiff.file <- gsub(pattern=".nc",replacement='.tiff',data.file)
  count.file <- gsub(pattern=paste0(varname,'_stereo'),replacement=paste0(varname,'_count'),data.file)

  print('TIFF')
  work <- paste0("gdalwarp -t_srs 'EPSG:4326 '  ",data.file," ",tiff.file)
  system(work)

  print('Translate')
  work <- paste0(" gdal_translate -of netcdf ",tiff.file," ",count.file)
  system(work)

      
  move.back <- paste("rsync -av --progress ",count.file," ",gcm.dir,sep='')
  print(move.back)
  system(move.back)

      clean.up <- paste0('rm ',tiff.file)
      system(clean.up)
      clean.up <- paste0('rm ',count.file)
      system(clean.up)
      clean.up <- paste0('rm ',data.file)
      system(clean.up)

}##Scenarios

