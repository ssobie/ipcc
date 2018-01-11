###Script to assemble a human readable inventory of CMIP5 data at PCIC

slice.by.time <- function(space.file,time.file,era) {
  if (era) {
    time.write <- gsub(pattern='[0-9]{8}-[0-9]{8}',replacement='19000101-20161231',time.file)
    work <- paste('cdo -O seldate,1900-01-01T00:00,2016-12-31T23:59 ',space.file,' ',time.write,sep='')
  } else {
    time.write <- gsub(pattern='[0-9]{8}-[0-9]{8}',replacement='20060101-21001231',time.file)
    work <- paste('cdo -O seldate,2006-01-01T00:00,2100-12-31T23:59 ',space.file,' ',time.write,sep='')
  }
  system(work)
}

make.new.files <- function(nc.files) {

      ist <- regexpr('[0-9]{8}-',head(nc.files,1))
      zst <- ist + attr(ist, "match.length")-2      
      ien <- regexpr('-[0-9]{8}',tail(nc.files,1))      
      zen <- ien  + attr(ien, "match.length")-1 

      yst <- substr(head(nc.files,1),ist,zst)
      yen <- substr(tail(nc.files,1),ien+1,zen)
      new.file <- gsub('[0-9]{8}-[0-9]{8}',paste(yst,'-',yen,sep=''),head(nc.files,1))
      return(new.file)
}

group.files <- function(file.list,write.dir,base.dir) {

   ##All scenarios including historical
     flen <- length(file.list)
     if (flen > 1) {
        file.paste <- paste(paste(base.dir,file.list,sep=''),collapse=' ')
        new.file <- make.new.files(file.list)
        new.full <- paste(write.dir,new.file,sep='')
        work <- paste('ncrcat -O ',file.paste,' ',new.full,sep='')      
        print(work)
        system(work)
     } else if (flen==1) {
        print('File exists')
        work <- paste('cp ',base.dir,file.list,' ',write.dir,file.list,sep='')
        print(work)
        system(work)
     } else if (flen==0) {
        print(paste0('No files to concatenate for ',scenario))

     }

}


base.dir <- '/storage/data/climate/downscale/CMIP5/incoming/'

##Make sure the GCM and Centre names are paired in the same order below
##gcms <- c('bcc-csm1-1','CanESM2','CCSM4','CESM1-CAM5','CNRM-CM5','CSIRO-Mk3-6-0',
##          'GFDL-CM3','GFDL-ESM2M','HadGEM2-ES')
##gcms <-  c('inmcm4','IPSL-CM5A-LR','IPSL-CM5A-MR','MIROC5','MIROC-ESM','MIROC-ESM-CHEM','MPI-ESM-LR','MRI-CGCM3','NorESM1-M')

gcms <- c('bcc-csm1-1','CanESM2','CSIRO-Mk3-6-0', 'GFDL-CM3','HadGEM2-ES','IPSL-CM5A-MR','MIROC-ESM','MPI-ESM-MR',
              'ACCESS1-3',  'bcc-csm1-1-m','CCSM4','GFDL-ESM2G','HadGEM2-AO', 'inmcm4','IPSL-CM5B-LR','MIROC-ESM-CHEM','MRI-CGCM3',
              'BNU-ESM','CNRM-CM5','GFDL-ESM2M','HadGEM2-CC','IPSL-CM5A-LR','MIROC5','MPI-ESM-LR')

gcms <- 'ACCESS1-0'
variable.list <- 'sic' ##c('pr','tasmax','tasmin')


for (gcm in gcms) {
  ##write.dir <- paste0('/storage/data/climate/downscale/BCCAQ2/CMIP5/global/control/',gcm,'/')
  write.dir <- paste0('/storage/data/projects/CanSISE/n2w/CMIP5/',gcm,'/')
  if (!file.exists(write.dir)) {
     dir.create(write.dir,recursive=T)
  }

  for (variable in variable.list) {
    ##PR files
    pr.files <- list.files(path=paste(base.dir,gcm,'/download',sep=''),pattern=paste(variable,'_OImon',sep=''),recursive=TRUE)
    print(pr.files)
    plen <- length(pr.files)
    if (plen!=0) {
       ##print(pr.files)      
       gcm.dir <- paste0(base.dir,gcm,'/download/')
       group.files(pr.files,write.dir,gcm.dir)    
    }##If statement
  }##Variable Loop
}##GCM Loop

