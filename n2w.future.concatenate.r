###Script to assemble a human readable inventory of CMIP5 data at PCIC

slice.by.time <- function(space.file,time.file,past) {
  if (past) {
    time.write <- gsub(pattern='[0-9]{8}-[0-9]{8}',replacement='18500101-20051231',time.file)
    work <- paste('cdo -O seldate,1850-01-01T00:00,2005-12-31T23:59 ',space.file,' ',time.write,sep='')
  } else {
    time.write <- gsub(pattern='[0-9]{8}-[0-9]{8}',replacement='20060101-21001231',time.file)
    work <- paste('cdo -O seldate,2006-01-01T00:00,2100-12-31T23:59 ',space.file,' ',time.write,sep='')
  }

##  if (past) {
##    time.write <- gsub(pattern='[0-9]{6}-[0-9]{6}',replacement='185001-200612',time.file)
##    work <- paste('cdo -O seldate,1850-01-01,2005-12-31 ',space.file,' ',time.write,sep='')
##  } else {
##    time.write <- gsub(pattern='[0-9]{6}-[0-9]{6}',replacement='200601-210012',time.file)
##    work <- paste('cdo -O seldate,2006-01-01,2100-12-31 ',space.file,' ',time.write,sep='')
##  }

  print(work)
  system(work)
}

make.subsets <- function(tmp.dir) {

  files <- list.files(path=paste(tmp.dir,'grouptmp',sep=''))
  len <- length(files)
  for (i in 1:len) {
    file <- files[i]
    past <- grepl('_historical_',file)
    input.file <- paste(tmp.dir,'grouptmp/',file,sep='')
    time.file <- paste(tmp.dir,'timetmp/time_subset_',file,sep='')
    ##time.file <- paste(tmp.dir,'timetmp/',file,sep='')
    slice.by.time(input.file,time.file,past)
  }
}



concat.files <- function(tmp.dir,gcm) {

  time.dir <- paste(tmp.dir,'timetmp/',sep='')
  
  file.list <- list.files(path=time.dir,pattern='*nc')
  hist.file <- file.list[grep('_historical_',file.list)]
  rcp.files <- file.list[grep('rcp',file.list)]
  
  for (i in seq_along(rcp.files)) {
    gcm.dir <- paste(tmp.dir,gcm,'/',sep='')
    if (!file.exists(gcm.dir))
      dir.create(gcm.dir,recursive=TRUE)
    rcp.file <- rcp.files[i]
    rst <- regexpr(pattern='rcp[0-9]{2}',rcp.file)
    ren <- rst + attr(rst, "match.length")-1      
    rcp <- substr(rcp.file,rst,ren)

    ist <- regexpr('[0-9]{8}-',hist.file)
    zst <- ist[1] + attr(ist, "match.length")-2
    ien <- regexpr('-[0-9]{8}',tail(rcp.files,1))
    zen <- ien[1]  + attr(ien, "match.length")-1

    yst <- substr(head(file.list,1),ist,zst)
    yen <- substr(tail(file.list,1),ien+1,zen)
   
    cat.file <- gsub(pattern='[0-9]{8}-[0-9]{8}',replacement=paste0(yst,'-',yen),hist.file)      
    cat.file <- gsub(pattern='historical',replacement=paste('historical+',rcp,sep=''),cat.file)
    cat.file <- gsub(pattern='time_subset_',replacement='',cat.file)

    work <- paste('ionice -c 3 ncrcat -O ',time.dir,hist.file,' ',time.dir,rcp.file,' ',gcm.dir,cat.file,sep='')
    print(work)
    system(work)

  }
}
  

make.new.files <- function(nc.files,h.ix) {

      ist <- regexpr('[0-9]{8}-',head(nc.files[h.ix],1))
      zst <- ist + attr(ist, "match.length")-2      
      ien <- regexpr('-[0-9]{8}',tail(nc.files[h.ix],1))      
      zen <- ien  + attr(ien, "match.length")-1 

      yst <- substr(head(nc.files[h.ix],1),ist,zst)
      yen <- substr(tail(nc.files[h.ix],1),ien+1,zen)
      new.file <- gsub('[0-9]{8}-[0-9]{8}',paste(yst,'-',yen,sep=''),head(nc.files[h.ix],1))
      return(new.file)
}

group.files <- function(file.list,matrix.subset,tmp.dir,base.dir,scenarios) {

   ##All scenarios including historical
   for (scenario in scenarios) {
     nc.files <- unlist(file.list)
     h.ix <- grep(scenario,nc.files)
     
     if (1==0) {##(length(h.ix) > 1) {
        file.paste <- paste(paste(base.dir,'/',file.list[h.ix],sep=''),collapse=' ')
        new.file <- make.new.files(nc.files,h.ix)
        new.full <- paste(tmp.dir,'grouptmp/',new.file,sep='')
        work <- paste('ncrcat -O ',file.paste,' ',new.full,sep='')      
        print(work)
        system(work)
     } else {##if (length(h.ix) ==1) {
        print('File exists')
        work <- paste('cp ',base.dir,'/',file.list[h.ix],' ',tmp.dir,'grouptmp/',nc.files[h.ix],sep='')
        print(work)
        system(work)
     }
   }

}

##******************************************************

##args <- commandArgs(trailingOnly=TRUE)
##for(i in 1:length(args)){
##    eval(parse(text=args[[i]]))
##}

tmp.dir <- '/local_temp/ssobie/' ##tmpdir

if (!file.exists(tmp.dir)) {
  dir.create(tmp.dir,recursive=TRUE)
  dir.create(paste0(tmp.dir,'timetmp/'),recursive=TRUE)
  dir.create(paste0(tmp.dir,'grouptmp/'),recursive=TRUE)
}

##Clean tmpdir
##rm.tmp <- paste0("rm ",tmp.dir,"*nc")
##print(rm.tmp)
##system(rm.tmp)

##system('ls /local_temp/ssobie/*/*')
##system('rm /local_temp/ssobie/prism/interpolated/*')
##system('rm /local_temp/ssobie/prism/*')

varname <- 'tasmin'
##'ACCESS1-0','ACCESS1-3','bcc-csm1-1','bcc-csm1-1-m','CanESM2',
gcm.list <- c('CCSM4','CNRM-CM5',
                'CSIRO-Mk3-6-0','HadGEM2-AO','HadGEM2-CC','HadGEM2-ES','inmcm4','IPSL-CM5A-LR','IPSL-CM5A-MR',
                'MIROC5','MIROC-ESM','MIROC-ESM-CHEM','MPI-ESM-LR','MPI-ESM-MR','MRI-CGCM3','NorESM1-M')
##gcm.list <- c('MIROC-ESM','MIROC-ESM-CHEM','MPI-ESM-LR','MPI-ESM-MR','MRI-CGCM3','NorESM1-M')
##'CNRM-CM5','CSIRO-Mk3-6-0','HadGEM2-CC','HadGEM2-ES','MPI-ESM-MR')
gcm.list <- c('GFDL-CM3','GFDL-ESM2G','GFDL-ESM2M')
##
##gcm.list <- c('bcc-csm1-1-m','CanESM2','CCSM4','CNRM-CM5','CSIRO-Mk3-6-0','GISS-E2-H','GISS-E2-R',
##               'MIROC5','MIROC-ESM','MIROC-ESM-CHEM','MPI-ESM-LR','MPI-ESM-MR','MRI-CGCM3','NorESM1-M')

gcm.list <- c('NorESM1-ME')

for (gcm in gcm.list) {

if (1==1) {
  scen.list <- c('historical','rcp26','rcp45','rcp85') ##,'rcp85') ##'rcp26','rcp45','rcp85')

  ##read.dir <- paste0('/storage/data/projects/CanSISE/n2w/CMIP5/incoming/',gcm,'/download')
  ##write.dir <- '/storage/data/projects/CanSISE/n2w/CMIP5/monthly/'
  read.dir <- paste0('/storage/data/climate/downscale/CMIP5/incoming/',gcm,'/download')
  ##read.dir <- paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/CMIP5/global/',gcm,'/')
  write.dir <- paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/CMIP5/global/',gcm,'/')

  full.files <- list.files(path=read.dir,pattern=paste(varname,'_day_',gcm,sep=''))
  split.apart <- unlist(strsplit(unlist(full.files),'_'))
  runs <- unique(split.apart[grep('r*i1p*',split.apart)])

  move.to <- paste("rsync -av ",read.dir,'/',varname,'_day_',gcm,'*nc ',tmp.dir,sep='')
  print(move.to)
  system(move.to)

  ##PR files
  all.files <- as.list(list.files(path=tmp.dir,pattern=paste(varname,'_day',sep='')))
  
  for (run in runs) {

     pr.files <- all.files[grep(run,all.files)]
     print(unlist(pr.files))

     ##print(pr.files)
     fsplit <- lapply(pr.files,function(x) {strsplit(x,'_')[[1]]})
     file.matrix <- matrix(NA,nrow=length(pr.files),ncol=length(fsplit[[1]])+1)
     for (i in 1:length(pr.files)) {
        file.matrix[i,] <- c(fsplit[[i]][1:(length(fsplit[[1]])-1)],substr(fsplit[[i]][6],1,4),substr(fsplit[[i]][6],10,13))
        ##file.matrix[i,] <- c(fsplit[[i]][1:(length(fsplit[[1]])-1)],substr(fsplit[[i]][6],1,4),substr(fsplit[[i]][6],8,11))
     }

     group.files(pr.files,file.matrix,tmp.dir,tmp.dir,scen.list)    
     make.subsets(tmp.dir)
     concat.files(tmp.dir,gcm)
     print('Done concat')

    move.cat <- paste("rsync -av ",tmp.dir,gcm,"/ ",write.dir,sep='')
    ##move.cat <- paste("rsync -av ",tmp.dir,"timetmp/ ",write.dir,sep='')
    print(move.cat)
    system(move.cat)

    clean.group <- paste('rm ', tmp.dir,'grouptmp/*',sep='')
    system(clean.group)
    clean.time  <- paste('rm ', tmp.dir,'timetmp/*',sep='')
    system(clean.time)
}
    clean.in <- paste('rm ', tmp.dir,'*',sep='')
    system(clean.in)
    clean.in <- paste('rm ', tmp.dir,gcm,'/*nc',sep='')
    system(clean.in)
    clean.in <- paste('rmdir ', tmp.dir,gcm,'/',sep='')
    system(clean.in)

}

}

clean.data  <- paste('rm /local_temp/ssobie/*/*',sep='')
system(clean.data)
clean.data  <- paste('rm ', tmp.dir,'*/*nc',sep='')
system(clean.data)



system('ls /local_temp/ssobie/*/*')