##Script to assemble downloaded CMIP5 GCM files in Earth System Grid Format
##into continuous time files in the format actually used

library(ncdf4)
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

##-----------------------------------------------------------------------------

slice_by_time <- function(space.file,time.file,dst,den) {
  time.write <- gsub(pattern='[0-9]{8}-[0-9]{8}',
                     replacement=paste0(format(dst,'%Y%m%d'),'-',format(den,'%Y%m%d')),time.file)
  work <- paste0('cdo -O seldate,',dst,'T00:00:00,',den,'T23:59:59 ',space.file,' ',time.write)
  print(work)
  system(work)
}

##-----------------------------------------------------------------------------

get_dates_from_file <- function(file) {
   nc <- nc_open(file)
   time.range <- as.Date(as.character(range(netcdf.calendar(nc))))
   nc_close(nc)
   return(time.range)
}

make_subsets <- function(tmp.dir) {
  print('Making subsets')
  files <- list.files(path=paste0(tmp.dir,'grouptmp'),full.name=T)

  dates <- lapply(files,get_dates_from_file)
  
  hist.ix <- grep('_historical_',files)

  hist.dates <- dates[[hist.ix]]
  ##file.copy(from=files[[hist.ix]],
  ##          to=paste0(tmp.dir,'timetmp/time_subset_',basename(files[[hist.ix]])),
  ##          overwrite=TRUE)
  slice_by_time(files[[hist.ix]],
                paste0(tmp.dir,'timetmp/time_subset_',basename(files[[hist.ix]])),
                hist.dates[1],hist.dates[2])

  fst <- hist.dates[2]+1 ##Subset the future files one day after
                         ##the end of the past (want to maximize past time)
  fut.files <- files[-hist.ix]
  fut.dates <- dates[-hist.ix]
  len <- length(fut.files)
  for (i in 1:len) {
     file <- fut.files[i]
     time.file <- paste(tmp.dir,'timetmp/time_subset_',basename(file),sep='')
     slice_by_time(file,time.file,fst,fut.dates[[i]][2])
  }
}

##-----------------------------------------------------------------------------

concat_files <- function(tmp.dir,gcm) {
  print('Concatenating files')
  time.dir <- paste(tmp.dir,'timetmp/',sep='')
  
  file.list <- list.files(path=time.dir,pattern='*nc')
  hist.file <- file.list[grep('_historical_',file.list)]
  ssp.files <- file.list[grep('ssp',file.list)]

  for (i in seq_along(ssp.files)) {
    gcm.dir <- paste(tmp.dir,gcm,'/',sep='')
    if (!file.exists(gcm.dir))
      dir.create(gcm.dir,recursive=TRUE)
    ssp.file <- ssp.files[i]
    rst <- regexpr(pattern='ssp[0-9]{3}',ssp.file)
    ren <- rst + attr(rst, "match.length")-1      
    ssp <- substr(ssp.file,rst,ren)

    ist <- regexpr('[0-9]{8}-',hist.file)
    zst <- ist[1] + attr(ist, "match.length")-2
    ien <- regexpr('-[0-9]{8}',ssp.file)
    zen <- ien[1]  + attr(ien, "match.length")-1

    yst <- substr(hist.file,ist,zst)
    yen <- substr(ssp.file,ien+1,zen)
   
    cat.file <- gsub(pattern='[0-9]{8}-[0-9]{8}',replacement=paste0(yst,'-',yen),hist.file)      
    cat.file <- gsub(pattern='historical',replacement=paste('historical+',ssp,sep=''),cat.file)
    cat.file <- gsub(pattern='time_subset_',replacement='',cat.file)
    work <- paste('ionice -c 3 ncrcat -O ',time.dir,hist.file,' ',time.dir,ssp.file,' ',gcm.dir,cat.file,sep='')
    print(work)
    system(work)
  }
}
  
##-----------------------------------------------------------------------------

make_new_files <- function(files) {

      dates <- lapply(files,get_dates_from_file)
      yst <- format(do.call(min,dates),'%Y%m%d')
      yen <- format(do.call(max,dates),'%Y%m%d')
      new.file <- gsub('[0-9]{8}-[0-9]{8}',paste(yst,'-',yen,sep=''),head(basename(files),1))
      return(new.file)
}

##-----------------------------------------------------------------------------
check_for_duplicate_files <- function(files) {

    dates <- lapply(files,get_dates_from_file)
    series <- unlist(lapply(dates,function(x){seq(x[1],by='day',x[2])}))

    if (any((range(diff(series)) < 1 | range(diff(series)) > 10 ))) {
       print(files)
       print(range(diff(series)))
       print(char.dates[which(diff(series) < 1)])
       print(char.dates[which(diff(series) > 10)])

       ##browser()
       stop('Duplicate years or missing files')
    }
}

##-----------------------------------------------------------------------------

group_files <- function(files,tmp.dir,base.dir,scenarios) {
   used.files <- c()
   ##All scenarios including historical
   for (scenario in scenarios) {
     h.ix <- grep(paste0('_',scenario,'_'),files)    

     if (length(h.ix) > 1) {
        ##Check for duplicates
        check_for_duplicate_files(files[h.ix])       
        ##file.copy(from=files[h.ix],to=paste0(tmp.dir,'grouptmp/'),overwrite=TRUE)
        file.paste <- paste(files[h.ix],collapse=' ')
        ##file.paste <- paste(paste(base.dir,'grouptmp/',basename(files[h.ix]),sep=''),collapse=' ')
        new.file <- make_new_files(files[h.ix])
        new.full <- paste0(tmp.dir,'grouptmp/',new.file)
        work <- paste('ncrcat -O ',file.paste,' ',new.full,sep='')      
        print(work)
        system(work)
        used.files <- c(used.files,files[h.ix])
     } else if (length(h.ix) == 1) {
        print('One file exists')
        file.copy(from=files[h.ix],to=paste0(tmp.dir,'grouptmp/'),overwrite=TRUE)
        used.files <- c(used.files,files[h.ix])
     } else {
        print(paste0('No files for: ',scenario))
     }
   }
   return(used.files)
}

##-----------------------------------------------------------------------------
##*****************************************************************************

testing <- FALSE

if (testing) {
   varname <- 'tasmin'
   centre <- 'NCC'
   scenario <- 'ssp585'
   tmpdir <- '/local_temp/ssobie/assembly'
} else {
   args <- commandArgs(trailingOnly=TRUE)
   for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
   }
}


tmp.dir <- paste0(tmpdir,'/',centre,'/',varname,'/',scenario,'/') ##'/local_temp/ssobie/cmip5/' ##
 
if (!file.exists(tmp.dir)) {
  dir.create(tmp.dir,recursive=TRUE)
  dir.create(paste0(tmp.dir,'timetmp/'),recursive=TRUE)
  dir.create(paste0(tmp.dir,'grouptmp/'),recursive=TRUE)
}

scen.list <- c('historical',scenario) ###'ssp245','ssp585') 

base.dir <- '/storage/data/climate/CMIP6/incoming/CMIP6/'
write.dir <- '/storage/data/climate/CMIP6/assembled/'

##----------------------------------------------------------------------
##Find the GCM Centre information

centre.info <- list()

   centre.dir <- paste0(base.dir,'ScenarioMIP/',centre,'/')
   centre.gcms <- gsub('./','',list.files(path=centre.dir,recursive=FALSE,full.names=FALSE))

   glen <- length(centre.gcms)
   gcm.list <- list()
   for (gcm in centre.gcms) {
  
      gcm.dir <- paste0(base.dir,'ScenarioMIP/',centre,'/',gcm)
      gcm.ssp.files <- list.files(path=gcm.dir,
                              pattern=paste0(varname,'_day_',gcm),
                              recursive=TRUE,full.name=TRUE)    

      ##Runs         
      ist <- regexpr('r[0-9]{1,2}i1p[0-9]{1}f[0-9]{1}',basename(gcm.ssp.files))
      zst <- ist + attr(ist, "match.length")-1
      runs <- mapply(function(x,s,e){substr(x,s,e)},as.list(basename(gcm.ssp.files)),as.list(ist),as.list(zst))         

      ##Grid Label - n for native and r for regridded. If more than one, append an integer less than 10
      ##distinguish between resolutions.
      ist.g <- regexpr('(gn[0-9]*|gr[0-9]*)',basename(gcm.ssp.files))
      zst.g <- ist.g + attr(ist.g, "match.length")-1
      grids <- mapply(function(x,s,e){substr(x,s,e)},as.list(basename(gcm.ssp.files)),as.list(ist.g),as.list(zst.g))         

      gcm.hist.dir <- paste0(base.dir,'CMIP/',centre,'/',gcm)
      gcm.hist.files <- list.files(path=gcm.hist.dir,
                                   pattern=paste0(varname,'_day_',gcm,'_historical'),
                                   recursive=TRUE,full.name=TRUE)
      gcm.files <- c(gcm.hist.files,gcm.ssp.files)

      rv <- list(list(files=gcm.files,runs=unique(runs),grids=unique(grids)))       

      names(rv) <- gcm
      gcm.list <- append(gcm.list,rv)      
   }
   add.list <- list(gcm.list)
   names(add.list) <- centre
   centre.info <- append(centre.info,add.list)

   
##----------------------------------------------------------------------

   gcm.list <- centre.info[[centre]]
   centre.gcms <-  "KIOST-ESM" ######names(gcm.list) ##

   for (gcm in centre.gcms) {
      print(gcm)
      gcm.info <- gcm.list[[gcm]]     
      gcm.info$runs <- 'r1i1p1f1'

      used.files <- c()      
      for (run in gcm.info$runs) {
         print(run)
         run.grid.files <- gcm.info$files[grep(run,gcm.info$files)]
         for (grid in gcm.info$grids) {
            print(grid)
            run.files <- run.grid.files[grep(grid,run.grid.files)]

            if (length(run.files) > 1) {
               used.run.files <- group_files(run.files,tmp.dir,tmp.dir,scen.list)    
               used.files <- c(used.files,used.run.files)
               grouped.files <- list.files(path=paste0(tmp.dir,'grouptmp'),full.name=T)

               if (length(grouped.files) > 1) {
                  print('Making subsets and concatenating')
                  make_subsets(tmp.dir)
                  concat_files(tmp.dir,gcm)
                  print('Done concat')
               } else {
                  print('Copying assembled file')
                  file.copy(from=grouped.files,to=paste0(tmp.dir,gcm),overwrite=TRUE)
               }

               file.remove(list.files(path=paste0(tmp.dir,'grouptmp'),pattern='*.nc',full.name=TRUE))
               file.remove(list.files(path=paste0(tmp.dir,'timetmp'),pattern='*.nc',full.name=TRUE))
            } else {
               print(paste0('Only one file for run: ',run))
               ##print(run.files)
               ##used.files <- c(used.files,run.files)
            }
         }    
      }

      ux <- gcm.info$files %in% used.files
      print('Step before copying to assembly directory')
      print('Copying back')
      ##print('Files left to copy')
      ##print(gcm.info$files[!ux])
      if (!file.exists(paste0(write.dir,gcm,'/'))) {       
         print(list.files(path=paste0(tmp.dir,gcm)))
         file.copy(from=paste0(tmp.dir,gcm),to=write.dir,recursive=TRUE,overwrite=TRUE)
      } else {
         cp.files <- list.files(path=paste0(tmp.dir,gcm),full.name=TRUE)
         print(cp.files)
         file.copy(from=cp.files,to=paste0(write.dir,gcm,'/'),overwrite=TRUE)
      }

      ##Unused files - copy to new directory    
##      file.copy(from=gcm.info$files[!ux],to=paste0(write.dir,gcm),overwrite=TRUE)
      ##Clean Up
      files.rm <- list.files(path=paste0(tmp.dir,gcm,'/'),full.name=TRUE)
      file.remove(files.rm)
      file.remove(paste0(tmp.dir,gcm))

  }

