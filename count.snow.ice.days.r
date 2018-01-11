##Script to convert the gcms to the same grid as observations.

##******************************************************

scenarios <- 'rcp85' ##c('rcp26','rcp45') ##,'rcp85')

args <- commandArgs(trailingOnly=TRUE)
for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
}

tmp.dir <- tmpdir ##'/local_temp/ssobie/count/'
##varname <- 'sic'

if (!file.exists(tmp.dir)) {
  dir.create(tmp.dir,recursive=TRUE)
}
##Clean tmpdir
##rm.tmp <- paste0("rm ",tmp.dir,"*nc")
##system(rm.tmp)
system('rm /local_temp/ssobie/BCCAQ2/*nc')

##------------------------------------

##gcm.list <- c('CanESM2','CNRM-CM5','CSIRO-Mk3-6-0',
##              'HadGEM2-CC','HadGEM2-ES','inmcm4','IPSL-CM5A-LR',
##              'MIROC5','MIROC-ESM','MIROC-ESM-CHEM','MPI-ESM-LR')

##              'bcc-csm1-1-m','MPI-ESM-MR','MRI-CGCM3','NorESM1-M',
##              'ACCESS1-3','FGOALS-g2')



##for (gcm in gcm.list) {

base.dir <- '/storage/data/projects/CanSISE/n2w/CMIP5/'
gcm.dir <- paste(base.dir,gcm,sep='')
scen.files <- list.files(path=gcm.dir,pattern=paste0(varname,'_day_',gcm))

  
for (scenario in scenarios) {
  ice.files <- scen.files[grep(scenario,scen.files)]
  print(ice.files)
  if (!length(ice.files)==0) {
  
    split.apart <- unlist(strsplit(unlist(ice.files),'_'))
    runs <- unique(split.apart[grep('r*i1p1',split.apart)])
    print(runs)
    for (run in runs) {
      print(run)
      run.file <- ice.files[grep(run,ice.files)]
      print(run.file)

      move.to <- paste("rsync -av ",gcm.dir,'/',run.file, " " ,tmp.dir,sep='')
      print(move.to)
      system(move.to)
      data.file <- paste0(tmp.dir,run.file)
      print(data.file)
      count.file <- gsub(pattern=paste0(varname,'_day'),replacement=paste0(varname,'_count'),data.file)
      tmp.file <- gsub(pattern=paste0(varname,'_day'),replacement=paste0(varname,'_TMP'),data.file)
      yr.file <- gsub(pattern=paste0(varname,'_day'),replacement=paste0(varname,'_ann'),data.file)
    
      work <- paste0(" cdo -expr,\'",varname,"=(",varname,"<0.15)?0.0:",varname,";\' ",data.file," ",tmp.file)
      print(work)
      system(work)
      work <- paste0(" cdo -expr,\'",varname,"=(",varname,">=0.15)?1.0:",varname,";\' ",tmp.file," ",yr.file)
      system(work)
      work <- paste0(" cdo -s -O yearsum ",yr.file," ",count.file)
      system(work)

      move.back <- paste("rsync -av ",count.file," ",gcm.dir,sep='')
      print(move.back)
      system(move.back)

      clean.up <- paste0('rm ',tmp.file)
      system(clean.up)
      clean.up <- paste0('rm ',yr.file)
      system(clean.up)
      clean.up <- paste0('rm ',count.file)
      system(clean.up)
      clean.up <- paste0('rm ',data.file)
      system(clean.up)

    }##Run loop
  }##If statement  
}##Scenarios

##}