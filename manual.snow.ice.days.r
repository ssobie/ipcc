##Script to convert the gcms to the same grid as observations.

##******************************************************

scenarios <- c('rcp26','rcp45','rcp85') ##,'rcp85')

tmp.dir <- '/local_temp/ssobie/count/'

if (!file.exists(tmp.dir)) {
  dir.create(tmp.dir,recursive=TRUE)
}
##Clean tmpdir
##rm.tmp <- paste0("rm ",tmp.dir,"*nc")
##system(rm.tmp)


gcm <- 'MIROC-ESM-CHEM'
varname <- 'sic'
##------------------------------------


base.dir <- '/storage/data/projects/CanSISE/n2w/CMIP5/monthly/'
gcm.dir <- paste(base.dir,gcm,sep='')
all.files <- list.files(path=gcm.dir,pattern=paste0(varname,'_OImon_',gcm))
scen.files <- all.files[grep('1850',all.files)]
 
for (scenario in scenarios) {
  ice.files <- scen.files[grep(scenario,scen.files)]
  print(ice.files)
  if (!length(ice.files)==0) {
  
    split.apart <- unlist(strsplit(unlist(ice.files),'_'))
    runs <- 'r1i1p1' ##unique(split.apart[grep('r*i1p1',split.apart)])

    print(runs)
    for (run in runs) {
      print(run)
      run.file <- ice.files[grep(run,ice.files)]
      print(run.file)

      move.to <- paste("rsync -av --progress ",gcm.dir,'/',run.file, " " ,tmp.dir,sep='')
      print(move.to)
      system(move.to)
      data.file <- paste0(tmp.dir,run.file)
      print(data.file)
      count.file <- gsub(pattern=paste0(varname,'_OImon'),replacement=paste0(varname,'_count'),data.file)
      tmp.file <- gsub(pattern=paste0(varname,'_OImon'),replacement=paste0(varname,'_TMP'),data.file)
      yr.file <- gsub(pattern=paste0(varname,'_OImon'),replacement=paste0(varname,'_ann'),data.file)
    
##      work <- paste0(" cdo -expr,\'",varname,"=(",varname,"<0.15)?0.0:",varname,";\' ",data.file," ",tmp.file)
##      print(work)
##      system(work)
##      work <- paste0(" cdo -expr,\'",varname,"=(",varname,">=0.15)?1.0:",varname,";\' ",tmp.file," ",yr.file)
##      system(work)

      print('Dividing')
      work <- paste0(" cdo -s -O divc,100 ",data.file," ",tmp.file)
      system(work)

      print('Monthly')
      work <- paste0(" cdo -s -O muldpm ",tmp.file," ",yr.file)
      system(work)

      print('Yearly')
      work <- paste0(" cdo -s -O yearsum ",yr.file," ",count.file)
      system(work)
      
      move.back <- paste("rsync -av --progress ",count.file," ",gcm.dir,sep='')
      print(move.back)
      system(move.back)

##      clean.up <- paste0('rm ',tmp.file)
##      system(clean.up)
##      clean.up <- paste0('rm ',yr.file)
##      system(clean.up)
##      clean.up <- paste0('rm ',count.file)
##      system(clean.up)
##      clean.up <- paste0('rm ',data.file)
##      system(clean.up)



    }##Run loop
  }##If statement  
}##Scenarios

##}