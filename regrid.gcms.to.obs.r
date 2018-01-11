##Script to convert the gcms to the same grid as observations.

scenario <- 'rcp85'
obs <- 'hadcru4'

base.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/CMIP5/global/annual/'
grid.file <- paste0('/storage/home/ssobie/grid_files/',obs,'.grid.txt')

gcm.list <- c('ACCESS1-0', 'bcc-csm1-1','CanESM2','CSIRO-Mk3-6-0', 'GFDL-CM3','HadGEM2-ES','IPSL-CM5A-MR','MIROC-ESM','MPI-ESM-MR',
              'ACCESS1-3',  'bcc-csm1-1-m','CCSM4','GFDL-ESM2G','HadGEM2-AO', 'inmcm4','IPSL-CM5B-LR','MIROC-ESM-CHEM','MRI-CGCM3',
              'CNRM-CM5','GFDL-ESM2M','HadGEM2-CC','IPSL-CM5A-LR','MIROC5','MPI-ESM-LR','NorESM1-M')

##gcm.list <- c('ACCESS1-0','ACCESS1-3','BNU-ESM','CanESM2','CSIRO-Mk3-6-0','GFDL-ESM2M','IPSL-CM5A-MR','IPSL-CM5B-LR',
##              'MPI-ESM-LR','MPI-ESM-MR','NorESM1-M','bcc-csm1-1','bcc-csm1-1-m','CNRM-CM5','inmcm4',
##              'HadGEM2-CC','HadGEM2-ES','IPSL-CM5A-LR','MIROC5','MIROC-ESM','MIROC-ESM-CHEM','MRI-CGCM3')

gcm.list <- 'IPSL-CM5B-LR'
for (gcm in gcm.list) {
  print(gcm)

  gcm.dir <- paste(base.dir,gcm,sep='')
  
  all.files <- list.files(path=gcm.dir,pattern=paste('tas_ann_',gcm,sep=''),full.name=TRUE)
  scen.files <- all.files[grep('18500101',all.files)]
  ##print(scen.files)
  tas.file <- scen.files[grep(scenario,scen.files)]
  ##print('Tas file')
  ##print(tas.file)
  interp.file <- gsub(pattern='tas_ann',replacement=paste0('tas_ann_',obs,'_grid'),tas.file)
  ##print('interp file')
  ##print(interp.file)
  work <- paste0('cdo -s -O remapbil,',grid.file,' ',tas.file,' ',interp.file)
  print(work)
  system(paste('cdo -s -O remapbil,',grid.file,' ',tas.file,' ',interp.file,sep=''))

}