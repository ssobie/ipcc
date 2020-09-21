###Script to assemble a human readable inventory of CMIP5 data at PCIC

##base.dir <- '/storage/data/climate/CMIP5/incoming/output1/'
base.dir <- '/storage/data/climate/CMIP5/output/'
##base.dir <- '/storage/data/climate/downscale/BCCAQ2/CMIP5/'

write.dir <- '/storage/data/projects/rci/stat.downscaling/inventories/'
variable <- 'tasmax'

header <- c('Maximum Temperature','Location: ',base.dir,'')
##header <- c('Minimum Temperature','Location: ',base.dir,'')
##header <- c('Precipitation','Location: ',base.dir,'')

centres <- list.dirs(path=base.dir,full.names=FALSE,recursive=FALSE)
omit <- 'tmp_src'
centres <- centres[!(centres==omit)]

if (1==1) {
full.names <- full.list <- c('Centre','Model','Scenario','Run','Start','End','Version')

scen.list <- c('historical','historicalGHG','historicalNat','historicalMisc','rcp26','rcp45','rcp60','rcp85')

for (c in seq_along(centres)) {
  centre <- centres[c]
  print(centre)
  ##PR files
  pr.files <- as.list(list.files(path=paste(base.dir,centre,sep=''),pattern=paste(variable,'_day',sep=''),recursive=TRUE))

  if (length(pr.files)!=0) {
    print(pr.files)
    fsplit <- lapply(pr.files,function(x) {strsplit(x,'/')[[1]]})
    file.matrix <- matrix(NA,nrow=length(pr.files),ncol=length(full.names))
    for (i in 1:length(pr.files)) {
      file.only <- fsplit[[i]][length(fsplit[[i]])]           
      ##Remove the file name from the end
      file.slct <- fsplit[[i]][-length(fsplit[[i]])]           
      ##Centre Name
      file.matrix[i,1] <- centre
      ##Model Name
      file.matrix[i,2] <- file.slct[1]
      ##Scenario
      file.matrix[i,3] <- file.slct[grepl('(historical|rcp)',file.slct)]
      ##Run                       
      file.matrix[i,4] <- file.slct[grep('r[0-9]{,2}i[0-9]{,2}p[0-9]{,2}',file.slct)]
      ##Start and End Years
      years <- unlist(regmatches(file.only,gregexpr('[0-9]{8}-[0-9]{8}',file.only)))
      file.matrix[i,5] <- substr(years,1,4)
      file.matrix[i,6] <- substr(years,10,13)
      ##Version (If it has one)
      version <- file.slct[grep('v[0-9]{8}',file.slct)]    
      if (length(version) > 0) {file.matrix[i,7] <- version}
      else {file.matrix[i,7] <- 'None'}
      
    }
    full.list <- rbind(full.list,file.matrix)

  }
}
}

##Sort the data into a more readable format
read.list <- c('Model','Scenario','Runs','Years')
models <- unique(full.list[-1,2])

for (m in seq_along(models)) {
  model <- models[m]
  sub.list <- full.list[full.list[,2]==model,]
  if (is.null(dim(sub.list))) {
  } else {
    scenarios <- unique(sub.list[,3])
    for (scen in scenarios) {
      scen.list <- sub.list[sub.list[,3]==scen,]
      if (is.null(dim(scen.list))) {
        runs <- scen.list[4]
        years <- paste(scen.list[5],scen.list[6],sep='-')
      } else {
        runs <- paste(unique(scen.list[,4]),collapse=' | ')
        bnds <- cbind(as.character(scen.list[,5]),as.character(scen.list[,6]))
        years <- apply(bnds,1,paste,collapse='-')
      }
      if (length(unique(years))>1) {
        ##warning(paste('Years in models runs are not the same for: ', model,'-',scen,': ',years,sep=''))
         ylen <- length(unique(years))
         rv <- cbind(rep(model,ylen),rep(scen,ylen),rep(runs,ylen),unique(years))
      } else {
         rv <- as.character(c(model,scen,runs,unique(years)))
      }
      read.list <- rbind(read.list,rv)
    }
  }
}

##write.table(rbind(header,read.list),file=paste0(write.dir,variable,'_CMIP5_model_inventory.csv'),quote=FALSE,row.name=F,col.name=F,sep=',')
##write.table(rbind(header,read.list),file=paste0(write.dir,variable,'_CMIP5_incoming_model_inventory.csv'),quote=FALSE,row.name=F,col.name=F,sep=',')
write.table(rbind(header,read.list),file=paste0(write.dir,variable,'_CMIP5_output_model_inventory.csv'),quote=FALSE,row.name=F,col.name=F,sep=',')
