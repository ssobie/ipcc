###Script to assemble a human readable inventory of CMIP5 data at PCIC

base.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/CMIP5/global/'

write.dir <- '/storage/data/projects/rci/stat.downscaling/inventories/'


variable <- 'pr'
header <- c('Precipitation','Location: ',base.dir,'')
write.file <- paste0(write.dir,variable,'_Global_CMIP5_model_inventory.csv')

models <- list.dirs(path=base.dir,full.names=FALSE,recursive=FALSE)
omit <- c('annual','grouptmp','spacetmp','timetmp','historicalNAT','control')

models <- models[!(models %in% omit)]

full.names <- full.list <- c('Model','Scenario','Run','Start','End')

scen.list <- c('historical','rcp26','rcp45','rcp60','rcp85','piControl')

for (m in seq_along(models)) {
  model <- models[m]
  print(model)

  model.files <- as.list(list.files(path=paste(base.dir,model,sep=''),pattern=paste(variable,'_day',sep=''),recursive=TRUE))

  if (length(model.files)!=0) {
    print(model.files)
    fsplit <- lapply(model.files,function(x) {strsplit(x,'_')[[1]]})
    file.matrix <- matrix(NA,nrow=length(model.files),ncol=length(full.names))
    for (i in 1:length(model.files)) {
      file.only <- model.files[[i]]
      file.slct <- fsplit[[i]]
      file.matrix[i,1] <- model
      ##Scenario
      file.matrix[i,2] <- file.slct[grepl('(historical|rcp|Control)',file.slct)]
      ##Run                       
      file.matrix[i,3] <- file.slct[grep('r[0-9]{,2}i[0-9]{,2}p[0-9]{,2}',file.slct)]
      ##Start and End Years
      years <- unlist(regmatches(file.only,gregexpr('[0-9]{8}-[0-9]{8}',file.only)))
      file.matrix[i,4] <- substr(years,1,4)
      file.matrix[i,5] <- substr(years,10,13)

    }
    full.list <- rbind(full.list,file.matrix)
  }
}

##Sort the data into a more readable format
read.list <- c('Model','Scenario','Runs','Years')

for (m in seq_along(models)) {
  model <- models[m]
  sub.list <- full.list[full.list[,1]==model,]
  if (is.null(dim(sub.list))) {
    rv <- c(model,sub.list[2],sub.list[3],paste0(sub.list[4],'-',sub.list[5]))      
    read.list <- rbind(read.list,rv)
  } else {
    scenarios <- unique(sub.list[,2])
    for (scen in scenarios) {
      scen.list <- sub.list[sub.list[,2]==scen,]
      if (is.null(dim(scen.list))) {
        runs <- scen.list[3]
        years <- paste(scen.list[4],scen.list[5],sep='-')
      } else {
        runs <- paste(unique(scen.list[,3]),collapse=' | ')
        bnds <- cbind(as.character(scen.list[,4]),as.character(scen.list[,5]))
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

write.table(rbind(header,read.list),file=write.file,quote=FALSE,row.name=F,col.name=F,sep=',')

