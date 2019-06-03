###Script to assemble a human readable inventory of CMIP5 data at PCIC

write.dir <- '/storage/data/projects/rci/stat.downscaling/inventories/'
variable <- 'pr'
header <- c('Precipitation','','','')
incoming.dir <- '/storage/data/climate/CMIP5/incoming/output1/'
output.dir <- '/storage/data/climate/CMIP5/output/'

incoming.centres <- list.dirs(path=incoming.dir,full.names=FALSE,recursive=FALSE)
output.centres <-  list.dirs(path=output.dir,full.names=FALSE,recursive=FALSE)

centres <- sort(unique(c(incoming.centres,output.centres)))

full.names <- full.list <- c('Centre','Location','File','Output File Sizes','Incoming File Sizes', 'File Size Difference')

for (c in seq_along(centres)) {
  centre <- centres[c]
  print(centre)
  ##PR files
  incoming.paths <- as.list(list.files(path=paste(incoming.dir,centre,sep=''),pattern=paste(variable,'_day',sep=''),recursive=TRUE))
  output.paths <- as.list(list.files(path=paste(output.dir,centre,sep=''),pattern=paste(variable,'_day',sep=''),recursive=TRUE))
  
  if (length(incoming.paths)!=0 | length(output.paths)!=0) {

    incoming.fsplit <- lapply(incoming.paths,function(x) {strsplit(x,'/')[[1]]})
    incoming.files <- lapply(incoming.fsplit,function(x,y) {return(x[grep(y,x)])},paste0(variable,'_day'))

    output.fsplit <- lapply(output.paths,function(x) {strsplit(x,'/')[[1]]})
    output.files <- lapply(output.fsplit,function(x,y) {return(x[grep(y,x)])},paste0(variable,'_day'))
    
    incoming.match <- incoming.files %in% output.files
    incoming.sub <- incoming.files[incoming.match]
    output.match <- output.files %in% incoming.sub    
    output.sub <- output.files[output.match]
    incoming.sub <- incoming.files[incoming.files %in% output.sub]

    if (sum(output.match)!=0) {        
        output.paths <- unlist(output.paths)[output.match]
        incoming.paths <- unlist(incoming.paths)[incoming.files %in% output.sub]

        output.fsplit <- output.fsplit[output.match]
        output.dirs <- lapply(output.fsplit,function(x){return(paste0(x[-length(x)],collapse='/'))})  
        output.sizes <- file.size(paste0(output.dir,centre,'/',output.paths))         
        incoming.sizes <- file.size(paste0(incoming.dir,centre,'/',incoming.paths))         

        print(length(rep(centre,length(output.dirs))))
        print(length(output.dirs))
        print(length(output.sub))
        print(length(output.sizes))
        print(length(incoming.sizes))
        file.info <- cbind(rep(centre,length(output.dirs)),output.dirs,output.sub,output.sizes,incoming.sizes,output.sizes-incoming.sizes)
        full.list <- rbind(full.list,file.info)
                   
    }
  }
}


write.table(full.list,file=paste0(write.dir,variable,'_CMIP5_exact_duplicates.csv'),quote=FALSE,row.name=F,col.name=F,sep=',')

