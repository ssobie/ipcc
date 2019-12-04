###Script to assemble a human readable inventory of CMIP5 data at PCIC

base.dir <- '/storage/data/climate/CMIP5/daily/'

write.dir <- '/storage/data/projects/rci/stat.downscaling/inventories/'

var.list <- c('pr','tasmax','tasmin')
header.list <- list(c('Precipitation','Location: ',base.dir,''),
                        c('Maximum Temperature','Location: ',base.dir,''),
                        c('Minimum Temperature','Location: ',base.dir,''))

gcms <- list.dirs(path=base.dir,full.names=FALSE,recursive=FALSE)

scen.list <- c('historical','historicalGHG','historicalNat','historicalMisc',
               'rcp26','rcp45','rcp60','rcp85')

for (v in seq_along(var.list)) {
   full.names <- full.list <- c('Model','Scenario','Run','Start','End')
   variable <- var.list[v]
   header <- header.list[[v]]

   for (g in seq_along(gcms)) {
     print(g)
     gcm <- gcms[g]
     print(gcm)

     ##GCM files
     gcm.files <- list.files(path=paste(base.dir,gcm,sep=''),pattern=paste(variable,'_day',sep=''),recursive=TRUE)
     file.matrix <- matrix(NA,nrow=length(gcm.files),ncol=length(full.names))
     for (i in 1:length(gcm.files)) {
        file <- gcm.files[i]
        file.split <- strsplit(file,'_')[[1]]
        run <- file.split[grepl('r[0-9]{1,2}i1p[0-9]{1}',file.split)]
        scenario <- file.split[grepl('historical',file.split)]
        future <- grepl('rcp',scenario)

        ##GCM Name
        file.matrix[i,1] <- gcm
        ##Scenario
        file.matrix[i,2] <- scenario
        ##Run                       
        file.matrix[i,3] <- run
        ##Start and End Years
        years <- unlist(regmatches(file,gregexpr('[0-9]{8}-[0-9]{8}',file)))
        file.matrix[i,4] <- substr(years,1,4)
        file.matrix[i,5] <- substr(years,10,13)
     }
     full.list <- rbind(full.list,file.matrix) 
  }

  ##Sort the data into a more readable format
  read.list <- c('Model','Scenario','Runs','Years')

  for (m in seq_along(gcms)) {
    model <- gcms[m]
    sub.list <- full.list[full.list[,1]==model,]
    if (is.null(dim(sub.list))) {
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

  write.table(rbind(header,read.list),file=paste0(write.dir,variable,'_assembled_daily_CMIP5_model_inventory.csv'),quote=FALSE,row.name=F,col.name=F,sep=',')
}
