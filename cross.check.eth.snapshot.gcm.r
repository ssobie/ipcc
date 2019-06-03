##Script to check whether the CMIP5 snapshot at ETH Zurich has the same GCMs
##stored as we do.


##--------------------------------------------------------------------------
hydro.list <- rbind(c('ACCESS1-0','rcp45','r1i1p1'),
                    c('ACCESS1-0','rcp85','r1i1p1'),
                c('CanESM2','rcp45','r1i1p1'),
                c('CanESM2','rcp45','r2i1p1'),
                c('CanESM2','rcp45','r3i1p1'),
                c('CanESM2','rcp45','r4i1p1'),
                c('CanESM2','rcp45','r5i1p1'),
                c('CanESM2','rcp85','r1i1p1'),
                c('CanESM2','rcp85','r2i1p1'),
                c('CanESM2','rcp85','r3i1p1'),
                c('CanESM2','rcp85','r4i1p1'),
                c('CanESM2','rcp85','r5i1p1'),
                c('CCSM4','rcp45','r1i1p1'),
                c('CCSM4','rcp45','r2i1p1'),
                c('CCSM4','rcp45','r6i1p1'),
                c('CCSM4','rcp85','r1i1p1'),
                c('CCSM4','rcp85','r2i1p1'),
                c('CCSM4','rcp85','r6i1p1'),
                c('CNRM-CM5','rcp45','r1i1p1'),
                c('CNRM-CM5','rcp85','r1i1p1'),
                c('HadGEM2-ES','rcp45','r1i1p1'),
                c('HadGEM2-ES','rcp45','r2i1p1'),
                c('HadGEM2-ES','rcp45','r3i1p1'),
                c('HadGEM2-ES','rcp45','r4i1p1'),
                c('HadGEM2-ES','rcp85','r1i1p1'),
                c('HadGEM2-ES','rcp85','r2i1p1'),
                c('HadGEM2-ES','rcp85','r3i1p1'),
                c('HadGEM2-ES','rcp85','r4i1p1'),
                c('MPI-ESM-LR','rcp45','r1i1p1'),
                c('MPI-ESM-LR','rcp45','r2i1p1'),
                c('MPI-ESM-LR','rcp45','r3i1p1'),
                c('MPI-ESM-LR','rcp85','r1i1p1'),
                c('MPI-ESM-LR','rcp85','r2i1p1'),
                c('MPI-ESM-LR','rcp85','r3i1p1'))

##--------------------------------------------------------------------------

extract_years <- function(file) {

    ist <- regexpr('[0-9]{8}-',file)
    zst <- ist[1] + attr(ist, "match.length")-6
    ien <- regexpr('-[0-9]{8}',file)
    zen <- ien[1]  + attr(ien, "match.length")-5
    yst <- substr(file,ist,zst)
    yen <- substr(file,ien+1,zen)

    rv <- as.numeric(yst):as.numeric(yen)
    return(rv)
}

##--------------------------------------------------------------------------


if (!file.exists('/storage/data/projects/rci/stat.downscaling/inventories/cmip5.snapshot.files.RData')) {
   snapshots <- read.csv('/storage/data/projects/rci/stat.downscaling/inventories/cmip5.snapshot.files2.csv',header=T)
   snapshot.files <- snapshots[,2]
   files <- sapply(snapshot.files,as.character)
   rm(snapshots)
   rm(snapshot.files)
   save(files,file='/storage/data/projects/rci/stat.downscaling/inventories/cmip5.snapshot.files.RData')
} else {
   load('/storage/data/projects/rci/stat.downscaling/inventories/cmip5.snapshot.files.RData')
}

hydro.match <- list(pr=list(historical=list(),rcp26=list(),rcp45=list(),rcp85=list()),
                    tasmax=list(historical=list(),rcp26=list(),rcp45=list(),rcp85=list()),
                    tasmin=list(historical=list(),rcp26=list(),rcp45=list(),rcp85=list()))

files.split <- sapply(files,function(x){strsplit(x,'_')[[1]]})

var.names <- c('pr','tasmax','tasmin')
scenarios <- c('historical','rcp26','rcp45','rcp85')

for (var.name in var.names) {
   print(var.name)
   for (scenario in scenarios) {
   print(scenario)
##var.name <- 'tasmax'
##scenario <- 'rcp26'

      ix <- grepl(paste0(var.name,'_day'),files)
      var.files <- files[ix]

      hist.files <- var.files[grepl(paste0('_',scenario,'_'),var.files)]
      hist.split <- sapply(hist.files,function(x){strsplit(x,'_')[[1]]})

      ##Assumes the same file name format for indexing components
      models <- unique(hist.split[3,])

      file.info <- c('Model','Scenario','Run','Start','End','Continous?')

      for (i in seq_along(models)) {
         gcm <- models[i]
         model.files <- hist.files[grepl(gcm,hist.files)]
         model.split <- lapply(model.files,function(x){strsplit(x,'_')[[1]]})
         runs <- unlist(lapply(model.split,function(x){return(x[grep('r[0-9]',x)])}))
         unique.runs <- unique(runs)
         for (j in seq_along(unique.runs)) {
            run <- unique.runs[j]
            run.files <- model.files[runs==run]

            years <- unique(sort(unlist(sapply(run.files,extract_years))))

            ##yst <- sort(years[1,])
            ##yen <- sort(years[2,])
            year.range <- range(years)
            year.check <- unique(diff(years)) ###unlist(mapply(seq,yst,yen))))

            flag <- 'Yes'
            if (length(year.check) !=1) {
               'Years in this file are not continuous'
               print(year.check)

               full.years <- year.range[1]:year.range[2]
               missing <- full.years[!(full.years %in% years)]
               flag <- paste0('Missing:',paste0(missing,collapse='|'))  
            }
            info.add <- c(gcm,scenario,run,year.range[1],year.range[2],flag)
            file.info <- rbind(file.info,info.add)

         }
      }

      ##Check hydro list 
      hlen <- nrow(hydro.list)
      hydro.fill <- c('Model','Scenario','Run','Start','End','Continous?')
      for (k in 1:hlen) {
         hi <- which(file.info[,1] %in% hydro.list[k,1] &  file.info[,3] %in% hydro.list[k,3])
         if (length(hi)==1) {
            hydro.add <- file.info[hi,]
         } else {
            hydro.add <- c(hydro.list[k,1],scenario,hydro.list[k,3],'MISSING','MISSING','MISSING')
         }
         hydro.fill <- rbind(hydro.fill,hydro.add)
      }
      hydro.match[[var.name]][[scenario]] <- hydro.fill
      hydro.file <- paste0('/storage/data/projects/rci/stat.downscaling/inventories/hydro.check.',scenario,'.',var.name,'.csv')
      write.table(hydro.fill,hydro.file,sep=',',quote=FALSE,row.name=FALSE,col.name=F)
      
      write.file <- paste0('/storage/data/projects/rci/stat.downscaling/inventories/cmip5.snapshot.',scenario,'.',var.name,'.csv')
      write.table(file.info,write.file,sep=',',quote=FALSE,row.name=FALSE,col.name=F)




   }
}

