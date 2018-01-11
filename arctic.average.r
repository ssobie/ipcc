

get.arctic.average <- function(nc,obs.mask,obs.st,obs.en,clim.st,clim.en) {
                      
   lon <- ncvar_get(nc,'lon')
   lat <- ncvar_get(nc,'lat')
   time <- netcdf.calendar(nc)
  
   var.name <- 'tas'
   tas.series <- ncvar_get(nc,var.name) - 273
   
   tas.subset <- tas.series[,,obs.st:obs.en]
   rm(tas.series)              
   tas.subset[obs.mask] <- NA
   
   print(dim(tas.subset))

   tas.clim <- apply(tas.subset[,,clim.st:clim.en],c(1,2),mean,na.rm=T)
   ##tas.clim <- apply(tas.subset,c(1,2),mean,na.rm=T)
   tas.sd.clim <- apply(tas.subset,c(1,2),sd,na.rm=T)
   tas.avg <- aperm(apply(tas.clim,c(1,2),function(x){rep(x,each=dim(tas.subset)[3])}),c(2,3,1))
   tas.sd <- aperm(apply(tas.sd.clim,c(1,2),function(x){rep(x,each=dim(tas.subset)[3])}),c(2,3,1))
   tas.anom <- (tas.subset - tas.avg)##/tas.sd

   rm(tas.subset)
   ar <- get.areas(lon,lat,dim(tas.anom))
   wts <- ar[,1]/sum(ar[,1])

   ac.ix <- lat > 60

   arctic.wts <- ar[ac.ix,1]/sum(ar[ac.ix,1])
   global.wts <- ar[,1]/sum(ar[,1])

   arctic.subset <- tas.anom[,ac.ix,]
   arctic.mid <- apply(arctic.subset,c(1,3),weighted.mean,arctic.wts,na.rm=T)
   arctic.mean <- apply(apply(arctic.subset,c(1,3),weighted.mean,arctic.wts,na.rm=T),2,mean,na.rm=T)

   global.subset <- tas.anom
   global.mean <- apply(apply(global.subset,c(1,3),weighted.mean,global.wts,na.rm=T),2,mean,na.rm=T)
   nc_close(nc)

   if (length(arctic.mean) == 151) {
     arctic.mean <- arctic.mean[-1]
     global.mean <- global.mean[-1]
   }
   rv <- list(arctic=arctic.mean,       
              global=global.mean)

   return(rv)
}
