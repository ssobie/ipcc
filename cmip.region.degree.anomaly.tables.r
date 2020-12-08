##Script to produce tables of projected temperature and precipitation changes
##The output format is the same as the MOTI tables 

library(raster)
library(ncdf4)
library(rgeos)
library(PCICt)
library(rgdal)

##-----------------------------------------------------------------------------------------------
##Area average functions

get_shape_file <- function(region,shape.dir) {
   reg.shp <- readOGR(shape.dir,region,stringsAsFactors=F)
   return(reg.shp)
}

cell_mask <- function(file.brick,shape) {
   ras <- rasterize(shape,file.brick,getCover=T)
   ras[ras==0] <- NA
   return(ras)
}

##--------------------------------------------------------------------------------------

fix_gcm_coords <- function(gcm.file,tmp.dir) {

   file.name <- basename(gcm.file)
   file.fix <- paste0("coord_fix_",basename(gcm.file))
   if (length(file.name) > 1 | length(file.name)==0) { browser()}
   
   if (grepl('(FGOALS-g2|GFDL-ESM2G|GFDL-ESM2M|FGOALS-g3)',file.name)) {
      nc <- nc_open(paste0(tmp.dir,file.name))
      nlon <- nc$dim$lon$len
      nlat <- nc$dim$lat$len
      nc_close(nc)
      file.reg <- paste0("reg_fix_",file.name)
      work <- paste0("cdo -s remapcon,r",nlon,"x",nlat," ",tmp.dir,file.name," ",tmp.dir,file.reg)
      system(work)
      Sys.sleep(1)

      work <- paste0("cdo -s sellonlatbox,-180,180,-90,90 ",tmp.dir,file.reg," ",tmp.dir,file.fix)
      system(work)
      Sys.sleep(1)
      file.remove(paste0(tmp.dir,file.name))
      file.remove(paste0(tmp.dir,file.reg))  
   } else {
  
      work <- paste0("cdo -s sellonlatbox,-180,180,-90,90 ",tmp.dir,file.name," ",tmp.dir,file.fix)
      system(work)
      Sys.sleep(1)
      file.remove(paste0(tmp.dir,file.name))
   }
   return(file.fix)
}
 

##--------------------------------------------------------------------------------------

##Split apart area masked to remove area weight calculation to happen once

new_area_weighted_average <- function(clim.file,area.shape,area.info,tmp.dir) {

   if (grepl('return_periods',clim.file)) {
      file.copy(from=clim.file,to=tmp.dir,overwrite=TRUE)
      fix.file <- fix_gcm_coords(clim.file,tmp.dir)
   }  else {
      file.copy(from=clim.file,to=tmp.dir,overwrite=TRUE)
      fix.file <- basename(clim.file) ###fix_gcm_coords(clim.file,tmp.dir)
   }

   file.brick <- brick(paste0(tmp.dir,fix.file))
   area.masked <- mask(file.brick,area.info$overlay)

   area.weighted <- area.masked * area.info$weights
   area.average <- cellStats(area.weighted,sum)

   rm(file.brick)
   file.remove(paste0(tmp.dir,fix.file))
   
   return(area.average)   

}



##--------------------------------------------------------------------------------------

get_rounding_value <- function(var.name) {
                   
  zero.list <- c('pr','cdd','gdd','fdd','hdd','pas',
                 'fdETCCDI','suETCCDI','idETCCDI','gslETCCDI',
                 'wsdiETCCDI','csdiETCCDI',
                 'rx1dayETCCDI','rx2dayETCCDI','rx5dayETCCDI',
                 'sdiiETCCDI','r1mmETCCDI','r10mmETCCDI','r20mmETCCDI',
                 'cddETCCDI','cdd90ETCCDI','cddmaxETCCDI','cwdETCCDI',
                 'r95pETCCDI','r95daysETCCDI','r95distETCCDI','r95sepETCCDI',
                 'r99pETCCDI','r99daysETCCDI',
                 'prcptotETCCDI')
 
  one.list <- c('tasmax','tasmin','trETCCDI',
               'txxETCCDI','tnxETCCDI','txnETCCDI','tnnETCCDI',
                'tn10pETCCDI','tx10pETCCDI','tn90pETCCDI','tx90pETCCDI',
                'dtrETCCDI',
                'pr_RP5','pr_RP20','pr_RP50',
                'tasmax_RP5','tasmax_RP20',
                'tasmin_RP5','tasmin_RP20')
  rd <- 1
  if (length(grep(var.name,zero.list))>0)
    rd <- 0
  if (length(grep(var.name,one.list))>0)
    rd <- 1

  return(rd)
}

##--------------------------------------------------------------------------------------

get_variable_title <- function(var.name) {
  rv <- switch(var.name,
                      pr='Total Precipitation',
                      tasmax='Maximum Temperature',
                      tasmin='Minimum Temperature',
                      tas='Average Temperature',
                      snowdepth='Average Snowpack (m)',
                      fdETCCDI='Frost Days',
                      suETCCDI='Summer Days',
                      su30ETCCDI='Summer Hot Days',
                      idETCCDI='Ice Days',
                      gslETCCDI='Growing Season Length',
                      wsdiETCCDI='Warm Spell Duration',
                      csdiETCCDI='Cold Spell Duration',
                      rx1dayETCCDI='One Day Precipitation',
                      rx5dayETCCDI='Five day precipitation',
                      sdiiETCCDI='Simple Daily Intensity Index',
                      r10mmETCCDI='Heavy Precipitation Days',
                      r20mmETCCDI='Very Heavy Precipitation Days',
                      cddETCCDI='Consecutive Dry Days',
                      cddmaxETCCDI='Maximum Consecutive Dry Days',
                      cdd90ETCCDI='90th % Consecutive Dry Days',
                      cwdETCCDI='Consecutive Wet Days',
                      r95pETCCDI='Very Wet Days',
                      r95daysETCCDI='Number of Very Wet Days',
                      r95distETCCDI='Time between Wet Day Events',
                      r95sepETCCDI='Wet Day Intervals',
                      r99pETCCDI='Extremely Wet Days',
                      r99daysETCCDI='Number of Extremely Wet Days',
                      prcptotETCCDI='Annual Total Precipitation',
                      trETCCDI='Tropical Nights',
                      txxETCCD='Hottest Days',
                      tnxETCCDI='Hottest Nights',
                      txnETCCDI='Coldest Days',
                      tnnETCCDI='Coldest Nights',
                      tn10pETCCDI='Cool Nights',
                      tx10pETCCDI='Cool Days',
                      tn90pETCCDI='Warm Nights',
                      tx90pETCCDI='Warm Days',
                      dtrETCCDI='Diurnal Temperature Range',
                      ffd='Frost Free Days',
                      fdd='Freezing Degree Days',
                      cdd='Cooling Degree Days',
                      hdd='Heating Degree Days',
                      gdd='Growing Degree Days',
                      pas='Precipitation as Snow',
                      pr_RP5=paste('Maximum Precipitation (5-Year)',sep='') ,
                      pr_RP20=paste('Maximum Precipitation (20-Year)',sep='') ,
                      pr_RP50=paste('Maximum Precipitation (50-Year)',sep='') ,
                      tasmax_RP5=paste('Maximum Temperature (5-Year)',sep='') ,
                      tasmax_RP20=paste('Maximum Temperature (20-Year)',sep='') ,
                      tasmin_RP5=paste('Minimum Temperature (5-Year)',sep=''),
                      tasmin_RP20=paste('Minimum Temperature (20-Year)',sep=''))



  return(rv)
}

##--------------------------------------------------------------------------------------
##Correct the table formatting
format_tables <- function(mon.vals,models,runs,rd,var.name,region.title,cmip,scenario,type) {

  vals.avg <- apply(mon.vals,2,mean,na.rm=T)
  table.vals <- rbind(mon.vals,vals.avg)       

    vals.10 <- apply(mon.vals,2,quantile,0.1,na.rm=T)
    vals.50 <- apply(mon.vals,2,quantile,0.5,na.rm=T)
    vals.90 <- apply(mon.vals,2,quantile,0.9,na.rm=T)
    table.vals <- rbind(table.vals,vals.10,vals.50,vals.90)

  new.table <- round(table.vals,rd)
  new.table <- cbind(c(models,'Ens_Avg.','10th_%ile','Median','90th_%ile'),c(runs,rep('',4)),new.table)

  if (type=="t_and_p") {
     new.table <- rbind(c('Model','Run',month.abb,'Winter','Spring','Summer','Fall','Annual'),new.table) 
     title <- c(paste('Table: ',get_variable_title(var.name),' for ',region.title,sep=''),cmip,scenario,rep(' ',16))
  }
  if (type=='annual') {
     new.table <- rbind(c('Model','Run','Annual'),new.table) 
     title <- c(paste('Table: ',get_variable_title(var.name),' for ',region.title,sep=''),cmip,scenario)
  }
  if (type=="seasonal") {
     new.table <- rbind(c('Model','Run','Winter','Spring','Summer','Fall','Annual'),new.table) 
     title <- c(paste('Table: ',get_variable_title(var.name),' for ',region.title,sep=''),cmip,scenario,rep(' ',4))
  }
  new.table <- rbind(title,new.table)
  return(new.table)
}

##--------------------------------------------------------------------------------------

##Read in the regional time series
t_and_p_regional_values <- function(var.name,area.info,
                                    read.dir,tmp.dir,
                                    clip.shp,
                                    past.int,deg) {
                                   

  ##GCM and Regional averages grouping monthly and seasonal values - needs pre-computed seasonal and monthly files
    seas.files <- list.files(path=paste0(read.dir,'seasonal/climatologies/'),pattern=paste0(var.name,'_seasonal'),full.name=TRUE)
    mon.files <- list.files(path=paste(read.dir,'monthly/climatologies/',sep=''),pattern=paste0(var.name,'_monthly'),full.name=TRUE)
    ann.files <- list.files(path=paste(read.dir,'annual/climatologies/',sep=''),pattern=paste0(var.name,'_annual'),full.name=TRUE)

    ##Add option to calculate TAS files

    ##-------------------------------------------------
    past.seas.file <- seas.files[grep(past.int,seas.files)]
    past.mon.file <- mon.files[grep(past.int,mon.files)]
    past.ann.file <- ann.files[grep(past.int,ann.files)]
    ##Extract subset of data for region
    past.seas.avg <- new_area_weighted_average(past.seas.file,clip.shp,area.info,tmp.dir)
    past.mon.avg <- new_area_weighted_average(past.mon.file,clip.shp,area.info,tmp.dir)
    past.ann.avg <- new_area_weighted_average(past.ann.file,clip.shp,area.info,tmp.dir)

    seas.files <- list.files(path=paste0(read.dir,'seasonal/degree_anomaly_climatologies/'),pattern=paste0(var.name,'_seasonal'),full.name=TRUE)
    mon.files <- list.files(path=paste(read.dir,'monthly/degree_anomaly_climatologies/',sep=''),pattern=paste0(var.name,'_monthly'),full.name=TRUE)
    ann.files <- list.files(path=paste(read.dir,'annual/degree_anomaly_climatologies/',sep=''),pattern=paste0(var.name,'_annual'),full.name=TRUE)

browser()

    proj.seas.file <- seas.files[grep(deg,seas.files)]
    proj.mon.file <- mon.files[grep(deg,mon.files)]
    proj.ann.file <- ann.files[grep(deg,ann.files)]    

    proj.seas.avg <- new_area_weighted_average(proj.seas.file,clip.shp,area.info,tmp.dir)
    proj.mon.avg <- new_area_weighted_average(proj.mon.file,clip.shp,area.info,tmp.dir)
    proj.ann.avg <- new_area_weighted_average(proj.ann.file,clip.shp,area.info,tmp.dir)

    ##Function to extract subset of data for moti region
    past.values <- c(past.mon.avg,past.seas.avg,past.ann.avg)
    proj.values <- c(proj.mon.avg,proj.seas.avg,proj.ann.avg)

    abs.anoms <- proj.values - past.values
    prc.anoms <- (proj.values - past.values)/past.values*100    
  rv <- rbind(past.values,proj.values,abs.anoms,prc.anoms)
  ##print(rv)

  return(rv)
}


##--------------------------------------------------------------------------------------

degree_day_regional_values <- function(var.name,area.info,
                                       read.dir,tmp.dir,
                                       clip.shp,
                                       past.int,deg) {                                  

  ##GCM and Regional averages grouping monthly and seasonal values - needs pre-computed seasonal and monthly files

  dd.files <- list.files(path=paste0(read.dir,'degree_days/climatologies/'),pattern=paste0(var.name,'_annual'),full.name=TRUE)

  ##-------------------------------------------------
  past.dd.file <- dd.files[grep(past.int,dd.files)]
  past.values <- new_area_weighted_average(past.dd.file,clip.shp,area.info,tmp.dir)

  dd.files <- list.files(path=paste0(read.dir,'degree_days/degree_anomaly_climatologies/'),pattern=paste0(var.name,'_annual'),full.name=TRUE)

  proj.dd.file <- dd.files[grep(deg,dd.files)]
  proj.values <- new_area_weighted_average(proj.dd.file,clip.shp,area.info,tmp.dir)

  abs.anoms <- proj.values - past.values
  prc.anoms <- (proj.values - past.values)/past.values*100    
  rv <- rbind(past.values,proj.values,abs.anoms,prc.anoms)
  return(rv)
}

##--------------------------------------------------------------------------------------

return_period_regional_values <- function(var.name,area.info,
                                          read.dir,tmp.dir,
                                          clip.shp,
                                          past.int,deg) {                                  

  ##GCM and Regional averages grouping monthly and seasonal values - needs pre-computed seasonal and monthly files

  rp.files <- list.files(path=paste0(read.dir,'return_periods'),pattern=paste0('^',var.name,'_'),full.name=TRUE)

  ##-------------------------------------------------
  past.rp.file <- rp.files[grep(past.int,rp.files)]
  past.values <- new_area_weighted_average(past.rp.file,clip.shp,area.info,tmp.dir)

  rp.files <- list.files(path=paste0(read.dir,'return_periods/degree_anomaly_rps'),pattern=paste0('^',var.name,'_'),full.name=TRUE)

  proj.rp.file <- rp.files[grep(deg,rp.files)]
  proj.values <- new_area_weighted_average(proj.rp.file,clip.shp,area.info,tmp.dir)

  abs.anoms <- proj.values - past.values
  prc.anoms <- (proj.values - past.values)/past.values*100    
  rv <- rbind(past.values,proj.values,abs.anoms,prc.anoms)
  return(rv)
}

##--------------------------------------------------------------------------------------

climdex_regional_values <- function(var.name,area.info,
                                    read.dir,tmp.dir,
                                    clip.shp,
                                    past.int,deg) {                                  

  ##GCM and Regional averages grouping monthly and seasonal values - needs pre-computed seasonal and monthly files
  clim.name <- strsplit(var.name,'_')[[1]][1]
  freq <- strsplit(var.name,'_')[[1]][2]
  my.writedir <- paste0(write.dir,clim.name,'/')

  cx.files <- list.files(path=paste0(read.dir,'climdex/climatologies/'),pattern=paste0('^',clim.name,'_annual_'),full.name=TRUE)
  if (length(cx.files)==0) { 
     print('Annual climdex')
     browser()
  }  

  ##-------------------------------------------------
  past.cx.file <- cx.files[grep(past.int,cx.files)]
  if (length(past.cx.file) > 1) {
     print('past climdex')
     browser()
  }
  ###past.ann.values <- area_weighted_average(past.cx.file,clip.shp,tmp.dir)
  past.ann.values <- new_area_weighted_average(past.cx.file,clip.shp,area.info,tmp.dir)

  cx.files <- list.files(path=paste0(read.dir,'climdex/degree_anomaly_climatologies/'),pattern=paste0('^',clim.name,'_annual_'),full.name=TRUE)

  proj.cx.file <- cx.files[grep(deg,cx.files)]
  if (length(proj.cx.file) > 1) {
     print('proj climdex') 
     browser()
  }
  ###proj.ann.values <- area_weighted_average(proj.cx.file,clip.shp,tmp.dir)
  proj.ann.values <- new_area_weighted_average(proj.cx.file,clip.shp,area.info,tmp.dir)

  past.values <- past.ann.values
  proj.values <- proj.ann.values

  if (freq=='seasonal') {

     seas.files <- list.files(path=paste0(read.dir,'climdex/climatologies/'),pattern=paste0('^',clim.name,'_seasonal_'),full.name=TRUE)
     if (length(cx.files)==0) { browser()}  

     ##-------------------------------------------------
     past.seas.file <- seas.files[grep(past.int,seas.files)]
     if (length(past.seas.file) > 1) {browser()}
     ###past.seas.avg <- area_weighted_average(past.seas.file,clip.shp,tmp.dir)
     past.seas.avg <- new_area_weighted_average(past.seas.file,clip.shp,area.info,tmp.dir)

     seas.files <- list.files(path=paste0(read.dir,'climdex/degree_anomaly_climatologies/'),pattern=paste0('^',clim.name,'_seasonal_'),full.name=TRUE)

     proj.seas.file <- seas.files[grep(deg,seas.files)]
     if (length(proj.seas.file) > 1) {browser()}
     ###proj.seas.avg <- area_weighted_average(proj.seas.file,clip.shp,tmp.dir)
     proj.seas.avg <- new_area_weighted_average(proj.seas.file,clip.shp,area.info,tmp.dir)

     past.values <- c(past.seas.avg,past.ann.values)
     proj.values <- c(proj.seas.avg,proj.ann.values)
  }

  abs.anoms <- proj.values - past.values
  prc.anoms <- (proj.values - past.values)/past.values*100    
  rv <- rbind(past.values,proj.values,abs.anoms,prc.anoms)
  return(rv)
}


##-------------------------------------------------

load_area_info <- function(region,model,var.name,load.dir) {

  pr.vars <- c('pr',
               'prcptotETCCDI','sdiiETCCDI',
               'r1mmETCCDI','r10mmETCCDI','r20mmETCCDI',
               'rx1dayETCCDI','rx2dayETCCDI','rx5dayETCCDI',
               'r95pETCCDI','r95daysETCCDI','r99pETCCDI','r99daysETCCDI',
               'pr_RP5','pr_RP20','pr_RP50',
               'cddETCCDI','cdd90ETCCDI','cddmaxETCCDI','cwdETCCDI')
   var.pr <- FALSE
   var.pr <- any(pr.vars %in% var.name)

   if (model=='KIOST-ESM' & var.pr) {
      load.file <- paste0(load.dir,region,'_KIOST-ESM_pr_area_info.RData')
   } else {
      load.file <- paste0(load.dir,region,'_',model,'_area_info.RData')
   }

   load(load.file)
   return(area.info)
}

##*********************************************************************
##*********************************************************************

make_seas_mon_tables <- function(model.list,var.list,regional_values_function,type,
                                 cmip,region,region.title,scenario,clip.shp,
                                 past.int,proj.int,
                                 read.dir,write.dir,area.dir,tmp.dir) {
   
  ##Climate parameters
  for (var.name in var.list) {
    print('------------------------------')
    print(var.name)    

    rd <- get_rounding_value(var.name)
       
    col.val <- switch(type,
                      t_and_p=17,
                      seasonal=5) 

    if (grepl('ETCCDI',var.name)) {
       clim.name <- strsplit(var.name,'_')[[1]][1]
       my.writedir <- paste0(write.dir,clim.name,'/')
    } else {
       my.writedir <- write.dir
    }

    region.avgs <- vector(mode='list',length=length(model.list))
    gcm.list <- run.list <- rep('',length(model.list))

    for (m in seq_along(model.list)) {

       model <- model.list[m]
       print(model)
       model.info <- strsplit(model,'_')[[1]]
       gcm.list[m] <- model.info[1]
       run.list[m] <- model.info[3]

       model.dir <- paste0(read.dir,model,'/')

       area.info <- load_area_info(region,model.info[1],var.name,area.dir)

       region.avgs[[m]] <- regional_values_function(var.name=var.name,area.info=area.info,
                                                  read.dir=model.dir,tmp.dir=tmp.dir,
                                                  clip.shp=clip.shp,
                                                  past.int,proj.int)                                                            
    }

    if (!file.exists(my.writedir))
      dir.create(my.writedir,recursive=TRUE)

    past.mon.vals <- matrix(unlist(lapply(region.avgs,function(x) {return(x[1,])})),nrow=length(region.avgs),ncol=col.val,byrow=TRUE)    
    past.table <- format_tables(past.mon.vals,gcm.list,run.list,rd,var.name,region.title,cmip,scenario,type) 
    write.table(past.table,file=paste(my.writedir,'past.',var.name,'.',scenario,'.',cmip,'.',past.int,'.csv',sep=''),sep=',',quote=FALSE,col.name=FALSE,row.name=FALSE)

    future.mon.vals <- matrix(unlist(lapply(region.avgs,function(x) {return(x[2,])})),nrow=length(region.avgs),ncol=col.val,byrow=TRUE)
    future.table <- format_tables(future.mon.vals,gcm.list,run.list,rd,var.name,region.title,cmip,scenario,type)
    write.table(future.table,file=paste(my.writedir,'future.',var.name,'.',scenario,'.',cmip,'.',proj.int,'.csv',sep=''),sep=',',quote=FALSE,col.name=FALSE,row.name=FALSE)
      
    abs.mon.vals <- matrix(unlist(lapply(region.avgs,function(x) {return(x[3,])})),nrow=length(region.avgs),ncol=col.val,byrow=TRUE)
    abs.table <- format_tables(abs.mon.vals,gcm.list,run.list,rd,var.name,region.title,cmip,scenario,type)
    write.table(abs.table,file=paste(my.writedir,'abs.anomalies.',var.name,'.',scenario,'.',cmip,'.',proj.int,'.csv',sep=''),sep=',',quote=FALSE,col.name=FALSE,row.name=FALSE)
      
    prc.mon.vals <- matrix(unlist(lapply(region.avgs,function(x) {return(x[4,])})),nrow=length(region.avgs),ncol=col.val,byrow=TRUE)
    prc.table <- format_tables(prc.mon.vals,gcm.list,run.list,1,var.name,region.title,cmip,scenario,type)
    write.table(prc.table,file=paste(my.writedir,'percent.anomalies.',var.name,'.',scenario,'.',cmip,'.',proj.int,'.csv',sep=''),
                  sep=',',quote=FALSE,col.name=FALSE,row.name=FALSE)

  }
}

##---------------------------------------------------------------------------

make_annual_tables <- function(model.list,var.list,regional_values_function,
                               cmip,region,region.title,scenario,clip.shp,
                               past.int,proj.int,
                               read.dir,write.dir,area.dir,tmp.dir) {
   
  ##Climate parameters
  for (var.name in var.list) {
    print('------------------------------')
    print(var.name)    

    rd <- get_rounding_value(var.name)    
    col.val <- 1 
    type <- 'annual'

    if (grepl('ETCCDI',var.name)) {
       clim.name <- strsplit(var.name,'_')[[1]][1]
       my.writedir <- paste0(write.dir,clim.name,'/')
    } else {
       my.writedir <- write.dir
    }

    region.avgs <- vector(mode='list',length=length(model.list))
    gcm.list <- run.list <- rep('',length(model.list))

    for (m in seq_along(model.list)) {

       model <- model.list[m]
       print(model)
       model.info <- strsplit(model,'_')[[1]]
       gcm.list[m] <- model.info[1]
       run.list[m] <- model.info[3]

       model.dir <- paste0(read.dir,model,'/')

       area.info <- load_area_info(region,model.info[1],var.name,area.dir)

       region.avgs[[m]] <- regional_values_function(var.name=var.name,area.info=area.info,
                                                    read.dir=model.dir,tmp.dir=tmp.dir,
                                                    clip.shp=clip.shp,
                                                    past.int,proj.int)                                                            
    }

    if (!file.exists(my.writedir))
      dir.create(my.writedir,recursive=TRUE)

    past.mon.vals <- matrix(unlist(lapply(region.avgs,function(x) {return(x[1])})),nrow=length(region.avgs),ncol=col.val,byrow=TRUE)    
    past.table <- format_tables(past.mon.vals,gcm.list,run.list,rd,var.name,region.title,cmip,scenario,type) 
    write.table(past.table,file=paste(my.writedir,'past.',var.name,'.',scenario,'.',cmip,'.',past.int,'.csv',sep=''),sep=',',quote=FALSE,col.name=FALSE,row.name=FALSE)

    future.mon.vals <- matrix(unlist(lapply(region.avgs,function(x) {return(x[2])})),nrow=length(region.avgs),ncol=col.val,byrow=TRUE)
    future.table <- format_tables(future.mon.vals,gcm.list,run.list,rd,var.name,region.title,cmip,scenario,type)
    write.table(future.table,file=paste(my.writedir,'future.',var.name,'.',scenario,'.',cmip,'.',proj.int,'.csv',sep=''),sep=',',quote=FALSE,col.name=FALSE,row.name=FALSE)
      
    abs.mon.vals <- matrix(unlist(lapply(region.avgs,function(x) {return(x[3])})),nrow=length(region.avgs),ncol=col.val,byrow=TRUE)
    abs.table <- format_tables(abs.mon.vals,gcm.list,run.list,rd,var.name,region.title,cmip,scenario,type)
    write.table(abs.table,file=paste(my.writedir,'abs.anomalies.',var.name,'.',scenario,'.',cmip,'.',proj.int,'.csv',sep=''),sep=',',quote=FALSE,col.name=FALSE,row.name=FALSE)
      
    prc.mon.vals <- matrix(unlist(lapply(region.avgs,function(x) {return(x[4])})),nrow=length(region.avgs),ncol=col.val,byrow=TRUE)
    prc.table <- format_tables(prc.mon.vals,gcm.list,run.list,1,var.name,region.title,cmip,scenario,type)
    write.table(prc.table,file=paste(my.writedir,'percent.anomalies.',var.name,'.',scenario,'.',cmip,'.',proj.int,'.csv',sep=''),
                  sep=',',quote=FALSE,col.name=FALSE,row.name=FALSE)

  }
}



##---------------------------------------------------------------------
##*********************************************************************

testing <- TRUE

if (testing) {

   ##------------------------------------
   ##Testing
   cmip <- "CMIP5"
   scenario <- "rcp85"

   regions <- "canada_boundary" ##c("canada_boundary","bc","prairies","ontario","atlantic_canada",
               ##"northern_canada") ##"quebec",

   tmpdir <- "/local_temp/ssobie/cmip_regions"
   var.class <- 't_and_p'
   ##------------------------------------
} else {

   args <- commandArgs(trailingOnly=TRUE)
   for(i in 1:length(args)){
       eval(parse(text=args[[i]]))
   }
   var.class <- type
}

for (region in regions) {

region.title <- toupper(region)

tmp.dir <- paste0(tmpdir,'/tables_degrees_',cmip,'_',scenario,'_',region,'/')


if (!file.exists(tmp.dir)) {
   dir.create(tmp.dir,recursive=TRUE)
}


if (cmip=='CMIP6') {
   read.dir <- "/storage/data/climate/CMIP6/Derived/"
}
if (cmip=='CMIP5') {
   read.dir <- "/storage/data/climate/CMIP5/daily/Derived/"
}

area.dir <- paste0("/storage/data/projects/rci/data/cas/",tolower(cmip),"/area_info/",region,"/")

shape.dir <- "/storage/data/projects/rci/data/assessments/shapefiles/canada_regions/"
clip.shp <- get_shape_file(region,shape.dir)

write.base <- paste0("/storage/data/projects/rci/data/cas/",tolower(cmip),"/degree_tables/",region,"/")

scenario.models <- list.files(read.dir,pattern=scenario)

cmip5.list <- list(c('ACCESS1-0','r1i1p1'),
                   c('bcc-csm1-1','r1i1p1'),
                   c('bcc-csm1-1-m','r1i1p1'),
                   c('BNU-ESM','r1i1p1'),
                   c('CanESM2','r1i1p1'),
                   c('CCSM4','r2i1p1'),
                   c('CNRM-CM5','r1i1p1'),
                   c('CSIRO-Mk3-6-0','r1i1p1'),
                   c('FGOALS-g2','r1i1p1'),
                   c('GFDL-CM3','r1i1p1'),
                   c('GFDL-ESM2G','r1i1p1'),
                   c('GFDL-ESM2M','r1i1p1'),
                   c('HadGEM2-AO','r1i1p1'),
                   c('HadGEM2-CC','r1i1p1'),
                   c('HadGEM2-ES','r1i1p1'),
                   c('inmcm4','r1i1p1'),
                   c('IPSL-CM5A-LR','r1i1p1'),
                   c('IPSL-CM5A-MR','r1i1p1'),
                   c('MIROC-ESM','r1i1p1'),
                   c('MIROC-ESM-CHEM','r1i1p1'),
                   c('MIROC5','r3i1p1'),
                   c('MPI-ESM-LR','r3i1p1'),
                   c('MPI-ESM-MR','r1i1p1'),
                   c('MRI-CGCM3','r1i1p1'),
                   c('NorESM1-M','r1i1p1'))
cmip5.list <- list(c('ACCESS1-0','r1i1p1'),
                   c('bcc-csm1-1','r1i1p1'),
                   c('bcc-csm1-1-m','r1i1p1'))


cmip6.list <- list(c('ACCESS-CM2','r1i1p1f1'),
                   c('ACCESS-ESM1-5','r1i1p1f1'),
                   c('BCC-CSM2-MR','r1i1p1f1'),
                   c('CanESM5','r1i1p2f1'),
                   c('CNRM-CM6-1','r1i1p1f2'),
                   c('CNRM-ESM2-1','r1i1p1f2'),
                   c('EC-Earth3','r4i1p1f1'),
                   c('EC-Earth3-Veg','r1i1p1f1'),
                   c('FGOALS-g3','r1i1p1f1'),
                   c('GFDL-ESM4','r1i1p1f1'),
                   c('HadGEM3-GC31-LL','r1i1p1f3'),
                   c('INM-CM4-8','r1i1p1f1'),
                   c('INM-CM5-0','r1i1p1f1'),
                   c('IPSL-CM6A-LR','r1i1p1f1'),
                   c('KACE-1-0-G','r2i1p1f1'),
                   c('KIOST-ESM','r1i1p1f1'),
                   c('MIROC6','r1i1p1f1'),
                   c('MIROC-ES2L','r1i1p1f2'),
                   c('MPI-ESM1-2-HR','r1i1p1f1'),
                   c('MPI-ESM1-2-LR','r1i1p1f1'),
                   c('MRI-ESM2-0','r1i1p1f1'),
                   c('NESM3','r1i1p1f1'),
                   c('NorESM2-LM','r1i1p1f1'),
                   c('NorESM2-MM','r1i1p1f1'),
                   c('UKESM1-0-LL','r1i1p1f2'))

##cmip6.list <- list(c('ACCESS-CM2','r1i1p1f1'),
##                   c('ACCESS-ESM1-5','r1i1p1f1'),
##                   c('BCC-CSM2-MR','r1i1p1f1'))



model.dirs <- switch(cmip,
                     CMIP6=cmip6.list,
                     CMIP5=cmip5.list) 


if (scenario=='rcp26') {
   omx <- c(grep('ACCESS1-0',model.dirs),
            grep('HadGEM2-CC',model.dirs),
            grep('inmcm4',model.dirs))
   model.dirs <- model.dirs[-omx]
}

model.list <- unlist(lapply(model.dirs,function(x,y){paste(x[1],y,x[2],sep='_')},scenario))

past.int <- "1971-2000"

degrees <- c('one_deg','two_deg','three_deg')

for (degree in degrees) {

##Degree Days
if (var.class=='degree_days') {
   var.list <- c('cdd','fdd','gdd','hdd')
   write.dir <- paste0(write.base,scenario,'/degree_days/')
   regional_values_function <- degree_day_regional_values
   make_annual_tables(model.list,var.list,regional_values_function,
                                 cmip,region,region.title,scenario,clip.shp,
                                 past.int,degree,
                                 read.dir,write.dir,area.dir,tmp.dir)
}

##Temperature and Precipitation
if (var.class=='t_and_p') {
   var.list <- c('pr','tasmax','tasmin')
   write.dir <- paste0(write.base,scenario,'/t_and_p/')
   regional_values_function <- t_and_p_regional_values
   make_seas_mon_tables(model.list,var.list,regional_values_function,type='t_and_p',
                    cmip,region,region.title,scenario,clip.shp,
                    past.int,degree,
                    read.dir,write.dir,area.dir,tmp.dir)
}


if (var.class=='return_periods') {
  ##Return Periods
  var.list <- c('pr_RP5','pr_RP20','pr_RP50',
             'tasmax_RP5','tasmax_RP20',
             'tasmin_RP5','tasmin_RP20')

  write.dir <- paste0(write.base,scenario,'/return_periods/')
  regional_values_function <- return_period_regional_values
  make_annual_tables(model.list,var.list,regional_values_function,
                               cmip,region,region.title,scenario,clip.shp,
                               past.int,degree,
                               read.dir,write.dir,area.dir,tmp.dir)
}

if (var.class=='climdex') {
##Climdex
  ann.var.list <- c('idETCCDI','trETCCDI','fdETCCDI',
                  'suETCCDI','su30ETCCDI','gslETCCDI',
                  'sdiiETCCDI','r1mmETCCDI','r10mmETCCDI','r20mmETCCDI','cddETCCDI','cwdETCCDI',
                  'r95pETCCDI','r99pETCCDI','r95daysETCCDI','r99daysETCCDI',
                  'prcptotETCCDI')

  write.dir <- paste0(write.base,scenario,'/climdex/')
  regional_values_function <- climdex_regional_values 
  make_annual_tables(model.list,paste0(ann.var.list,'_annual'),regional_values_function,
                               cmip,region,region.title,scenario,clip.shp,
                               past.int,degree,
                               read.dir,write.dir,area.dir,tmp.dir)

  seas.var.list <- c('txxETCCDI','tnxETCCDI','txnETCCDI','tnnETCCDI','dtrETCCDI',
                   'rx1dayETCCDI','rx2dayETCCDI','rx5dayETCCDI')
  regional_values_function <- climdex_regional_values
  make_seas_mon_tables(model.list,paste0(seas.var.list,'_seasonal'),regional_values_function,
                     type='seasonal',
                     cmip,region,region.title,scenario,clip.shp,
                     past.int,degree,
                     read.dir,write.dir,area.dir,tmp.dir)

}


}

}