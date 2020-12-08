##Script to sort out the data requested for Okanagan

source('/storage/home/ssobie/code/repos/assessments/summary.table.comments.r',chdir=T)

library(openxlsx)

##------------------------------------------------------------- 
##------------------------------------------------------------- 

get_units_pane <- function(var.name) {
  leg.label <- NA
  if (grepl("(tas|txx|tnn|txn|tnx|tmax|tmin|dtr)", var.name))
    leg.label <- c(rep('degC',7),rep('%',6))
  if (grepl("(pr|rx|r10|r20|r9|RP|sdii|prcptot)", var.name))
    leg.label <- c(rep('mm',7),rep('%',6))
  if (grepl("(dd)", var.name))
    leg.label <- c(rep('degree days',7),rep('%',6))
  if (grepl("(fdE|cddE|cdd90|cddmax|cwd|su|gsl|id|trE|su30|r95daysE|r99daysE)", var.name))
    leg.label <- c(rep('days',7),rep('%',6))
  if (grepl("(dtr)", var.name))
    leg.label <- c(rep('degC',7),rep('%',6))
  if (grepl("(r95sep|r95dist|r95days|r99days)", var.name))
    leg.label <- c(rep('days',7),rep('%',6))
  return(leg.label)
} 

##------------------------------------------------------------- 
get_round_val <- function(var.name) {
  rd <- 0
  if (grepl("(dd)", var.name))
    rd <- 0    
  if (grepl("(tas|txx|tnn|tnx|txn|tmax|tmin)", var.name))
    rd <- 1
  if (grepl("(pr|rx|r9|RP|rp|tx90|tn10|trE|cddE|cdd90|cddmax|cwdE|idE|dtrE|wsdiE|csdiE|r95sep)", var.name))
    rd <- 0
  if (grepl("(pas|snowdepth)", var.name))
    rd <- 0
  return(rd)
} 

##------------------------------------------------------------- 
##Separate variables into precip and temperature
## 'E' denotes a climdex index
##Set the order that the variables should be arranged here.

filter_input_variables <- function(table.vars) {
  all.vars <- unlist(lapply(table.vars,function(x){return(x[1])}))
  pr.vars <- c('pr',
               'prcptotETCCDI','sdiiETCCDI','r10mmETCCDI','r20mmETCCDI',
               'rx1dayETCCDI','rx2dayETCCDI','rx5dayETCCDI',
               'r95pETCCDI','r95daysETCCDI','r99pETCCDI','r99daysETCCDI',
               'pr.maximum','pr.minimum','pr.standard_deviation',
               'pr_rp5','pr_rp20','pr_rp50',
               'cddETCCDI','cdd90ETCCDI','cddmaxETCCDI','cwdETCCDI') 
  pr.ix <- match(pr.vars,all.vars)
  pr.selected <- table.vars[pr.ix[!is.na(pr.ix)]]

  tas.vars <- c('tasmax','tas','tasmin',
                'txxETCCDI','tnnETCCDI','txnETCCDI','tnxETCCDI',
                'dtrETCCDI','suETCCDI','su30ETCCDI','trETCCDI',
                'idETCCDI','fdETCCDI','gslETCCDI',
                'cdd','gdd','hdd','fdd',
                'tasmax.annual_quantile_975',
                'tasmax.annual_quantile_990',       
                'tasmax.annual_quantile_996',                        
                'tasmin.annual_quantile_004',
                'tasmin.annual_quantile_010',
                'tasmin.annual_quantile_025',
                'tasmax_rp5','tasmax_rp20',
                'tasmin_rp5','tasmin_rp20')                          
  tas.ix <- match(tas.vars,all.vars)
  tas.selected <- table.vars[tas.ix[!is.na(tas.ix)]]
  rv <- list(pr=pr.selected,tas=tas.selected)
  return(rv)
}

##------------------------------------------------------------- 

find_row_locations <- function(sorted.vars) {
                   
   pr.vars <- sorted.vars$pr
   pr.rows <- vector(length=length(pr.vars)+1,mode='list')

   pr.rows[[1]] <- 3:5 ##Header Rows
   for (i in 1:length(pr.vars)) {
      row.len <- switch(pr.vars[[i]][[2]],annual=2,seasonal=6)
      pr.rows[[i+1]] <- seq(tail(pr.rows[[i]],1)+1,length.out=row.len)
   }

   tas.start <- tail(pr.rows[[length(pr.vars)+1]],1)
   tas.vars <- sorted.vars$tas
   tas.rows <- vector(length=length(tas.vars)+1,mode='list')
   tas.rows[[1]] <- (1:3)+tas.start
   for (i in 1:length(tas.vars)) {
      row.len <- switch(tas.vars[[i]][[2]],annual=2,seasonal=6)
      tas.rows[[i+1]] <- seq(tail(tas.rows[[i]],1)+1,length.out=row.len)
   }

   rv <- list(pr=pr.rows,tas=tas.rows)
   return(rv)

}
##---------------------------------------------

get_scenarios <- function(var.name,seas,interval,scenario,gcm.list,cmip,type,read.dir,rp=NULL) {
  

  file.name <- paste(type,'.',var.name,'.',scenario,'.',toupper(cmip),'.',interval,'.csv',sep='')

  ##Return Periods  
  if (!is.null(rp)) {
    file.name <- paste(type,'.',var.name,'_RP',rp,'.',scenario,".",toupper(cmip),".",interval,'.csv',sep='')
    file.scen <- as.matrix(read.csv(paste(read.dir,'return_periods/',file.name,sep=''),
                                    header=TRUE,as.is=TRUE))

  } else if (grepl('ETCCDI',var.name)) {

    type.files <- list.files(path=paste0(read.dir,'climdex/',var.name,'/'),pattern=type)
    file.name <- type.files[grep(interval,type.files)]

    ##seas.name <- switch(seas,Winter='seasonal',Spring='seasonal',Summer='seasonal',Fall='seasonal',Annual='annual')
    ##file.name <- paste(type,'.',var.name,'_',tolower(seas.name),'.',scenario,'.',toupper(cmip),'.',interval,'.csv',sep='')
    file.scen <- as.matrix(read.csv(paste(read.dir,'climdex/',var.name,'/',file.name,sep=''),header=TRUE,as.is=TRUE))

  } else if (grepl('(^cdd$|hdd|gdd|fdd|pas)',var.name)) {
    file.scen <- as.matrix(read.csv(paste(read.dir,'degree_days/',file.name,sep=''),header=TRUE,as.is=TRUE))

  } else if (grepl('^pr$|^tasmax$|^tasmin$',var.name)) {   
    file.scen <- as.matrix(read.csv(paste(read.dir,'t_and_p/',file.name,sep=''),header=TRUE,as.is=TRUE))
  } else {
    file.scen <- as.matrix(read.csv(paste(read.dir,file.name,sep=''),header=TRUE,as.is=TRUE))
  }

  models <- file.scen[,1]
  mod.sub <- models %in% gcm.list
  data.names <- file.scen[1,]
  seas.sub <- data.names %in% seas
  data.scen <- file.scen[which(mod.sub),which(seas.sub)]
  return(data.scen)
}
 
get_data <- function(var.name,seas,interval,scenario,gcm.list,cmip,type,read.dir,rp=NULL) {
  data.rcp <- get_scenarios(var.name,seas,interval,scenario,gcm.list,cmip,type,read.dir,rp)
  data.comb <- as.numeric(data.rcp)
  data.stats <- c(mean(data.comb,na.rm=T,finite=T),
                  quantile(data.comb,0.1,na.rm=T,finite=T),
                  quantile(data.comb,0.9,na.rm=T,finite=T))
  names(data.stats) <- c('avg','10%','90%')
  return(data.stats)
}
 
get_seasonal_data <- function(var.name,scenario,gcm.list,cmip,read.dir) {
  seasons <- c('Winter','Spring','Summer','Fall','Annual')
  seas.values <- matrix(0,nrow=5,ncol=13)
  rd <- get_round_val(var.name)
  for (s in seq_along(seasons)) {
    seas <- seasons[s]
    vals <- c(round(get_data(var.name,seas,'1971-2000',scenario,gcm.list,cmip,'past',read.dir)[1],rd),
              round(get_data(var.name,seas,'2011-2040',scenario,gcm.list,cmip,'abs.anomalies',read.dir),rd),
              round(get_data(var.name,seas,'2041-2070',scenario,gcm.list,cmip,'abs.anomalies',read.dir),rd),
              round(get_data(var.name,seas,'2071-2100',scenario,gcm.list,cmip,'abs.anomalies',read.dir),rd),              
              round(get_data(var.name,seas,'2011-2040',scenario,gcm.list,cmip,'percent.anomalies',read.dir),rd),
              round(get_data(var.name,seas,'2041-2070',scenario,gcm.list,cmip,'percent.anomalies',read.dir),rd),
              round(get_data(var.name,seas,'2071-2100',scenario,gcm.list,cmip,'percent.anomalies',read.dir),rd))
    lower <- vals[names(vals)=='10%']
    upper <- vals[names(vals)=='90%']
    brackets <- paste('(',paste(lower,upper,sep=' to '),')',sep='')
    avgs <- vals[names(vals)=='avg']
    proj <- avgs[-1]
    len  <- length(c(proj,brackets))
    result <- rep(0,length=len)
    result[seq(1,len,2)] <- proj
    result[seq(2,len,2)] <- brackets
    result <- c(avgs[1],result)
    
    seas.values[s,] <- result
  }
  seas.data <- cbind(seasons,seas.values)  
  return(seas.data)
}


get_annual_data <- function(var.name,scenario,gcm.list,cmip,read.dir,rp) {

  seas <- 'Annual'
  rd <- get_round_val(var.name)
  vals <- c(round(get_data(var.name,seas,'1971-2000',scenario,gcm.list,cmip,'past',read.dir,rp)[1],rd),
            round(get_data(var.name,seas,'2011-2040',scenario,gcm.list,cmip,'abs.anomalies',read.dir,rp),rd),
            round(get_data(var.name,seas,'2041-2070',scenario,gcm.list,cmip,'abs.anomalies',read.dir,rp),rd),
            round(get_data(var.name,seas,'2071-2100',scenario,gcm.list,cmip,'abs.anomalies',read.dir,rp),rd),              
            round(get_data(var.name,seas,'2011-2040',scenario,gcm.list,cmip,'percent.anomalies',read.dir,rp),rd),
            round(get_data(var.name,seas,'2041-2070',scenario,gcm.list,cmip,'percent.anomalies',read.dir,rp),rd),
            round(get_data(var.name,seas,'2071-2100',scenario,gcm.list,cmip,'percent.anomalies',read.dir,rp),rd))
  lower <- vals[names(vals)=='10%']
  upper <- vals[names(vals)=='90%']
  brackets <- paste('(',paste(lower,upper,sep=' to '),')',sep='')
  avgs <- vals[names(vals)=='avg']
  proj <- avgs[-1]
  len  <- length(c(proj,brackets))
  result <- rep(0,length=len)
  result[seq(1,len,2)] <- proj
  result[seq(2,len,2)] <- brackets
  result <- c('Annual',avgs[1],result)
  return(result)
}

seasonal_table <- function(var.name,scenario,gcm.list,cmip,read.dir,rp=NULL) {

  no.percent <- '(tas|tasmax|tasmin|txxETCCDI|tnnETCCD|tnxETCCDI|txnETCCDI|trETCCDI|suETCCDI|su30ETCCDI)'
  result <- get_seasonal_data(var.name,scenario,gcm.list,cmip,read.dir)
  if (grepl(no.percent,var.name))
    result[,9:14] <- 'NA'
  return(as.data.frame(result))
}


annual_table <- function(var.name,scenario,gcm.list,cmip,read.dir,rp=NULL) {
  
  no.percent <- '(tas|tasmax|tasmin|txxETCCDI|tnnETCCD|trETCCDI|r95sep|r99days|r95days)'
  result <- get_annual_data(var.name,scenario,gcm.list,cmip,read.dir,rp) 

  if (grepl(no.percent,var.name))
    result[9:14] <- 'NA'
  return(as.data.frame(t(result)))
}


##----------------------------------------------------------------------------------------
write_variables <- function(wb,sorted.vars,row.locs,type,gcm.list,scenario,cmip,read.dir) {

  len <- length(sorted.vars)
  for (i in 1:len) {               
    current.var <- sorted.vars[[i]]
    var.name <- current.var[1]
    print(var.name)
    season <- current.var[2]
    var.title <- current.var[3]
    print(var.title)
    current.row <- row.locs[[i+1]]

    s1 <- createStyle(fontSize = 12, fontColour = "black", textDecoration = c("BOLD"))
    s2 <- createStyle(fontSize = 12, fontColour = "black")
    c1 <- createComment(comment = variable.comment(var.name),style=c(s1,s2),visible=FALSE)
    writeComment(wb, 1, col = 1, row = current.row[1], comment = c1)
    seas.fx <- switch(season,   
                      annual=annual_table,   
                      seasonal=seasonal_table)
    if (grepl('rp',var.name)) {
      rp <- as.numeric(gsub('rp','',strsplit(var.name,'_')[[1]][2]))
      rpvar <- gsub('rp','',strsplit(var.name,'_')[[1]][1])
      var.entry <- seas.fx(rpvar,scenario,gcm.list,cmip,read.dir,rp=rp)
    } else {                      
      var.entry <- seas.fx(var.name,scenario,gcm.list,cmip,read.dir)                      
    }

    pane.colour <- switch(type,pr='lightblue',tas='tan1')             

    hdstyle <- createStyle(fgFill = pane.colour, halign = "CENTER", textDecoration = "Bold",
                             border = "TopBottomLeftRight", fontColour = "black")  
    units.pane <- c(var.title,get_units_pane(var.name))                             
    writeData(wb, sheet=1, as.data.frame(t(units.pane)), startRow = current.row[1], startCol = 1, headerStyle = hdstyle,
              borders = "rows", borderStyle = "medium",colNames=FALSE)
    addStyle(wb,sheet=1,hdstyle,rows=current.row[1],cols=1:14,gridExpand=FALSE,stack=FALSE)

    datastyle <- createStyle(fgFill = 'white', halign = "RIGHT",
                             border = "TopBottomLeftRight", fontColour = "black")                              
    writeData(wb, sheet=1, var.entry, startRow = current.row[2], startCol = 1, headerStyle = hdstyle,
              borders = "rows", borderStyle = "medium",colNames=FALSE)
    addStyle(wb,sheet=1,datastyle,rows=current.row[-1],cols=1:14,gridExpand=TRUE,stack=FALSE)

    highlight <- createStyle(fgFill = 'lightyellow', halign = "CENTER", 
                             border = "TopBottomLeftRight", fontColour = "black")  
    addStyle(wb,sheet=1,highlight,rows=current.row[-1],cols=2,gridExpand=FALSE,stack=FALSE)
    addStyle(wb,sheet=1,highlight,rows=current.row[-1],cols=5,gridExpand=FALSE,stack=FALSE)
    addStyle(wb,sheet=1,highlight,rows=current.row[-1],cols=6,gridExpand=FALSE,stack=FALSE)
    addStyle(wb,sheet=1,highlight,rows=current.row[-1],cols=11,gridExpand=FALSE,stack=FALSE)
    addStyle(wb,sheet=1,highlight,rows=current.row[-1],cols=12,gridExpand=FALSE,stack=FALSE)

    if (grepl('(suETCCDI|su30ETCCDI)',var.name) | var.name =='cdd') {
        prctstyle <- createStyle(fgFill = 'lightgray', halign = "RIGHT",
                                 border = "TopBottomLeftRight", fontColour = "black")                              
        writeData(wb, sheet=1, var.entry[9:14], startRow = current.row[2], startCol = 9, headerStyle = prctstyle,
                  borders = "rows", borderStyle = "medium",colNames=FALSE)
        addStyle(wb,sheet=1,prctstyle,rows=current.row[-1],cols=9:14,gridExpand=TRUE,stack=FALSE)
        s1 <- createStyle(fontSize = 12, fontColour = "black", textDecoration = c("BOLD"))
        s2 <- createStyle(fontSize = 12, fontColour = "black")
        f1 <- createComment(comment = c('CAUTION\n','Percent changes from a low baseline value can result in deceptively large percent values'),
                            style=c(s1,s2),visible=FALSE)
        for (j in 9:14) {                     
          writeComment(wb, 1, col = j, row = current.row[1]+1, comment = f1)
        }

    }


  }
     
}
##----------------------------------------------------------------------------------------
      ##Top Frozen Pane

create_frozen_top_pane <- function(wb) {

      pane.titles <- list(' ','Past','2020s Change',' ',
                              '2050s Change',' ',
                              '2080s Change',' ',
                              '2020s Percent Change',' ',
                              '2050s Percent Change',' ',
                              '2080s Percent Change',' ')

      fz1 <- createStyle(fgFill = "gray94", halign = "CENTER", textDecoration = "Bold",
                         border = "TopBottomLeftRight", fontColour = "black")
      writeData(wb, sheet=1, pane.titles, startRow = 1, startCol = 1, headerStyle = fz1,
                borders = "rows", borderStyle = "medium")
      addStyle(wb,sheet=1,fz1,rows=1,cols=1:14,gridExpand=FALSE,stack=FALSE)
      mergeCells(wb,sheet=1,cols=c(3,4),rows=1)
      mergeCells(wb,sheet=1,cols=c(5,6),rows=1)
      mergeCells(wb,sheet=1,cols=c(7,8),rows=1)
      mergeCells(wb,sheet=1,cols=c(9,10),rows=1)
      mergeCells(wb,sheet=1,cols=c(11,12),rows=1)
      mergeCells(wb,sheet=1,cols=c(13,14),rows=1)
      ##freezePane(wb,sheet=1,firstRow=TRUE)

      prct.header <- list(' ',' ','Average','10th%-90th%','Average','10th%-90th%','Average','10th%-90th%',
                                  'Average','10th%-90th%','Average','10th%-90th%','Average','10th%-90th%')          
      writeData(wb, sheet=1, prct.header, startRow = 2, startCol = 1, headerStyle = fz1,
                borders = "rows", borderStyle = "medium")
      addStyle(wb,sheet=1,fz1,rows=2,cols=1:14,gridExpand=FALSE,stack=FALSE)
      ##freezePane(wb,sheet=1,firstActiveRow=3)
}

##-----------------------------------------------------------------------------------------
      ##Precipitation Header Rows
create_title_panes <- function(wb,var.name,start.row) {

      pane.titles <- list(' ','Past','2020s Change',' ',
                               '2050s Change',' ',
                              '2080s Change',' ',
                              '2020s Percent Change',' ',
                              '2050s Percent Change',' ',
                              '2080s Percent Change',' ')

      pane.colour <- switch(var.name,pr='lightblue',tas='tan1')             
      hdstyle <- createStyle(fgFill = pane.colour, halign = "CENTER", textDecoration = "Bold",
                             border = "TopBottomLeftRight", fontColour = "black") 
      writeData(wb, sheet=1, pane.titles, startRow = start.row, startCol = 1, headerStyle = hdstyle,
                borders = "rows", borderStyle = "medium")
      addStyle(wb,sheet=1,hdstyle,rows=start.row,cols=1:14,gridExpand=FALSE,stack=FALSE)
      mergeCells(wb,sheet=1,cols=c(3,4),rows=start.row)
      mergeCells(wb,sheet=1,cols=c(5,6),rows=start.row)
      mergeCells(wb,sheet=1,cols=c(7,8),rows=start.row)
      mergeCells(wb,sheet=1,cols=c(9,10),rows=start.row)
      mergeCells(wb,sheet=1,cols=c(11,12),rows=start.row)
      mergeCells(wb,sheet=1,cols=c(13,14),rows=start.row)
      prct.header <- list(' ',' ','Average','10th%-90th%','Average','10th%-90th%','Average','10th%-90th%',
                                  'Average','10th%-90th%','Average','10th%-90th%','Average','10th%-90th%')          
      writeData(wb, sheet=1, prct.header, startRow = start.row+1, startCol = 1, headerStyle = hdstyle,
                borders = "rows", borderStyle = "medium")
      addStyle(wb,sheet=1,hdstyle,rows=start.row+1,cols=1:14,gridExpand=FALSE,stack=FALSE)

      mergeCells(wb,sheet=1,cols=1:14,rows=start.row+2)
      pane.header <- list(switch(var.name,pr='Precipitation',tas='Temperature'))
      writeData(wb, sheet=1, pane.header, startRow = start.row+2, startCol = 1, headerStyle = hdstyle,
                borders = "rows", borderStyle = "medium")
      addStyle(wb,sheet=1,hdstyle,rows=start.row+2,cols=1:14,gridExpand=FALSE,stack=FALSE)
}


##*******************************************************************
##
##*******************************************************************

cmip6.models <- c("ACCESS-CM2","ACCESS-ESM1-5",
                  "BCC-CSM2-MR","CanESM5","CNRM-CM6-1","CNRM-ESM2-1",
                  "EC-Earth3","EC-Earth3-Veg","FGOALS-g3","GFDL-ESM4","HadGEM3-GC31-LL",
                  "INM-CM4-8","INM-CM5-0","IPSL-CM6A-LR","KACE-1-0-G","KIOST-ESM",
                  "MIROC6","MIROC-ES2L","MPI-ESM1-2-HR","MPI-ESM1-2-LR","MRI-ESM2-0",
                  "NESM3","NorESM2-LM","NorESM2-MM",
                  "UKESM1-0-LL")


cmip <- 'cmip6'
scenario <- 'ssp245'

gcm.list <- switch(cmip,cmip6=cmip6.models,cmip5=cmip5.models)

base.dir <- paste0("/storage/data/projects/rci/data/cas/",cmip,"/tables/")

reg.list <- 'prairies'

table.vars <- list(c('pr','seasonal','PR'),
                   c('rx1dayETCCDI','seasonal','RX1DAY'),
                   c('rx5dayETCCDI','seasonal','RX5DAY'),
                   c('r95pETCCDI','annual','R95P'),
                   c('r95daysETCCDI','annual','R95DAYS'),
                   c('r99pETCCDI','annual','R99P'),
                   c('r99daysETCCDI','annual','R99DAYS'),
                   c('pr_rp20','annual','RP20 PR'),
                   c('pr_rp5','annual','RP5 PR'),
                   c('pr_rp50','annual','RP50 PR'),
                   c('cddETCCDI','annual','CDD'),
                   c('cwdETCCDI','annual','CWD'),
                   c('tasmax','seasonal','TASMAX'),
                   c('tasmin','seasonal','TASMIN'),
                   c('txxETCCDI','seasonal','TXX'),
                   c('txnETCCDI','seasonal','TXN'),
                   c('tnnETCCDI','seasonal','TNN'),
                   c('tnxETCCDI','seasonal','TNX'),
                   c('dtrETCCDI','seasonal','DTR'),
                   c('suETCCDI','annual','SU'),
                   c('su30ETCCDI','annual','SU30'),
                   c('trETCCDI','annual','TR'),
                   c('idETCCDI','annual','ID'),
                   c('fdETCCDI','annual','FD'),
                   c('gslETCCDI','annual','GSL'),
                   c('cdd','annual','CDD'),
                   c('gdd','annual','GDD'),
                   c('hdd','annual','HDD'),
                   c('fdd','annual','FDD'),
                   c('tasmax_rp20','annual','RP20 TX'),
                   c('tasmin_rp20','annual','RP20 TN'),
                   c('tasmax_rp5','annual','RP5 TX'),
                   c('tasmin_rp5','annual','RP5 TN'))

##                   c('tas','seasonal','TAS'),

#table.vars <- list(c('pr','seasonal','PR'),
#                   c('pr_rp20','annual','RP20 PR'),
#                   c('pr_rp5','annual','RP5 PR'),
#                   c('pr_rp50','annual','RP50 PR'),
#                   c('tasmax','seasonal','TASMAX'),
#                   c('tasmin','seasonal','TASMIN'),
#                   c('tasmax_rp20','annual','RP20 TX'),
#                   c('tasmin_rp20','annual','RP20 TN'),
#                   c('tasmax_rp5','annual','RP5 TX'),
#                   c('tasmin_rp5','annual','RP5 TN'))



sorted.vars <- filter_input_variables(table.vars) 

row.locs <- find_row_locations(sorted.vars)

##------------------------------------------------------------------------

for (i in seq_along(reg.list)) {
  region.info <- reg.list[[i]]
  region <- region.info[1]
  readloc <- region.info[2]
  writeloc <- readloc
  read.dir <- paste0(base.dir,region,"/",scenario,"/")  


  ##Formatted Table
  wb <- createWorkbook()
  addWorksheet(wb, "Regional Averages")
  setColWidths(wb, sheet = 1, cols = 1:14, widths = 14) ##Set fixed width
  create_frozen_top_pane(wb)
  create_title_panes(wb,var.name='pr',start.row=row.locs$pr[[1]][1])
  write_variables(wb,sorted.vars$pr,row.locs$pr,'pr',gcm.list,scenario,cmip,read.dir)

  create_title_panes(wb,var.name='tas',start.row=row.locs$tas[[1]][1])
  write_variables(wb,sorted.vars$tas,row.locs$tas,'tas',gcm.list,scenario,cmip,read.dir)
  freezePane(wb,sheet=1,firstActiveCol=2,firstActiveRow=3)

  write.dir <- paste0(base.dir,region,"/")  
  saveWorkbook(wb, paste0(write.dir,region,"_",scenario,"_",cmip,"_variable_table_rcp85.xlsx"), overwrite = TRUE)

}

