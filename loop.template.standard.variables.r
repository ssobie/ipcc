
##---------------------------------------------------------------------------

base.dir <- '/storage/data/climate/CMIP5/daily/'

##---------------------------------------------------------------------------

var.name <- 'tasmin'
type <- 'annual_extremes'

gcm.list <- list.files(base.dir)
scenario.list <- c('rcp26','rcp45','rcp85')

##gcm.list <- 'EC-Earth'
##scenario.list <- 'rcp45'

##---------------------------------------------------------------------------


derived.template <- "/storage/home/ssobie/code/repos/ipcc/template.standard.variables.pbs"

for (gcm in gcm.list) {
   gcm.dir <- paste0(base.dir,gcm,'/')

   for (scenario in scenario.list) {
   
         scen.files <- list.files(path=gcm.dir,pattern=paste0("historical\\+",scenario,"_"))
         var.files <- scen.files[grep(paste0(var.name,"_day_",gcm),scen.files)]
         for (file in var.files) { 
            print('Submitting')
            print(file)

            file.split <- strsplit(file,'_')[[1]]
            run <- file.split[grepl('r[0-9]{1,2}i1p[0-9]{1}',file.split)]
            job.name <- paste0(substr(gcm,0,3),'.',substr(run,0,2),'.',substr(scenario,4,5))
            work <- paste0("qsub -N ",job.name," -v gcm='",gcm,"',run='",run,"',scenario='",scenario,"',type='",type,"',varname='",var.name,"',pctl='000' run.standard.variables.pbs")
            print(work)
            system(work)
            Sys.sleep(1)

         }
   }
}




