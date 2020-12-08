#!/bin/bash                                                                                           
                 
gcm="EC-Earth3-Veg"
run='r1i1p1f1'
scenario='ssp126'

#qsub -N "ann.pr" -v gcm=$gcm,run=$run,scenario=$scenario,type='annual',varname="pr",pctl='000' run.standard.variables.pbs
#qsub -N "ann.tx" -v gcm=$gcm,run=$run,scenario=$scenario,type='annual',varname="tasmax",pctl='000' run.standard.variables.pbs
#qsub -N "ann.tn" -v gcm=$gcm,run=$run,scenario=$scenario,type='annual',varname="tasmin",pctl='000' run.standard.variables.pbs

#qsub -N "seas.pr" -v gcm=$gcm,run=$run,scenario=$scenario,type='seasonal',varname="pr",pctl='000' run.standard.variables.pbs
#qsub -N "seas.tx" -v gcm=$gcm,run=$run,scenario=$scenario,type='seasonal',varname="tasmax",pctl='000' run.standard.variables.pbs
#qsub -N "seas.tn" -v gcm=$gcm,run=$run,scenario=$scenario,type='seasonal',varname="tasmin",pctl='000' run.standard.variables.pbs

qsub -N "mon.pr" -v gcm=$gcm,run=$run,scenario=$scenario,type='monthly',varname="pr",pctl='000' run.standard.variables.pbs
qsub -N "mon.tx" -v gcm=$gcm,run=$run,scenario=$scenario,type='monthly',varname="tasmax",pctl='000' run.standard.variables.pbs
qsub -N "mon.tn" -v gcm=$gcm,run=$run,scenario=$scenario,type='monthly',varname="tasmin",pctl='000' run.standard.variables.pbs

#qsub -N "ext.pr" -v gcm=$gcm,run=$run,scenario=$scenario,type='annual_extremes',varname="pr",pctl='000' run.standard.variables.pbs
#qsub -N "ext.tx" -v gcm=$gcm,run=$run,scenario=$scenario,type='annual_extremes',varname="tasmax",pctl='000' run.standard.variables.pbs
#qsub -N "ext.tn" -v gcm=$gcm,run=$run,scenario=$scenario,type='annual_extremes',varname="tasmin",pctl='000' run.standard.variables.pbs

#qsub -N "qntl.tx" -v gcm=$gcm,run=$run,scenario=$scenario,type='annual_quantiles',varname="tasmax",pctl="975" run.standard.variables.pbs
#qsub -N "qntl.tx" -v gcm=$gcm,run=$run,scenario=$scenario,type='annual_quantiles',varname="tasmax",pctl="990" run.standard.variables.pbs
#qsub -N "qntl.tx" -v gcm=$gcm,run=$run,scenario=$scenario,type='annual_quantiles',varname="tasmax",pctl="996" run.standard.variables.pbs

#qsub -N "qntl.tn" -v gcm=$gcm,run=$run,scenario=$scenario,type='annual_quantiles',varname="tasmin",pctl="004" run.standard.variables.pbs
#qsub -N "qntl.tn" -v gcm=$gcm,run=$run,scenario=$scenario,type='annual_quantiles',varname="tasmin",pctl="010" run.standard.variables.pbs
#qsub -N "qntl.tn" -v gcm=$gcm,run=$run,scenario=$scenario,type='annual_quantiles',varname="tasmin",pctl="025" run.standard.variables.pbs
