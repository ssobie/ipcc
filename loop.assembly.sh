#!/bin/bash

centre='CanESM5'
echo $centre

qsub -N "assemble.pr" -v centre=$centre,varname='pr' run.assembly.pbs 
qsub -N "assemble.tx" -v centre=$centre,varname='tasmax' run.assembly.pbs 
qsub -N "assemble.tn" -v centre=$centre,varname='tasmin' run.assembly.pbs


