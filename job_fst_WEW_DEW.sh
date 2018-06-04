#! /bin/bash
#PBS -d /storage/cozzip/Projects/danara
#PBS -l nodes=1:ppn=64
#PBS -N Fst_WEW_DEW
#PBS -q general

##########################################
#                                        #
#   Output some useful job information.  #
#                                        #
##########################################

echo ------------------------------------------------------
echo -n 'Job is running on node '; cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH
echo ------------------------------------------------------

# Starting job
echo "job $PBS_JOBNAME started: `whoami` `hostname` `pwd` `date`."

##########################################
#                                        #
#   Environment settings                 #
#                                        #
##########################################

# get reserved CPUs
export CORES=$(wc -l $PBS_NODEFILE | cut -f1 --delim=" ")

# activate conda environment
source activate R-3.4

##########################################
#                                        #
#   Calling R script                     #
#                                        #
##########################################

# executing Rscript
Rscript --slave --vanilla Fst_WEW_DEW.R

##########################################
#                                        #
#   Ending job                           #
#                                        #
##########################################

# finishing job
echo "job $PBS_JOBNAME finished `whoami` `hostname` `pwd` `date`."
