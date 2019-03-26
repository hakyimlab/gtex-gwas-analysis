#!/bin/bash
#PBS -N gtex_v8_count_hits
#PBS -S /bin/bash
#PBS -l walltime=1:00:00
#PBS -l mem=4gb
#PBS -l nodes=1:ppn=1
#PBS -t 1-92

#PBS -o logs_hits/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e logs_hits/${PBS_JOBNAME}.e${PBS_JOBID}.err

module load gcc/6.2.0
module load R/3.4.1

cd $PBS_O_WORKDIR
echo $PBS_O_WORKDIR

Rscript count_hits.R ${PBS_ARRAYID}