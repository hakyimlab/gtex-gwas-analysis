#!/bin/bash
#PBS -N predixcan_expression_summary
#PBS -S /bin/bash
#PBS -l walltime=8:00:00
#PBS -l mem=8gb
#PBS -l nodes=1:ppn=1
#PBS -t 1-87

#PBS -o logs_summary/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e logs_summary/${PBS_JOBNAME}.e${PBS_JOBID}.err

module load gcc/6.2.0
module load R/3.4.1

cd $PBS_O_WORKDIR
echo $PBS_O_WORKDIR

Rscript gtex_paper_enloc_multi_stats.R ${PBS_ARRAYID}
