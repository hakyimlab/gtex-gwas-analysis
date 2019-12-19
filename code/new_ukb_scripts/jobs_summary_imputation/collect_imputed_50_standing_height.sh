#!/bin/bash
#PBS -N collect_imputed_50_standing_height
#PBS -S /bin/bash
#PBS -l walltime=1:00:00
#PBS -l mem=16gb
#PBS -l nodes=1:ppn=1
#PBS -o logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e logs/${PBS_JOBNAME}.e${PBS_JOBID}.err
########################################################################################################################
#CRI submission dandruff

module load gcc/6.2.0
module load python/3.5.3

cd $PBS_O_WORKDIR

python3 /gpfs/data/im-lab/nas40t2/abarbeira/software/genomic_tools/genomic_tools/src/gwas_summary_imputation_postprocess.py \
-gwas_file /scratch/abarbeira3/v8_process/ubk_sign/results_gwas_ukbn/50_standing_height.txt.gz \
-folder /scratch/abarbeira3/v8_process/ubk_sign/results_summary_imputation \
-pattern 50_standing_height.* \
-parsimony 7 \
-output processed_summary_imputation/imputed_50_standing_height.txt.gz