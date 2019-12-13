#!/bin/bash
#PBS -N 50_standing_height_chr14_sb4_by_region
#PBS -S /bin/bash
#PBS -l walltime=4:00:00
#PBS -l mem=4gb
#PBS -l nodes=1:ppn=1
#PBS -o logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e logs/${PBS_JOBNAME}.e${PBS_JOBID}.err
########################################################################################################################
#CRI submission dandruff

module load gcc/6.2.0
module load python/3.5.3

export MKL_NUM_THREADS=1
export OPEN_BLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

cd $PBS_O_WORKDIR

python3 /gpfs/data/im-lab/nas40t2/abarbeira/software/genomic_tools/genomic_tools/src/gwas_summary_imputation.py \
-by_region_file /gpfs/data/im-lab/nas40t2/abarbeira/projects/gtex_v8/data/eur_ld.bed.gz \
-gwas_file /scratch/abarbeira3/v8_process/ubk_sign/results_gwas_ukbn/50_standing_height.txt.gz \
-parquet_genotype /gpfs/data/im-lab/nas40t2/abarbeira/projects/gtex_v8/data_formatting/model_training_to_parquet/results/parquet_eur_maf0.01_biallelic/gtex_v8_eur_itm.chr14.variants.parquet \
-parquet_genotype_metadata /gpfs/data/im-lab/nas40t2/abarbeira/projects/gtex_v8/data_formatting/model_training_to_parquet/results/parquet_eur_maf0.01_biallelic/gtex_v8_eur_itm.variants_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 14 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 4 \
--standardise_dosages \
-output results_summary_imputation/50_standing_height_chr14_sb4_reg0.1_ff0.01_by_region.txt.gz