#!/bin/bash
#PBS -N gwas_parsing_ukb_sample
#PBS -S /bin/bash
#PBS -l walltime=1:00:00
#PBS -l mem=32gb
#PBS -l nodes=1:ppn=1

#PBS -o logs_kk/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e logs_kk/${PBS_JOBNAME}.e${PBS_JOBID}.err

module load gcc/6.2.0
module load python/3.5.3

cd $PBS_O_WORKDIR 

python3 /gpfs/data/im-lab/nas40t2/abarbeira/software/genomic_tools/genomic_tools/src/gwas_parsing.py \
-gwas_file /scratch/abarbeira3/v8_process/ubk_sign/gwas/50_standing_height.txt.gz \
--insert_value variant_id NA \
-output_column_map variant var  -split_column var ':' chromosome position non_effect_allele effect_allele  \
-output_column_map beta effect_size -output_column_map se standard_error -output_column_map n_complete_samples sample_size --chromosome_format  \
-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size n_cases \
-liftover /gpfs/data/im-lab/nas40t2/abarbeira/data/liftover/hg19ToHg38.over.chain.gz \
-snp_reference_metadata /gpfs/data/im-lab/nas40t2/abarbeira/projects/gtex_v8/data_formatting/vcf_process/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.txt.gz METADATA \
-output results_gwas_ukbn/50_standing_height.txt.gz

#