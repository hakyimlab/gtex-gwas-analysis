GTEX_DIR="/gpfs/data/gtex-group/v8"
SUBDIR="59349/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_all_associations/"
EQTL_DIR=${GTEX_DIR}/${SUBDIR}

for FILE in `ls $EQTL_DIR`; do
  qsub -v DIR=${DIR},FILE=${FILE} extract_beta_and_se.pbs
done
