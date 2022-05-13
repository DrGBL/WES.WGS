#!/bin/bash

set -eu
################################################################################
# Please fill in the below variables
################################################################################
# Metadata
STUDY_NAME="BQC-19-Sweden"
ANALYST_LAST_NAME="Butler-Laporte"
DATE="$(date +'%Y%m%d')"
OUTNAME=/project/richards/guillaume.butler-laporte/WGS/bqc.individual.wgs.20200908/finalScripts/PC_PROJECTION_REVIEWS/results/"${STUDY_NAME}.${ANALYST_LAST_NAME}.${DATE}"
################################################################################
# Location of downloaded input files
PCA_LOADINGS="/project/richards/guillaume.butler-laporte/WGS/bqc.individual.wgs.20200908/finalScripts/PC_PROJECTION_REVIEWS/loadings_freq/grch38_loading.tsv"
PCA_AF="/project/richards/guillaume.butler-laporte/WGS/bqc.individual.wgs.20200908/finalScripts/PC_PROJECTION_REVIEWS/loadings_freq/grch38_freq.tsv"
################################################################################
# Location of imputed genotype files
# [Recommended]
# PLINK 2 binary format: a prefix (with directories) of .pgen/.pvar/.psam files
PFILE=""
# [Acceptable]
# PLINK 1 binary format: a prefix of .bed/.bim/.fam files
BFILE="/scratch/richards/guillaume.butler-laporte/WGS/regenieInputs/BQC19.hc.080921.swedes.exome.Eur.normID.GTflt.AB.noChrM.vqsr.flt"
################################################################################

function error_exit() {
  echo "${1:-"Unknown Error"}" 1>&2
  exit 1
}

# Input checks
if [[ -z "${STUDY_NAME}" ]]; then
  error_exit "Please specify \$STUDY_NAME."
fi

if [[ -z "${ANALYST_LAST_NAME}" ]]; then
  error_exit "Please specify \$ANALYST_LAST_NAME."
fi

if [[ -z "${PCA_LOADINGS}" ]]; then
  error_exit "Please specify \$PCA_LOADINGS."
fi

if [[ -z "${PCA_AF}" ]]; then
  error_exit "Please specify \$PCA_AF."
fi

if [[ -n "${PFILE}" ]]; then
  input_command="--pfile ${PFILE}"
elif [[ -n "${BFILE}" ]]; then
  input_command="--bfile ${BFILE}"
else
  error_exit "Either \$PFILE or \$BFILE should be specified"
fi

# Run plink2 --score
plink2 \
  ${input_command} \
  --score ${PCA_LOADINGS} \
  variance-standardize \
  cols=-scoreavgs,+scoresums \
  list-variants \
  header-read \
  --score-col-nums 3-12 \
  --read-freq ${PCA_AF} \
  --out ${OUTNAME}


  
