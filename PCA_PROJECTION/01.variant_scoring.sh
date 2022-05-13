#!/bin/bash

set -eu
################################################################################
# Please fill in the below variables
################################################################################
# Metadata
STUDY_NAME="your cohort name"
ANALYST_LAST_NAME="" #your name
DATE="$(date +'%Y%m%d')"
OUTNAME=/path/to/results/folder/"${STUDY_NAME}.${ANALYST_LAST_NAME}.${DATE}"
################################################################################
# Location of downloaded input files
PCA_LOADINGS="" #either grch38_loading.tsv or grch37_loading.tsv, depending on your reference genome
PCA_AF="" #either grch38_freq.tsv or grch37_freq.tsv, depending on your reference genome
################################################################################
# Location of your plink file. Can use the results of $pathPlink in 08.WGS.vcf.to.plink.sh
# [Recommended]
# PLINK 2 binary format: a prefix (with directories) of .pgen/.pvar/.psam files
PFILE=""
# [Acceptable]
# PLINK 1 binary format: a prefix of .bed/.bim/.fam files
BFILE=""
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


  
