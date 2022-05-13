sscore_file=/path/to/sscore/file #the output from the previous step
pheno_file=/path/to/phenotype/file #the same phenotype file you provided as input to regenie for the association analyses
pheno_col="string" #A, B, or C, or whatever else you called them
study_name="string" #your study name
ref_score=/path/to/1000G_snps_scores.txt.gz
path_out=/path/out/folder/ #only the folder, the name of the files are auto-filled

#if only one ancestry, also need to add the following (see below if multiple ancestries)
anc="string" #AFR, AMR, EAS, EUR, MID, or SAS

Rscript plot_projected_pc.R \
  --sscore ${sscore_file} \
  --phenotype-file ${pheno_file} \
  --phenotype-col ${pheno_col} \
  --ancestry ${anc} \
  --plot-pc-num 10 \
  --study ${study_name} \
  --reference-score-file ${ref_score} \
  --out ${path_out}/${study_name}_${pheno_col}_projection


#if multiple ancestries
