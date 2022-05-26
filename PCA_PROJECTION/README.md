# PCA_PROJECTION

***Addendum: I created a docker image for the R part of the analysis, since it was causing problems with version and dependencies, see at the end (before the example), for details.***

This code will allow all cohorts in the consortium to project their participant's exome on the same reference panel, using the same variant loadings. That is, it'll make our principal components comparable. As usual file paths need to be adjusted.

Step 01 can be ran once, but step 02 should be ran once for each phenotype. In step 02, you will have the option of running one ancestry at a time, or all of them together (preferred if your cohort has multiple ancestries).

Once the code has run, please compress all .rds outputs, all .png outputs, and the .projected.pca.tsv.gz output in one file, and send it to me. These should not contain any identifiable information, and will allow me to arrange the plots in one larger plot.

The code was heavily inspired by the one provided by the COVID-19 Host Genetics Initiative, available here for reference: https://github.com/covid19-hg/pca_projection. My thanks to their coding team. However, note that the reference used is different, and if you provided plots for the COVID-19 HGI GWAS, then outputs from their code and the one here cannot be directly compared (although clustering with the 1000G reference ancestries should be preserved).

Here's a summary of each file:

Files `grch3X_freq.tsv`, `grch3X_loading.tsv`, and `1000G_snps_scores.txt.gz` are reference files obtained from the exonic regions of the 1000G reference panel. They include around 25,000 variants pruned for LD and with high allele frequency (MAF>10% in the 1000G panel). These files are used in the two steps below. Please use the one corresponding to your reference genome (grch37 or 38). ***The variants were named using the following convention: chrN:pos:ref:alt i.e. each variant is contained in both grch37 and grch38 files, but with different variant IDs. For the functions to work, the variants in your plink files (in step 01) must be named the same way, otherwise you will need to rename either your variants, or the files here.***

`01.variant_scoring.sh`: this scores each variant and adds up their contribution to each principal components. The file that ends in sscore will be used for step 2.

`02.pca_projection.sh`: this calls the `plot_projected_pc.R` R script, which plots the principal components obtained above in the same space as the 1000G reference panel. This outputs both png files and rds objects. Please use the part that applies to the number of ancestries in your cohort (and adjust phenotype accordingly).

## Docker image
A docker image for the R component of the code is available for easier control of dependencies. To pull the docker image simply run the following:
```
docker pull monsieurbl/pc_proj_r
```

Then run the following in the directory which contains the expected inputs for `02.pca_projection.sh`.
```
docker run -ti --rm -v `pwd`:`pwd` -w `pwd` monsieurbl/pc_proj_r bash
```
And you will be able to directly run the correct Rscript from `02.pca_projection.sh` afterwards. For ease, I have added the reference file (`ref_score`) and the R script (`plot_projected_pc.R`) in the image. Once you run the command above (`docker run`), you can refer to these as `/app/loadings_freq/1000G_snps_scores.txt.gz` and `/app/plot_projected_pc.R`, respectively. Hence you would run the following script (example for one ancestry only):
```
Rscript /app/plot_projected_pc.R \
  --sscore ${sscore_file} \
  --phenotype-file ${pheno_file} \
  --phenotype-col ${pheno_col} \
  --ancestry ${anc} \
  --plot-pc-num 10 \
  --study ${study_name} \
  --reference-score-file /app/loadings_freq/1000G_snps_scores.txt.gz \
  --out ${path_out}/${study_name}_${pheno_col}_projection
```

## Example with 1000G:
In the end, the resulting file ending with PC1-10.png should overlap with the principal component plots below, which were obtained with the same variant loadings and code as provided here, but using participants from the 1000G reference panel (simply disregard the "cases" and "controls" in the legend, they aren't plotted here). ***There might be differences in your results and the example below due to varying number of variants that will be used, these will cause the scales to be a bit different. I'll fix this in the final combined plots.***
![1000G_plots projected pca ancestry all PC1-10](https://user-images.githubusercontent.com/25112827/168510862-9bb2ed43-3489-4e1c-9a24-09ae98e18bf5.png)
