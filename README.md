# Covid-19 HGI WES/WGS burden test BQC-19 pipeline

***June 12, 2021: Added hard filters and QC plots in the hail qc step (step2), and updated step 8 to recode the X chromosome as diploid (i.e. 0/2 for males) by making all participants "females" in plink.***

This is code made available to other members of the consortium, in an effort to help with their local pipelines. The full analysis protocol (v5) can be found here: https://docs.google.com/document/d/1Ouii904IqUArMECXHWynZjBiZ8QO-r2-3rsymB815Bo/edit

Note that all code is written assuming that the GRCh38 genome is used. It also assumes that the initial input is a ***joint-called*** vcf file with all samples to be used for the analysis.

Every effort was made to prioritize code clarity over efficiency, and there are multiple ways this code can be made more efficient e.g. using PBS arrays if your local computing cluster supports it.

The code is numbered to give the order which it was ran at McGill, however these aren't always necessary, and you can sometimes jump ahead.

Dependencies: R (with tidyverse and randomForest), bcftools, plink (v1.9), regenie (v2), VEP (with dbNSFP plugin), and hail (python 3.7 module)

***Again, the most important tweak to be done locally is to modify the path names in each of the functions, depending on how your cluster works. However, I've made efforts to make these easy to spot and modify. That is, the functions below will not work "out-of-the-box", you still need to tell them where to save or load some of the files they output or require as input.***

Here's a summary of what each function does, with more comments in each specific file:

`exclusion.list.sh`: this is the suggested code to obtain the list of variants to be used to build our list of exclusion variants.

`01.norm.ID.sh`: normalizes and left aligns variants using the reference genome. Also removes the Y chromosome and the mitochondrial chromosomes.

`02.hail.py`: QC using hail. At the BQC sex imputation and QC had been done by our genome center, as well as the variant recalibration, so they are omitted here. For the full hail QC code can be found here: https://github.com/mkveerapen/covid19_sequencing

`03.ancestryPCA.sh` and `03.ancestryPCA.R`: use 1000G to train a random forest classifier to infer continental ancestry in your cohort. Note that this may not be fully necessary if homogeneous ancestry is expected, and may not be sufficient if the cohort has a very admixed ancestry either.

`04.selectAncestry.sh`: subsets samples from a specific ancestry as obtained in the previous function (here I used europeans as an example)

`05.PCA.sh`: calculates common and rare variant PCAs for your selected samples, to be used in the covar input file in regenie.

`06.splitChrom.sh`: here I split the vcf file by chromosome, as the next few functions were too slow otherwise.

`07.finalAnnot.sh`: uses VEP to annotate variants. This requires the dbNSFP plugins.

`08.WGS.vcf.to.plink.sh`: this transforms the full VCF file in plink format, and removes the rare variants. This will be used for step 1 of regenie (which does not use Firth regression)

`09.variantPrep.sh`: this uses VEP annotations to build the --set-list and --anno-file regenie step 2 inputs.

`10.make.aaf.file.sh`: this uses VEP annotations and your cohort's QCed vcf file to assign the correct allele frequencies to each variants used in the burden tests.

`11.regenieAnalysis.sh`: this is where we use regenie to perform the burden tests.

`transpose.awk`: a vector transpose function, useful in steps 09 and 10.
