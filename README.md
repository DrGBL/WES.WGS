# Covid-19 HGI WES/WGS burden test BQC-19 pipeline

This is code made available to other members of the consortium, in an effort to help with their local pipelines. The full analysis protocol (v4) can be found here: https://docs.google.com/document/d/1Rr1mhvY6JViHpkGg_QtoP0VfB0-YNqb2vxscNY7PbVw/edit

Note that all code is written assuming that the GRCh38 genome is used.

Every effort was made to prioritize code clarity over efficiency, and there are multiple ways this code can be made more efficient e.g. using PBS arrays if your local computing cluster supports it.

The code is numbered to give the order which it was ran at McGill, however these aren't always necessary, and you can sometimes jump ahead.

Dependencies: R (with tidyverse), bcftools, plink, regenie, VEP (with LOFTEE and dbNSFP plugins), and hail (python 3.7 module)

***The most important tweak to be done locally is to modify the path names in each of the functions, depending on how your cluster works. That is, the functions below will not work "out-of-the-box", you still need to tell them where to save or load some of the files they output or require as input.***

Here's a summary of what each function does, with more comments in each specific file:

`exclusionList.sh`: this is the suggested code to obtain the list of variants to be used to build our list of exclusion variants list.

`01.autosomalXY.sh`: removes mitochondrial DNA (chrM).

`02.norm.ID.sh`: normalizes and left aligns variants using the reference genome.

`03.hail.py`: QC using hail. At the BQC sex imputation and QC had been done by our genome center, as well as the variant recalibration, so they are omitted here. For the full hail QC code can be found here: https://github.com/mkveerapen/covid19_sequencing

`04.reheader.sh`: advanced hail functions did not work on our local cluster, so the rest of our analysis was done on regular vcf. This function reheaders the meta-data the vcf obtained from the hail QC. It will require some local tweaking.

`05.ancestryPCA.sh`: uses 1000G to train a random forest classifier to infer continental ancestry in your cohort. Note that this may not be fully necessary if homogeneous ancestry is expected.

`06.selectAncestry.sh`: subsets samples from a specific ancestry as obtained in the previous function (here I used europeans as an example)

`07.PCA.sh`: calculates common and rare variant PCAs for your selected samples

`08.splitChrom.sh`: here I split the vcf file by chromosome, as the next few functions were too slow otherwise.

`09.finalAnnot.sh`: uses VEP to annotate variants. This requires LOFTEE and dbNSFP plugins.

`10.WGS.vcf.to.plink.sh`: this transforms the full VCF file in plink format, and removes the rare variants. This will be used for step 1 of regenie (which does not use Firth regression)

`11.makeRegenieVCF.sh`: this is where we build the different masks, and code genes as though they were SNPs, to be used in step 2 of regenie. This function depends on 4 R files: `M1.GT.R`, `M2.GT.R`, `M3.GT.R`, and `M4.GT.R`.

`12.regenieAnalysis.sh`: this is where the magic happens with regenie.





