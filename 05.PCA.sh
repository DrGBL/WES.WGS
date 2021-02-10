
#this function outputs PCA results from common and rare variants from the selected ancestry obtained in the previous step.
#these PCs will be used in the covar file input in regenie (code to build that file is not shown in the git, since it's too dependent on how your data is stored

#paths

#path to QCed vcf restricted to ancestry of interest
pathAncestry=/scratch/richards/guillaume.butler-laporte/WGS/allSamples.Eur.normID.GTflt.AB.noChrM.vqsr.flt.vcf.gz

#path to ouput directory
pathPCA=/scratch/richards/guillaume.butler-laporte/WGS/PCA/

#common
plink --vcf $pathAncestry \
--biallelic-only strict \
--chr 1-22 \
--geno 0.05 \
--snps-only 'just-acgt' \
--hwe 1E-6 midp \
--indep-pairwise 50 5 0.05 \
--keep-allele-order \
--mac 5 \
--maf 0.01 \
--out "${pathPCA}"commonAllelesPruned


plink --vcf $pathAncestry \
--extract "${pathPCA}"commonAllelesPruned.prune.in \
--pca 10 \
--out "${pathPCA}"commonPCA.txt


#rare

plink --vcf $pathAncestry \
--biallelic-only strict \
--chr 1-22 \
--geno 0.05 \
--snps-only 'just-acgt' \
--hwe 1E-6 midp \
--indep-pairwise 50 5 0.05 \
--keep-allele-order \
--max-maf 0.01 \
--mac 2 \
--out "${pathPCA}"rareAllelesPruned


plink --vcf $pathAncestry \
--extract "${pathPCA}"rareAllelesPruned.prune.in \
--pca 20 \
--out "${pathPCA}"rarePCA.txt
