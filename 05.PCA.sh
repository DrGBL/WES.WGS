#this function outputs PCA results from common and rare variants from the selected ancestry obtained in the previous step.
#again I assume european ancestry here, as per previous step.

mkdir -p PCA

#common
plink --vcf /path/to/sequence.file.Eur.normID.GTflt.AB.noChrM.vcf.gz \
--biallelic-only strict \
--chr 1-22 \
--geno 0.05 \
--snps-only 'just-acgt' \
--hwe 1E-6 midp \
--indep-pairwise 50 5 0.05 \
--keep-allele-order \
--mac 5 \
--maf 0.01 \
--out /PCA/commonAllelesPruned


plink --vcf /path/to/sequence.file.Eur.normID.GTflt.AB.noChrM.vcf.gz \
--extract /PCA/commonAllelesPruned.in \
--pca 10 \
--out /PCA/commonPCA.txt


#rare, given the relatively small size of my cohort, I set a MAC threshold of 2. This can be adjusted by each cohort.

plink --vcf /path/to/sequence.file.Eur.normID.GTflt.AB.noChrM.vcf.gz \
--biallelic-only strict \
--chr 1-22 \
--geno 0.05 \
--snps-only 'just-acgt' \
--hwe 1E-6 midp \
--indep-pairwise 50 5 0.05 \
--keep-allele-order \
--max-maf 0.01 \
--mac 2 \
--out /PCA/rareAllelesPruned


plink --vcf /path/to/sequence.file.Eur.normID.GTflt.AB.noChrM.vcf.gz \
--extract /PCA/rareAllelesPruned.in \
--pca 20 \
--out /PCA/rarePCA.txt
