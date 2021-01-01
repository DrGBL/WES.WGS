#this takes the filtered vcf file from the WGS sequencing, removes variants with MAF < 1%, MAC < 5, and for whom HWE is not reached, and outputs the plink binary file.
#genotype call rate and MAF threshold is set high to avoid quasi-separation, since step 1 does not use Firth
#this file is to be used in step 1 of regenie. regenie does not accept vcf as inputs.
#warning: may take a bit of time.

mkdir -p regenieInputs

plink --vcf /path/to/sequence.file.Eur.normID.rehead.GTflt.AB.noChrM.vcf.gz \
 --hwe 1E-6 midp \
 --make-bed \
 --maf 0.01 \
 --mac 5 \
 --geno 0.1 \
 --out /regenieInputs/sequence.file.Eur.normID.rehead.GTflt.AB.noChrM 


#suggest to also prune for LD, to speed up regenie step 1, if this is done, then the following code should be ran (instead of the plink call above).

plink --vcf /path/to/sequence.file.Eur.normID.rehead.GTflt.AB.noChrM.vcf.gz \
 --hwe 1E-6 midp \
 --make-bed \
 --maf 0.01 \
 --mac 5 \
 --geno 0.1 \
 --indep-pairwise 50 5 0.05 \
 --out /regenieInputs/sequence.file.Eur.normID.rehead.GTflt.AB.noChrM 
 
