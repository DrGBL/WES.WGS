#this function does two things necessary for regenie step 1 and 2 outputs
# 1) it outputs the QCed vcf restricted to the chosen ancestry to a plink binary file (used in step 1 and step 2)
# 2) it assigns a female sex to all samples, in order for the chromosome X to be coded as 0/2 for males. This is just a trick to make plink readjust the phenotypes
# and sex is added separately as a covariate in the analysis (step 11).
# 3) from this plink binary file, it further filters variants to retain only those in linkage equilibrium with MAF>1% and that attained HWE (used in step 1 only)
# the filtering of variants above is used to build the polygenic score in regenie's step 1.

#paths

#path to given ancestry QCed vcf
pathAncestry=/scratch/richards/guillaume.butler-laporte/WGS/allSamples.swedes.exome.Eur.normID.GTflt.AB.noChrM.vqsr.flt.vcf.gz

#path to chrX adjusted plink
pathPlinkTmp=/scratch/richards/guillaume.butler-laporte/WGS/regenieInputs/allSamples.swedes.exome.Eur.normID.GTflt.AB.noChrM.vqsr.flt.tmp

#path to FID IID and new sex (all female) file (input for plink)
allFemale=/scratch/richards/guillaume.butler-laporte/WGS/regenieInputs/all_female.txt

#path to given ancestry QCed plink binary
pathPlink=/scratch/richards/guillaume.butler-laporte/WGS/regenieInputs/allSamples.swedes.exome.Eur.normID.GTflt.AB.noChrM.vqsr.flt

#path to pruned variants
pathPruned=/scratch/richards/guillaume.butler-laporte/WGS/regenieInputs/allSamples.swedes.exome.regenie.LD.prune.maf0.01.geno0.1.Eur.normID.GTflt.AB.noChrM.vqsr.flt


plink --vcf $pathAncestry  \
  --make-bed \
  --out $pathPlinkTmp

awk '{print $1, $2, "2"}' $pathPlinkTmp.fam > $allFemale  

plink --bfile $pathPlinkTmp \
  --update-sex $allFemale \
  --make-bed \
  --out $pathPlink

# LD pruned variants for regenie step 1
plink --bfile $pathPlink \
 --hwe 1E-15 midp \
 --maf 0.01 \
 --geno 0.1 \
 --indep-pairwise 50 5 0.05 \
 --out $pathPruned
