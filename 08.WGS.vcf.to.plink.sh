#this function does two things necessary for regenie step 1 and 2 outputs
# 1) it outputs the QCed vcf restricted to the chosen ancestry to a plink binary file (used in regenie step 2)
# 2) from this plink binary file, it further filters variants to retain only those in linkage equilibrium with MAF>1% and that attained HWE (used in step 1 only)
# the filtering of variants above is used to build the polygenic score in regenie's step 1.

#paths

#path to given ancestry QCed vcf
pathAncestry=/scratch/richards/guillaume.butler-laporte/WGS/allSamples.Eur.normID.GTflt.AB.noChrM.vqsr.flt.vcf.gz

#path to given ancestry QCed plink binary
pathPlink=/scratch/richards/guillaume.butler-laporte/WGS/regenieInputs/allSamples.Eur.normID.GTflt.AB.noChrM.vqsr.flt

#path to pruned variants
pathPruned=/scratch/richards/guillaume.butler-laporte/WGS/regenieInputs/allSamples.regenie.LD.prune.maf0.01.geno0.1.Eur.normID.GTflt.AB.noChrM.vqsr.flt


plink --vcf $pathAncestry  \
  --make-bed \
  --out $pathPlink

# LD pruned variants for regenie step 1
plink --bfile $pathPlink \
 --hwe 1E-15 midp \
 --maf 0.01 \
 --geno 0.1 \
 --indep-pairwise 50 5 0.05 \
 --out $pathPruned
