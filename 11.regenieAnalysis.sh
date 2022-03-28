#this performs the actual regenie analyses (steps 1 and 2, with and without the common variant exclusion file)
#note that the --htp option in regenie step 2 is important for the meta-analysis, so please make sure to use it

#paths

#path to QCed plink binaries
pathPlink=/scratch/richards/guillaume.butler-laporte/WGS/regenieInputs/allSamples.Eur.normID.GTflt.AB.noChrM.vqsr.flt
#path to covariates file
pathCov=/scratch/richards/guillaume.butler-laporte/WGS/regenieInputs/covar.txt
#path to phenotype file
pathPheno=/scratch/richards/guillaume.butler-laporte/WGS/regenieInputs/pheno.txt
#path to regenie tmp folder
pathTmpReg=/scratch/richards/guillaume.butler-laporte/WGS/regenieRes/tmp/
#path to pruned variants
pathPruned=/scratch/richards/guillaume.butler-laporte/WGS/regenieInputs/allSamples.regenie.LD.prune.maf0.01.geno0.1.Eur.normID.GTflt.AB.noChrM.vqsr.flt
#path to regenie inputs, where the mask.def, set.list, anno.file, and aaf.file are located. mask.def can be found on the git
pathReg=/scratch/richards/guillaume.butler-laporte/WGS/regenieInputs/
#path to regenie results folder
pathOut=/scratch/richards/guillaume.butler-laporte/WGS/regenieRes/

#name of your cohort
name=name.of.your.cohort



#step 1 is common to all

regenie \
  --step 1 \
  --bed $pathPlink \
  --covarFile $pathCov \
  --phenoFile $pathPheno \
  --bt \
  --lowmem \
  --lowmem-prefix "${pathTmpReg}"tmp_rg \
  --extract "${pathPruned}".prune.in \
  --bsize 1000 \
  --out "${pathOut}"step1AllPhenoLD


#with the pooled exclusion list (named grch38.maf.X.id.regenie.txt below), to be given to the participating cohorts
#note that since two exclusion lists are done (one for MAF>1% and one for MAR>0.1%), we need to do this step twice.

#for MAF>1%
plink \
  --bfile $pathPlink \
  --exclude /scratch/richards/guillaume.butler-laporte/WGS/regenieInputs/to.give.final/MAF.1.percent/for.use.with.grch38/grch38.maf.1.id.regenie.txt \
  --make-bed \
  --out ${pathPlink}.maf1perc

regenie \
  --step 2 \
  --minMAC 1 \
  --covarFile $pathCov \
  --phenoFile $pathPheno \
  --bed ${pathPlink}.maf1perc \
  --aaf-bins 0.01 \
  --build-mask 'max' \
  --write-mask \
  --mask-def "${pathReg}"regenie.mask.def.txt \
  --set-list "${pathReg}"regenie.set.list.txt \
  --anno-file "${pathReg}"regenie.anno.file.txt \
  --aaf-file "${pathReg}"regenie.aaf.file.txt \
  --bt \
  --htp $name \
  --firth --approx \
  --firth-se \
  --pred "${pathOut}"step1AllPhenoLD_pred.list \
  --bsize 200 \
  --out "${pathOut}"burden.res.common.1.perc
 
#For MAF>0.1%
plink \
  --bfile $pathPlink \
  --exclude /scratch/richards/guillaume.butler-laporte/WGS/regenieInputs/to.give.final/MAF.0.1.percent/for.use.with.grch38/grch38.maf.0.1.id.regenie.txt \
  --make-bed \
  --out ${pathPlink}.maf0.1perc

regenie \
  --step 2 \
  --minMAC 1 \
  --covarFile $pathCov \
  --phenoFile $pathPheno \
  --bed ${pathPlink}.maf0.1perc \
  --aaf-bins 0.001 \
  --build-mask 'max' \
  --write-mask \
  --mask-def "${pathReg}"regenie.mask.def.txt \
  --set-list "${pathReg}"regenie.set.list.txt \
  --anno-file "${pathReg}"regenie.anno.file.txt \
  --aaf-file "${pathReg}"regenie.aaf.file.txt \
  --bt \
  --htp $name \
  --firth --approx \
  --firth-se \
  --pred "${pathOut}"step1AllPhenoLD_pred.list \
  --bsize 200 \
  --out "${pathOut}"burden.res.common.0.1.perc
