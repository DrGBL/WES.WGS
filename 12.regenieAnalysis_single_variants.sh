#paths

#path to QCed plink binaries
pathPlink=/scratch/richards/guillaume.butler-laporte/WGS/regenieInputs/allSamples.Eur.normID.GTflt.AB.noChrM.vqsr.flt
#path to covariates file (same as for burden test, but without the rare variant PCs)
pathCov=/scratch/richards/guillaume.butler-laporte/WGS/regenieInputs/covar_single_variants.txt
#path to phenotype file (same as for burden test)
pathPheno=/scratch/richards/guillaume.butler-laporte/WGS/regenieInputs/pheno.txt
#path to regenie tmp folder
pathTmpReg=/scratch/richards/guillaume.butler-laporte/WGS/regenieRes/tmp/
#path to pruned variants
pathPruned=/scratch/richards/guillaume.butler-laporte/WGS/regenieInputs/allSamples.regenie.LD.prune.maf0.01.geno0.1.Eur.normID.GTflt.AB.noChrM.vqsr.flt
#path to regenie results folder
pathOut=/scratch/richards/guillaume.butler-laporte/WGS/regenieRes/

#name of your cohort
name=name.of.your.cohort


#step 1

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
  --out "${pathOut}"BQC19.hc.080921.eur.step1.single_variant
 
#step 2

regenie \
  --step 2 \
  --minMAC 6 \
  --covarFile $pathCov \
  --phenoFile $pathPheno \
  --bed ${pathPlink} \
  --bt \
  --htp $name \
  --firth --approx \
  --firth-se \
  --maxstep-null 5 \
  --maxiter-null 10000 \
  --pred "${pathOut}"BQC19.hc.080921.eur.step1.single_variant_pred.list \
  --bsize 200 \
  --out "${pathOut}"BQC19.hc.080921.eur.step2.single_variant
  
