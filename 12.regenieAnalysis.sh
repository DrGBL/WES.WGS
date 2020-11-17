#the following is the regenie code. Step 1 needs only to be done once. I encourage everyone to therefore put all phenotypes in the phenotype input file (A2,B2,C2).
#the phenotype and covariates inputs from regenie are described on their github and will be referred to as pheno.txt and covar.txt below
#regenie github: https://rgcgithub.github.io/regenie/

mkdir -p regenieRes
mkdir -r regenieRes/M1
mkdir -r regenieRes/M2
mkdir -r regenieRes/M3
mkdir -r regenieRes/M4


#step 1
regenie \
  --step 1 \
  --bed /path/to/regenieInputs/sequence.file.Eur.normID.rehead.GTflt.AB.noChrM \
  --covarFile /path/to/covar.txt \
  --phenoFile /path/to/pheno.txt \
  --bt --lowmem \
  --bsize 100 \
  --lowmem-prefix /regenieRes/ \
  --out /regenieRes/


 #step 2 for M1, outcome number 1
regenie \
  --step 2 \
  --bed /path/to/M1/finalMask/MergeM1 \
  --covarFile /path/to/covar.txt \
  --phenoFile /path/to/pheno.txt \
  --bt \
  --bsize 100 \
  --firth --approx \
  --pThresh 0.05 \
  --pred /regenieRes/regenieRes_1.loco \
  --split \
  --out /regenieRes/M1/M1.step2
  
   #step 2 for M2, outcome number 1
regenie \
  --step 2 \
  --bed /path/to/M2/finalMask/MergeM2 \
  --covarFile /path/to/covar.txt \
  --phenoFile /path/to/pheno.txt \
  --bt \
  --bsize 100 \
  --firth --approx \
  --pThresh 0.05 \
  --pred /regenieRes/regenieRes_1.loco \
  --split \
  --out /regenieRes/M2/M2.step2
  
   #step 2 for M3, outcome number 1
regenie \
  --step 2 \
  --bed /path/to/M3/finalMask/MergeM3 \
  --covarFile /path/to/covar.txt \
  --phenoFile /path/to/pheno.txt \
  --bt \
  --bsize 100 \
  --firth --approx \
  --pThresh 0.05 \
  --pred /regenieRes/regenieRes_1.loco \
  --split \
  --out /regenieRes/M3/M3.step2
  
   #step 2 for M4, outcome number 1
regenie \
  --step 2 \
  --bed /path/to/M4/finalMask/MergeM4 \
  --covarFile /path/to/covar.txt \
  --phenoFile /path/to/pheno.txt \
  --bt \
  --bsize 100 \
  --firth --approx \
  --pThresh 0.05 \
  --pred /regenieRes/regenieRes_1.loco \
  --split \
  --out /regenieRes/M4/M4.step2
