#some of the next few steps were too slow to do on the whole genome, and so I split my file in chromosomes.
#this is one of the steps that should 100% be performed using PBS arrays, if available, to speed up computing time

#path to QCed vcf restricted to given ancestry
pathAncestry=/scratch/richards/guillaume.butler-laporte/WGS/allSamples.Eur.normID.GTflt.AB.noChrM.vqsr.flt.vcf.gz

#path to split chromosomes output
pathSplit=/scratch/richards/guillaume.butler-laporte/WGS/splitChrom/allSamples.chr${x}.Eur.normID.GTflt.AB.noChrM.vqsr.flt.vcf.gz

for chr in {1..22} X; do
  bcftools filter -r chr${x}  $pathAncestry -Oz > ${pathSplit}
  tabix -p vcf ${pathSplit};
done
