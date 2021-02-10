#this subsets the vcf to only the ancestry of interest. I took european ancestry as an example. This will need to be changed based on what you're working on.

#paths

#path to ancestry IDs
pathID=/scratch/richards/guillaume.butler-laporte/WGS/1000G.PCA/eurIDsPCA.txt

#path to QCed vcf
pathQC=/scratch/richards/guillaume.butler-laporte/WGS/QCed/allSamples.normID.GTflt.AB.noChrM.vqsr.flt.vcf.bgz

#path to QCed vcf restricted to the ancestry of interest
pathAncestry=/scratch/richards/guillaume.butler-laporte/WGS/allSamples.Eur.normID.GTflt.AB.noChrM.vqsr.flt.vcf.gz

#this splits the vcf to only the ancestry of interest. To change based on what you're working on.

bcftools view -S $pathID $pathQC -Oz > $pathAncestry
  
tabix -p vcf $pathAncestry

