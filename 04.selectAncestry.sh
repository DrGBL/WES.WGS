#this subsets the vcf to only the ancestry of interest. I took european ancestry as an example. This will need to be changed based on what you're working on.

bcftools view -S /path/to/eurIDsPCA.txt \
  /path/to/sequence.file.normID.GTflt.AB.noChrM.vcf.bgz -Oz > /path/to/sequence.file.Eur.normID.GTflt.AB.noChrM.vcf.gz
  
tabix -p vcf /path/to/sequence.file.Eur.normID.GTflt.AB.noChrM.vcf.gz
