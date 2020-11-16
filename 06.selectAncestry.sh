#this subsets the vcf to only the ancestry of interest. I took european ancestry as an example. This will need to be changeg based on what you're working on.

bcftools view -S /path/toeurIDsPCA.txt \
  /path/to/sequence.file.normID.rehead.GTflt.AB.noChrXYM.vcf.gz -Oz > /path/to/sequence.file.Eur.normID.rehead.GTflt.AB.noChrXYM.vcf.gz
  
tabix -p vcf /path/to/sequence.file.Eur.normID.rehead.GTflt.AB.noChrXYM.vcf.gz
