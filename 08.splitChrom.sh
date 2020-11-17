#some of the next few steps were too slow to do on the whole genome, and so I split my file in chromosomes.
#this is one of the steps that should 100% be performed using PBS arrays, if available, to speed up computing time

for chr in {1..22}; do
  bcftools filter -r chr"${chr}" /path/to/sequence.file.Eur.normID.rehead.GTflt.AB.noChrM.vcf.gz -Oz | \
    > /path/to/sequence.file.chr"${chr}".Eur.normID.rehead.GTflt.AB.noChrM.vcf.gz
    tabix -p vcf /path/to/sequence.file.chr"${chr}".Eur.normID.rehead.GTflt.AB.noChrM.vcf.gz;
done
