#The "sequencing.file.vcf.gz" file below is the variant call vcf. Note that at our center it comes with variant recalibration already performed.

#This removes all non autosomal chromosomes from further analysis.

bcftools filter -r chr1,chr2,chr3,ch4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22 \
  /path/to/sequencing.file.vcf.gz \
  -Oz > /path/to/sequencing.file.noChrXYM.vcf.gz
tabix -p vcf /scratch/richards/guillaume.butler-laporte/WGS/rawVCF/allSamples.noChrXYM.vqsr.flt.vcf.gz
