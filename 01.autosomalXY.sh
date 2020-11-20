#The "sequencing.file.vcf.gz" file below is the variant call vcf. Note that at our center it comes with variant recalibration already performed.

#This removes mitochondrial DNA

bcftools filter -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY \
  /path/to/sequencing.file.vcf.gz \
  -Oz > /path/to/sequencing.file.noChrM.vcf.gz
tabix -p vcf /path/to/sequencing.file.noChrM.vcf.gz
