#This function normalizes variants and left-aligns them. Will be useful for further downstream analysis, including ancestry ascertainment and variant annotation

#it requires a fasta reference genome file (here the GRCh38), this can be downloaded with the following command (which then needs unzipping):
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gzip -d GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

#path to your vcf file
pathVCF=/scratch/richards/guillaume.butler-laporte/WGS/rawVCF/allSamples.hc.011421.vqsr.flt.vcf.gz

#path to the reference genome .fna file (downloaded above for GRCh38)
pathRef=/scratch/richards/guillaume.butler-laporte/WGS/1000G.PCA/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna

#path to your normalized, left-aligned, annotated with variants IDs vcf output
pathOut=/scratch/richards/guillaume.butler-laporte/WGS/normalizedVCF/allSamples.normID.noChrM.vqsr.flt.vcf.gz


#the first line below lists the chromosome you want to keep, I removed chromosome Y and mitochondrial chromosome for the rest of the analysis.
bcftools filter -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX \
  "${pathVCF}" -Ou | \
  bcftools norm -m -any --check-ref w -f "${pathRef}" -Ou | \
  bcftools annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -Oz > "${pathOut}"
tabix -p vcf "${pathOut}"

