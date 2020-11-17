#This function normalizes variants and left-aligns them. Will be useful for further downstream analysis, including ancestry ascertainment and variant annotation

#it requires a fasta reference genome file (here the GRCh38), this can be downloaded with the following command (which then needs unzipping):
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gzip -d GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

bcftools norm -m -any --check-ref w -f /path/to/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna /path/to/sequence.file.rehead.GTflt.AB.noChrXY.vcf.gz -Ou | \
  bcftools annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -Oz > /path/to/sequence.file.normID.rehead.GTflt.AB.noChrXY.vcf.gz
tabix -p vcf /path/to/sequence.file.normID.rehead.GTflt.AB.noChrXY.vcf.gz


