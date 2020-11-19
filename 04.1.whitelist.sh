
# The function is in 2 parts. First, it repeats what was done in function 04.normID.sh, in case you've already ran it.
# Second, it outputs chromosome, position, ID, alleles, and MAF.

# Part 1 (again, can be skipped if 04.normID.sh was done already)
# It should still be done only after variant, genotype, and sample QC was performed
# It requires a fasta reference genome file (here the GRCh38), this can be downloaded with the following command (which then needs unzipping):
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gzip -d GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

bcftools norm -m -any --check-ref w -f /path/to/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna /path/to/sequence.file.rehead.GTflt.AB.noChrM.vcf.gz -Ou | \
  bcftools annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -Oz > /path/to/sequence.file.normID.rehead.GTflt.AB.noChrM.vcf.gz
tabix -p vcf /path/to/sequence.file.normID.rehead.GTflt.AB.noChrM.vcf.gz

# Part 2, using the output from part 1 above (or from 04.normID.sh)
bcftools query -H -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF\n' /path/to/sequence.file.normID.rehead.GTflt.AB.noChrM.vcf.gz > fullVariantList.txt
