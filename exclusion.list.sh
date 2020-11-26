
# This function requires normalized and left-aligned variants that passed QC. So in the pipeline in this github, this is done after the hail QC.
# It filters for MAF>1% (the --min-af 0.01:minor part) i.e. AF between 1% and 99%
# Then if pipes it to another bcftools function to query only the information that is required. You can add more if needed.
# Please remember to do this for each ancestry you plan on doing the analysis on. 
# Also provide separately a sample size for each ancestry used (in a separate document).

bcftools view --min-af 0.01:minor /path/to/sequence.file.QCed.normID.noChrM.vcf.gz | \
  bcftools query -H -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' > fullVariantList.txt
