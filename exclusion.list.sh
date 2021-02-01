# This function requires normalized and left-aligned variants that passed QC. In the pipeline on this github, this is done after the hail QC.

# In the first part, it filters for MAF>1% (the --min-af 0.01:minor part) i.e. AF between 1% and 99%
# Then if pipes it to another bcftools function to query only the information that is required. You can add more if needed.
# The result of this list will be used for the exclusion list in the MAF<1% analysis

# Note that in order to make it a smaller file, and since we will be using whole exome anyway, if you are using WGS you can restrict the analysis to the exome sequences in the exomePos.txt (bed format)
# I provide an example of such a file (the one I used for the BQC19) on the git, based on GENCODE reference: https://www.gencodegenes.org/

bcftools view -R /path/to/exomePos.txt --min-af 0.01:minor \
  /path/to/sequence.file.Eur.normID.GTflt.AB.noChrM.vcf.gz | \
  bcftools query -H -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' > common.one.perc.ExomeEur.txt

# In the second part, it filters for MAF>0.1% AND MAC>6
# This will be used in the MAF<0.1% analysis

bcftools view -R /path/to/exomePos.txt --min-ac 6:minor --min-af 0.01:minor \
  /path/to/sequence.file.Eur.normID.GTflt.AB.noChrM.vcf.gz | \
  bcftools query -H -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' > common.zero.one.perc.ExomeEur.txt
  
# Please remember to do both parts for each ancestry you plan on doing the analysis on. Here I used european as an example.
# Also provide separately a sample size for each ancestry used (in a separate document).
