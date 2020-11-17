#now the problem with the output of the new QC'ed file (post Hail) is that the vcf header metadata needs fixing, or you will get downstream problems.
#unfortunately part of this needs to be done by hand. This is slightly annoying, apologies.

#first obtain the header from pre-hail vcf
bcftools view -h /path/to/sequence.file.noChrM.vcf.gz > /path/to/header.txt

#obtain the header from the hail output vcf
bcftools view -h /path/to/sequence.file.GTflt.AB.noChrM.vcf.gz > /path/to/postHailHeader.txt

#the postHailHeader.txt header will need to be modified manually so that the metadata works
#e.g. in the postHailHeader you will have: ##FORMAT=<ID=AD,Number=.,Type=Integer,Description="">
#this needs to be fixed to #FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed"> 

#Once this is done for all relevant lines proceed with the following:
bcftools reheader -h /path/to/postHailHeader.txt /path/to/sequence.file.GTflt.AB.noChrM.vcf.gz | \
  bgzip > /path/to/sequence.file.rehead.GTflt.AB.noChrM.vcf.gz
tabix -p vcf /path/to/sequence.file.rehead.GTflt.AB.noChrM.vcf.gz
