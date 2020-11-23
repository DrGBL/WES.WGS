#now the problem with the output of the new QC'ed file (post Hail) is that the vcf header metadata needs fixing, or you will get downstream problems.
#unfortunately part of this needs to be done by hand. This is slightly annoying, apologies.

#first obtain the header from pre-hail vcf
bcftools view -h /path/to/sequence.file.noChrM.vcf.gz > /path/to/header.txt

#obtain the header from the hail output vcf
#Here I used zcat and grep because bcftools misbehaved on the post-hail vcf file. I suspect it's because the header was not formatted well. Another reason to rehead it.
#the "head" call is simply because I know that the header is the first few lines, so there's no reason for grep to go through the whole genome, which would be very slow.
zcat /path/to/sequence.file.GTflt.AB.noChrM.vcf.gz | head -n 5000 | grep "#" > /path/to/postHailHeader.txt

#the postHailHeader.txt header will need to be modified manually so that the metadata works
#e.g. in the postHailHeader you will have: ##FORMAT=<ID=AD,Number=.,Type=Integer,Description="">
#this needs to be fixed to #FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed"> 

#Once this is done for all relevant lines proceed with the following:
zcat /path/to/sequence.file.GTflt.AB.noChrM.vcf.gz | \
  bcftools reheader -h /path/to/postHailHeader.txt | \
  bgzip > /path/to/sequence.file.rehead.GTflt.AB.noChrM.vcf.gz
tabix -p vcf /path/to/sequence.file.rehead.GTflt.AB.noChrM.vcf.gz
