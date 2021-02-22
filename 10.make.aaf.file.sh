# this is the idea behind the code:
# 1) only keep variants that are either high or moderate VEP impacts, as these will not be used for the analysis anyway. This shrinks the list down.
# 2) find all remaining variants with a gnomAD/ESP MAF more than 1%, or a local MAF more than 1%, and give them an allele frequency of 10% for our purposes
# 3) of the remaining variants (i.e. singleton or MAF<1% without gnomAD/ESP annotation), check if they have a MAF<0.1%. If they do, assign them a MAF of 0.1%, otherwise a MAF of 0.5% for our purposes.
# 4) of the remaining variants, find the ones with a gnomAD/ESP MAF more than 0.1% (and hence less than 1%). Give them an allele frequency of 0.5% for our purposes.


#paths
#path to where the many temporary files will be stored
pathTmp=/scratch/richards/guillaume.butler-laporte/WGS/tmp/

#path to the regenie inputs
pathReg=/scratch/richards/guillaume.butler-laporte/WGS/regenieInputs/

for x in {1..22} X; do

#path to the VEP annotation outputs
 pathAnnot=/scratch/richards/guillaume.butler-laporte/WGS/annotation/finalAnnot.chr${x}.txt

#path to split chromosomes output
 pathSplit=/scratch/richards/guillaume.butler-laporte/WGS/splitChrom/allSamples.chr${x}.Eur.normID.GTflt.AB.noChrM.vqsr.flt.vcf.gz

# Hence in the same order as above

# 1) only keep variants that are either high or moderate VEP impacts, as these will not be used for the analysis anyway. This shrinks the list down.

 awk '/CANONICAL=YES/ && (/IMPACT=HIGH/ || /IMPACT=MODERATE/)' \
   $pathAnnot > "${pathTmp}"deleterious.var.chr${x}.txt

 awk '/CANONICAL=YES/ && (/IMPACT=HIGH/ || /IMPACT=MODERATE/)' \
   $pathAnnot | \
   awk '{print $1}' > "${pathTmp}"deleterious.id.chr${x}.txt


# 2) find all remaining variants with a gnomAD/ESP MAF more than 1%, or a local MAF more than 1%, and give them an allele frequency of 10% for our purposes

# this outputs variant MAF>1% in gnomAD
 awk '/AFR_AF=0\.0[1-9]|AFR_AF=0\.[1-8][0-9][0-9]|AFR_AF=0\.9[0-8]/ ||
   /AMR_AF=0\.0[1-9]|AMR_AF=0\.[1-8][0-9][0-9]|AMR_AF=0\.9[0-8]/ ||
   /EAS_AF=0\.0[1-9]|EAS_AF=0\.[1-8][0-9][0-9]|EAS_AF=0\.9[0-8]/ ||
   /EUR_AF=0\.0[1-9]|EUR_AF=0\.[1-8][0-9][0-9]|EUR_AF=0\.9[0-8]/ ||
   /SAS_AF=0\.0[1-9]|SAS_AF=0\.[1-8][0-9][0-9]|SAS_AF=0\.9[0-8]/ ||
   /AA_AF=0\.0[1-9]|AA_AF=0\.[1-8][0-9][0-9]|AA_AF=0\.9[0-8]/ ||
   /EA_AF=0\.0[1-9]|EA_AF=0\.[1-8][0-9][0-9]|EA_AF=0\.9[0-8]/ ||
   /gnomAD_AFR_AF=0\.0[1-9]|gnomAD_AFR_AF=0\.[1-8][0-9][0-9]|gnomAD_AFR_AF=0\.9[0-8]/ ||
   /gnomAD_AMR_AF=0\.0[1-9]|gnomAD_AMR_AF=0\.[1-8][0-9][0-9]|gnomAD_AMR_AF=0\.9[0-8]/ ||
   /gnomAD_ASJ_AF=0\.0[1-9]|gnomAD_ASJ_AF=0\.[1-8][0-9][0-9]|gnomAD_ASJ_AF=0\.9[0-8]/ ||
   /gnomAD_EAS_AF=0\.0[1-9]|gnomAD_EAS_AF=0\.[1-8][0-9][0-9]|gnomAD_EAS_AF=0\.9[0-8]/ ||
   /gnomAD_FIN_AF=0\.0[1-9]|gnomAD_FIN_AF=0\.[1-8][0-9][0-9]|gnomAD_FIN_AF=0\.9[0-8]/ ||
   /gnomAD_NFE_AF=0\.0[1-9]|gnomAD_NFE_AF=0\.[1-8][0-9][0-9]|gnomAD_NFE_AF=0\.9[0-8]/ ||
   /gnomAD_OTH_AF=0\.0[1-9]|gnomAD_OTH_AF=0\.[1-8][0-9][0-9]|gnomAD_OTH_AF=0\.9[0-8]/ ||
   /gnomAD_SAS_AF=0\.0[1-9]|gnomAD_SAS_AF=0\.[1-8][0-9][0-9]|gnomAD_SAS_AF=0\.9[0-8]/' \
   "${pathTmp}"deleterious.var.chr${x}.txt | \
   awk '!/#/' | \
   awk '{ print $1 }' | \
   sort -u -k 1.6 > "${pathTmp}"gnomAD.above.1perc.chr${x}.txt

# this outputs variants with MAF>1% in the cohort
 bcftools filter -i "%ID=@"${pathTmp}"deleterious.id.chr${x}.txt" \
   $pathSplit -Ou | \
   bcftools filter -i 'INFO/MAF>0.01' -Ou | \
   bcftools query -f '%ID \n' | sort -u -k 1.6 | > "${pathTmp}"cohort.above.1perc.chr${x}.txt

# combine the locally common and gnomAD common (MAF>1%) variants in one file
 cat "${pathTmp}"gnomAD.above.1perc.chr${x}.txt "${pathTmp}"cohort.above.1perc.chr${x}.txt | \
   sort -u -k 1.6 | \
   awk '$(NF+1) = "0.1"' > "${pathTmp}"above.1perc.chr${x}.txt

# this outputs variant not in the list above (i.e. those with MAF<1%), with their annotations, for further analysis

 awk 'NR==FNR{a[$1]=$1;next} !($1 in a)' \
   "${pathTmp}"above.1perc.chr${x}.txt \
   "${pathTmp}"deleterious.var.chr${x}.txt > "${pathTmp}"annot.deleterious.below.1perc.chr${x}.txt

 awk '{ print $1 }' "${pathTmp}"annot.deleterious.below.1perc.chr${x}.txt > "${pathTmp}"id.deleterious.below.1perc.chr${x}.txt

# remove garbage
 rm "${pathTmp}"deleterious.var.chr${x}.txt "${pathTmp}"deleterious.id.chr${x}.txt


# 3) of the remaining variants (i.e. singleton or MAF<1% without gnomAD/ESP annotation), check if they have a MAF<0.1%. If they do, assign them a MAF of 0.1%, otherwise a MAF of 0.5% for our purposes.

#find those with gnomAD annot more than 0.001 MAF
 awk '/AFR_AF=0\.00[1-9]|AFR_AF=0\.0[1-9]|AFR_AF=0\.[1-8][0-9][0-9]|AFR_AF=0\.9[0-8]|AFR_AF=0\.99[0-8]/ ||
   /AMR_AF=0\.00[1-9]|AMR_AF=0\.0[1-9]|AMR_AF=0\.[1-8][0-9][0-9]|AMR_AF=0\.9[0-8]|AMR_AF=0\.99[0-8]/ ||
   /EAS_AF=0\.00[1-9]|EAS_AF=0\.0[1-9]|EAS_AF=0\.[1-8][0-9][0-9]|EAS_AF=0\.9[0-8]|EAS_AF=0\.99[0-8]/ ||
   /EUR_AF=0\.00[1-9]|EUR_AF=0\.0[1-9]|EUR_AF=0\.[1-8][0-9][0-9]|EUR_AF=0\.9[0-8]|EUR_AF=0\.99[0-8]/ ||
   /SAS_AF=0\.00[1-9]|SAS_AF=0\.0[1-9]|SAS_AF=0\.[1-8][0-9][0-9]|SAS_AF=0\.9[0-8]|SAS_AF=0\.99[0-8]/ ||
   /AA_AF=0\.00[1-9]|AA_AF=0\.0[1-9]|AA_AF=0\.[1-8][0-9][0-9]|AA_AF=0\.9[0-8]|AA_AF=0\.99[0-8]/ ||
   /EA_AF=0\.00[1-9]|EA_AF=0\.0[1-9]|EA_AF=0\.[1-8][0-9][0-9]|EA_AF=0\.9[0-8]|EA_AF=0\.99[0-8]/ ||
   /gnomAD_AFR_AF=0\.00[1-9]|gnomAD_AFR_AF=0\.0[1-9]|gnomAD_AFR_AF=0\.[1-8][0-9][0-9]|gnomAD_AFR_AF=0\.9[0-8]|gnomAD_AFR_AF=0\.99[0-8]/ ||
   /gnomAD_AMR_AF=0\.00[1-9]|gnomAD_AMR_AF=0\.0[1-9]|gnomAD_AMR_AF=0\.[1-8][0-9][0-9]|gnomAD_AMR_AF=0\.9[0-8]|gnomAD_AMR_AF=0\.99[0-8]/ ||
   /gnomAD_ASJ_AF=0\.00[1-9]|gnomAD_ASJ_AF=0\.0[1-9]|gnomAD_ASJ_AF=0\.[1-8][0-9][0-9]|gnomAD_ASJ_AF=0\.9[0-8]|gnomAD_ASJ_AF=0\.99[0-8]/ ||
   /gnomAD_EAS_AF=0\.00[1-9]|gnomAD_FIN_AF=0\.0[1-9]|gnomAD_FIN_AF=0\.[1-8][0-9][0-9]|gnomAD_FIN_AF=0\.9[0-8]|gnomAD_FIN_AF=0\.99[0-8]/ ||
   /gnomAD_FIN_AF=0\.00[1-9]|gnomAD_FIN_AF=0\.0[1-9]|gnomAD_FIN_AF=0\.[1-8][0-9][0-9]|gnomAD_FIN_AF=0\.9[0-8]|gnomAD_FIN_AF=0\.99[0-8]/ ||
   /gnomAD_NFE_AF=0\.00[1-9]|gnomAD_NFE_AF=0\.0[1-9]|gnomAD_NFE_AF=0\.[1-8][0-9][0-9]|gnomAD_NFE_AF=0\.9[0-8]|gnomAD_NFE_AF=0\.99[0-8]/ ||
   /gnomAD_OTH_AF=0\.00[1-9]|gnomAD_OTH_AF=0\.0[1-9]|gnomAD_OTH_AF=0\.[1-8][0-9][0-9]|gnomAD_OTH_AF=0\.9[0-8]|gnomAD_OTH_AF=0\.99[0-8]/ ||
   /gnomAD_SAS_AF=0\.00[1-9]|gnomAD_SAS_AF=0\.0[1-9]|gnomAD_SAS_AF=0\.[1-8][0-9][0-9]|gnomAD_SAS_AF=0\.9[0-8]|gnomAD_SAS_AF=0\.99[0-8]/' \
   "${pathTmp}"annot.deleterious.below.1perc.chr${x}.txt | \
   awk '!/#/' | \
   awk '{ print $1 }' | \
   sort -u -k 1.6 > "${pathTmp}"gnomad.01.to.1perc.chr${x}.txt

#find the singleton or those with MAF<0.1% in the cohort, and assign them 0.001 MAF
 bcftools filter -i "%ID=@"${pathTmp}"id.deleterious.below.1perc.chr${x}.txt" \
   $pathSplit -Ou | \
   bcftools filter -i 'INFO/MAF<=0.001 || MAC==1' -Ou | \
   bcftools filter -e "%ID=@"${pathTmp}"gnomad.01.to.1perc.chr${x}.txt" -Ou | \
   bcftools query -f '%ID \n' | \
   sort -u -k 1.6 | \
   awk '$(NF+1) = "0.0005"' > "${pathTmp}"deleterious.below.0.1perc.chr${x}.txt

# 4) of the remaining variants, find the ones with a gnomAD/ESP MAF more than 0.1% (and hence less than 1%). Give them an allele frequency of 0.5% for our purposes.

#find the non singleton with either gnomAD MAF>0.1% or cohort MAF>0.1% and assign them 0.5%
 awk 'NR==FNR{a[$1]=$1;next} !($1 in a) {print $1}' \
   "${pathTmp}"deleterious.below.0.1perc.chr${x}.txt \
   "${pathTmp}"annot.deleterious.below.1perc.chr${x}.txt |
   awk '$(NF+1) = "0.005"' > "${pathTmp}"deleterious.0.1.to.1perc.chr${x}.txt


 rm "${pathTmp}"annot.deleterious.below.1perc.chr${x}.txt
 rm "${pathTmp}"cohort.above.1perc.chr${x}.txt
 rm "${pathTmp}"gnomad.01.to.1perc.chr${x}.*
 rm "${pathTmp}"gnomAD.above.1perc.chr${x}.*
 rm "${pathTmp}"id.deleterious.below.1perc.chr${x}.*;  
 done
 
#now build the actual regenie inputs
cat above.1perc.chr* deleterious.below.0.1perc* deleterious.0.1.to.1perc* > "${pathReg}"regenieInputs/regenie.aaf.file.txt
rm above.1perc.chr* deleterious.below.0.1perc* deleterious.0.1.to.1perc*
 
