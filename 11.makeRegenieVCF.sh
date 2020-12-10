#this is somewhat of a long file, because it is split in four masks.
#it uses the M1.GT.R M2.GT.R M3.GT.R and M4.GT.R functions

#first filter out variants with any of 1000G, ESP, or gnomAD populations MAF>1%
for i in {1..22} X Y; do awk '!/AFR_AF=0\.0[1-9]|AFR_AF=0\.[1-8][0-9][0-9]|AFR_AF=0\.9[0-8]/ &&
    !/AMR_AF=0\.0[1-9]|AMR_AF=0\.[1-8][0-9][0-9]|AMR_AF=0\.9[0-8]/ &&
    !/EAS_AF=0\.0[1-9]|EAS_AF=0\.[1-8][0-9][0-9]|EAS_AF=0\.9[0-8]/ &&
    !/EUR_AF=0\.0[1-9]|EUR_AF=0\.[1-8][0-9][0-9]|EUR_AF=0\.9[0-8]/ &&
    !/SAS_AF=0\.0[1-9]|SAS_AF=0\.[1-8][0-9][0-9]|SAS_AF=0\.9[0-8]/ &&
    !/AA_AF=0\.0[1-9]|AA_AF=0\.[1-8][0-9][0-9]|AA_AF=0\.9[0-8]/ &&
    !/EA_AF=0\.0[1-9]|EA_AF=0\.[1-8][0-9][0-9]|EA_AF=0\.9[0-8]/ &&
    !/gnomAD_AFR_AF=0\.0[1-9]|gnomAD_AFR_AF=0\.[1-8][0-9][0-9]|gnomAD_AFR_AF=0\.9[0-8]/ &&
    !/gnomAD_AMR_AF=0\.0[1-9]|gnomAD_AMR_AF=0\.[1-8][0-9][0-9]|gnomAD_AMR_AF=0\.9[0-8]/ &&
    !/gnomAD_ASJ_AF=0\.0[1-9]|gnomAD_ASJ_AF=0\.[1-8][0-9][0-9]|gnomAD_ASJ_AF=0\.9[0-8]/ &&
    !/gnomAD_EAS_AF=0\.0[1-9]|gnomAD_EAS_AF=0\.[1-8][0-9][0-9]|gnomAD_EAS_AF=0\.9[0-8]/ &&
    !/gnomAD_FIN_AF=0\.0[1-9]|gnomAD_FIN_AF=0\.[1-8][0-9][0-9]|gnomAD_FIN_AF=0\.9[0-8]/ &&
    !/gnomAD_NFE_AF=0\.0[1-9]|gnomAD_NFE_AF=0\.[1-8][0-9][0-9]|gnomAD_NFE_AF=0\.9[0-8]/ &&
    !/gnomAD_OTH_AF=0\.0[1-9]|gnomAD_OTH_AF=0\.[1-8][0-9][0-9]|gnomAD_OTH_AF=0\.9[0-8]/ &&
    !/gnomAD_SAS_AF=0\.0[1-9]|gnomAD_SAS_AF=0\.[1-8][0-9][0-9]|gnomAD_SAS_AF=0\.9[0-8]/' | 
    /path/to/annotation/finalAnnot.chr${x}.txt > /path/to/annotation/rareAnnot.rare.chr${x}.txt;
done


#builds the header for the vcf file I'm building
zcat /path/to/sequence.file.Eur.normID.rehead.GTflt.AB.noChrM.vcf.gz | head -n 5000 | grep '#CHROM'  > /path/to/columns.txt

###### M1 note that this uses max AF among 1000 genomes, ESP, and gnomAD populations, combined. ######
mkdir -p M1
mkdir -p M1/finalMask
#first filter for the correct variants
for i in {1..22}; do awk '/LoF=HC/ && /CANONICAL=YES/' /path/to/annotation/rareAnnot.chr${i}.txt > /M1/M1.annot.chr${i}.txt; done

#now obtain files listing the variants and their corresponding genes, ordered by chromosomal position
for i in {1..22}; do awk '{ print $4, $1 }' /M1/M1.annot.chr${i}.txt | sort -k 2 | uniq > /M1/M1.long.chr${i}.txt; done
for i in {1..22}; do awk '{ print $2 }' /M1/M1.long.chr${i}.txt > /M1/M1.variants.chr${i}.txt; done

#now use bcftools to view only those variants and obtain each sample's genotype
for i in {1..22}; do bcftools view -i "%ID=@/M1/M1.variants.chr${i}.txt" /path/to/sequence.file.chr${i}.Eur.normID.rehead.GTflt.AB.noChrM.vcf.gz -Ou | bcftools filter -i 'INFO/MAF<=0.01' -Ou | bcftools query -f '%ID [ %GT] \n' > /M1/M1.GT.chr${i}.txt; done

#now use R to build the genotype
Rscript /path/to/M1.GT.R
for i in {1..22}; do cat headerStep1.txt /M1/finalMask/M1.chr${i}.txt > /M1/finalMask/M1.chr${i}.vcf; done

#now use plink to build 
for i in {1..22}; do plink --vcf /M1/finalMask/M1.chr${i}.vcf --make-bed --out /M1/finalMask/M1.chr${i}; done

#get a list of those plink files
find /M1/finalMask/ -name "*.bim" | grep -e "chr" > /M1/finalMask/ForMerge.list ;

sed -i 's/.bim//g' /M1/finalMask/ForMerge.list ;

#merge all files in one
plink --merge-list /M1/finalMask/ForMerge.list --out /M1/finalMask/MergeM1 ;



###### M2 note that this uses max AF among 1000 genomes, ESP, and gnomAD populations, combined. ######
mkdir -p M2
mkdir -p M2/finalMask
#first filter for the correct variants
for i in {1..22}; do awk '(/LoF=HC/ || /missense_variant/) && /CANONICAL=YES/' /path/to/annotation/rareAnnot.chr${i}.txt > /M2/M2.annot.chr${i}.txt; done

#now obtain files listing the variants and their corresponding genes, ordered by chromosomal position
for i in {1..22}; do awk '{ print $4, $1 }' /M2/M2.annot.chr${i}.txt | sort -k 2 | uniq > /M2/M2.long.chr${i}.txt; done
for i in {1..22}; do awk '{ print $2 }' /M2/M2.long.chr${i}.txt > /M2/M2.variants.chr${i}.txt; done

#now use bcftools to view only those variants and obtain each sample's genotype
for i in {1..22}; do bcftools view -i "%ID=@/M2/M2.variants.chr${i}.txt" /path/to/sequence.file.chr${i}.Eur.normID.rehead.GTflt.AB.noChrM.vcf.gz -Ou | bcftools filter -i 'INFO/MAF<=0.01' -Ou | bcftools query -f '%ID [ %GT] \n' > /M2/M2.GT.chr${i}.txt; done

#now use R to build the genotype
Rscript /path/to/M2.GT.R
for i in {1..22}; do cat headerStep1.txt /M2/finalMask/M2.chr${i}.txt > /M2/finalMask/M2.chr${i}.vcf; done

#now use plink to build 
for i in {1..22}; do plink --vcf /M2/finalMask/M2.chr${i}.vcf --make-bed --out /M2/finalMask/M2.chr${i}; done

#get a list of those plink files
find /M2/finalMask/ -name "*.bim" | grep -e "chr" > /M2/finalMask/ForMerge.list ;

sed -i 's/.bim//g' /M2/finalMask/ForMerge.list ;

#merge all files in one
plink --merge-list /M2/finalMask/ForMerge.list --out /M2/finalMask/MergeM2 ;


#for M3 and M4, need to play around with the dbNSFP vep annotation a bit to match the canonical variant
#I will create a tmp file here due to this part creating a lot of temporary files

mkdir -p path/to/tmp

for i in {1..22} X Y; 
  do awk '/missense_variant/ && /CANONICAL=YES;/ && /VEP_canonical=/' /path/to/annotation/rareAnnot.rare.chr${x} > /path/to/tmp/Missense.annot.chr${i}.txt
  
#sift  
  grep "VEP_canonical=[,\.]*Y" /path/to/tmp/Missense.annot.chr${i}.txt | grep -o -P "(SIFT_pred=[,\.DT]*D[,\.DT]*;)" | grep -o -P "(=[,\.DT]*D[,\.DT]*;)" > /path/to/tmp/sift.count.chr${i}.txt
  grep "SIFT_pred=[,\.DT]*D[,\.DT]*;" /path/to/tmp/Missense.annot.chr${i}.txt | grep -o -P "(VEP_canonical=.*)" | grep -o -P "(=[,\.]*Y)" | awk '{ print length }' > /path/to/tmp/siftCanonical.count.chr${i}.txt
  grep "VEP_canonical=[,\.]*Y" /path/to/tmp/Missense.annot.chr${i}.txt | grep "SIFT_pred=[,\.DT]*D[,\.DT]*;" | awk '{ print $1 }' > /path/to/tmp/siftVariant.chr${i}.txt
  paste /path/to/tmp/siftVariant.chr${i}.txt /path/to/tmp/sift.count.chr${i}.txt /path/to/tmp/siftCanonical.count.chr${i}.txt > /path/to/tmp/sift.chr${i}.txt
  awk '{ if (substr($2,$3,1) ~ /D/) { print $1 } }' /path/to/tmp/sift.chr${i}.txt | sort | uniq > /path/to/tmp/sift.del.canon.chr${i}.txt
  
#polyphen HVAR  
  grep "VEP_canonical=[,\.]*Y" /path/to/tmp/Missense.annot.chr${i}.txt | grep -o -P "(Polyphen2_HVAR_pred=[,DPB\.]*D[,DPB\.]*;)" | grep -o -P "(=[,DPB\.]*D[,DPB\.]*;)"  > /path/to/tmp/pphenHVAR.count.chr${i}.txt
  grep "Polyphen2_HVAR_pred=[,DPB\.]*D[,DPB\.]*;" /path/to/tmp/Missense.annot.chr${i}.txt | grep -o -P "(VEP_canonical=.*)" | grep -o -P "(=[,\.]*Y)" | awk '{ print length }' > /path/to/tmp/pphenHVARCanonical.count.chr${i}.txt
  grep "VEP_canonical=[,\.]*Y" /path/to/tmp/Missense.annot.chr${i}.txt | grep "Polyphen2_HVAR_pred=[,DPB\.]*D[,DPB\.]*;" | awk '{ print $1 }' > /path/to/tmp/pphenHVARVariant.chr${i}.txt
  paste //path/to/tmp/pphenHVARVariant.chr${i}.txt /path/to/tmp/pphenHVAR.count.chr${i}.txt /path/to/tmp/pphenHVARCanonical.count.chr${i}.txt > /path/to/tmp/pphenHVAR.chr${i}.txt
  awk '{ if (substr($2,$3,1) ~ /D/) { print $1 } }' /path/to/tmp/pphenHVAR.chr${i}.txt | sort | uniq > /path/to/tmp/pphenHVAR.del.canon.chr${i}.txt

#polyphen HDIV
  grep "VEP_canonical=[,\.]*Y" /path/to/tmp/Missense.annot.chr${i}.txt | grep -o -P "(Polyphen2_HDIV_pred=[,DPB\.]*D[,DPB\.]*;)" | grep -o -P "(=[,DPB\.]*D[,DPB\.]*;)"  > /path/to/tmp/pphenHDIV.count.chr${i}.txt
  grep "Polyphen2_HDIV_pred=[,DPB\.]*D[,DPB\.]*;" /path/to/tmp/Missense.annot.chr${i}.txt | grep -o -P "(VEP_canonical=.*)" | grep -o -P "(=[,\.]*Y)" | awk '{ print length }' > /path/to/tmp/pphenHDIVCanonical.count.chr${i}.txt
  grep "VEP_canonical=[,\.]*Y" /path/to/tmp/Missense.annot.chr${i}.txt | grep "Polyphen2_HDIV_pred=[,DPB\.]*D[,DPB\.]*;" | awk '{ print $1 }' > /path/to/tmp/pphenHDIVVariant.chr${i}.txt
  paste /path/to/tmp/pphenHDIVVariant.chr${i}.txt /path/to/tmp/pphenHDIV.count.chr${i}.txt /path/to/tmp/pphenHDIVCanonical.count.chr${i}.txt > /path/to/tmp/pphenHDIV.chr${i}.txt
  awk '{ if (substr($2,$3,1) ~ /D/) { print $1 } }' /path/to/tmp/pphenHDIV.chr${i}.txt | sort | uniq > /path/to/tmp/pphenHDIV.del.canon.chr${i}.txt

#mutation taster
  grep "VEP_canonical=[,\.]*Y" /path/to/tmp/Missense.annot.chr${i}.txt | grep -o -P "(MutationTaster_pred=[,ADNP\.]*[AD][,ADNP\.]*;)" | grep -o -P "(=[,ADNP\.]*[AD][,ADNP\.]*;)"  > /path/to/tmp/mutationTaster.count.chr${i}.txt
  grep "MutationTaster_pred=[,ADNP\.]*[AD][,ADNP\.]*;" /path/to/tmp/Missense.annot.chr${i}.txt | grep -o -P "(VEP_canonical=.*)" | grep -o -P "(=[,\.]*Y)" | awk '{ print length }' > /path/to/tmp/mutationTasterCanonical.count.chr${i}.txt
  grep "VEP_canonical=[,\.]*Y" /path/to/tmp/Missense.annot.chr${i}.txt | grep "MutationTaster_pred=[,ADNP\.]*[AD][,ADNP\.]*;" | awk '{ print $1 }' > /path/to/tmp/mutationTasterVariant.chr${i}.txt
  paste /path/to/tmp/mutationTasterVariant.chr${i}.txt /path/to/tmp/mutationTaster.count.chr${i}.txt /path/to/tmp/mutationTasterCanonical.count.chr${i}.txt > /path/to/tmp/mutationTaster.chr${i}.txt
  awk '{ if (substr($2,$3,1) ~ /[AD]/) { print $1 } }' /path/to/tmp/mutationTaster.chr${i}.txt | sort | uniq > /path/to/tmp/mutationTaster.del.canon.chr${i}.txt

#LRT
  awk '(/LRT_pred=D;/ && /CANONICAL=YES;/)' /path/to/tmp/Missense.annot.chr${i}.txt | awk '{ print $1 }' | sort | uniq  > /path/to/tmp/LRT.del.canon.chr${i}.txt

#loftee 
  awk '{ print $1 }' /scratch/richards/guillaume.butler-laporte/WGS/Masks/M1/M1.annot.chr${i}.txt > /path/to/tmp/LOF.chr${i}.txt
 
#combine for M3 mask 
  comm -12 /path/to/tmp/sift.del.canon.chr${i}.txt /path/to/tmp/pphenHVAR.del.canon.chr${i}.txt | sort | uniq > /path/to/tmp/tmp1.chr${i}.txt
  comm -12 /path/to/tmp/tmp1.chr${i}.txt /path/to/tmp/pphenHDIV.del.canon.chr${i}.txt | sort | uniq > /path/to/tmp/tmp2.chr${i}.txt
  comm -12 /path/to/tmp/tmp2.chr${i}.txt /path/to/tmp/mutationTaster.del.canon.chr${i}.txt | sort | uniq > /path/to/tmp/tmp3.chr${i}.txt
  comm -12 /path/to/tmp/tmp3.chr${i}.txt /path/to/tmp/LRT.del.canon.chr${i}.txt | sort | uniq > /path/to/tmp/tmp4.chr${i}.txt
  cat /path/to/tmp/tmp4.chr${i}.txt /path/to/tmp/LOF.chr${i}.txt | sort | uniq > /path/to/tmp/M3.ID.chr${i}.txt

#combine for M4 mask
  cat /path/to/tmp/sift.del.canon.chr${i}.txt /path/to/tmp/pphenHVAR.del.canon.chr${i}.txt /path/to/tmp/pphenHDIV.del.canon.chr${i}.txt /path/to/tmp/mutationTaster.del.canon.chr${i}.txt /path/to/tmp/LRT.del.canon.chr${i}.txt /path/to/tmp/LOF.chr${i}.txt | sort | uniq > /path/to/tmp/M4.ID.chr${i}.txt
  
done



###### M3 note that this uses max AF among 1000 genomes, ESP, and gnomAD populations, combined. ######
mkdir -p M3
mkdir -p M3/finalMask
#first filter for the correct variants
for i in {1..22} X Y; 
  do awk '(/LoF=HC/ || (/missense_variant/ && /SIFT_pred=[,DT\.]*D[,DT\.]*;/ && /Polyphen2_HVAR_pred=[,DPB\.]*D[,DPB\.]*;/ && /Polyphen2_HDIV_pred=[,DPB\.]*D[,DPB\.]*;/ && (/MutationTaster_pred=[,ADNP\.]*D[,ADNP\.]*;/ || /MutationTaster_pred=[,ADNP\.]*A[,ADNP\.]*;/) && /LRT_pred=D;/)) && /CANONICAL=YES;/' /path/to/annotation/rareAnnot.chr${i}.txt | sort > /path/to/tmp/sorted.M3.Annot.chr${i}.txt 

  join -j 1 /path/to/tmp/M3.ID.chr${i}.txt /path/to/tmp/sorted.M3.Annot.chr${i}.txt > /M3/M3.annot.chr${i}.txt
done

#now obtain files listing the variants and their corresponding genes, ordered by chromosomal position
for i in {1..22}; do awk '{ print $4, $1 }' /M3/M3.annot.chr${i}.txt | sort -k 2 | uniq > /M3/M3.long.chr${i}.txt; done
for i in {1..22}; do awk '{ print $2 }' /M3/M3.long.chr${i}.txt > /M3/M3.variants.chr${i}.txt; done

#now use bcftools to view only those variants and obtain each sample's genotype
for i in {1..22}; do bcftools view -i "%ID=@/M3/M3.variants.chr${i}.txt" /path/to/sequence.file.chr${i}.Eur.normID.rehead.GTflt.AB.noChrM.vcf.gz -Ou | bcftools filter -i 'INFO/MAF<=0.01' -Ou | bcftools query -f '%ID [ %GT] \n' > /M3/M3.GT.chr${i}.txt; done

#now use R to build the genotype
Rscript  /path/to/M3.GT.R
for i in {1..22}; do cat headerStep1.txt /M3/finalMask/M3.chr${i}.txt > /M3/finalMask/M3.chr${i}.vcf; done

#now use plink to build 
for i in {1..22}; do plink --vcf /M3/finalMask/M3.chr${i}.vcf --make-bed --out /M3/finalMask/M3.chr${i}; done

#get a list of those plink files
find /M3/finalMask/ -name "*.bim" | grep -e "chr" > /M3/finalMask/ForMerge.list ;

sed -i 's/.bim//g' /M3/finalMask/ForMerge.list ;

#merge all files in one
plink --merge-list /M3/finalMask/ForMerge.list --out /M3/finalMask/MergeM3 ;




###### M4 note that this uses max AF among 1000 genomes, ESP, and gnomAD populations, combined. ######
mkdir -p M4
mkdir -p M4/finalMask
#first filter for the correct variants
for i in {1..22}; do awk '(/loF=HC/ || (/missense_variant/ && (/SIFT_pred=[,DT]*D/ || /Polyphen2_HVAR_pred=[,DPB]*D[,DPB]*;/ || /Polyphen2_HDIV_pred=[,DPB]*D[,DPB]*;/ || /MutationTaster_pred=[,ADNP]*D[,ADNP]*;/ || /MutationTaster_pred=[,ADNP]*A[,ADNP]*;/ || /LRT_pred=[;DNU]*D[;DNU]*;/))) && /CANONICAL=YES/' /path/to/annotation/rareAnnot.chr${i}.txt > /M4/M4.annot.chr${i}.txt; done

#now obtain files listing the variants and their corresponding genes, ordered by chromosomal position
for i in {1..22}; do awk '{ print $4, $1 }' /M4/M4.annot.chr${i}.txt | sort -k 2 | uniq > /M4/M4.long.chr${i}.txt; done
for i in {1..22}; do awk '{ print $2 }' /M4/M4.long.chr${i}.txt > /M4/M4.variants.chr${i}.txt; done

#now use bcftools to view only those variants and obtain each sample's genotype
for i in {1..22}; do bcftools view -i "%ID=@/M4/M4.variants.chr${i}.txt" /path/to/sequence.file.chr${i}.Eur.normID.rehead.GTflt.AB.noChrM.vcf.gz -Ou | bcftools filter -i 'INFO/MAF<=0.01' -Ou | bcftools query -f '%ID [ %GT] \n' > /M4/M4.GT.chr${i}.txt; done

#now use R to build the genotype
Rscript  /path/to/M4.GT.R
for i in {1..22}; do cat headerStep1.txt /M4/finalMask/M4.chr${i}.txt > /M4/finalMask/M4.chr${i}.vcf; done

#now use plink to build 
for i in {1..22}; do plink --vcf /M4/finalMask/M4.chr${i}.vcf --make-bed --out /M4/finalMask/M4.chr${i}; done

#get a list of those plink files
find /M4/finalMask/ -name "*.bim" | ggrep -e "chr" > /M4/finalMask/ForMerge.list ;

sed -i 's/.bim//g' /M4/finalMask/ForMerge.list ;

#merge all files in one
plink --merge-list /M4/finalMask/ForMerge.list --out /M4/finalMask/MergeM4 ;
