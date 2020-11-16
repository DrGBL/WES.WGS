#this is somewhat of a long file, because it is split in four masks.
#it uses the M1.GT.R M2.GT.R M3.GT.R and M4.GT.R functions

#builds the header for the vcf file I'm building
printf "##fileformat=VCFv4.2\n##fileDate=20201110\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">" > headerStep1Final.txt
zcat /scratch/richards/guillaume.butler-laporte/WGS/allSamples.Eur.normID.rehead.GTflt.AB.noChrXYM.vqsr.flt.vcf.gz | head -n 5000 | grep '#CHROM'  > /project/richards/guillaume.butler-laporte/WGS/bqc.individual.wgs.20200908/VCF.sct.regenie/columns.txt

###### M1 note that this uses max AF among 1000 genomes, ESP, and gnomAD populations, combined. ######
mkdir -p M1
mkdir -p M1/finalMask
#first filter for the correct variants
for i in {1..22}; do awk '/LoF=HC/ && /CANONICAL=YES/ && !/MAX_AF=0.[0-9][1-9]*;/' /path/to/annotation/finalAnnot.chr${i}.txt > /M1/M1.annot.chr${i}.txt; done

#now obtain files listing the variants and their corresponding genes, ordered by chromosomal position
for i in {1..22}; do awk '{ print $4, $1 }' /M1/M1.annot.chr${i}.txt | sort -k 2 | uniq > /M1/M1.long.chr${i}.txt; done
for i in {1..22}; do awk '{ print $2 }' /M1/M1.long.chr${i}.txt > /M1/M1.variants.chr${i}.txt; done

#now use bcftools to view only those variants and obtain each sample's genotype
for i in {1..22}; do bcftools view -i "%ID=@/M1/M1.variants.chr${i}.txt" /path/to/sequence.file.chr${i}.Eur.normID.rehead.GTflt.AB.noChrXYM.vcf.gz -Ou | bcftools filter -i 'INFO/MAF<=0.01' -Ou | bcftools query -f '%ID [ %GT] \n' > /M1/M1.GT.chr${i}.txt; done

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
for i in {1..22}; do awk '(/LoF=HC/ || /missense_variant/) && /CANONICAL=YES/ && !/MAX_AF=0.[0-9][1-9]*;/' /path/to/annotation/finalAnnot.chr${i}.txt > /M2/M2.annot.chr${i}.txt; done

#now obtain files listing the variants and their corresponding genes, ordered by chromosomal position
for i in {1..22}; do awk '{ print $4, $1 }' /M2/M2.annot.chr${i}.txt | sort -k 2 | uniq > /M2/M2.long.chr${i}.txt; done
for i in {1..22}; do awk '{ print $2 }' /M2/M2.long.chr${i}.txt > /M2/M2.variants.chr${i}.txt; done

#now use bcftools to view only those variants and obtain each sample's genotype
for i in {1..22}; do bcftools view -i "%ID=@/M2/M2.variants.chr${i}.txt" /path/to/sequence.file.chr${i}.Eur.normID.rehead.GTflt.AB.noChrXYM.vcf.gz -Ou | bcftools filter -i 'INFO/MAF<=0.01' -Ou | bcftools query -f '%ID [ %GT] \n' > /M2/M2.GT.chr${i}.txt; done

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



###### M3 note that this uses max AF among 1000 genomes, ESP, and gnomAD populations, combined. ######
mkdir -p M3
mkdir -p M3/finalMask
#first filter for the correct variants
for i in {1..22}; do awk '(/loF=HC/ || (/missense_variant/ && /SIFT_pred=[,DT]*D/ && /Polyphen2_HVAR_pred=[,DPB]*D[,DPB]*;/ && /Polyphen2_HDIV_pred=[,DPB]*D[,DPB]*;/ && (/MutationTaster_pred=[,ADNP]*D[,ADNP]*;/ || /MutationTaster_pred=[,ADNP]*A[,ADNP]*;/) && /LRT_pred=[;DNU]*D[;DNU]*;/)) && /CANONICAL=YES/ && !/MAX_AF=0.[0-9][1-9]*;/' /path/to/annotation/finalAnnot.chr${i}.txt > /M3/M3.annot.chr${i}.txt; done

#now obtain files listing the variants and their corresponding genes, ordered by chromosomal position
for i in {1..22}; do awk '{ print $4, $1 }' /M3/M3.annot.chr${i}.txt | sort -k 2 | uniq > /M3/M3.long.chr${i}.txt; done
for i in {1..22}; do awk '{ print $2 }' /M3/M3.long.chr${i}.txt > /M3/M3.variants.chr${i}.txt; done

#now use bcftools to view only those variants and obtain each sample's genotype
for i in {1..22}; do bcftools view -i "%ID=@/M3/M3.variants.chr${i}.txt" /path/to/sequence.file.chr${i}.Eur.normID.rehead.GTflt.AB.noChrXYM.vcf.gz -Ou | bcftools filter -i 'INFO/MAF<=0.01' -Ou | bcftools query -f '%ID [ %GT] \n' > /M3/M3.GT.chr${i}.txt; done

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
for i in {1..22}; do awk '(/loF=HC/ || (/missense_variant/ && (/SIFT_pred=[,DT]*D/ || /Polyphen2_HVAR_pred=[,DPB]*D[,DPB]*;/ || /Polyphen2_HDIV_pred=[,DPB]*D[,DPB]*;/ || /MutationTaster_pred=[,ADNP]*D[,ADNP]*;/ || /MutationTaster_pred=[,ADNP]*A[,ADNP]*;/ || /LRT_pred=[;DNU]*D[;DNU]*;/))) && /CANONICAL=YES/ && !/MAX_AF=0.[0-9][1-9]*;/' /path/to/annotation/finalAnnot.chr${i}.txt > /M4/M4.annot.chr${i}.txt; done

#now obtain files listing the variants and their corresponding genes, ordered by chromosomal position
for i in {1..22}; do awk '{ print $4, $1 }' /M4/M4.annot.chr${i}.txt | sort -k 2 | uniq > /M4/M4.long.chr${i}.txt; done
for i in {1..22}; do awk '{ print $2 }' /M4/M4.long.chr${i}.txt > /M4/M4.variants.chr${i}.txt; done

#now use bcftools to view only those variants and obtain each sample's genotype
for i in {1..22}; do bcftools view -i "%ID=@/M4/M4.variants.chr${i}.txt" /path/to/sequence.file.chr${i}.Eur.normID.rehead.GTflt.AB.noChrXYM.vcf.gz -Ou | bcftools filter -i 'INFO/MAF<=0.01' -Ou | bcftools query -f '%ID [ %GT] \n' > /M4/M4.GT.chr${i}.txt; done

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
