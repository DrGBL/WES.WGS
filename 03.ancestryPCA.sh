#this may take some time to run in its entirety

#based on this https://www.biostars.org/p/335605/
#requires plink 1.9 and bcftools >1.3

#the idea here is as follows:
#1) obtain the 1000G reference files (in GRCh38 always, only SNPs, and autosomal chromosomes), normalize and left-align them with the reference genome
#2) from this file, only use the variants also found in your specific cohort, and further prune them to MAF over 10% and in linkage equilibrium
#3) from your cohort variant file, only select those pruned variants from step 2
#4) merge the 1000G file and you cohort's file, with only those pruned variants
#5) perform PCA on 1000G using this variant set, and project your cohort's genotype on the resulting PCs
#6) train a random forest with 6 principal components on the 1000G dataset
#7) predict the ancestry in your cohort using that trained random forest and your cohort's projection on 1000G's PCs
#8) output files named like "afrIDsPCA.txt", for the study IDs of individuals predicted to be of a certain ancestry (here african, for example).


# paths

# path to folder where computations will be done and results saved
pathWork=/scratch/richards/guillaume.butler-laporte/WGS/1000G.PCA/

# path to QCed vcf file (e.g. from hail)
pathQC=/scratch/richards/guillaume.butler-laporte/WGS/QCed/allSamples.normID.GTflt.AB.noChrM.vqsr.flt.vcf.bgz


cd $pathWork
mkdir -p cohortSample
mkdir -p Pruned

# download the GRCh38 1000G data, if not done already

prefix="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr"
suffix=".shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz"

for chr in {1..22}; do
  wget "$prefix""$chr""$suffix" "$prefix""$chr""$suffix".tbi
done

#now will need to download ancestry information from the 1000G
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped

#download reference genome fasta file (GRCh38)
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gzip -d GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

#convert 1000G vcf to bcf, normalizing and left-aligning variants
for chr in {1..22}; do
    bcftools norm -m-any --check-ref w -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
      ALL.chr"${chr}".shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz | \
      bcftools annotate -x ID -I +'%CHROM:%POS:%REF:%ALT' | \
        bcftools norm -Ob --rm-dup both \
          > ALL.chr"${chr}".shapeit2_integrated_v1a.GRCh38.20181129.phased.bcf ;

    bcftools index ALL.chr"${chr}".shapeit2_integrated_v1a.GRCh38.20181129.phased.bcf ;
done

#convert bcf to plink format
for chr in {1..22}; do
    plink --noweb \
      --bcf ALL.chr"${chr}".shapeit2_integrated_v1a.GRCh38.20181129.phased.bcf \
      --keep-allele-order \
      --vcf-idspace-to _ \
      --allow-extra-chr 0 \
      --split-x b38 no-fail \
      --make-bed \
      --out ALL.chr"${chr}".shapeit2_integrated_v1a.GRCh38.20181129.phased ;
done

#Obtain the list of variants for ancestry inference PCA
for chr in {1..22}; do
  bcftools view -r chr${chr} $pathQC | \
  bcftools query -f "%ID\n" > /cohortSample/variants.chr${chr}.txt;
done


#prune for common variants in linkage equilibrium, and only keeping variants also in your cohort
for chr in {1..22}; do
    plink --noweb \
      --bfile ALL.chr"${chr}".shapeit2_integrated_v1a.GRCh38.20181129.phased \
      --extract /cohortSample/variants.chr"${chr}".txt \
      --maf 0.10 --indep 50 5 1.5 \
      --out Pruned/ALL.chr"${chr}".shapeit2_integrated_v1a.GRCh38.20181129.phased ;

    plink --noweb \
      --bfile ALL.chr"${chr}".shapeit2_integrated_v1a.GRCh38.20181129.phased \
      --extract Pruned/ALL.chr"${chr}".shapeit2_integrated_v1a.GRCh38.20181129.phased.prune.in \
      --make-bed \
      --out Pruned/ALL.chr"${chr}".shapeit2_integrated_v1a.GRCh38.20181129.phased ;
done

#get a list of plink files
find . -name "*.bim" | grep -e "Pruned" > ForMerge.list ;

sed -i 's/.bim//g' ForMerge.list ;

#merge all biles in one
plink --merge-list ForMerge.list --out Merge ;

#make plink file from the full all sample VCF from your cohort, only keeping the pruned variants, then merge with the 1000G plink files
awk '{ print $2 }' Merge.bim > MergeVariants.txt

plink --vcf $pathQC \
 --noweb \
 --extract MergeVariants.txt \
 --make-bed \
 --out /cohortSample/cohortSample

printf "Merge\n./cohortSample/cohortSample" > ForMergeFull.list

plink --merge-list ForMergeFull.list --out MergeFullForPCA 
 
#now divide between 1000G and cohort samples
awk '{ print $1,$2 }' Merge.fam | awk '$(NF+1) = "1000G"' > 1000G.cluster.txt
awk '{ print $1,$2 }' cohortSample/cohortSample.fam | awk '$(NF+1) = "Cohort"' > cohort.cluster.txt
cat 1000G.cluster.txt cohort.cluster.txt > clusters.txt

#now do PCA on 1000G dataset, and project the cohort genotype on these PCs
plink --bfile MergeFullForPCA \
 --pca-cluster-names 1000G \
 --pca \
 --within clusters.txt
