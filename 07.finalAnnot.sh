#this annotates each chromosome using VEP and dbNSFP
#again, this could be sped up using PBS arrays if your cluster permits
#it assumes that vep was installed to work offline, and with a cache directory
#the --fork 10 makes it faster, but can give trouble sometimes
#here I did it on the european population, but it can be done once on the whole sample (then used in each ancestry stratified analysis).

#paths
#path to output folder
pathAnnot=/scratch/richards/guillaume.butler-laporte/WGS/annotation/

#path to VEP cache directory
pathCache=/project/richards/guillaume.butler-laporte/bin/.vep

for x in {1..22} X; do
  #set the path to splitchromosomes
  pathSplit=/scratch/richards/guillaume.butler-laporte/WGS/splitChrom/allSamples.chr"${x}".Eur.normID.GTflt.AB.noChrM.vqsr.flt.vcf.gz

  vep -i $pathSplit \
    --plugin dbNSFP,/project/richards/guillaume.butler-laporte/bin/dbSNFP4.0a_files/dbNSFP4.0a.gz,Ensembl_transcriptid,Uniprot_acc,VEP_canonical,LRT_pred,SIFT_pred,MutationTaster_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred \
    -everything \
    --buffer_size 100000 \
    --force_overwrite \
    --offline \
    --fork 10 \
    --dir_cache $pathCache \
    --cache -o "${pathAnnot}"finalAnnot.chr${x}.txt;
done

