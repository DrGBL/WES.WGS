#this annotates each chromosome using VEP and dbNSFP
#again, this could be sped up using PBS arrays if your cluster permits
#it assumes that vep was installed to work offline, and with a cache directory
#the --fork 10 makes it faster, but can give trouble sometimes
#here I did it on the european population, but it can be done once on the whole sample (then used in each ancestry stratified analysis).

mkdir -p annotation

for chr in {1..22} X; do
  vep -i /path/to/sequence.file.chr"${chr}".Eur.normID.GTflt.AB.noChrM.vcf.gz \
    --plugin dbNSFP,/project/richards/guillaume.butler-laporte/bin/dbSNFP4.0a_files/dbNSFP4.0a.gz,Ensembl_transcriptid,Uniprot_acc,VEP_canonical,LRT_pred,SIFT_pred,MutationTaster_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred \
    -everything \
    --buffer_size 10000 \
    --fork 10 \
    --offline \
    --dir_cache /path/to/cache/directory \
    --cache -o /annotation/finalAnnot.chr"${chr}".txt ;
done

