#this annotates each chromosome using VEP, loftee, and dbNSFP
#again, this could be sped up using PBS arrays if your cluster permits
#it assumes that vep was installed to work offline, and with a cache directory
#it also assumes so files useful for loftee (e.g. human_ancestor.fa.gz), please see the very useful loftee github for more info: https://github.com/konradjk/loftee

mkdir -p annotation

for chr in {1..22}; do
  vep -i /path/to/sequence.file.chr"${chr}".Eur.normID.rehead.GTflt.AB.noChrM.vcf.gz \
    --plugin LoF,loftee_path:/path/to/.vep/Plugins/loftee,human_ancestor_fa:/path/to/human_ancestor.fa.gz,conservation_file:/path/to/phylocsf_gerp.sql \
    --plugin dbNSFP,/project/richards/guillaume.butler-laporte/bin/dbSNFP4.0a_files/dbNSFP4.0a.gz,Ensembl_transcriptid,Uniprot_acc,VEP_canonical,LRT_pred,SIFT_pred,MutationTaster_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred \
    -everything \
    --buffer_size 10000 \
    --fork 10 \
    --offline \
    --dir_plugins /path/to/.vep/Plugins/loftee \
    --dir_cache /path/to/cache/directory \
    --cache -o /annotation/finalAnnot.chr"${chr}".txt ;
done

