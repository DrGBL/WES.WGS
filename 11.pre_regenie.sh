cd /scratch/richards/guillaume.butler-laporte/WGS/tmp/

for chr in {1..22} X; do
  cat pLoF.chr"${chr}".* moderate.non.missense.chr"${chr}".* missense.5in5.chr"${chr}".* missense.1in5.chr"${chr}".* missense.0in5.chr"${chr}".* > regenie.anno.file.chr"${chr}".txt
  cat above.1perc.chr"${chr}".* deleterious.below.0.1perc.chr"${chr}".* deleterious.0.1.to.1perc.chr"${chr}".* | sort -u -k 1.6 > regenie.aaf.file.tmp.chr"${chr}".txt
  awk 'NR==FNR{a[$1]=$1;next} ($1 in a)' \
    regenie.anno.file.chr"${chr}".txt \
    regenie.aaf.file.tmp.chr"${chr}".txt > regenie.aaf.file.chr"${chr}".txt
done

cat regenie.aaf.file.chr* > /scratch/richards/guillaume.butler-laporte/WGS/regenieInputs/regenie.aaf.file.txt

cat pLoF.chr* moderate.non.missense.chr* missense.5in5.chr* missense.1in5.chr* missense.0in5.chr* > /scratch/richards/guillaume.butler-laporte/WGS/regenieInputs/regenie.anno.file.txt

set list!
cat regenie.set.list.chr* > regenie.set.list.tmp.txt
awk 'NR==FNR{a[$2]=$1;next} ($1 in a)' \
  /scratch/richards/guillaume.butler-laporte/WGS/regenieInputs/regenie.anno.file.txt \
  regenie.set.list.tmp.txt > /scratch/richards/guillaume.butler-laporte/WGS/regenieInputs/regenie.set.list.txt





rm pLoF.chr* moderate.non.missense.chr* missense.5in5.chr* missense.1in5.chr* missense.0in5.chr*
rm above.1perc.chr* deleterious.below.0.1perc* deleterious.0.1.to.1perc* regenie.aaf.file.chr*
rm regenie.set.list.chr* regenie.anno.file.chr* regenie.aaf.file.tmp.chr*



