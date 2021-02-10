#this is a python script loosely based on Kumar and Konrad's effort here: https://github.com/mkveerapen/covid19_sequencing
#again, some of the QC at our institution was done by our genome center, and therefore you should refer to the above link for more thorough QC
#specifically, variant recalibration should still be done, even if not shown here, can discuss with me on how to do it using gatk.

import hail as hl
import hail.expr.aggregators as agg

#the different paths below

# this is where some of the temporary computations are done. I would make sure to assign it to a folder that does not have a strict data cap.
pathTMP = '/scratch/richards/guillaume.butler-laporte/tmpSort/'
hl.init(spark_conf=None, tmp_dir= pathTMP)

# this is the path to your vcf from step 1
pathNorm = '/scratch/richards/guillaume.butler-laporte/WGS/normalizedVCF/allSamples.normID.noChrM.vqsr.flt.vcf.gz'

# this is the path of the matrix table hail will create
pathMT = '/scratch/richards/guillaume.butler-laporte/WGS/hailFiles/hail.full.normID.noChrM.mt'

# this is path to the QCed vcf file output. Important that the end of the file be ".vcf.bgz" otherwise hail will not provide correct format
pathQC = '/scratch/richards/guillaume.butler-laporte/WGS/QCed/allSamples.normID.GTflt.AB.noChrM.vqsr.flt.vcf.bgz'



#import the data and sample QC
hl.import_vcf(pathNorm, min_partitions=4, reference_genome='GRCh38', force_bgz=True).write(pathMT, overwrite=True)
metaData = hl.get_vcf_metadata(pathNorm)

mtAll = hl.read_matrix_table(pathMT)
mtAll= mtAll.annotate_entries(AB = (mtAll.AD[1] / hl.sum(mtAll.AD) ))
mtAll=hl.sample_qc(mtAll)
mtAll = mtAll.filter_cols((mtAll.sample_qc.call_rate >= 0.97) & (mtAll.sample_qc.dp_stats.mean >= 20))
mtAll = mtAll.filter_entries( (mtAll.GQ>=20) &
                 (mtAll.DP >= 10) &
                 ((mtAll.GT.is_hom_ref() & (mtAll.AB <= 0.1)) |
                        (mtAll.GT.is_het() & (mtAll.AB >= 0.25) & (mtAll.AB <= 0.75)) |
                        (mtAll.GT.is_hom_var() & (mtAll.AB >= 0.9))))
                        
hl.export_vcf(mtAll, pathQC, metadata=metaData, tabix=True)
