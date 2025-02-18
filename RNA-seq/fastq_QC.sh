/media/gsadmin/vd1/heguangliang/bcbio3/thirdTools/bin/fastqc  -t 50 -o /media/gsadmin/vd3/xudongliang/project/20190906/data_analysis/data_eval/fastqc/45 --extract /media/gsadmin/vd3/xudongliang/project/20190906/raw_data/fastq/45_R1.fastq.gz  &>> /media/gsadmin/vd3/xudongliang/project/20190906/data_analysis/data_eval/log/45.log
cp /media/gsadmin/vd3/xudongliang/project/20190906/data_analysis/data_eval/fastqc/45/45_R1_fastqc/Images/per_base_quality.png /media/gsadmin/vd3/xudongliang/project/20190906/report/01_Quality_Statistic/Raw_Reads_Quality_Statistic/45/45.R1.qual.png
/media/gsadmin/vd1/heguangliang/bcbio3/thirdTools/bin/seqtk fqchk /media/gsadmin/vd3/xudongliang/project/20190906/raw_data/fastq/45_R1.fastq.gz > /media/gsadmin/vd3/xudongliang/project/20190906/data_analysis/data_eval/fastx/45.R1.quality.stat.txt
grep -v "ALL\|POS\|min_len" /media/gsadmin/vd3/xudongliang/project/20190906/data_analysis/data_eval/fastx/45.R1.quality.stat.txt > /media/gsadmin/vd3/xudongliang/project/20190906/data_analysis/data_eval/fastx/45.quality.stat
/media/gsadmin/vd1/heguangliang/bcbio3/anaconda/bin/Rscript /media/gsadmin/vd3/xudongliang/pipeline/RNA_seq/util/base_composition_se.R /media/gsadmin/vd3/xudongliang/project/20190906/data_analysis/data_eval/fastx/45.quality.stat /media/gsadmin/vd3/xudongliang/project/20190906/report/01_Quality_Statistic/Raw_Reads_Quality_Statistic/45/45.base.composition.pdf 45 &>> /media/gsadmin/vd3/xudongliang/project/20190906/data_analysis/data_eval/log/45.log
/media/gsadmin/vd1/heguangliang/bcbio3/anaconda/bin/Rscript /media/gsadmin/vd3/xudongliang/pipeline/RNA_seq/util/error_rate_se.R /media/gsadmin/vd3/xudongliang/project/20190906/data_analysis/data_eval/fastx/45.quality.stat /media/gsadmin/vd3/xudongliang/project/20190906/report/01_Quality_Statistic/Raw_Reads_Quality_Statistic/45/45.error.rate.pdf 45 &>> /media/gsadmin/vd3/xudongliang/project/20190906/data_analysis/data_eval/log/45.log
touch /media/gsadmin/vd3/xudongliang/project/20190906/data_analysis/data_eval/fastx/45.finish
