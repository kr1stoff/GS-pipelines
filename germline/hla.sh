mkdir -p /media/gsadmin/vd2/data/wes/GST0001/WES_hg19_work_dir/1.process/GST0001/HLA /media/gsadmin/vd2/data/wes/GST0001/WES_hg19_work_dir/2.result/GST0001/result_variation/HLA
razers3 -i 95 -m 1 -dr 0 -o /media/gsadmin/vd2/data/wes/GST0001/WES_hg19_work_dir/1.process/GST0001/HLA/GST0001.1.bam /media/gsadmin/vd1/heguangliang/pipeline/WGS/software/OptiType/data/hla_reference_dna.fasta /media/gsadmin/vd2/data/wes/GST0001/Rawdata/TEST20200612_1.fq.gz && /media/gsadmin/vd1/heguangliang/pipeline/WGS/bin/samtools fastq  /media/gsadmin/vd2/data/wes/GST0001/WES_hg19_work_dir/1.process/GST0001/HLA/GST0001.1.bam >/media/gsadmin/vd2/data/wes/GST0001/WES_hg19_work_dir/1.process/GST0001/HLA/GST0001.1.hla.fastq && rm /media/gsadmin/vd2/data/wes/GST0001/WES_hg19_work_dir/1.process/GST0001/HLA/GST0001.1.bam &razers3 -i 95 -m 1 -dr 0 -o /media/gsadmin/vd2/data/wes/GST0001/WES_hg19_work_dir/1.process/GST0001/HLA/GST0001.2.bam /media/gsadmin/vd1/heguangliang/pipeline/WGS/software/OptiType/data/hla_reference_dna.fasta /media/gsadmin/vd2/data/wes/GST0001/Rawdata/TEST20200612_2.fq.gz && /media/gsadmin/vd1/heguangliang/pipeline/WGS/bin/samtools fastq  /media/gsadmin/vd2/data/wes/GST0001/WES_hg19_work_dir/1.process/GST0001/HLA/GST0001.2.bam >/media/gsadmin/vd2/data/wes/GST0001/WES_hg19_work_dir/1.process/GST0001/HLA/GST0001.2.hla.fastq && rm /media/gsadmin/vd2/data/wes/GST0001/WES_hg19_work_dir/1.process/GST0001/HLA/GST0001.2.bam && wait
OptiTypePipeline.py -i /media/gsadmin/vd2/data/wes/GST0001/WES_hg19_work_dir/1.process/GST0001/HLA/GST0001.1.hla.fastq /media/gsadmin/vd2/data/wes/GST0001/WES_hg19_work_dir/1.process/GST0001/HLA/GST0001.2.hla.fastq --dna -v -o /media/gsadmin/vd2/data/wes/GST0001/WES_hg19_work_dir/1.process/GST0001/HLA -p GST0001
cp /media/gsadmin/vd2/data/wes/GST0001/WES_hg19_work_dir/1.process/GST0001/HLA/GST0001* /media/gsadmin/vd2/data/wes/GST0001/WES_hg19_work_dir/2.result/GST0001/result_variation/HLA
echo HLA GST0001 finished >>/media/gsadmin/vd2/data/wes/GST0001/WES_hg19_work_dir/finishedSampleAnalysis.log 
