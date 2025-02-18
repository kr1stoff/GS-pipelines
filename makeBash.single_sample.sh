mkdir -p Target_hg19_work_dir/0.shell Target_hg19_work_dir/1.process Target_hg19_work_dir/2.result
mkdir -p Target_hg19_work_dir/1.process/splitBed Target_hg19_work_dir/1.process/XRL/tmp_alignment Target_hg19_work_dir/1.process/XRL/snpindel/mutect2 Target_hg19_work_dir/1.process/XRL/snpindel/vardict Target_hg19_work_dir/1.process/XRL/HLA Target_hg19_work_dir/1.process/XRL/cnvsv/cnvkit /media/gsadmin/vd2/data/zl/20190306XRL/Target_hg19_work_dir/1.process/XRL/cnvsv/lumpy/  /media/gsadmin/vd2/data/zl/20190306XRL/Target_hg19_work_dir/1.process/XRL/cnvsv/delly /media/gsadmin/vd2/data/zl/20190306XRL/Target_hg19_work_dir/1.process/XRL/cnvsv/manta /media/gsadmin/vd2/data/zl/20190306XRL/Target_hg19_work_dir/1.process/XRL/cnvsv/whamg /media/gsadmin/vd2/data/zl/20190306XRL/Target_hg19_work_dir/1.process/XRL/cnvsv/metasv/ 
mkdir -p Target_hg19_work_dir/2.result/report Target_hg19_work_dir/2.result/XRL/clean_data Target_hg19_work_dir/2.result/XRL/result_alignment Target_hg19_work_dir/2.result/XRL/result_variation/snp_indel Target_hg19_work_dir/2.result/XRL/result_variation/HLA Target_hg19_work_dir/2.result/XRL/result_variation/cnv_sv 
cp /media/gsadmin/vd2/data/zl/20181229LYF/Target_hg19_work_dir/1.process/splitBed/* Target_hg19_work_dir/1.process/splitBed/
cp /media/gsadmin/vd2/data/zl/20190306LBH/Target_hg19_work_dir/final.sh Target_hg19_work_dir/
cp /media/gsadmin/vd2/data/zl/20190306LBH/Target_hg19_work_dir/0.shell/*sh Target_hg19_work_dir/0.shell
ls Target_hg19_work_dir/0.shell/*sh |xargs -n1 sed -ie 's/20190306LBH/20190306XRL/g'
ls Target_hg19_work_dir/0.shell/*sh |xargs -n1 sed -ie 's/LBH/XRL/g'
sed -ie 's/20190306LBH/20190306XRL/g' Target_hg19_work_dir/final.sh && sed -ie 's/LBH/XRL/g' Target_hg19_work_dir/final.sh
sed -ie 's/宫颈癌/??/' Target_hg19_work_dir/0.shell/align.sh
