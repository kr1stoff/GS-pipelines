mkdir -p Target_hg19_work_dir/0.shell Target_hg19_work_dir/1.process Target_hg19_work_dir/2.result
mkdir -p Target_hg19_work_dir/1.process/splitBed Target_hg19_work_dir/1.process/tumor/tmp_alignment Target_hg19_work_dir/1.process/tumor/snpindel/mutect2 Target_hg19_work_dir/1.process/tumor/snpindel/strelka Target_hg19_work_dir/1.process/tumor/snpindel/vardict Target_hg19_work_dir/1.process/tumor/snpindel/HLA Target_hg19_work_dir/1.process/tumor/cnvsv/cnvkit
mkdir -p Target_hg19_work_dir/1.process/normal/tmp_alignment Target_hg19_work_dir/1.process/normal/snpindel/mutect2 Target_hg19_work_dir/1.process/normal/snpindel/strelka Target_hg19_work_dir/1.process/normal/snpindel/vardict Target_hg19_work_dir/1.process/normal/snpindel/HLA Target_hg19_work_dir/1.process/normal/cnvsv/cnvkit /media/gsadmin/vd2/data/zl/20181229LYF/Target_hg19_work_dir/1.process/normal/HLA
mkdir -p Target_hg19_work_dir/2.result/report Target_hg19_work_dir/2.result/tumor/clean_data Target_hg19_work_dir/2.result/tumor/result_alignment Target_hg19_work_dir/2.result/tumor/result_variation/snp_indel Target_hg19_work_dir/2.result/tumor/result_variation/HLA Target_hg19_work_dir/2.result/tumor/result_variation/cnv_sv /media/gsadmin/vd2/data/zl/20181229LYF/Target_hg19_work_dir/1.process/tumor/HLA
mkdir -p Target_hg19_work_dir/2.result/normal/clean_data Target_hg19_work_dir/2.result/normal/result_alignment Target_hg19_work_dir/2.result/normal/result_variation/snp_indel Target_hg19_work_dir/2.result/normal/result_variation/HLA Target_hg19_work_dir/2.result/normal/result_variation/cnv_sv
cp /media/gsadmin/vd2/data/zl/20181229LYF/Target_hg19_work_dir/1.process/splitBed/* Target_hg19_work_dir/1.process/splitBed/
cp /media/gsadmin/vd2/data/zl/20181229LYF/Target_hg19_work_dir/final.sh Target_hg19_work_dir/
cp /media/gsadmin/vd2/data/zl/20181229LYF/Target_hg19_work_dir/0.shell/*sh Target_hg19_work_dir/0.shell
ls Target_hg19_work_dir/0.shell/*sh |xargs -n1 sed -ie 's/20181229LYF/20190329CSW/g'
ls Target_hg19_work_dir/0.shell/*sh |xargs -n1 sed -ie 's/LYF/CSW/g'
sed -ie 's/20181229LYF/20190329CSW/g' Target_hg19_work_dir/final.sh && sed -ie 's/LYF/CSW/g' Target_hg19_work_dir/final.sh
sed -ie 's/宫颈癌/??/' Target_hg19_work_dir/0.shell/align.sh #ALL, 宫颈癌, 肝癌, 肝内胆管癌, 肺癌, 胃癌, 结直肠癌, 食管和食管胃结合部癌, 胰腺癌, 乳腺癌, 卵巢癌, 甲状腺癌, 肾癌, 鼻咽癌, 膀胱癌及尿路肿瘤
