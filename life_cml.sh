samtools view -h basecaller_results/IonXpress_071_rawlib.basecaller.bam |perl -ne 'if(/^@/){print $_;}else{my @arr = split /\t/;next if length($arr[9])<=50;print qq{$_};}' | samtools view -b -  > IonXpress_071_gt50.bam
tmap mapall -n 4 -f /results/analysis/output/tmap-f3/pathogen_reference_accurate_20200312/pathogen_reference_accurate_20200312.fasta -r IonXpress_071_gt50.bam -i bam -s IonXpress_071_rawlib.realigned.bam -v -Y -u --prefix-exclude 5 -o 2 --do-repeat-clip --end-repair 15 --context stage1 map4
samtools sort IonXpress_071_rawlib.realigned.bam IonXpress_071_rawlib.realigned.sort
samtools index IonXpress_071_rawlib.realigned.sort.bam
bash /results/plugins/coverageAnalysis/run_coverage_analysis.sh -ag -O IonXpress_071_rawlib.realigned.sort -L pathogen -B /results/uploads/BED/67/pathogen_reference_accurate_20200312/unmerged/detail/pathogen_20200305.bed /results/analysis/output/tmap-f3/pathogen_reference_accurate_20200312/pathogen_reference_accurate_20200312.fasta IonXpress_071_rawlib.realigned.sort.bam

grep AMPL4310705 IonXpress_071_rawlib.realigned.sort.amplicon.cov.xls
samtools tview IonXpress_050_rawlib.bam /results/analysis/output/tmap-f3/pathogen_20200305reference/pathogen_20200305reference.fasta -p COVID19:22659

# tmap mapall -n 4 -f /results/analysis/output/tmap-f3/pathogen_20200305reference/pathogen_20200305reference.fasta -r basecaller_results/IonXpress_071_rawlib.basecaller.bam -i bam -s IonXpress_071_rawlib.realigned.bam -v -Y -u --prefix-exclude 5 -o 2 -J 25 --end-repair 15 --do-repeat-clip --context stage1 map4
