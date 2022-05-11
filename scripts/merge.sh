#!/usr/bin/zsh
s1=$1
s2=$2
s3=$3
out_dir=$4
if [ ! -d "$out_dir" ]; then
  mkdir $out_dir
fi
#s1
bcftools head $s1 > $out_dir/header_s1.sv.vcf
cp $out_dir/header_s1.sv.vcf $out_dir/header_s1.snp.vcf
bcftools view $s1 -H | grep "RNAMES=" >$out_dir/s1.sv.vcf
bcftools view $s1 -H | grep -v "RNAMES=" >$out_dir/s1.snp.vcf
cat $out_dir/s1.sv.vcf >> $out_dir/header_s1.sv.vcf
cat $out_dir/s1.snp.vcf >> $out_dir/header_s1.snp.vcf
#bgzip $out_dir/header_s1.sv.vcf
bgzip $out_dir/header_s1.snp.vcf
#tabix $out_dir/header_s1.sv.vcf.gz
tabix $out_dir/header_s1.snp.vcf.gz

#s2
bcftools head $s2 > $out_dir/header_s2.sv.vcf
cp $out_dir/header_s2.sv.vcf $out_dir/header_s2.snp.vcf
bcftools view $s2 -H | grep "RNAMES=" >$out_dir/s2.sv.vcf
bcftools view $s2 -H | grep -v "RNAMES=" >$out_dir/s2.snp.vcf
cat $out_dir/s2.sv.vcf >> $out_dir/header_s2.sv.vcf
cat $out_dir/s2.snp.vcf >> $out_dir/header_s2.snp.vcf
#bgzip $out_dir/header_s2.sv.vcf
bgzip $out_dir/header_s2.snp.vcf
#tabix $out_dir/header_s2.sv.vcf.gz
tabix $out_dir/header_s2.snp.vcf.gz


#s3
bcftools head $s3 > $out_dir/header_s3.sv.vcf
cp $out_dir/header_s3.sv.vcf $out_dir/header_s3.snp.vcf
bcftools view $s3 -H | grep "RNAMES=" >$out_dir/s3.sv.vcf
bcftools view $s3 -H | grep -v "RNAMES=" >$out_dir/s3.snp.vcf
cat $out_dir/s3.sv.vcf >> $out_dir/header_s3.sv.vcf
cat $out_dir/s3.snp.vcf >> $out_dir/header_s3.snp.vcf
#bgzip $out_dir/header_s3.sv.vcf
bgzip $out_dir/header_s3.snp.vcf
#tabix $out_dir/header_s3.sv.vcf.gz
tabix $out_dir/header_s3.snp.vcf.gz

touch $out_dir/sur.input
echo $out_dir/header_s1.sv.vcf >> $out_dir/sur.input
echo $out_dir/header_s2.sv.vcf >> $out_dir/sur.input
echo $out_dir/header_s3.sv.vcf >> $out_dir/sur.input

bcftools merge -Oz $out_dir/header_s1.snp.vcf.gz $out_dir/header_s2.snp.vcf.gz $out_dir/header_s3.snp.vcf.gz -o $out_dir/snps.vcf.gz
/home/caronkey/docs/cityu/SURVIVOR/Debug/SURVIVOR merge $out_dir/sur.input 0.2 1 1 1 0 10 $out_dir/svs.vcf
tabix $out_dir/snps.vcf.gz
bcftools sort $out_dir/svs.vcf -Oz -o $out_dir/svs.sorted.vcf.gz
#bgzip $out_dir/svs.vcf
tabix $out_dir/svs.sorted.vcf.gz


bcftools concat --threads 4 $out_dir/snps.vcf.gz $out_dir/svs.sorted.vcf.gz -Oz -a -o $out_dir/all.vcf.gz
