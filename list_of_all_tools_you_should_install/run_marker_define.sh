#!/bin/bash
# run_marker_define.sh
# find snp markers using Illumina short reads
# 01 Oct 2021
# Hequan Sun

# usage

if [ $# -ne 5 ]
  then
echo  
echo "Usage: run_marker_define.sh <1.genome.fa> <2.R1.reads.fq.gz> <3.R2.reads.fq.gz> <4.nthreads> <5.output_label>"
echo "   "
echo "   1-5. self-explained"
echo "   "
exit 1
fi
#
#
wd=$(pwd)
cd ${wd}
#
echo
echo 
echo "##############################################################################################################" 
echo "Defining markers starting..."
current_date=`date`
echo ${current_date}

echo "##############################################################################################################" 
echo "collect reference genome: $1"
hap_genome=$1
echo "collect R1 reads: $2" 
R1=$2
echo "collect R2 reads: $3" 
R2=$3
echo "set up threads for bowtie2 alignment: $4"
nthreads=$4
echo "get some output label"
olabel=$5
echo "##############################################################################################################" 
echo
#
echo -e "*\t*\t*\t*\t2" > ploidy
echo "aligning reads using bowtie2 with ${nthreads} threads..."
bowtie2 -p ${nthreads} -x ${hap_genome} -1 ${R1} -2 ${R2} 2> bowtie2.err | samtools view -@ ${nthreads} -bS - | samtools sort -@ ${nthreads} -o marker.bam -
echo "aligning reads done."
ls -l marker.bam
echo "getting bam depth..."
samtools depth -Q 1 -a marker.bam > marker_depth.txt 
samtoolsDepthHisto marker_depth.txt >samtoolsDepthHisto.log
echo "getting bam depth done."
echo
#
# -A necessary 
echo "mpileuping alignments and variant calling..."
bcftools mpileup -A -Oz -o marker_mpileup.gz -f ${hap_genome} marker.bam
samtools index marker.bam
echo "bam indexed."
bcftools call -m -v -Ov --ploidy-file ploidy marker_mpileup.gz > marker.vcf
ls -l marker.vcf
rm marker_mpileup.gz
echo "mpileuping alignments and variant calling done."
echo
#
echo "formatting marker.vcf into plain text..."
allele_counter convert --marker marker.vcf --folder allele_counter_converted -runid ${olabel} -no-r -no-c >convert.log
ls -l allele_counter_converted
cd allele_counter_converted/
awk '$6>=50' ${olabel}_converted_variant.txt > ${olabel}_converted_variant_iMQ50.txt
echo "selected variants in ${olabel}_converted_variant_iMQ50.txt."
echo "formatting done."
#
echo
echo "Defining markers done."
current_date=`date`
echo ${current_date}
echo "##############################################################################################################" 
