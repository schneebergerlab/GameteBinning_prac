#!/bin/bash
# gamete_consensus_calling.sh
# call variations against reference genome
# 25 JAN 2022
# Hequan Sun, MPIPZ
# sun@mpipz.mpg.de/sunhequan@gmail.com
#
#
# usage
if [ $# -ne 8 ]
  then
echo  
echo "Usage: gamete_consensus_calling.sh <1.haploid_genome.fasta> <2.gamete_id> <3.coverage_value> <4.threads> <5.snp_marker> <6.contigs.sizes> <7.R1.fq.gz> <8.R2.fq.gz>"
echo "   "
echo "   1. haploid genome assembled"
echo "   2. the gamete identifier, e.g., 1,2,..."
echo "   3. coverage value regarding the sequencing, e.g., 0.1,0.5 etc"
echo "   4. number of threads for read alignment"
echo "   5. SNP markers to count alleles at each gamete genome"
echo "   6. sizes of contigs in the initial assembly, format: \"chr size\" with tab-separated"
echo "   7. R1 read file"
echo "   8. R2 read file"
echo "   "
echo "Dependency: bowtie2, samtools, bcftools, allele_counter."
echo "   "
exit 1
fi
#
# sub-tools: 
# bowtie2 version 2.2.8
# samtools 1.9
# bcftools 1.9
# allele_counter
#
echo 
echo
echo "##############################################################################################################" 
echo "Variant calling for gamete $3 starting ..."
current_date=`date`
echo ${current_date}
#
# step 0. get variables
echo "##############################################################################################################" 
echo "get haploid genome as reference: $1"
hap_genome=$1
echo "set up output gamete label: $2"
gi=$2
echo "get coverage info: $3"
cov=$3
echo "set up threads: $4"
nthreads=$4
echo "get SNP markers: $5"
markers=$5
echo "get sizes of contigs: $6"
contig_sizes=$6
echo "get read files: $7,$8"
R1=$7
R2=$8
#
date_string=20211105 # to give as parameter in future
echo "##############################################################################################################" 
echo 
#
##############################
# echo "creating folders for gamete_${gi} at coverage ${cov}x.."
# cd ${cellpath}
# cd gamete_${gi}
# mkdir ${cov}
# cd ${cov}
# ll -rt ../${R1} ../${R2}
# echo "creating folders done."
#
echo "aligning reads of for gamete_${gi} at coverage ${cov}x.."
bowtie2 -p ${nthreads} -N 1 -x ${hap_genome} -1 ${R1} -2 ${R2} 2> bowtie2.err | samtools view -@ ${nthreads} -bS - | samtools sort -@ ${nthreads} -o gamete${gi}_${cov}x.bam -
echo "aligning reads done."
#
echo "calling consensus for gamete_${gi} at coverage ${cov}x.."
bcftools mpileup -A -Oz -o gamete${gi}_${cov}x_PE_mpileup.gz -f ${hap_genome} gamete${gi}_${cov}x.bam
echo -e "*\t*\t*\t*\t2" > ploidy
bcftools call -A -m -Ov --ploidy-file ./ploidy gamete${gi}_${cov}x_PE_mpileup.gz > gamete${gi}_${cov}x_PE.vcf
echo "calling consensus done. "
#
echo "reformating variant calling for gamete_${gi} at coverage ${cov}x.."
allele_counter convert --marker gamete${gi}_${cov}x_PE.vcf --folder allele_counter_converted -runid 20211105 -no-r >convert.log
rm *_PE_mpileup.gz *.vcf
echo "reformatting done. "
#
echo "counting alleles at given snp markers for gamete_${gi} at coverage ${cov}x.."
cd allele_counter_converted
allele_counter extract --marker ${markers} --chrsizes ${contig_sizes} --folder . --consen 20211105_converted_consen.txt 
rm 20211105_converted_consen.txt
echo "counting allelels done."
#
##############################
#
# 
echo "Variant calling for gamete ${gamete} done."
current_date=`date`
echo ${current_date}
echo 
echo
