#!/bin/bash
# genome_simulation.sh
# simulate a heterozygous genome with (2) linkage groups
# 01 Oct 2021
# Hequan Sun

# usage

if [ $# -eq 0 ]
  then
echo  
echo "Usage: genome_simulation.sh <1.work path> <2.haploid_genome.fasta> <3.Chr1,Chr4> <4.mutation rate>"
echo "   "
echo "   1. work path"
echo "   2. haploid genome, for example, we can use Col-0 reference genome here"
echo "   3. chrs to keep in the simulation, here we only use Chr1,Chr4 to save space and time"
echo "   4. rate of single-nucleotide polymorphisms to be created between homologous chromosomes"
echo "   "
exit 1
fi

# sub-tools: 
# 	fasta_name_selecter
# 	fasta_length
# 	fasta_separator
#       mutate_dna

echo 
echo
echo "##############################################################################################################" 
echo "Genome (or, two haplotype-specific genomes) simulation starting ..."
current_date=`date`
echo ${current_date}

# step 0. get variables

echo "##############################################################################################################" 
echo "set up work path: $1" #  wd=$(pwd)
wd=$1
echo "get initialized haploid genome: $2" # here we use Col-0 genome: /biodata/dep_mercier/grp_schneeberger/data/Athal/AMPRIL/Col.chr.all.v2.0.fasta
ref_genome=$2
echo "get chrs to extract: $3" #: Chr1,Chr4
chrs_selected=$3
echo "set up expected mutation rate: $4" # 0.005
mutation_rate=$4
echo "##############################################################################################################" 

# step 1. prepare HapA (Col-0) and HapB ("Ler") genomes, where B is simply mutated from A genomes.
# ....... select chrs for creating "F1" heterozygous genome 

echo
echo "##############################################################################################################" 
cd ${wd}
mkdir ss1_chr_selection
cd ss1_chr_selection
cleaned_fa=`ls *.fa *.fasta`
if [ -z "$cleaned_fa" ]; then 
        echo "no cleaning of existing fasta needed."
    else 
        echo "warning: existing fasta files under folder ss1_chr_selection would be cleaned: "
        echo ${cleaned_fa}
        rm ${cleaned_fa}
fi
echo "##############################################################################################################" 

echo
echo "##############################################################################################################" 
echo "check chr lengths in ${ref_genome} "
# fasta_length ${ref_genome} | grep '>' | sed 's/>//g' > col.sizes
fasta_length ${ref_genome} > col.sizes
cat col.sizes
echo "##############################################################################################################" 

echo
echo "##############################################################################################################" 
echo "select chrs to simulate: " # Chr1,Chr4
#
>col_chr_selected1_for_unique
for i in $(echo $chrs_selected | tr "," "\n")
do
  echo -e $i >> col_chr_selected1_for_unique
done
cat col_chr_selected1_for_unique
#
fasta_name_selecter ${ref_genome} col_chr_selected1_for_unique
targt_fasta=$(ls *.fasta)
mv ${targt_fasta} hapA_${targt_fasta}
#
echo "size of selected chrs: "
#fasta_length hapA_${targt_fasta} | grep '>' | sed 's/>//g'
fasta_length hapA_${targt_fasta} 
echo "##############################################################################################################" 
#
echo
echo "##############################################################################################################" 
echo "get hapA chrs in fa: "
#
fasta_separator hapA_${targt_fasta}
hapA_chr_fasta=`ls chr*.fa`
#
echo "get hapA chrs in fa done. "
echo "##############################################################################################################" 
#
echo 
echo "##############################################################################################################" 
echo "creat hapB chrs in fasta by mutating hapA chrs..."
#
for fa in ${hapA_chr_fasta}; do 
   mutate_dna ${fa} ${mutation_rate} hapB_${fa}
done
hapB_chr_fasta=`ls hapB_*.rate${mutation_rate}.fasta`
echo "creat hapB chrs in fasta by mutating hapA chrs done."
echo "##############################################################################################################" 

# step 2. merge chrs of haps A and B as haploid genomes
#
cd ${wd}
echo
echo "##############################################################################################################" 
mkdir ss2_F1_genome_creation
cd ss2_F1_genome_creation
echo
echo "merge chr-wise sequences from the same hap..."
>F1_Col_hapA.fa
>F1_Ler_hapB.fa
for i in $(echo $chrs_selected | tr "," "\n")
do
  cat ../ss1_chr_selection/chr${i}.fa >> F1_Col_hapA.fa
  cat ../ss1_chr_selection/hapB_chr${i}.fa.rate${mutation_rate}.fasta >> F1_Ler_hapB.fa
done
echo "merge chr-wise sequences from the same hap done."
echo "Info: haplotype A (Col-0) sequences given by F1_Col_hapA.fa with length: "
# fasta_length F1_Col_hapA.fa | grep '>' | sed 's/>//g'
fasta_length F1_Col_hapA.fa
#
echo "Info: haplotype B (\"Ler\") sequences given by F1_Ler_hapB.fa with length: "
# fasta_length F1_Ler_hapB.fa | grep '>' | sed 's/>//g' 
fasta_length F1_Ler_hapB.fa
echo "##############################################################################################################" 
echo
#
# 
echo "Genome (or, two haplotype-specific genomes) simulation done."
current_date=`date`
echo ${current_date}
echo 
echo







