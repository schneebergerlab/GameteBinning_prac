#!/bin/bash
# contig_marker_phasing.sh
# correct MSTmap-generated genetic map, output phase info of contigs - need inputs from MSTmap and gamete_binning_dip
# 19 Nov 2021
# Hequan Sun

# usage

if [ $# -eq 0 ]
  then
echo "contig_marker_phasing.sh <1.path to gamete binning folder> <2.common_header> <3.marker genotype> <4.max_recomb_freq> <5.population_size>"
echo "   "
echo "   1. work path / gamete binning result folder "
echo "   2. common header required by MSTmap "
echo "   3. s2_genotype_contig_seq.txt from gamete binning process"
echo "   4. percentage of difference between two contig marker genotypes"
echo "   5. size of population: here it is the number of gamete genomes."
echo "   "
exit 1
fi

# step 0. get variables

# work path
GB_DIR=$1
# Common header required by MSTmap 
MST_Header=$2
# s2_genotype_contig_seq.txt from gamete binning process
marker_gt=$3
# maximum differences between patterns of contig end markers - percentage as [0, 100]
max_recomb_freq=$4

# step 1. reformat genotype matrix to fit MSTmap

cd ${GB_DIR}
mkdir MSTmap_correction
cd MSTmap_correction
marker_gt=../${marker_gt}

# common header: common_header.txt
cat ${MST_Header}               > s2_genotype_contig_seq_reformat.txt
# number of markers
M=$(grep 'sequence' ${marker_gt} | wc -l)
echo "number_of_loci "$M       >> s2_genotype_contig_seq_reformat.txt
# number of gametes
N=$5
echo "number_of_individual "$N >> s2_genotype_contig_seq_reformat.txt
echo                           >> s2_genotype_contig_seq_reformat.txt
# marker column indication and gamate ids
echo -n -e "locus_name\t" >> s2_genotype_contig_seq_reformat.txt
for i in $(seq 1 $N); do echo -n -e "i"$i"\t" >> s2_genotype_contig_seq_reformat.txt; done
echo >> s2_genotype_contig_seq_reformat.txt
# marker id with sequence: P->A, M->B, U->U
grep 'sequence' ${marker_gt} |sed 's/\tleft/_le/g' |sed 's/\tright/_ri/g' | sed 's/-sequence//g' | awk '{printf("%s\t%s\n", $2, $1)}' | sed 's/P/A\t/g' | sed 's/M/B\t/g' | sed 's/U/U\t/g' >> s2_genotype_contig_seq_reformat.txt
# remove redundant info 
sed -i ':a;N;$!ba;s/\t\n/\n/g' s2_genotype_contig_seq_reformat.txt


# step 2. linkage group with MSTmap

MSTmap s2_genotype_contig_seq_reformat.txt s2_genotype_contig_seq_reformat_lg_res.txt > MSTmap.log

# step 3. correct MSTmap

MSTmap=s2_genotype_contig_seq_reformat_lg_res.txt
num_lg_to_outpu=2
MSTmap_corrector ${MSTmap} ${num_lg_to_outpu} ${marker_gt} ${max_recomb_freq} > MSTmap_corrector.log

# step 4. list result 

echo 
echo "Info: results in folder ${GB_DIR}/MSTmap_correction: "
echo
echo "Result 1: linkage group wise genetic map: "
echo
ls -lrt final_corrected_lg_*_map.txt
echo 
echo "Result 2: markers with gentyopes at individuals, corresponding to genetic map: "
echo
ls -lrt final_corrected_lg_*_marker_gt_phase.txt
echo 

# step 5. prepare files for later usage: these files are not used for long_read_genotyper anymore, because the markers they have are limited. But the following files are required by asCaffolder. 20211128 to generate updated map and phasing info of markers.

wd=${GB_DIR}/MSTmap_correction
echo ${wd}/final_corrected_lg_*_map.txt | sed 's/ /\n/g' > final_meta_genetic_map_file_list.txt
echo ${wd}/final_corrected_lg_*_marker_gt_phase.txt | sed 's/ /\n/g' > final_meta_genetic_map_marker_phase_file_list.txt

echo 
echo "Result 3: files prepared for asCaffolder: "
echo
ls -lrt final_meta_genetic_map*.txt
echo 






