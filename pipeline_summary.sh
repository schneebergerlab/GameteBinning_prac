#!/bin/bash
# gamete binning practice - one run to the end
# 28 JAN 2022
# Hequan Sun

# Here is the practice of the gamete binning pipeline including simulation of data.

# Note 1: update overall path before running everything.
wd=/netscratch/dep_mercier/grp_schneeberger/projects/Method/gamete_binning_book_chapter_v9_ms_torm_test_zenodo/
cd ${wd}
# Note 2: update "--sample 10 --barcode 11" required by "del_marker_genotyper" at line 298
#         Example: path "/netscratch/dep_mercier/grp_schneeberger/projects/Method/gamete_binning_book_chapter_v9_ms/phased_assembly/s4_gamete_sequencing_simulation/gam_read_snpr0.004/gamete_0/0.1/" have 11 levels of subfolders, where 1st: netscratch, 2nd: dep_mercier, until sample info at 10th level, i.e., "gamete_0", while coverage/barcode info at 11st level, i.e., "0.1"

########################################################################################################################

#### create the working folder structure
mkdir data            # some dependency data should be collected here
mkdir -p software/bin # where all binaries should be installed
mkdir phased_assembly # where the assembly would be performed
phased_assembly=${wd}"/phased_assembly"

#### set up environment variable
export PATH=${wd}/software/bin:${PATH}
#### Binaries of all tools should be collected to the same folder. To conveniently call each tool, we set up the system variable,
export PATH=/software/bin/:${PATH}

#### global parameters
mutation_rate=0.004
cov=0.1
N=500

#### All necessary data would be simulated as below - HiFi, initial assembly and single-gamete sequencings
####     if you want to simulate different cases with the global parameters, comment lines 35 and 153 to turn on simulation.

#### <<simulaton_code

#### 3.1 Simulation of heterozygous diploid genome. 
####     1. Create working directory for genome simulation
cd ${phased_assembly}
mkdir s1_genome_simulation
cd s1_genome_simulation
####     2. Set up parameters, and create a sub-folder to collect result. Here, the two chromosomes 1 (30.4 Mb) and 4 (18.6 Mb) of Arabidopsis genome are used as haplotype A, and the sequences of haplotype B are simulated by randomly altering haplotype A with a base-level mutation rate of 0.4%.
mutation_rate=0.004
mkdir snpr${mutation_rate}
cd snpr${mutation_rate}
wd1=${phased_assembly}/s1_genome_simulation/snpr${mutation_rate}
ref_genome=${phased_assembly}/../data/Col.chr.all.v2.0.fasta
chrs_selected=Chr1,Chr4
####     3. Simulate a heterozygous diploid genome: (~1 minute)
genome_simulation.sh ${wd1} ${ref_genome} ${chrs_selected} ${mutation_rate} > genome_simulation.log
####     4. Update names of the chromosomes of haplotype A and B (30 seconds)
cd ${wd1}/ss2_F1_genome_creation
sed -i 's/>/>new_/g' F1_Col_hapA.fa
sed -i 's/>/>new_/g' F1_Ler_hapB.fa


#### 3.2 Simulation of long reads from the heterozygous genome. 
####     1. Create a directory where result is collected
cd ${phased_assembly}
mkdir s1p_1_pacbio_simulation
cd s1p_1_pacbio_simulation
####     2. Create a director to collect result for the respective mutation rate,
mutation_rate=0.004
mkdir pb_snpr${mutation_rate}
cd pb_snpr${mutation_rate}
####     3. Simulate long reads (with a mean length of 5 kb) from two haplotypes, 30x per haplotype.
refAchr=${phased_assembly}/s1_genome_simulation/snpr${mutation_rate}/ss2_F1_genome_creation/F1_Col_hapA.fa
refBchr=${phased_assembly}/s1_genome_simulation/snpr${mutation_rate}/ss2_F1_genome_creation/F1_Ler_hapB.fa
refA=Col
refB=Ler
pbsim --prefix ${refA}_HiFi_sim --data-type CLR --length-min 3000 --length-max 7000 --accuracy-min 0.98 --accuracy-max 1.00 --model_qc ${phased_assembly}/../data/model_qc_clr --length-mean 5000 --length-sd 2000 --accuracy-mean 0.99 --accuracy-sd 0.01 --depth 30 ${refAchr}
pbsim --prefix ${refB}_HiFi_sim --data-type CLR --length-min 3000 --length-max 7000 --accuracy-min 0.98 --accuracy-max 1.00 --model_qc ${phased_assembly}/../data/model_qc_clr --length-mean 5000 --length-sd 2000 --accuracy-mean 0.99 --accuracy-sd 0.01 --depth 30 ${refBchr}
####     4. Update names of the reads (such that from a read name we can judge which haplotype a read is from, for the purpose of haplotyping validation later), and compress the files.
sed 's/@S/@Col1_S/g' Col_HiFi_sim_0001.fastq | sed 's/+S/+Col1_S/g' | gzip > upd_Col_HiFi_sim_0001.fastq.gz
sed 's/@S/@Col2_S/g' Col_HiFi_sim_0002.fastq | sed 's/+S/+Col2_S/g' | gzip > upd_Col_HiFi_sim_0002.fastq.gz
sed 's/@S/@Ler1_S/g' Ler_HiFi_sim_0001.fastq | sed 's/+S/+Ler1_S/g' | gzip > upd_Ler_HiFi_sim_0001.fastq.gz
sed 's/@S/@Ler2_S/g' Ler_HiFi_sim_0002.fastq | sed 's/+S/+Ler2_S/g' | gzip > upd_Ler_HiFi_sim_0002.fastq.gz

#### 3.3 Simulation of a consensus haploid genome assembly (Note 1)
####     1. Create a directory where the result would be collected
cd ${phased_assembly}
mkdir s1p_1_assembly_simulation
cd s1p_1_assembly_simulation
####     2. Create a directory to collect result for the respective mutation rate,
mutation_rate=0.004
mkdir asm_snpr${mutation_rate}
cd asm_snpr${mutation_rate}
####     3. Create a directory to results of shuffling variants
mkdir shuffled_variants
cd shuffled_variants
####     4. Select 50% alleles from haplotype B randomly and mix them with alleles in haplotype A (~1.5 mins)
for chr in 1 4; do 
seq_path=../../../s1_genome_simulation/snpr${mutation_rate}/ss1_chr_selection
col_seq=chrChr${chr}
marker=hapB_chrChr${chr}.fa.rate${mutation_rate}.marker.txt
N=$(($(wc -l ${seq_path}/${marker} | cut -d' ' -f1) / 2))
shuf -n ${N} ${seq_path}/${marker} --random-source=${seq_path}/${marker} | sort -k2,2 -k3,3n > half_${marker}
insertMarkerToFasta ${seq_path}/${col_seq}.fa half_${marker} ${col_seq}_mixed
done
####     5. Update names of chromosomes from “ChrI” to “new_ChrI” (a few seconds)
sed -i 's/>Chr1/>new_Chr1/g' chrChr1_mixed.fasta
sed -i 's/>Chr4/>new_Chr4/g' chrChr4_mixed.fasta
####     6. Combine the sequences of Chr1 and Chr4 into one file
cat chrChr1_mixed.fasta chrChr4_mixed.fasta > ../mixed_CL_hap.fa
####     7. Simulate sizes of contigs (mean length of 500 kb and a standard deviation of 200 kb - Note 2) (in seconds)
cd ${phased_assembly}/s1p_1_assembly_simulation/asm_snpr${mutation_rate}
Rscript ${phased_assembly}/../software/bin/simulate_contig_sizes.R ${mutation_rate} ./
####     8. Extract contigs from the simulated haploid genome (a few seconds)
genome=mixed_CL_hap.fa
bedtools getfasta -fi ${genome} -bed chr1_contig_interval.bed > chr1_contigs_cl_mixed.fa
bedtools getfasta -fi ${genome} -bed chr2_contig_interval.bed > chr2_contigs_cl_mixed.fa
####     9. Combine the contigs into one file (in seconds)
cat chr1_contigs_cl_mixed.fa chr2_contigs_cl_mixed.fa > simulated_contigs.fa
####     10. Update names of contigs (in seconds)
sed -i 's/:/_/g' simulated_contigs.fa
sed -i 's/-/_/g' simulated_contigs.fa 
####     11. Get length of contigs (in seconds)
fasta_length simulated_contigs.fa > simulated_contigs.sizes
####     12. Index the fasta file of simulated assembly for read alignment later (~ 2 minutes)
bowtie2-build simulated_contigs.fa simulated_contigs.fa


#### 3.4 Simulation of crossovers (between haplotype-specific chromosomes) and gamete genomes
####     1. Create a directory where the result of crossover landscape data would be collected
cd ${phased_assembly}
mutation_rate=0.004
mkdir s2_coland_snpr${mutation_rate}
cd s2_coland_snpr${mutation_rate}
####     2. Simulate a limited number of crossovers but showing a consistent landscape to the real one (a few seconds)
Rscript ${phased_assembly}/../software/bin/calculate_CO_landscape_s1.R ./ ${phased_assembly}/../data/FileS2r1.csv
####     3. Create a directory where the (simulated) gamete genomes would be collected
cd ${phased_assembly}
mkdir s3_gamete_simulation
cd s3_gamete_simulation
mutation_rate=0.004
mkdir gam_snpr${mutation_rate}
cd gam_snpr${mutation_rate}
####     4. Simulate 500 haploid recombinant gamete genomes (~30 minutes)
hapA=../../s1_genome_simulation/snpr${mutation_rate}/ss2_F1_genome_creation/F1_Col_hapA.fa
hapB=../../s1_genome_simulation/snpr${mutation_rate}/ss2_F1_genome_creation/F1_Ler_hapB.fa
COland=../../s2_coland_snpr${mutation_rate}/Chr1_4_CO_counts_winsize100kb_step100kb.txt
N_gamete=500
prefix="sim"
gamete_genome_simulator ${hapA} ${hapB} ${COland} ${N_gamete} ${prefix} ${mutation_rate} 1 > sim.log
####     5. Create a directory where the (simulated) short reads of each gamete genome would be collected
cd ${phased_assembly}
mkdir s4_gamete_sequencing_simulation
cd s4_gamete_sequencing_simulation
mutation_rate=0.004
mkdir gam_read_snpr${mutation_rate}
cd gam_read_snpr${mutation_rate}

#### simulaton_code

####     6. Simulate short reads from each gamete genome at a coverage of 0.1x (~5 seconds per gamete; Note 3)
####        bsub -o pirs.log -e pirs.err -q short -R "span[hosts=1] rusage[mem=100]" -M 1000 "pirs simulate -i ${genome} -s ${baseprofile} -m 500 -v 25 -l 100 -x ${cov} -e 0.001 -a 0 -g 0 -o gamete${gi}_${cov}x > gamete${gi}_${cov}x_sim.log"
for gi in {0..499}; do
cd ${phased_assembly}/s4_gamete_sequencing_simulation/gam_read_snpr${mutation_rate}
mkdir gamete_${gi}
cd gamete_${gi}
genome=../../../s3_gamete_simulation/gam_snpr${mutation_rate}/sim_recombinant_gamete_${gi}.fa
baseprofile=../../../../data/pirs_Profiles/Base-Calling_Profiles/humNew.PE100.matrix.gz
cov=0.1
pirs simulate -i ${genome} -s ${baseprofile} -m 500 -v 25 -l 100 -x ${cov} -e 0.001 -a 0 -g 0 -o gamete${gi}_${cov}x > gamete${gi}_${cov}x_sim.log
done
####


#### 3.5 Definition of SNP/deletion-like markers using short read alignments (pooled from single-gamete sequencings) to the haploid genome assembly
####     1. Create a folder where the SNP markers would be created
mutation_rate=0.004
cd ${phased_assembly}/s1p_1_assembly_simulation/asm_snpr${mutation_rate}
mkdir marker_definition
cd marker_definition
####     2. Build the read set (leading to 50x per haplotype) 
cat ../../../s4_gamete_sequencing_simulation/gam_read_snpr${mutation_rate}/gamete_{0..499}/gamete*_0.1x_100_500_1.fq.gz > gamete0to499_0.1x_100_500_1.fq.gz
cat ../../../s4_gamete_sequencing_simulation/gam_read_snpr${mutation_rate}/gamete_{0..499}/gamete*_0.1x_100_500_2.fq.gz > gamete0to499_0.1x_100_500_2.fq.gz
####     3. Define SNP markers (Note 3)
####        bsub -q multicore20 -n ${nthread} -R "span[hosts=1] rusage[mem=10000]" -M 12000 -o snp_marker.log -e snp_marker.err "run_marker_define.sh ${genome} ${R1} ${R2} ${nthread} ${outlabel} > run_marker_define.log"
R1=gamete0to499_0.1x_100_500_1.fq.gz
R2=gamete0to499_0.1x_100_500_2.fq.gz
genome=../simulated_contigs.fa
nthread=8
outlabel=20211105
run_marker_define.sh ${genome} ${R1} ${R2} ${nthread} ${outlabel} > run_marker_define.log
####     4. Create a folder where the deletion-like markers would be created (Note 4)
mutation_rate=0.004
cd ${phased_assembly}/s1p_1_assembly_simulation/asm_snpr${mutation_rate}
mkdir marker_definition_del
cd marker_definition_del
####     5. Define deletion-like markers (a few seconds)
snp_marker=../marker_definition/allele_counter_converted/20211105_converted_variant_iMQ50.txt
contig_sizes=../simulated_contigs.sizes
del_marker_finder --marker ${snp_marker} --min-del-size 2000 --chrsizes ${contig_sizes} -o del_like > define_del_isize2000.log
####     6. Update the suffix to the name of the file of deletion-like markers
mv del_like_isize2000_20211105_converted_variant_iMQ50.txt del_like_isize2000_20211105_converted_variant_iMQ50.bed
####     7. Check sequencing depth at deletion-like markers (~1 minute)
bam=../marker_definition/marker.bam
marker_bed=del_like_isize2000_20211105_converted_variant_iMQ50.bed
samtools depth -a -b ${marker_bed} ${bam} > ${marker_bed}.depth
####     8. Define deletion-like markers as “hap” or “hom” (a few seconds; Note 5)
marker_bed=del_like_isize2000_20211105_converted_variant_iMQ50.bed
cutoff=20
del_depth_finder --region ${marker_bed} --depth ${marker_bed}.depth --hap-hom-cutoff ${cutoff} -o del_like_isize2000_converted_variant > del_depth_finder


#### 3.6 Consensus calling for each gamete genome with short read alignments to the haploid genome assembly
####     1. Create a folder where read alignment and variant calling would be performed
cd ${phased_assembly}
mkdir s5_read_alignment_var_calling
cd s5_read_alignment_var_calling
mutation_rate=0.004
mkdir var_call_snpr${mutation_rate}
cd var_call_snpr${mutation_rate}
####     2. Read alignment and variant calling for each gamete genome (3~5 minutes per gamete; Note 3)
####        bsub -o gamete_consensus_calling.log -e gamete_consensus_calling.err -q short -R "span[hosts=1] rusage[mem=1000]" -M 2000 "gamete_consensus_calling.sh ${genome} ${gi} ${cov} ${nthreads} ${marker} ${contig_sizes} ${R1} ${R2} > gamete${gi}_consensus_calling.log"
cellpath=${phased_assembly}/s4_gamete_sequencing_simulation/gam_read_snpr${mutation_rate}/
for gi in {0..499}; do
    cd ${cellpath}/gamete_${gi}
    cov=0.1
    mkdir ${cov}    
    cd ${cov}
    R1=../gamete${gi}_${cov}x_100_500_1.fq.gz
    R2=../gamete${gi}_${cov}x_100_500_2.fq.gz
    genome=${phased_assembly}/s1p_1_assembly_simulation/asm_snpr${mutation_rate}/simulated_contigs.fa
    marker=${phased_assembly}/s1p_1_assembly_simulation/asm_snpr${mutation_rate}/marker_definition/allele_counter_converted/20211105_converted_variant_iMQ50.txt
    contig_sizes=${phased_assembly}/s1p_1_assembly_simulation/asm_snpr${mutation_rate}/simulated_contigs.sizes
    nthreads=1
    gamete_consensus_calling.sh ${genome} ${gi} ${cov} ${nthreads} ${marker} ${contig_sizes} ${R1} ${R2} > gamete${gi}_consensus_calling.log
    cd ..
done
####     3. Collect the consensus calling information from all gamete genomes
mutation_rate=0.004
cd ${phased_assembly}/s5_read_alignment_var_calling/var_call_snpr${mutation_rate}/
cov=0.1
N=500
mkdir final_cell_config_cov${cov}x_rep
cd final_cell_config_cov${cov}x_rep
echo ${phased_assembly}/s4_gamete_sequencing_simulation/gam_read_snpr${mutation_rate}/gamete_*/${cov}/allele_counter_converted/extracted_consensus_0.txt | sed 's/ /\n/g' > cov${cov}x_consen_${N}cells.txt


#### 3.7 Gamete binning 
####     1. Create the folder where gamete binning would carried out
mutation_rate=0.004
cov=0.1
N=500
cd ${phased_assembly}/s5_read_alignment_var_calling/var_call_snpr${mutation_rate}/final_cell_config_cov${cov}x_rep
mkdir gb_cov${cov}x_consen_${N}cells_rep
cd  gb_cov${cov}x_consen_${N}cells_rep
####     2. Phase the SNPs within each contig and generate contig markers (15 minutes)
####        bsub -q normal -o gb.log -e gb.err -R "span[hosts=1] rusage[mem=10000]" -M 12000 "gamete_binning_dip --marker ${marker} --pollen ${cells} --size ${contig_sizes} --corr --ims 0.81 -o ${outflag}_res > ${outflag}.log"
marker=../../../../s1p_1_assembly_simulation/asm_snpr${mutation_rate}/marker_definition/allele_counter_converted/20211105_converted_variant_iMQ50.txt
contig_sizes=../../../../s1p_1_assembly_simulation/asm_snpr${mutation_rate}/simulated_contigs.sizes
cells=../cov${cov}x_consen_${N}cells.txt
outflag=gb_cov${cov}x_consen_${N}cells_rep
gamete_binning_dip --marker ${marker} --pollen ${cells} --size ${contig_sizes} --corr --ims 0.81 -o ${outflag}_res > ${outflag}.log


#### 3.8 Genetic mapping (Note 6)
####     1. Enter the folder and set up parameters
mutation_rate=0.004
cov=0.1
N=500
workpath=${phased_assembly}/s5_read_alignment_var_calling/
cd ${workpath}/var_call_snpr${mutation_rate}/final_cell_config_cov${cov}x_rep/gb_cov${cov}x_consen_${N}cells_rep
GB_DIR=/${workpath}/var_call_snpr${mutation_rate}/final_cell_config_cov${cov}x_rep/gb_cov${cov}x_consen_${N}cells_rep
marker_gt=./gb_cov${cov}x_consen_${N}cells_rep_res_tmp_pollen_genotypes/s2_genotype_contig_seq.txt
MST_header=${phased_assembly}/../data/gm_common_header.txt
max_recomb_freq=40
####     2. Perform genetic mapping using genotype patterns at contig markers (a few seconds; note: cut_off_p_value in header file should be set up properly to get "expected" linkage groups)
contig_marker_phasing.sh ${GB_DIR} ${MST_header} ${marker_gt} ${max_recomb_freq} ${N} > MST_correction.log

#### 3.9 Completing the genetic map
####     1. Create the folder where the analysis is performed
cd ${phased_assembly}
mkdir s6_phasing_del_like
cd s6_phasing_del_like
####     2. Get the read coverage at each deletion-like marker for each gamete genome (a few seconds for all gamete genomes)
####        bsub -o bedtools.log -e bedtools.err -q short -R "rusage[mem=100]" -M 1000 "bedtools coverage -counts -a ${del_marker} -b gamete${gi}_${cov}x.bam -bed > ${cov}_del_like_read_count.bed"
cellpath=${phased_assembly}/s4_gamete_sequencing_simulation/gam_read_snpr${mutation_rate}/
cov=0.1
del_marker=../../../../s1p_1_assembly_simulation/asm_snpr${mutation_rate}/marker_definition_del/del_like_isize2000_20211105_converted_variant_iMQ50.bed
for gi in {0..499}; do
    cd ${cellpath}/gamete_${gi}/${cov}
    bedtools coverage -counts -a ${del_marker} -b gamete${gi}_${cov}x.bam -bed > ${cov}_del_like_read_count.bed
done
####     3. Create the folder to collect genotype patterns at deletion-like markers (--sample and --barcode give the sub-folder level where bowtie2.err for the gamete can be found)
cd ${phased_assembly}/s6_phasing_del_like
mutation_rate=0.004
cov=0.1
mkdir del_snpr${mutation_rate}_rep
cd del_snpr${mutation_rate}_rep
mkdir del_phasing_cov${cov}
cd del_phasing_cov${cov}
N=500
mkdir gb_cov${cov}x_consen_${N}cells_rep
cd gb_cov${cov}x_consen_${N}cells_rep
####     4. Build up genotype patterns across all gamete genomes at deletion-like markers according to coverage (1 minute; note: --sample gamete_i--folder-count, --barcode cov-folder-count)
leaf_depth=${phased_assembly}/s1p_1_assembly_simulation/asm_snpr${mutation_rate}/marker_definition_del/del_like_isize2000_converted_variant_del_like_interval_avg_depth_sorted.bed
cells=${phased_assembly}/s5_read_alignment_var_calling/var_call_snpr${mutation_rate}/final_cell_config_cov${cov}x_rep/cov${cov}x_consen_${N}cells.txt
out_label=20211223
del_marker_genotyper --pollen ${cells} --leaf-depth ${leaf_depth} --sample 10 --barcode 11 -o ${out_label} > ${out_label}_del_marker_genotyper.log


####     5. Incorporate markers to the genetic map comprehensively (1 second)
N=500
gmap=${phased_assembly}/s5_read_alignment_var_calling/var_call_snpr${mutation_rate}/final_cell_config_cov${cov}x_rep/gb_cov${cov}x_consen_${N}cells_rep/MSTmap_correction/final_meta_genetic_map_file_list.txt
phase=${phased_assembly}/s5_read_alignment_var_calling/var_call_snpr${mutation_rate}/final_cell_config_cov${cov}x_rep/gb_cov${cov}x_consen_${N}cells_rep/MSTmap_correction/final_meta_genetic_map_marker_phase_file_list.txt
marker_all=${phased_assembly}/s5_read_alignment_var_calling/var_call_snpr${mutation_rate}/final_cell_config_cov${cov}x_rep/gb_cov${cov}x_consen_${N}cells_rep/gb_cov${cov}x_consen_${N}cells_rep_res_tmp_pollen_genotypes/s2_genotype_contig_seq.txt
del_like=./20211223_tmp_pollen_del_like_genotypes_addi/s2_genotype_contig_seq_del_like.txt
cellcut=0.55
asCaffolder --map ${gmap} --phase ${phase} --marker ${marker_all} --marker-del-like ${del_like} --hap2hom ${cellcut} -o phased_0p55
####     6. Set up meta information of files with phase information for later references
mutation_rate=0.004
N=500
cov=0.1
wd=${phased_assembly}/s6_phasing_del_like/del_snpr${mutation_rate}_rep/del_phasing_cov${cov}/gb_cov${cov}x_consen_${N}cells_rep/z_genetic_maps_updated_with_PMsimilarity_of_snp_plus_del_like_contigs_${cellcut}
cd ${wd}
echo ${wd}/upd_final_corrected_lg*_map.txt | sed 's/ /\n/g' > final_meta_genetic_map_file_list.txt
echo ${wd}/PM_phase_updfinal_corrected_lg_*_map.txt | sed 's/ /\n/g' >final_meta_genetic_map_marker_phase_file_list.txt


#### 3.10 Separation of long reads to haplotype-specific chromosomes
####     1. Align long reads to the haploid assembly (half an hour with 8 threads)
####        bsub -o hifi_align.log -e hifi_align.err -q multicore20 -n ${nthread} "minimap2 -ax map-pb -t ${nthread} -N 1 --secondary=no ${genome} ${simHiFi} > long_read_alignment.sam"
mutation_rate=0.004
cd ${phased_assembly}/s1p_1_assembly_simulation/asm_snpr${mutation_rate}
mkdir HiFi_alignment
cd HiFi_alignment
simHiFi="../../../s1p_1_pacbio_simulation/pb_snpr${mutation_rate}/*.fastq.gz"
nthread=8
genome=../simulated_contigs.fa
minimap2 -ax map-pb -t ${nthread} -N 1 --secondary=no ${genome} ${simHiFi} > long_read_alignment.sam
####     2. Create the folder where the separation of long reads would be carried out
mutation_rate=0.004
N=500
cov=0.1
wd=${phased_assembly}/s6_phasing_del_like/del_snpr${mutation_rate}_rep/del_phasing_cov${cov}/gb_cov${cov}x_consen_${N}cells_rep
cd ${wd}
mkdir pbread_separation_cov${cov}x
cd pbread_separation_cov${cov}x	
####     3. Separate long reads to haplotypes using phased variations (15 minutes)
sam=${phased_assembly}/s1p_1_assembly_simulation/asm_snpr${mutation_rate}/HiFi_alignment/long_read_alignment.sam
snp_phased=${phased_assembly}/s5_read_alignment_var_calling/var_call_snpr${mutation_rate}/final_cell_config_cov${cov}x_rep/gb_cov${cov}x_consen_${N}cells_rep/gb_cov${cov}x_consen_${N}cells_rep_res_tmp_pollen_genotypes/s4_phased_markers.txt
del_phased=../phased_0p55_s2_genotype_contig_seq_del_like.txt
cellcut=0.55
MSTphase=../z_genetic_maps_updated_with_PMsimilarity_of_snp_plus_del_like_contigs_${cellcut}/final_meta_genetic_map_marker_phase_file_list.txt
long_read_genotyper --sam ${sam} --marker ${snp_phased} --marker2 ${del_phased} --phase ${MSTphase} --ims 0.9 -o sim_pb_separation_with_dels_0p55 > sim_pb_separation_with_dels_op55.log 
####     4. Check how reads have been separated (Note 8)
####        1-(966+20+999+19+138+17+14+82)/(966+20+999+19+138+17+14+82+187464+187569+116904+117007)=0.9963
cd sim_pb_separation_with_dels_0p55_snp_marker_separated_pbreads/
echo "Haplotyping accuracy checking using reads information: "
for gt in MMM PPP; do 
for lg in 1 2; do
   echo "# lg "${lg}" of "${gt}
   grep '>' PM_phase_updfinal_corrected_lg_${lg}_map.txt_${gt}_pbreads.fa | grep 'snp-phased' | sed 's/_/ /g' | cut -d' ' -f1,2 | sort | uniq -c
done
done


#### 3.11 Independent assembly with each group of reads (Notes 3 & 9). Add: careful on "Illegal instruction", if hifiasm compiled on dell, run on dell; compiled on hpc, run on hpc.
mutation_rate=0.004
N=500
cov=0.1
wd=${phased_assembly}/s6_phasing_del_like/del_snpr${mutation_rate}_rep/del_phasing_cov${cov}/gb_cov${cov}x_consen_${N}cells_rep/pbread_separation_cov${cov}x
cd ${wd}/sim_pb_separation_with_dels_0p55_snp_marker_separated_pbreads/
for gt in MMM PPP; do 
for lg in 1 2; do
   mkdir hifiasm_lg${lg}_${gt}
   cd hifiasm_lg${lg}_${gt}
   hifiasm -t 4 -l2 -o lg_wise ../PM_phase_updfinal_corrected_lg_${lg}_map.txt_${gt}_pbreads.fa > hifiasm_${fa}.log
   cd ..
done
done

