# reformat genotype matrix to fit MSTmap

cov=0.05

cp s2_genotype_contig_seq.txt s2_genotype_contig_seq_${cov}x.txt


# common header
cat common_header.txt          > s2_genotype_contig_seq_${cov}x_reformat.txt
# number of gametes
M=$(grep 'sequence' s2_genotype_contig_seq_${cov}x.txt | wc -l)
echo "number_of_loci "$M       >> s2_genotype_contig_seq_${cov}x_reformat.txt
# number of markers
N=300
echo "number_of_individual "$N >> s2_genotype_contig_seq_${cov}x_reformat.txt
echo                           >> s2_genotype_contig_seq_${cov}x_reformat.txt
# marker column indication and gamate ids
echo -n -e "locus_name\t" >> s2_genotype_contig_seq_${cov}x_reformat.txt
for i in $(seq 1 $N); do echo -n -e "i"$i"\t" >> s2_genotype_contig_seq_${cov}x_reformat.txt; done
echo >> s2_genotype_contig_seq_${cov}x_reformat.txt
# marker id with sequence: P->A, M->B, U->U
grep 'sequence' s2_genotype_contig_seq_${cov}x.txt |sed 's/\tleft/_le/g' |sed 's/\tright/_ri/g' | sed 's/-sequence//g' | awk '{printf("%s\t%s\n", $2, $1)}' | sed 's/P/A\t/g' | sed 's/M/B\t/g' | sed 's/U/U\t/g' >> s2_genotype_contig_seq_${cov}x_reformat.txt
# remove redundant info 
sed -i ':a;N;$!ba;s/\t\n/\n/g' s2_genotype_contig_seq_${cov}x_reformat.txt


# linkage grouping 

../MSTmap s2_genotype_contig_seq_${cov}x_reformat.txt s2_genotype_contig_seq_${cov}x_reformat_lg_res.txt > MSTmap_${cov}x.log

# correction of MSTmap

MSTmap=s2_genotype_contig_seq_${cov}x_reformat_lg_res.txt
num_lg_to_outpu=2
marker_gt=s2_genotype_contig_seq_${cov}x.txt
max_recomb_freq=10

MSTmap_corrector ${MSTmap} ${num_lg_to_outpu} ${marker_gt} ${max_recomb_freq} > MSTmap_corrector

