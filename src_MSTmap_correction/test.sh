
MSTmap=s2_genotype_contig_seq_0.50x_reformat_lg_res.txt
num_lg_to_outpu=2
marker_gt=s2_genotype_contig_seq_0.50x.txt
max_recomb_freq=20

MSTmap_corrector ${MSTmap} ${num_lg_to_outpu} ${marker_gt} ${max_recomb_freq} > torm
