gmap=./data/final_meta_genetic_map_file_list.txt
phase=./data/final_meta_genetic_map_marker_phase_file_list.txt
marker_all=./data/s2_genotype_contig_seq.txt
del_like=./data/s2_genotype_contig_seq_del_like.txt
cellcut=0.55

asCaffolder_dip --map ${gmap} --phase ${phase} --marker ${marker_all} --marker-del-like ${del_like} --hap2hom ${cellcut} -o phased_0p55
