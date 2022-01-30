#
gamete_genome_simulator ln_Col ln_Ler Chr1_4_CO_counts_winsize100kb_step100kb.txt 300 torm 0 > torm
grep 'bp' torm | sed 's/new_Chr//g' > 300co_smiulated.list
#
gamete_genome_simulator ln_Col ln_Ler Chr1_4_CO_counts_winsize100kb_step100kb.txt 1000 torm 0 > torm
grep 'bp' torm | sed 's/new_Chr//g' > 1000co_smiulated.list
