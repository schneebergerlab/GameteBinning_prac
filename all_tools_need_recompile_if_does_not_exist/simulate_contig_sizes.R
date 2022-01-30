#!/bin/env Rscript

# Here we simulate sizes of contigs along chromosomes, where a bed file is generated for each chr for fasta extraction.
## get options from cmd
###############################################################################
args<-commandArgs(TRUE)
if(length(args) < 2) {
  cat(length(args))
  cat("\n   Info: simulate sizes of contigs along selected chromosomes, where a bed file is generated for each chr for fasta extraction.")
  cat("   USAGE: Rscript simulate_contig_sizes.R random_seed work_path\n")
  cat("          random_seed: to randomnize generated size of contigs. \n\n")
}else
{
  #
  rseed=as.numeric(args[[1]])
  wd=args[[2]]
  setwd(wd)
  #
  generate_bed<-function(chrid, chrsize, contig_sizes)
  {
    contig_numb <- length(contig_sizes)
    contig_size <- matrix(data = NA, nrow = contig_numb, ncol = 3) # chr (sta end]
    winsta      <- 1
    for(ci in 1:contig_numb)
    {
      contig_size[ci, 1] <- chrid
      contig_size[ci, 2] <- winsta-1
      contig_size[ci, 3] <- winsta+contig_sizes[ci]-1
      winsta             <- winsta+contig_sizes[ci]-1+1
    }
    return(contig_size)
  }
  #
  set.seed(rseed*1000/2)
  # Expected chr sizes
  chrsizes <- c(30427671,18585056)
  chr1_cnt <- round(30427671/sum(chrsizes)*100)
  chr2_cnt <- round(18585056/sum(chrsizes)*100)
  # Geneated contig sizes
  chr1_contig_size <- round(rnorm(chr1_cnt, mean=500000, sd=200000))
  chr2_contig_size <- round(rnorm(chr2_cnt, mean=500000, sd=200000))
  # 
  sum(chr1_contig_size) # 32763871 > 30427671
  sum(chr2_contig_size) # 19413874 > 18585056
  # Correct contig size
  correct_chr1        <- round((sum(chr1_contig_size) - chrsizes[1])/chr1_cnt)
  correct_chr2        <- round((sum(chr2_contig_size) - chrsizes[2])/chr2_cnt)
  chr1_contig_size    <- abs(chr1_contig_size - correct_chr1)
  chr2_contig_size    <- abs(chr2_contig_size - correct_chr2)
  chr1_contig_size[1] <- chr1_contig_size[1] - (sum(chr1_contig_size) - chrsizes[1])
  chr2_contig_size[1] <- chr2_contig_size[1] - (sum(chr2_contig_size) - chrsizes[2])
  # Now as expected
  sum(chr1_contig_size) # 30427671 vs. 30427671
  sum(chr2_contig_size) # 18585056 vs. 18585056
  #
  # generate 
  chr1_contig_interval<- generate_bed("new_Chr1", chrsize[1], chr1_contig_size)
  chr2_contig_interval<- generate_bed("new_Chr4", chrsize[1], chr2_contig_size)
  # check
  ax<-as.numeric(chr1_contig_interval[, 3])-as.numeric(chr1_contig_interval[,2])
  bx<-as.numeric(chr2_contig_interval[, 3])-as.numeric(chr2_contig_interval[,2])
  sum(ax)
  sum(bx)
  # output
  write.table(x = chr1_contig_interval, file = paste(wd, "chr1_contig_interval.bed", sep=""), append = F, quote = F, row.names = F, col.names = F, sep = "\t")
  write.table(x = chr2_contig_interval, file = paste(wd, "chr2_contig_interval.bed", sep=""), append = F, quote = F, row.names = F, col.names = F, sep = "\t")
}
