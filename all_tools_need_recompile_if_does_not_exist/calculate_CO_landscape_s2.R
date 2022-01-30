#!/bin/env Rscript

# Here we create the CO frequency along chromosomes.
## get options from cmd
###############################################################################
args<-commandArgs(TRUE)
if(length(args) < 4) {
  cat(length(args))
  cat("\n   Info: create the CO frequency along chromosomes.")
  cat("   USAGE: Rscript calculate_CO_landscape_s2.R work_path CO_list.csv 500co_smiulated.list 1000co_smiulated.list\n\n")
}else
{
  options(scipen=9990)
  # parameters
  wd=args[[1]]
  setwd(wd)
  #
  co_file=args[[2]]  
  file_500co=args[[3]]
  file_1000co=args[[4]]
  # Here we create the CO frequency along chromosomes.
  #
  window_counter<-function(raw_coordinate, chr_size, win_size, step_size)
  {
    # raw_coordinate: raw coordinates of potential breakpoints along a certain chromosome
    # chr_size: in bp
    # win_size: in bp
    n_chrco <- length(raw_coordinate)
    len <- floor((chr_size-win_size)/step_size)
    rem <- (chr_size-win_size)%%step_size
    if(rem > 0)
    {
      len <- len + 1
    }
    count <- matrix(0, len+1, 4)
    
    for(istep in c(1:(len+1)))
    {
      winsta       <- (istep-1)*step_size + 1
      winend       <- min(winsta+win_size-1, chr_size)
      count[istep, 1] <- winsta
      count[istep, 2] <- winend
      count[istep, 3] <- round(length(raw_coordinate[raw_coordinate>=winsta & raw_coordinate<=winend])/n_chrco*100, digits = 3)
      count[istep, 4] <- length(raw_coordinate[raw_coordinate>=winsta & raw_coordinate<=winend])
    }
    return(count)
  }
  #
  pdf(paste("Suppl_Fig_Rowan_et_al_2019_slidwin_s2_real_sim.pdf", sep=""), family="Helvetica", height=2.5, width=8.26772/2)
  par(ps=8)
  #
  m <- rbind(c(1,1,1),
             c(2,2,2) )
  layout(m, heights = c(1, 1))
  par(mai = c(0.3, 0.5, 0.1, 0.2)) # margin: bottom, left, top, right
  # general info
  chrsize  = c(30427671,19698289,23459830,18585056,26975502)
  chrtick  = list(c(0,5000000,10000000,15000000,20000000,25000000,30427671),
                  c(0,5000000,10000000,15000000,19698289),
                  c(0,5000000,10000000,15000000,20000000,23459830),
                  c(0,5000000,10000000,15000000,18585056),
                  c(0,5000000,10000000,15000000,20000000,26975502) )
  chrlabel = list(c("0.0 (Mb)","5.0","10.0","15.0","20.0","25.0","30.4"),
                  c("0.0 (Mb)","5.0","10.0","15.0","19.7"),
                  c("0.0 (Mb)","5.0","10.0","15.0","20.0","23.5"),
                  c("0.0 (Mb)","5.0","10.0","15.0","18.6"),
                  c("0.0 (Mb)","5.0","10.0","15.0","20.0","27.0")   )
  # non-overlapping 100 kb windows
  windowsize <- 100000
  stepsize   <- 100000
  # CO list from Rowan et al 2019
  coInfo         <- read.csv(co_file)
  coInfo_sim500  <- read.table(file_500co)
  coInfo_sim1000 <- read.table(file_1000co)
  #
  coInfo_sim500  <- data.frame(coInfo_sim500[coInfo_sim500$V1!="bp_breakpoint_noCO", c(1,3,4,5,6,7)])
  coInfo_sim1000 <- data.frame(coInfo_sim1000[coInfo_sim1000$V1!="bp_breakpoint_noCO", c(1,3,4,5,6,7)])
  #
  colnames(coInfo_sim500) <- paste("V", 1:6, sep="")
  colnames(coInfo_sim1000) <- paste("V", 1:6, sep="")
  #
  all_count <- c()
  for (chr in c(1,4))
  {
    #### real
    chrcoinfo  <- coInfo[coInfo$chr==chr, ]
    xxx        <- chrcoinfo$breakpoint.pos
    wincount   <- window_counter(sort(xxx), chrsize[chr], windowsize, stepsize) # 100,000bp window
    all_count  <- rbind(all_count, cbind(rep(paste("new_Chr", chr, sep=""), length(wincount[,1])), wincount))  # chr name changed 
    ####
    ylimt      <- 4
    plot( (wincount[, 1]+wincount[,2])/2, 
          wincount[, 3], 
          type="l",
          lwd=1.5,
          xlim = c(1, chrsize[1]), 
          ylim=c(0, ylimt),
          #xlab=paste(" chr ", chr, " pos", sep=""),
          ylab="",
          frame.plot=FALSE,
          xaxt = "n",
          yaxt = "n",
          col="orangered")  
    #### sim 300
    chrcoinfo_sim  <- coInfo_sim500[coInfo_sim500$V3==chr, ]
    xxx_sim        <- chrcoinfo_sim$V6
    wincount_sim   <- window_counter(sort(xxx_sim), chrsize[chr], windowsize, stepsize) # 100,000bp window
    lines(x = (wincount_sim[, 1]+wincount_sim[,2])/2, y = wincount_sim[, 3], type = "l", col="deepskyblue", lwd=0.5)
    #### sim 1000
    chrcoinfo_sim  <- coInfo_sim1000[coInfo_sim1000$V3==chr, ]
    xxx_sim        <- chrcoinfo_sim$V6
    wincount_sim   <- window_counter(sort(xxx_sim), chrsize[chr], windowsize, stepsize) # 100,000bp window
    lines(x = (wincount_sim[, 1]+wincount_sim[,2])/2, y = wincount_sim[, 3], type = "l", col="yellowgreen", lwd=0.5)  
    ####
    axis(1, at=c(seq(0,chrsize[chr], 5000000), chrsize[chr]), labels = NA, cex.axis=1.5, tck=-0.05)
    axis(2, at=seq(0, ylimt, ylimt/2), labels = NA, cex.axis=1.5, tck=-0.05)
    tickxx <- seq(0,chrsize[chr], 1000000)
    axis(1, tickxx, rep("", length(tickxx)), tck=-0.03)
    axis(1, at=chrtick[[chr]], labels = chrlabel[[chr]], cex.axis=1.5, lwd = 0, line = -.7)
    axis(2, at=seq(0, ylimt, ylimt/2), labels = seq(0, ylimt, ylimt/2), cex.axis=1.5, lwd = 0, line = -.6)
    ####
    legend("topright", pch = c(15, 15), col=c("orangered", "deepskyblue", "yellowgreen"), legend = c("Real", "Sim 500 CO", "Sim 1000 CO"), border = NA, bty="n")
  }
  # close pdf
  dev.off()
  # output counts
  options(scipen=9990)
  write.table(x = all_count, file = "Chr1_4_CO_counts_winsize100kb_step100kb.txt", append = F, quote = F, sep = "\t", row.names = F, col.names = F)
}




