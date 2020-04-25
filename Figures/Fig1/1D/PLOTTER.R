setwd("/home/shriyas/cifs-lab/RIP_manuscript/Revised_Figures/Fig1/Fig1D/temp_files/")
#load functions
Col_CHIP <- read.table("CHIP_input1.csv", sep = ",", header=FALSE)
Nrpe_CHIP <- read.table("CHIP_input2.csv", sep = ",", header=FALSE)
Col_RIP <- read.table("RIP_input1.csv", sep = ",", header=FALSE)
Nrpe_RIP <- read.table("RIP_input2.csv", sep = ",", header=FALSE)

result_chip <- (((Col_CHIP / 28.86) + (0.5 / 28.86)) / ((Nrpe_CHIP / 18.1) + (0.5 / 28.86)))
result_rip <- (((Col_RIP / 74.14) + (0.5 / 74.14)) / ((Nrpe_RIP / 60.14) + (0.5 / 74.14)))

# result_chip <- (Col_CHIP+0.01)/((Nrpe_CHIP+0.01)*(28.86/18.1))
# result_rip <- (Col_RIP+0.01)/((Nrpe_RIP+0.01)*(71.14/60.14))

CHIP_avg <- colMeans(result_chip[,2:1202])
RIP_avg <- colMeans(result_rip[,2:1202])

plot_rip <- RIP_avg[seq(1, length(RIP_avg), 10)]
plot_chip <- CHIP_avg[seq(1, length(CHIP_avg), 10)]

  pdf("../Fig1D.pdf", width=8, height=5)
  matplot(cbind(plot_chip, plot_rip), type="l", bty="n", lwd=2, log="", axes=FALSE, col=c("blue","maroon"), main=NULL, ylab=NULL)
  axis(2, tck=0.015, ylab=NULL, lwd=0.5)
  #p.  locations of all ticks
  axis(1, c(0, 60, 120), tck=-0.01, lwd=0.5)
  dev.off()
