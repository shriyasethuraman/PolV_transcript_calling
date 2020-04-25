#load functions
setwd("/home/shriyas/cifs-lab/RIP_manuscript/Revised_Figures/Fig5/Fig5C/")
source('heatmap_tools_v07.r')
#Calculate scaled profile from CSV files
#If need to plot more than one profile, combine them by cbind
profile_scaled_ago_chip <- calculate_scaled_data("temp/Col_CHIP_input.csv", "temp/Nrpe_CHIP_input.csv", "temp/CHIP_LENGTHS.csv", 19.2,19.33)
profile_scaled_ago_rip <- calculate_scaled_data("temp/Col_RIP_input.csv", "temp/Nrpe_RIP_input.csv", "temp/RIP_LENGTHS.csv", 20.57, 17.2)
profile_scaled_rip <- calculate_scaled_data("~/cifs-lab/RIP_manuscript/Revised_Figures/Fig1/Fig1B/input1.csv", "~/cifs-lab/RIP_manuscript/Revised_Figures/Fig1/Fig1B/input2.csv", "~/cifs-lab/RIP_manuscript/Revised_Figures/Fig1/Fig1B/LENGTHS.csv", 74.14, 60.14)

x=colMeans(profile_scaled_ago_chip, na.rm=TRUE)
y=colMeans(profile_scaled_ago_rip, na.rm=TRUE)
z=colMeans(profile_scaled_rip, na.rm=TRUE)

#Plot the profile in the PDF format
pdf("Chip_profile_scaled.pdf", width=8, height=5)
#       line      box      line   scale   no axes     colors      title      y labels   scale
matplot((x-min(x))/(max(x)-min(x)), type="l", bty="n", lwd=2, log="", axes=FALSE, col="maroon", main=NULL, ylab=NULL)
#matplot(cbind(colMeans(profile_scaled_nrpe/max(colMeans(profile_scaled_nrpe, na.rm=TRUE)), na.rm=TRUE),colMeans(profile_scaled_nrpd/max(colMeans(profile_scaled_nrpd, na.rm=TRUE)), na.rm=TRUE),colMeans(profile_scaled_rip/max(colMeans(profile_scaled_rip, na.rm=TRUE)), na.rm=TRUE), type="l", bty="n", lwd=1, log="", axes=FALSE, col=c(2,4,6), main=NULL, ylab=NULL))
#p. tick length (p: 1-bottom, 2-left, 4-right)
axis(2, tck=0.015, ylab=NULL, lwd=0.5)
#p.  locations of all ticks
axis(1, c(0, 300, 800, 1100), tck=-0.01, lwd=0.5)
#abline(h=0, col="grey")
dev.off()

pdf("Rip_profile_scaled.pdf", width=8, height=5)
#       line      box      line   scale   no axes     colors      title      y labels   scale
matplot((y-min(y))/(max(y)-min(y)), type="l", bty="n", lwd=2, log="", axes=FALSE, col="green", main=NULL, ylab=NULL)
#matplot(cbind(colMeans(profile_scaled_nrpe/max(colMeans(profile_scaled_nrpe, na.rm=TRUE)), na.rm=TRUE),colMeans(profile_scaled_nrpd/max(colMeans(profile_scaled_nrpd, na.rm=TRUE)), na.rm=TRUE),colMeans(profile_scaled_rip/max(colMeans(profile_scaled_rip, na.rm=TRUE)), na.rm=TRUE), type="l", bty="n", lwd=1, log="", axes=FALSE, col=c(2,4,6), main=NULL, ylab=NULL))

#p. tick length (p: 1-bottom, 2-left, 4-right)
axis(2, tck=0.015, ylab=NULL, lwd=0.5)
#p.  locations of all ticks
axis(1, c(0, 300, 800, 1100), tck=-0.01, lwd=0.5)
#abline(h=0, col="grey")
dev.off()

pdf("PVRip_profile_scaled.pdf", width=8, height=5)
#       line      box      line   scale   no axes     colors      title      y labels   scale
matplot((z-min(z))/(max(z)-min(z)), type="l", bty="n", lwd=2, log="", axes=FALSE, col="blue", main=NULL, ylab=NULL)
#matplot(cbind((colMeans(profile_scaled_ago_chip, na.rm=TRUE)/max(colMeans(profile_scaled_ago_chip, na.rm=TRUE))),(colMeans(profile_scaled_ago_rip, na.rm=TRUE)/max(colMeans(profile_scaled_ago_rip, na.rm=TRUE))),(colMeans(profile_scaled_rip, na.rm=TRUE)/max(colMeans(profile_scaled_rip, na.rm=TRUE)))), type="l", bty="n", lwd=1, log="", axes=FALSE, col=c("maroon","green","blue"), main=NULL, ylab=NULL)

#p. tick length (p: 1-bottom, 2-left, 4-right)
axis(2, tck=0.015, ylab=NULL, lwd=0.5)
#p.  locations of all ticks
axis(1, c(0, 300, 800, 1100), tck=-0.01, lwd=0.5)
#abline(h=0, col="grey")
dev.off()




#Calculate significance of differences using unscaled profile
profile_unscaled <- calculate_unscaled_data("temp/Col_sense_input.csv", "temp/Nrpd_sense_input.csv", "temp/LENGTHS.csv", 3.35,1.49)
sink("sense_stats.txt")
print(significance(profile_unscaled, 100, 350, 100, 1000))	#Options: data, start of region 1, start of region 2, length of regions, number of permutations
print(significance(profile_unscaled, 100, 350, 100, 10000))	#Options: data, start of region 1, start of region 2, length of regions, number of permutations
sink()

