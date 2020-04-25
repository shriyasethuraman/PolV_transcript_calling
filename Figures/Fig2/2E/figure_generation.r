#load functions
source('heatmap_tools_v07.r')
setwd("/home/shriyas/cifs-lab/RIP_manuscript/Revised_Figures/Fig2/GUDRUN_2E/")
#Calculate scaled profile from CSV files
#If need to plot more than one profile, combine them by cbind
profile_scaled <- calculate_scaled_data("299_TSS-PolV-tsct_8995_H3K4me2.csv", "299_TSS-PolV-tsct_8995_H3.csv", "300_TSS-PolV-tsct_8994_sizes.csv", 699.60, 8167.81)
profile_scaled2 <- calculate_scaled_data("299_TSS-PolV-tsct_8995_H3K9me2.csv", "299_TSS-PolV-tsct_8995_H3.csv", "300_TSS-PolV-tsct_8994_sizes.csv", 8004.31, 8167.81)

#Plot the profile in the PDF format
pdf("profile_scaled.pdf", width=8, height=5)
#       line      box      line   scale   no axes     colors      title      y labels   scale
x=colMeans(profile_scaled, na.rm=TRUE)
y=colMeans(profile_scaled2, na.rm=TRUE)
matplot(cbind((x-min(x))/(max(x)-min(x)),(y-min(y))/(max(y)-min(y))), type="l", bty="n", lwd=1, log="", axes=FALSE, col=c("red","blue"), main=NULL, ylab=NULL)
#p. tick length (p: 1-bottom, 2-left, 4-right)
axis(2, tck=0.015, ylab=NULL, lwd=0.5)
#p.  locations of all ticks
axis(1, c(0, 300, 800, 1100), tck=-0.01, lwd=0.5)
#abline(h=0, col="grey")
dev.off()

pdf("H3K4me2_profile_scaled.pdf", width=8, height=5)
#       line      box      line   scale   no axes     colors      title      y labels   scale
matplot((x-min(x))/(max(x)-min(x)), type="l", bty="n", lwd=1, log="", axes=FALSE, col=c(2,4), main=NULL, ylab=NULL)
#p. tick length (p: 1-bottom, 2-left, 4-right)
axis(2, tck=0.015, ylab=NULL, lwd=0.5)
#p.  locations of all ticks
axis(1, c(0, 300, 800, 1100), tck=-0.01, lwd=0.5)
#abline(h=0, col="grey")
dev.off()

pdf("H3K9me2_profile_scaled.pdf", width=8, height=5)
#       line      box      line   scale   no axes     colors      title      y labels   scale
matplot((y-min(y))/(max(y)-min(y)), type="l", bty="n", lwd=1, log="", axes=FALSE, col="blue", main=NULL, ylab=NULL)
#p. tick length (p: 1-bottom, 2-left, 4-right)
axis(2, tck=0.015, ylab=NULL, lwd=0.5)
#p.  locations of all ticks
axis(1, c(0, 300, 800, 1100), tck=-0.01, lwd=0.5)
#abline(h=0, col="grey")
dev.off()

#Calculate significance of differences using unscaled profile
profile_unscaled <- calculate_unscaled_data("299_TSS-PolV-tsct_8995_H3K4me2.csv", "299_TSS-PolV-tsct_8995_H3.csv", "300_TSS-PolV-tsct_8994_sizes.csv", 699.60, 8167.81)
sink("H3K4me2_stats.txt")
print(significance(profile_unscaled, 350, 100, 100, 1000))	#Options: data, start of region 1, start of region 2, length of regions, number of permutations
print(significance(profile_unscaled, 350, 100, 100, 10000))	#Options: data, start of region 1, start of region 2, length of regions, number of permutations
sink()

profile_unscaled_anti <- calculate_unscaled_data("299_TSS-PolV-tsct_8995_H3K9me2.csv", "299_TSS-PolV-tsct_8995_H3.csv", "300_TSS-PolV-tsct_8994_sizes.csv", 8004.31, 8167.81)
sink("H3K9me2_stats.txt")
print(significance(profile_unscaled_anti, 100, 350, 100, 1000))		#Options: data, start of region 1, start of region 2, length of regions, number of permutations
print(significance(profile_unscaled_anti, 100, 350, 100, 10000))		#Options: data, start of region 1, start of region 2, length of regions, number of permutations
sink()
