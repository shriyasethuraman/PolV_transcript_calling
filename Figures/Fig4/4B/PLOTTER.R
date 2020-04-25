setwd("/home/shriyas/cifs-lab/RIP_manuscript/Revised_Figures/Fig4/Fig4B/")
system("Fig4B_BASH") #needs correction

profile_scaled <- calculate_scaled_data("sense_temp_files/input1.csv", "sense_temp_files/input2.csv", "sense_temp_files/LENGTHS.csv", 39.33, 27.61)
profile_scaled2 <- calculate_scaled_data("sense_temp_files/input3.csv", "sense_temp_files/input2.csv", "sense_temp_files/LENGTHS.csv", 33, 27.61)
profile_scaled3 <- calculate_scaled_data("sense_temp_files/input4.csv", "sense_temp_files/input2.csv", "sense_temp_files/LENGTHS.csv", 41.44, 27.61)

pdf("col0_profile_scaled.pdf", width=8, height=5)
#       line      box      line   scale   no axes     colors      title      y labels   scale
matplot(colMeans(profile_scaled, na.rm=TRUE), type="l", bty="n", lwd=1, log="", axes=FALSE, col=c(2,4), main=NULL, ylab=NULL)
#p. tick length (p: 1-bottom, 2-left, 4-right)
axis(2, tck=0.015, ylab=NULL, lwd=0.5)
#p.  locations of all ticks
axis(1, c(0, 300, 800, 1100), tck=-0.01, lwd=0.5)
#abline(h=0, col="grey")
dev.off()

pdf("ago4_profile_scaled.pdf", width=8, height=5)
#       line      box      line   scale   no axes     colors      title      y labels   scale
matplot(colMeans(profile_scaled2, na.rm=TRUE), type="l", bty="n", lwd=1, log="", axes=FALSE, col=c(2,4), main=NULL, ylab=NULL)
#p. tick length (p: 1-bottom, 2-left, 4-right)
axis(2, tck=0.015, ylab=NULL, lwd=0.5)
#p.  locations of all ticks
axis(1, c(0, 300, 800, 1100), tck=-0.01, lwd=0.5)
#abline(h=0, col="grey")
dev.off()

pdf("idn2_profile_scaled.pdf", width=8, height=5)
#       line      box      line   scale   no axes     colors      title      y labels   scale
matplot(colMeans(profile_scaled3, na.rm=TRUE), type="l", bty="n", lwd=1, log="", axes=FALSE, col=c(2,4), main=NULL, ylab=NULL)
#p. tick length (p: 1-bottom, 2-left, 4-right)
axis(2, tck=0.015, ylab=NULL, lwd=0.5)
#p.  locations of all ticks
axis(1, c(0, 300, 800, 1100), tck=-0.01, lwd=0.5)
#abline(h=0, col="grey")
dev.off()

pdf("ALL_profile_scaled.pdf", width=8, height=5)
matplot(cbind(colMeans(profile_scaled, na.rm=TRUE),colMeans(profile_scaled2, na.rm=TRUE),colMeans(profile_scaled3, na.rm=TRUE)), type="l", bty="n", lwd=1, log="", axes=FALSE, col=c(2,4), main=NULL, ylab=NULL)
axis(2, tck=0.015, ylab=NULL, lwd=0.5)
#p.  locations of all ticks
axis(1, c(0, 300, 800, 1100), tck=-0.01, lwd=0.5)
#abline(h=0, col="grey")
dev.off()

profile_unscaled1 <- calculate_unscaled_data("sense_temp_files/input1.csv", "sense_temp_files/input2.csv", "sense_temp_files/LENGTHS.csv", 39.33, 27.61)
sink("Col-Ago_stats.txt")
print(significance_compare(profile_unscaled2, profile_unscaled1, 350, 100, 1000))	#Options: data, start of region 1, start of region 2, length of regions, number of permutations
print(significance_compare(profile_unscaled2, profile_unscaled1, 350, 100, 10000))	#Options: data, start of region 1, start of region 2, length of regions, number of permutations
sink()

profile_unscaled2 <- calculate_unscaled_data("sense_temp_files/input3.csv", "sense_temp_files/input2.csv", "sense_temp_files/LENGTHS.csv", 33, 27.61)
sink("Col-Idn_stats.txt")
print(significance_compare(profile_unscaled3, profile_unscaled1, 350,100, 1000))	#Options: data, start of region 1, start of region 2, length of regions, number of permutations
print(significance_compare(profile_unscaled3, profile_unscaled1, 350,100, 10000))	#Options: data, start of region 1, start of region 2, length of regions, number of permutations
sink()

profile_unscaled3 <- calculate_unscaled_data("sense_temp_files/input4.csv", "sense_temp_files/input2.csv", "sense_temp_files/LENGTHS.csv", 41.44, 27.61)
sink("Ago-Idn_stats.txt")
print(significance_compare(profile_unscaled2, profile_unscaled3, 350, 100, 1000))	#Options: data, start of region 1, start of region 2, length of regions, number of permutations
print(significance_compare(profile_unscaled2, profile_unscaled3, 350, 100, 10000))	#Options: data, start of region 1, start of region 2, length of regions, number of permutations
sink()

sink("Compare_stats.txt")
print(significance_compare(profile_scaled, profile_scaled2, 550, 100, 100))		#Options: data, start of region 1, start of region 2, length of regions, number of permutations
print(significance_compare(profile_scaled, profile_scaled3, 550, 100, 100))		#Options: data, start of region 1, start of region 2, length of regions, number of permutations
print(significance_compare(profile_scaled2, profile_scaled3, 550, 100, 100))		#Options: data, start of region 1, start of region 2, length of regions, number of permutations
sink()
