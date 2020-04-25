#load functions
setwd("/home/shriyas/cifs-lab/RIP_manuscript/Revised_Figures/Fig5/Fig5B/")
source('heatmap_tools_v07.r')
#Calculate scaled profile from CSV files
#If need to plot more than one profile, combine them by cbind
profile_scaled_nrpd <- calculate_scaled_data("temp/Col_sense_input.csv", "temp/Nrpd_sense_input.csv", "temp/LENGTHS.csv", 3.35,1.49)
profile_scaled_nrpe <- calculate_scaled_data("temp/Col_nrpe_control_input.csv", "temp/Nrpe_sense_input.csv", "temp/LENGTHS.csv", 1.28,0.59)
profile_scaled_rip <- calculate_scaled_data("~/cifs-lab/RIP_manuscript/Revised_Figures/Fig1/Fig1B/input1.csv", "~/cifs-lab/RIP_manuscript/Revised_Figures/Fig1/Fig1B/input2.csv", "~/cifs-lab/RIP_manuscript/Revised_Figures/Fig1/Fig1B/LENGTHS.csv", 74.14, 60.14)


#Plot the profile in the PDF format FOR FIG S5B
pdf("../../FigS5/FigS5B/Nrpd_profile_scaled.pdf", width=8, height=5)
#       line      box      line   scale   no axes     colors      title      y labels   scale
matplot(colMeans(profile_scaled_nrpd, na.rm=TRUE), type="l", bty="n", lwd=2, log="", axes=FALSE, col="maroon", main=NULL, ylab=NULL)
#matplot(cbind(colMeans(profile_scaled_nrpe/max(colMeans(profile_scaled_nrpe, na.rm=TRUE)), na.rm=TRUE),colMeans(profile_scaled_nrpd/max(colMeans(profile_scaled_nrpd, na.rm=TRUE)), na.rm=TRUE),colMeans(profile_scaled_rip/max(colMeans(profile_scaled_rip, na.rm=TRUE)), na.rm=TRUE), type="l", bty="n", lwd=1, log="", axes=FALSE, col=c(2,4,6), main=NULL, ylab=NULL))

#p. tick length (p: 1-bottom, 2-left, 4-right)
axis(2, tck=0.015, ylab=NULL, lwd=0.5)
#p.  locations of all ticks
axis(1, c(0, 300, 800, 1100), tck=-0.01, lwd=0.5)
#abline(h=0, col="grey")
dev.off()

#Plot the heatmap in the PNG format
reduced <- profile_scaled_nrpd[seq(1, nrow(profile_scaled_nrpd), 2),]						#keeps every N-th sequence	
png("../../FigS5/FigS5B/Nrpd_scaled_heatmap.png", width=ncol(reduced), height=nrow(reduced), units="px")	
#par(mar=c(0,0,0,0))
breaks <- matrix(NA, nrow=nrow(reduced), ncol=2)
breaks[,1] <- 300																				#create matrix with column breaks
breaks[,2] <- 800		
my_palette <- colorRampPalette(c("#fee8c8", "#fdbb84", "#e34a33"))(n=299)				#define color palette - shades of red
col_breaks = c(seq(quantile(reduced, 0.10, na.rm=TRUE),quantile(reduced, 0.90, na.rm=TRUE),length=300))	#define color breaks using percentiles
heatmap.2(reduced, col=my_palette, trace="none", breaks=col_breaks,dendrogram="none", Colv="NA", Rowv="NA", 
          key=FALSE, labRow=FALSE, labCol=FALSE,  colsep=breaks, sepcolor="grey", lwid=c(0,1), lhei=c(0,1), margins=c(0,0))
dev.off()

pdf("../../FigS5/FigS5C/Nrpe_profile_scaled.pdf", width=8, height=5)
#       line      box      line   scale   no axes     colors      title      y labels   scale
matplot(colMeans(profile_scaled_nrpe, na.rm=TRUE), type="l", bty="n", lwd=2, log="", axes=FALSE, col="maroon", main=NULL, ylab=NULL)
#matplot(cbind(colMeans(profile_scaled_nrpe/max(colMeans(profile_scaled_nrpe, na.rm=TRUE)), na.rm=TRUE),colMeans(profile_scaled_nrpd/max(colMeans(profile_scaled_nrpd, na.rm=TRUE)), na.rm=TRUE),colMeans(profile_scaled_rip/max(colMeans(profile_scaled_rip, na.rm=TRUE)), na.rm=TRUE), type="l", bty="n", lwd=1, log="", axes=FALSE, col=c(2,4,6), main=NULL, ylab=NULL))

#p. tick length (p: 1-bottom, 2-left, 4-right)
axis(2, tck=0.015, ylab=NULL, lwd=0.5)
#p.  locations of all ticks
axis(1, c(0, 300, 800, 1100), tck=-0.01, lwd=0.5)
#abline(h=0, col="grey")
dev.off()

#Plot the heatmap in the PNG format
reduced <- profile_scaled_nrpe[seq(1, nrow(profile_scaled_nrpe), 2),]						#keeps every N-th sequence	
png("../../FigS5/FigS5C/Nrpe_scaled_heatmap.png", width=ncol(reduced), height=nrow(reduced), units="px")	
#par(mar=c(0,0,0,0))
breaks <- matrix(NA, nrow=nrow(reduced), ncol=2)
breaks[,1] <- 300																				#create matrix with column breaks
breaks[,2] <- 800		
my_palette <- colorRampPalette(c("#fee8c8", "#fdbb84", "#e34a33"))(n=299)				#define color palette - shades of red
col_breaks = c(seq(quantile(reduced, 0.10, na.rm=TRUE),quantile(reduced, 0.90, na.rm=TRUE),length=300))	#define color breaks using percentiles
heatmap.2(reduced, col=my_palette, trace="none", breaks=col_breaks,dendrogram="none", Colv="NA", Rowv="NA", 
          key=FALSE, labRow=FALSE, labCol=FALSE,  colsep=breaks, sepcolor="grey", lwid=c(0,1), lhei=c(0,1), margins=c(0,0))
dev.off()


pdf("Nrpd_profile_scaled.pdf", width=8, height=5)
#       line      box      line   scale   no axes     colors      title      y labels   scale
x=colMeans(profile_scaled_nrpd, na.rm=TRUE)
matplot((x-min(x))/(max(x)-min(x)), type="l", bty="n", lwd=2, log="", axes=FALSE, col="maroon", main=NULL, ylab=NULL)
#matplot(cbind(colMeans(profile_scaled_nrpe/max(colMeans(profile_scaled_nrpe, na.rm=TRUE)), na.rm=TRUE),colMeans(profile_scaled_nrpd/max(colMeans(profile_scaled_nrpd, na.rm=TRUE)), na.rm=TRUE),colMeans(profile_scaled_rip/max(colMeans(profile_scaled_rip, na.rm=TRUE)), na.rm=TRUE), type="l", bty="n", lwd=1, log="", axes=FALSE, col=c(2,4,6), main=NULL, ylab=NULL))

#p. tick length (p: 1-bottom, 2-left, 4-right)
axis(2, tck=0.015, ylab=NULL, lwd=0.5)
#p.  locations of all ticks
axis(1, c(0, 300, 800, 1100), tck=-0.01, lwd=0.5)
#abline(h=0, col="grey")
dev.off()

pdf("Nrpe_profile_scaled.pdf", width=8, height=5)
#       line      box      line   scale   no axes     colors      title      y labels   scale
x=colMeans(profile_scaled_nrpe, na.rm=TRUE)
matplot((x-min(x))/(max(x)-min(x)), type="l", bty="n", lwd=2, log="", axes=FALSE, col="green", main=NULL, ylab=NULL)
#matplot(cbind(colMeans(profile_scaled_nrpe/max(colMeans(profile_scaled_nrpe, na.rm=TRUE)), na.rm=TRUE),colMeans(profile_scaled_nrpd/max(colMeans(profile_scaled_nrpd, na.rm=TRUE)), na.rm=TRUE),colMeans(profile_scaled_rip/max(colMeans(profile_scaled_rip, na.rm=TRUE)), na.rm=TRUE), type="l", bty="n", lwd=1, log="", axes=FALSE, col=c(2,4,6), main=NULL, ylab=NULL))

#p. tick length (p: 1-bottom, 2-left, 4-right)
axis(2, tck=0.015, ylab=NULL, lwd=0.5)
#p.  locations of all ticks
axis(1, c(0, 300, 800, 1100), tck=-0.01, lwd=0.5)
#abline(h=0, col="grey")
dev.off()

pdf("Rip_profile_scaled.pdf", width=8, height=5)
#       line      box      line   scale   no axes     colors      title      y labels   scale
x=colMeans(profile_scaled_rip, na.rm=TRUE)
matplot((x-min(x))/(max(x)-min(x)), type="l", bty="n", lwd=2, log="", axes=FALSE, col="blue", main=NULL, ylab=NULL)
#matplot(cbind(colMeans(profile_scaled_nrpe/max(colMeans(profile_scaled_nrpe, na.rm=TRUE)), na.rm=TRUE),colMeans(profile_scaled_nrpd/max(colMeans(profile_scaled_nrpd, na.rm=TRUE)), na.rm=TRUE),colMeans(profile_scaled_rip/max(colMeans(profile_scaled_rip, na.rm=TRUE)), na.rm=TRUE), type="l", bty="n", lwd=1, log="", axes=FALSE, col=c(2,4,6), main=NULL, ylab=NULL))

#p. tick length (p: 1-bottom, 2-left, 4-right)
axis(2, tck=0.015, ylab=NULL, lwd=0.5)
#p.  locations of all ticks
axis(1, c(0, 300, 800, 1100), tck=-0.01, lwd=0.5)
#abline(h=0, col="grey")
dev.off()




#Calculate significance of differences using unscaled profile
profile_unscaled <- calculate_unscaled_data("temp/Col_sense_input.csv", "temp/Nrpd_sense_input.csv", "temp/LENGTHS.csv", 3.35,1.49)
sink("../../FigS5/FigS5B/sense_stats.txt")
print(significance(profile_unscaled, 100, 350, 100, 1000))	#Options: data, start of region 1, start of region 2, length of regions, number of permutations
print(significance(profile_unscaled, 100, 350, 100, 10000))	#Options: data, start of region 1, start of region 2, length of regions, number of permutations
sink()

profile_unscaled <- calculate_unscaled_data("temp/Col_nrpe_control_input.csv", "temp/Nrpe_sense_input.csv", "temp/LENGTHS.csv", 1.28,0.59)
sink("../../FigS5/FigS5C/Nrpe_stats.txt")
print(significance(profile_unscaled, 100, 350, 100, 1000))	#Options: data, start of region 1, start of region 2, length of regions, number of permutations
print(significance(profile_unscaled, 100, 350, 100, 10000))	#Options: data, start of region 1, start of region 2, length of regions, number of permutations
sink()