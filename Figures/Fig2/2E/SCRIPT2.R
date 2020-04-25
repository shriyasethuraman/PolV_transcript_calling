setwd("/home/shriyas/cifs-lab/RIP_manuscript/Revised_Figures/Fig2/Fig2E/")
system("RUN_INPUT_CREATOR_BASH") #needs correction

# H3 <- read.table("H3_input.csv", header=FALSE)
# H3K4me3 <- read.table("H3K4me3_input.csv", header=FALSE)
# H3K9Ac <- read.table("H3K9ac_input.csv", header=FALSE)
# H3K36me3 <- read.table("H3K36me3_input.csv", header=FALSE)
# H3K27me3 <- read.table("H3K27me3_input.csv", header=FALSE)
# LENGTHS <- read.table("H3_LENGTHS.csv", header=FALSE)

profile_scaled <- calculate_scaled_data("H3K4me2_input.csv", "H3_input.csv", "H3_LENGTHS.csv", 2.31, 5.8)
profile_scaled2 <- calculate_scaled_data("H3K9me2_input.csv", "H3_input.csv", "H3_LENGTHS.csv", 0.94, 5.8)

pdf("H3K4me2_profile_scaled.pdf", width=8, height=5)
#       line      box      line   scale   no axes     colors      title      y labels   scale
matplot(colMeans(profile_scaled, na.rm=TRUE), type="l", bty="n", lwd=1, log="", axes=FALSE, col=c(2,4), main=NULL, ylab=NULL)
#p. tick length (p: 1-bottom, 2-left, 4-right)
axis(2, tck=0.015, ylab=NULL, lwd=0.5)
#p.  locations of all ticks
axis(1, c(0, 300, 800, 1100), tck=-0.01, lwd=0.5)
#abline(h=0, col="grey")
dev.off()

# #Plot the heatmap in the PNG format
# reduced <- profile_scaled[seq(1, nrow(profile_scaled), 2),]						#keeps every N-th sequence	
# png("H3K4me2_scaled_heatmap.png", width=ncol(reduced), height=nrow(reduced), units="px")	
# #par(mar=c(0,0,0,0))
# breaks <- matrix(NA, nrow=nrow(reduced), ncol=2)
# breaks[,1] <- 300																				#create matrix with column breaks
# breaks[,2] <- 800		
# my_palette <- colorRampPalette(c("#fee8c8", "#fdbb84", "#e34a33"))(n=299)				#define color palette - shades of red
# col_breaks = c(seq(quantile(reduced, 0.10, na.rm=TRUE),quantile(reduced, 0.90, na.rm=TRUE),length=300))	#define color breaks using percentiles
# heatmap.2(reduced, col=my_palette, trace="none", breaks=col_breaks,dendrogram="none", Colv="NA", Rowv="NA", 
#           key=FALSE, labRow=FALSE, labCol=FALSE,  colsep=breaks, sepcolor="grey", lwid=c(0,1), lhei=c(0,1), margins=c(0,0))
# dev.off()


#Plot the profile in the PDF format H3K9AC
pdf("H3K9me2_profile_scaled.pdf", width=8, height=5)
#       line      box      line   scale   no axes     colors      title      y labels   scale
matplot(colMeans(profile_scaled2, na.rm=TRUE), type="l", bty="n", lwd=1, log="", axes=FALSE, col=c(2,4), main=NULL, ylab=NULL)
#p. tick length (p: 1-bottom, 2-left, 4-right)
axis(2, tck=0.015, ylab=NULL, lwd=0.5)
#p.  locations of all ticks
axis(1, c(0, 300, 800, 1100), tck=-0.01, lwd=0.5)
#abline(h=0, col="grey")
dev.off()

# #Plot the heatmap in the PNG format
# reduced <- profile_scaled[seq(1, nrow(profile_scaled), 2),]						#keeps every N-th sequence	
# png("H3K9me2_scaled_heatmap.png", width=ncol(reduced), height=nrow(reduced), units="px")	
# #par(mar=c(0,0,0,0))
# breaks <- matrix(NA, nrow=nrow(reduced), ncol=2)
# breaks[,1] <- 300																				#create matrix with column breaks
# breaks[,2] <- 800		
# my_palette <- colorRampPalette(c("#fee8c8", "#fdbb84", "#e34a33"))(n=299)				#define color palette - shades of red
# col_breaks = c(seq(quantile(reduced, 0.10, na.rm=TRUE),quantile(reduced, 0.90, na.rm=TRUE),length=300))	#define color breaks using percentiles
# heatmap.2(reduced, col=my_palette, trace="none", breaks=col_breaks,dendrogram="none", Colv="NA", Rowv="NA", 
#           key=FALSE, labRow=FALSE, labCol=FALSE,  colsep=breaks, sepcolor="grey", lwid=c(0,1), lhei=c(0,1), margins=c(0,0))
# dev.off()

pdf("ALL_profile_scaled.pdf", width=8, height=5)
matplot(cbind(colMeans(profile_scaled, na.rm=TRUE),colMeans(profile_scaled2, na.rm=TRUE)), type="l", bty="n", lwd=1, log="", axes=FALSE, col=c(2,4), main=NULL, ylab=NULL)
axis(2, tck=0.015, ylab=NULL, lwd=0.5)
#p.  locations of all ticks
axis(1, c(0, 300, 800, 1100), tck=-0.01, lwd=0.5)
#abline(h=0, col="grey")
dev.off()
