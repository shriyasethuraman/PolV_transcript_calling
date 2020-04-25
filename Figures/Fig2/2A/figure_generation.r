#load functions
source('heatmap_tools_v07.r')

#Calculate scaled profile from CSV files
#If need to plot more than one profile, combine them by cbind
profile_scaled <- calculate_scaled_data("RIP1_PolV_C_same.csv", "RIP1_PolV_N_same.csv", "lengths.csv", 20, 20)

#Plot the profile in the PDF format
pdf("profile_scaled.pdf", width=8, height=5)
#                                             line      box      line   scale   no axes     colors      title      y labels   scale
matplot(colMeans(profile_scaled, na.rm=TRUE), type="l", bty="n", lwd=1, log="", axes=FALSE, col=c(2,4), main=NULL, ylab=NULL)
#p. tick length (p: 1-bottom, 2-left, 4-right)
axis(2, tck=0.015, ylab=NULL, lwd=0.5)
#p.  locations of all ticks
axis(1, c(0, 300, 800, 1100), tck=-0.01, lwd=0.5)
#abline(h=0, col="grey")
dev.off()

#Plot the heatmap in the PNG format
reduced <- profile_scaled[seq(1, nrow(profile_scaled), 2),]						#keeps every N-th sequence	
png("scaled_heatmap.png", width=ncol(reduced), height=nrow(reduced), units="px")	
#par(mar=c(0,0,0,0))
breaks <- matrix(NA, nrow=nrow(reduced), ncol=2)
breaks[,1] <- 300																				#create matrix with column breaks
breaks[,2] <- 800		
my_palette <- colorRampPalette(c("#fee8c8", "#fdbb84", "#e34a33"))(n=299)				#define color palette - shades of red
col_breaks = c(seq(quantile(reduced, 0.10, na.rm=TRUE),quantile(reduced, 0.90, na.rm=TRUE),length=300))	#define color breaks using percentiles
heatmap.2(reduced, col=my_palette, trace="none", breaks=col_breaks,dendrogram="none", Colv="NA", Rowv="NA", 
          key=FALSE, labRow=FALSE, labCol=FALSE,  colsep=breaks, sepcolor="grey", lwid=c(0,1), lhei=c(0,1), margins=c(0,0))
dev.off()

#Calculate significance of differences using unscaled profile
profile_unscaled <- calculate_unscaled_data("RIP1_PolV_C_same.csv", "RIP1_PolV_N_same.csv", "lengths.csv", 20, 20)
significance(profile_unscaled, 100, 350, 100, 1000)		#Options: data, start of region 1, start of region 2, length of regions, number of permutations