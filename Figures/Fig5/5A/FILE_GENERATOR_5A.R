setwd("/home/shriyas/cifs-lab/RIP_manuscript/Revised_Figures/Fig5/Fig5A/")

# source('hp-profile-heatmap_modGB.r')
# figureME(d,300,300)

# full2<- read.table("input_MET.csv", header=FALSE, sep=",", na.strings = "NA")
# length(which(is.na(full)))
# full<- read.table("input_WT.csv", header=FALSE, sep=",", na.strings = "NA")

a <- calculate_values("input_WT.csv","input_NRPE.csv",1,1,1)

d <- DMRs(a, 0.012, 10, 3)	#cutoff, distance, minimum number of bases
sizes_2 <- load_sizes("LENGTHS.csv", a)
values_trim2 <- percent_trim(d, 5)
sizes_trim2 <- percent_trim(sizes_2, 5)
profile_scaled_2 <- scaler(values_trim2, sizes_trim2)
c <- min(colMeans(profile_scaled_2, na.rm=TRUE))
z <- max(colMeans(profile_scaled_2, na.rm=TRUE))

pdf("scaled_profile_final.pdf", width=8, height=5)
#  line      box      line   scale   no axes     colors      title      y labels
x=colMeans(profile_scaled_2, na.rm=TRUE)
matplot((x-min(x))/(max(x)-min(x)), type="l", bty="n", lwd=1, log="", axes=FALSE, col=c(2,4), main=NULL, ylab=NULL)
#p. tick length (p: 1-bottom, 2-left, 4-right)
axis(2, tck=0.015, ylab=NULL, lwd=0.5)
#p.  locations of all ticks
axis(1, c(0, 300, 800, 1100), tck=-0.01, lwd=0.5)
dev.off()

pdf("merged_profile_final.pdf", width=8, height=5)
#  line      box      line   scale   no axes     colors      title      y labels
y=colMeans(profile_scaled_rip, na.rm=TRUE)
matplot(cbind((x-min(x))/(max(x)-min(x)),(y-min(y))/(max(y)-min(y))), type="l", bty="n", lwd=1, log="", axes=FALSE, col=c(2,4), main=NULL, ylab=NULL)
#p. tick length (p: 1-bottom, 2-left, 4-right)
axis(2, tck=0.015, ylab=NULL, lwd=0.5)
#p.  locations of all ticks
axis(1, c(0, 300, 800, 1100), tck=-0.01, lwd=0.5)
dev.off()

pdf("../../FigS5/FigS5A/scaled_profile_final.pdf", width=8, height=5)
#                                         line      box      line   scale   no axes     colors      title      y labels
matplot(colMeans(profile_scaled_2, na.rm=TRUE), type="l", bty="n", lwd=1, log="", axes=FALSE, col=c(2,4), main=NULL, ylab=NULL)
#p. tick length (p: 1-bottom, 2-left, 4-right)
axis(2, tck=0.015, ylab=NULL, lwd=0.5)
#p.  locations of all ticks
axis(1, c(0, 300, 800, 1100), tck=-0.01, lwd=0.5)
dev.off()

reduced <- profile_scaled_2[seq(1, nrow(profile_scaled_2), 5),]	#keeps every N-th sequence	
png("../../FigS5/FigS5A/scaled5_heatmap_final.png", width=ncol(reduced), height=nrow(reduced), units="px")	
#par(mar=c(0,0,0,0))
breaks <- matrix(NA, nrow=nrow(reduced), ncol=2)
breaks[,1] <- 300				#create matrix with column breaks
breaks[,2] <- 800		
my_palette <- colorRampPalette(c("#fee8c8", "#fdbb84", "#e34a33"))(n=299)						#define color palette
col_breaks = c(seq(0,0.43,length=300))      #cutoffs used for CG methylation: in this case same cutoffs as for DMR calling was used
heatmap.2(reduced, col=my_palette, trace="none", breaks=col_breaks,dendrogram="none", Colv="NA", Rowv="NA", 
          key=FALSE, labRow=FALSE, labCol=FALSE,  colsep=breaks, sepcolor="grey", lwid=c(0,1), lhei=c(0,1), margins=c(0,0)) 
dev.off()

sink("stats.txt")
print(significance(values_trim2, 100, 350, 100, 1000))	#Options: data, start of region 1, start of region 2, length of regions, number of permutations
print(significance(values_trim2, 100, 350, 100, 10000))	#Options: data, start of region 1, start of region 2, length of regions, number of permutations
sink()
