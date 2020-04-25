setwd('~/cifs-lab/RIP_manuscript/Revised_Figures/Fig2/Fig2F_new/')

rand <- read.table("norm_perm_ovl_PolV-ChIP_met.bed", header=F, sep="\t")

tab1 <- read.table("col_met1_norm_PVT_overlap.bed", header=F, sep="\t") #All PV-transcripts
#tab1 <- read.table("PVT_100reads_2Xenrich.bed", header=F, sep="\t") # only PVTs with Col0>4reads and enrichment>2

pdf("boxplot_PolV-ChIP_met_ovl.pdf")
#pdf("boxplot_PolV-ChIP_100reads_2Xenrich_met_ovl.pdf")
boxplot(tab1$V10,tab1$V11,rand$V4,rand$V5, outline=F,
        ylab="Pol V ChIP [RPM]",
        xlab=c("col0, met1-3, random1000_col0, random1000_met1-3"))

#abline(h=5.1, lty=2, lwd=2, col="darkgrey")
dev.off()

library(MASS)
sink("STATS.txt")
print("stats_Fig2F_new")
wilcox.test(tab1$V10,tab1$V11)
