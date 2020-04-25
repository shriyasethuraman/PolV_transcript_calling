setwd('~/cifs-lab/RIP_manuscript/Revised_Figures/Fig3/Fig3B/')
rip <- read.table("/home/gudrunb/RIP/resubmission/Fig_3B/norm_SHORT_TRANSCRIPTS_12_1_15_norm-sense-AGO4-CNA", header=F, sep="\t")
pdf("boxplot_sense-AGO4-RIP_on_PolV-transcripts_12-13-15.pdf")
boxplot(rip$V11/rip$V13, rip$V11/rip$V12, rip$V12/rip$V13, outline=F,
        ylab="sense AGO4-RIP [RPM]", main="Pol V-transcripts",
        xlab="Col-0/ago4   Col-0/nrpe1   nrpe1/ago4")
abline(h=1, lty=2, lwd=2, col="darkgrey")
dev.off()

library(MASS)

wilcox.test(rip$V11,rip$V13)
wilcox.test(rip$V11,rip$V12)
wilcox.test(rip$V12,rip$V13)

t.test(rip$V11,rip$V13)
t.test(rip$V11,rip$V12)
t.test(rip$V12,rip$V13)

#wilcox.test(rip$V11/rip$V13,rip$V11/rip$V12)
sink("Stats_Fig3B.txt")
print("Col/Ago vs Nrpe/Ago")
print(wilcox.test(rip$V11/rip$V13,rip$V12/rip$V13))
print("Col/Nrpe vs Nrpe/Ago")
print(wilcox.test(rip$V11/rip$V12,rip$V12/rip$V13))
sink()
