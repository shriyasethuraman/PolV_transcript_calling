i= seq(1,30000001,500000)
j= seq(500000,30500001,500000)

setwd("~/cifs-lab/RIP_manuscript/Revised_Figures/Fig2/Fig2A/")
window= cbind(1,i,j)
dir.create("~/cifs-lab/RIP_manuscript/Revised_Figures/Fig2/Fig2A/temp_docs/", showWarnings = TRUE, recursive = FALSE)
write.table(window,"~/cifs-lab/RIP_manuscript/Revised_Figures/Fig2/Fig2A/temp_docs/Fig2a.bed", sep="\t")

## Run the command Fig2a_plot_bash in Linux
system("sh -e Fig2a_plot_bash")

counts=read.table("~/cifs-lab/RIP_manuscript/Revised_Figures/Fig2/Fig2A/temp_docs/Fig2a_counts.txt", header=FALSE)
final_counts=counts[,c(4,7)]
colnames(final_counts)=c("500kb windows on Chr1","Number of PolV Transcripts")

pdf("~/cifs-lab/RIP_manuscript/Revised_Figures/Fig2/Fig2A/Fig2a.pdf", width=11, height=8)
plot(final_counts, type="l", col="maroon",lwd=2)
title(main="Distribution of PolV-transcripts(black) on chromosome1")
dev.off()

png("~/cifs-lab/RIP_manuscript/Revised_Figures/Fig2/Fig2A/Fig2a.png", width = 480, height = 360, units = "px")
plot(final_counts, type="l", col="maroon",lwd=2)
title(main="Distribution of PolV-transcripts(black) on chromosome1")
dev.off()

mrna_counts=read.table("~/cifs-lab/RIP_manuscript/Revised_Figures/Fig2/Fig2A/temp_docs/mRNA_counts.txt", header=FALSE)
TE_counts=read.table("~/cifs-lab/RIP_manuscript/Revised_Figures/Fig2/Fig2A/temp_docs/TE_counts.txt", header=FALSE)
count1=mrna_counts[,c(4,7)]
count2=TE_counts[,c(7)]
figb_count=cbind(count1,count2)

pdf("~/cifs-lab/RIP_manuscript/Revised_Figures/Fig2/Fig2A/Figb.pdf", width=11, height=5)
plot(figb_count$V4, figb_count$count2, type="l", col="cyan",lwd=2)
lines(figb_count$V4, figb_count$V7, type="l", col="brown",lwd=2)
title(main="Distribution of PolV-transcripts(black) on chromosome1")
dev.off()
