setwd("/home/shriyas/cifs-lab/RIP_manuscript/Revised_figures/Fig1/Fig1B/")

system("sh -e CALLING1.sh")
#setwd("/home/shriyas/cifs-lab/Shriya/Transcript_Call/Correct_Transcripts/")

source("http://bioconductor.org/biocLite.R")
biocLite("qvalue")
# install.packages("NBPSeq", dependencies = TRUE)
library(NBPSeq)
full_table <- read.table("head_rep2All.txt", header=TRUE)
#pvtranscripts <- full_table[,c(7,8,9, 10)]
#pvtranscripts <- full_table[,c(9,10,11, 12)]
pvtranscripts <- full_table[,c(11, 12, 13, 14)]
rownames(pvtranscripts) <- full_table$Name
pvt <- data.matrix(pvtranscripts)
lib.sizes <- c(34810333, 32531662, 39327986, 27607999)
names(lib.sizes) <- colnames(pvt)
grp.ids <- c(1, 2, 1, 2)
grp1 <- 1
grp2 <- 2
set.seed(999)
res <- nbp.test(pvt, grp.ids, grp1, grp2, lib.sizes=lib.sizes)
#result <- cbind(pvt, res$p.value)
#result <- cbind(result, res$q.value)
result <- cbind(res$pseudo.counts, res$p.value, res$q.value)

tmp <- as.data.frame(result)

full_table$pvalue <- tmp$V5
full_table$qvalue <- tmp$V6

altered_reads <- cbind(full_table$Rep1_Col,((full_table$Rep1_Nrpe)*(34810333/32531662)),full_table$Rep2_Col,((full_table$Rep2_Nrpe)*(39327986/27607999)))

altered_reads <- cbind(altered_reads,full_table$pvalue)
altered_reads <- cbind(altered_reads,full_table$qvalue)
colnames(altered_reads) <- c("Rep1_Col","Rep1_Nrpe","Rep2_Col","Rep2_Nrpe","pvalue","qvalue")


library(ggplot2)
alt_head <- as.data.frame(altered_reads)
alt_head <- cbind(alt_head, c())
#alt_head[which(alt_head[(which(alt_head$pvalue > 0.05)),6] > 0.05),7] <- 0
#alt_head[which(alt_head[(which(alt_head$pvalue <= 0.05)),6] > 0.05),7] <- 1
#alt_head[which(alt_head[(which(alt_head$pvalue > 0.05)),6] <= 0.05),7] <- 2
#alt_head[which(alt_head[(which(alt_head$pvalue <= 0.05)),6] <= 0.05),7] <- 3

alt_head$V7[alt_head$pvalue > 0.05] <- "p>0.05"
alt_head$V7[alt_head$pvalue <= 0.05] <-"p<=0.05"
alt_head$V7[alt_head$pvalue <= 0.01] <-"p<=0.01"
alt_head$V7[alt_head$pvalue <= 0.001] <-"p<=0.001"

reduced <- profile_scaled[seq(1, nrow(profile_scaled), 2),]						#keeps every N-th sequence	
png("plot_pvalue.png", width=nrow(reduced), height=nrow(reduced), units="px")
ggplot() + 
  geom_point(data=alt_head, aes(x = log2(Rep1_Col - Rep1_Nrpe), y = log2(Rep2_Col - Rep2_Nrpe), color=factor(V7))) #+ geom_abline(lm(((Rep1_Col - Rep1_Nrpe)) ~ ((Rep2_Col - Rep2_Nrpe))))#+ geom_abline(intercept = log2(18.0987), slope = log2(2.37))
# ggplot() + 
#   geom_point(data=alt_head, aes(x = log2((Rep1_Col+0.01)/(Rep1_Nrpe+0.01)), y = log2((Rep2_Col+0.01)/(Rep2_Nrpe+0.01)), color=factor(V7))) #+ geom_abline(lm(((Rep1_Col - Rep1_Nrpe)) ~ ((Rep2_Col - Rep2_Nrpe))))#+ geom_abline(intercept = log2(18.0987), slope = log2(2.37))
dev.off()

alt_head <- cbind(alt_head, c())

alt_head$V8[alt_head$qvalue > 0.05] <- "q>0.05"
alt_head$V8[alt_head$qvalue <= 0.05] <-"q<=0.05"
alt_head$V8[alt_head$qvalue <= 0.01] <-"q<=0.01"
alt_head$V8[alt_head$qvalue <= 0.001] <-"q<=0.001"

png("plot_qvalue.png", width=nrow(reduced), height=nrow(reduced), units="px")

ggplot() + 
  geom_point(data=alt_head, aes(x = log2(Rep1_Col - Rep1_Nrpe), y = log2(Rep2_Col - Rep2_Nrpe), color=factor(V8))) #+ geom_abline(lm(((Rep1_Col - Rep1_Nrpe)) ~ ((Rep2_Col - Rep2_Nrpe))))#+ geom_abline(intercept = log2(18.0987), slope = log2(2.37))
dev.off()

sign <- full_table[which(full_table$pvalue <= 0.05),]
more_sign <- sign[which(sign$qvalue <= 0.05),]

write.table(more_sign, file="Very_Significant_q0.05_transcripts.txt", sep="\t")

write.table(sign, file="Significant_transcripts.txt", sep="\t")

system("sh -e CALLING2.sh")

system("sh -e NUCLEO_SENSE_INPUT.sh")
