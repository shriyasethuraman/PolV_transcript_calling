setwd("/home/shriyas/cifs-lab/RIP_manuscript/Revised_Figures/Fig4/Fig4A/")

system("RUN_Fig4A_BASH.txt") #needs correction
reads <- read.table("temp4.bed", header=FALSE)

result1 <- (((reads$V8 / 39.33) + (0.5 / 39.33)) / ((reads$V9 / 27.61) + (0.5 / 39.33)))
result2 <- (((reads$V10 / 33) + (0.5 / 33)) / ((reads$V9 / 27.61) + (0.5 / 33)))
result3 <- (((reads$V11 / 41.44) + (0.5 / 41.44)) / ((reads$V9 / 27.61) + (0.5 / 41.44)))
final_result <- cbind(result1,result2,result3)
pdf("Fig4A.pdf", width=8, height=5)
boxplot(final_result, ylim=c(0,200),outline=FALSE)
dev.off()

sink("Stats_Fig4A.txt")
print("Ago/Nrpe vs Col/Nrpe")
print(wilcox.test(result2,result1))
print("idn2/Nrpe vs Col/Nrpe")
print(wilcox.test(result3,result1))
sink()
