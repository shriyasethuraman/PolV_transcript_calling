setwd("~/cifs-lab/RIP_manuscript/Revised_Figures/Fig4/Fig4C/")
data = read.csv("Read_count_final.bed", sep="\t", header = F, stringsAsFactors = F)
colnames(data) = c("Chr","Start","End","ID","Score","Strand","Length","PolV1_Col0","PolV1_nrpe1","PolV2_Col0","PolV2_nrpe1","PolV2_ago4","PolV2_idn2","AGO4_Col0","AGO4_nrpe1", "AGO4_ago4")

data$PolV1_Col0 <- data$PolV1_Col0 / 34.81
data$PolV1_nrpe1 <- data$PolV1_nrpe1 / 32.53
data$PolV2_Col0 <- data$PolV2_Col0 / 39.32
data$PolV2_nrpe1 <- data$PolV2_nrpe1 / 27.60
data$PolV2_ago4 <- data$PolV2_ago4 / 33.00
data$PolV2_idn2 <- data$PolV2_idn2 / 41.43
data$AGO4_Col0 <- data$AGO4_Col0 / 20.57
data$AGO4_nrpe1 <- data$AGO4_nrpe1 / 17.20
data$AGO4_ago4 <- data$AGO4_ago4 / 15.86

library(ggplot2)

#Generate Fig. 4C ago4 vs. idn2
Fig4 <- c()

Fig4$IC_ratio <- data$PolV2_idn2-data$PolV2_Col0
Fig4$AC_ratio <- data$PolV2_ago4-data$PolV2_Col0

Fig4 <- data.frame(Fig4)
Fig4 <- Fig4[Fig4$IC_ratio > -4,]
Fig4 <- Fig4[Fig4$IC_ratio < 4,]
Fig4 <- Fig4[Fig4$AC_ratio > -4,]
Fig4 <- Fig4[Fig4$AC_ratio < 4,]

ggplot(data=Fig4, aes(x=IC_ratio, y=AC_ratio)) +
  geom_point(size=1.2) +
  geom_smooth(method=lm, se=FALSE) +
  xlab("idn2 - Col0 [RPM]") +
  ylab("ago4 - Col0 [RPM]") +
  scale_color_distiller(palette="Spectral") +
  geom_density2d(aes(color=..level..)) +
  #ylim(-5, 5) +
  #xlim(-5, 5) +
  theme_classic() 

sink("R-test.txt")
library(MASS)
ia.cor <- lm(IC_ratio ~ AC_ratio, data=Fig4)
print(summary(ia.cor))
print(cor.test(Fig4$IC_ratio, Fig4$AC_ratio))
sink()
