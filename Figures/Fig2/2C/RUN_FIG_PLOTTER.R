setwd("/home/shriyas/cifs-lab/RIP_manuscript/Revised_Figures/Fig2/Fig2C/")
 
# bins1=seq(1,61)
# bins2=seq(1,40)
# bins3=seq(1,47)
# bins4=seq(1,38)
# bins5=seq(1,54)

i= seq(1,30000001,500000)
j= seq(500000,30500001,500000)
window1= cbind(1,i,j)

i= seq(1,19500001,500000)
j= seq(500000,20000001,500000)
window2= cbind(2,i,j)

i= seq(1,23000001,500000)
j= seq(500000,23500001,500000)
window3= cbind(3,i,j)

i= seq(1,18500001,500000)
j= seq(500000,19000001,500000)
window4= cbind(4,i,j)

i= seq(1,26500001,500000)
j= seq(500000,27000001,500000)
window5= cbind(5,i,j)
window= rbind(window1,window2,window3,window4,window5)
write.table(window,"/home/shriyas/cifs-lab/RIP_manuscript/Revised_Figures/Fig2/Fig2C/window.bed", sep="\t")

system("sh -e INTERSECT_COUNT_BASH")

final <- read.table("/home/shriyas/cifs-lab/RIP_manuscript/Revised_Figures/Fig2/Fig2C/TE_RNA_PV_count.bed", header=FALSE)
colnames(final) = c("col1","col2","col3","col4","col5","col6")
final2=as.data.frame(final)

pdf("Fig2C_gg.pdf", width=7, height=6)
#plot((final$V5/final$V6),(final$V4/final$V6))
ggplot((data=final2), aes(x = col5/col6, y = col4/col6)) +
  geom_point() + geom_smooth(method=lm, se=FALSE) + theme_classic()#+ geom_abline(lm(((Rep1_Col - Rep1_Nrpe)) ~ ((Rep2_Col - Rep2_Nrpe))))#+ geom_abline(intercept = log2(18.0987), slope = log2(2.37))
dev.off()

sink("R-test.txt")
print(cor.test((final2$col5/final2$col6),(final2$col4/final2$col6)))
sink()
