#!/bin/sh

intersectBed -c -a ~/cifs-lab/RIP_manuscript/Revised_Figures/SHORT_TRANSCRIPTS_12_1_15.bed -b ~/cifs-lab/RIP_manuscript/Revised_Figures/Fig2/Fig2F_new/col0/sort_aligned.bed | intersectBed -c -a ~/cifs-lab/RIP_manuscript/Revised_Figures/SHORT_TRANSCRIPTS_12_1_15.bed -b ~/cifs-lab/RIP_manuscript/Revised_Figures/Fig2/Fig2F_new/met1-3/sort_aligned.bed | awk '{OFS="\t"}{print $0,$8/20.54+(0.5/20.54), $9/28.01 + (0.5/20.54)}' > col_met1_norm_PVT_overlap.bed

#awk '$8>100' col_met1_norm_PVT_overlap.bed | awk '$10/$11 > 2' > PVT_100reads_2Xenrich.bed
#awk '$8>4' col_met1_norm_PVT_overlap.bed | awk '$10/$11 > 2' > PVT_4reads_2Xenrich.bed


## PERMUTATIONS

i=1
while [ $i -le 1000 ]
do
shuffleBed -i ~/cifs-lab/RIP_manuscript/Revised_Figures/SHORT_TRANSCRIPTS_12_1_15.bed -g ../chrom_sizes.txt | sortBed -i stdin > temp1
intersectBed -c -a temp1 -b ~/cifs-lab/RIP_manuscript/Revised_Figures/Fig2/Fig2F_new/col0/sort_aligned.bed > temp3
intersectBed -c -a temp3 -b ~/cifs-lab/RIP_manuscript/Revised_Figures/Fig2/Fig2F_new/met1-3/sort_aligned.bed | awk '{print "1\t"$0}' | groupBy -i stdin -g 1 -c 9,10 -o median,median >> perm_ovl_PolV-ChIP_met.bed


((i++))
done

awk '{OFS="\t"}{print $0,$2/20.54+(0.5/20.54), $3/28.01 + (0.5/20.54)}' perm_ovl_PolV-ChIP_met.bed > norm_perm_ovl_PolV-ChIP_met.bed 


