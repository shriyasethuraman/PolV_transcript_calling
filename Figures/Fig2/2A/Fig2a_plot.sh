#!/bin/sh
cd ~/cifs-lab/RIP_manuscript/Revised_Figures/Fig2/Fig2A/temp_docs/
sed 1d Fig2a.bed > fig_alt1
cut -f 2- fig_alt1 > fig_alt2
awk '{OFS="\t"}{s+=1; print $0,s, "255", "+"}' fig_alt2 > Fig2a_alt.bed
/home/shriyas/Downloads/bedtools2/bin/intersectBed -c -a Fig2a_alt.bed -b ../../../FINAL_TRANSCRIPTS_12_1_15.bed > Fig2a_counts.txt
/home/shriyas/Downloads/bedtools2/bin/intersectBed -c -a Fig2a_alt.bed -b ../mRNA.bed > mRNA_counts.txt
/home/shriyas/Downloads/bedtools2/bin/intersectBed -c -a Fig2a_alt.bed -b ../Trans.bed > TE_counts.txt

