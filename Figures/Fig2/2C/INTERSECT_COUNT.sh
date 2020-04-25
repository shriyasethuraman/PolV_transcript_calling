#!/bin/sh
sed 1d window.bed > w1
awk '{OFS="\t"}{print $2,$3,$4}' w1 > final_window.bed
intersectBed -c -a final_window.bed -b ~/cifs-lab/RIP_manuscript/Revised_Figures/FINAL_TRANSCRIPTS_12_1_15.bed > PV_count.bed
intersectBed -c -a PV_count.bed -b genes_list.bed > RNA_PV_count.bed
intersectBed -c -a RNA_PV_count.bed -b final_TE.bed > TE_RNA_PV_count.bed

