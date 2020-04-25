#!/bin/sh
awk '{OFS="\t"}{if($6=="+") print $1, $2-300, $2+8995, $4, $5, $6, $7; else print $1, $3-8995, $3+300, $4, $5, $6, $7}' ../../SHORT_TRANSCRIPTS_12_1_15.bed > tss300.bed
sort -k7 -k1 -k2 -n tss300.bed > sort_tss300.bed

coverageBed -d -a sort_tss300.bed -b CHH_col0_rep2_GSM980986.bed > temp1.bed
awk '{OFS="\t"}{print $1,$2+$8-1,$2+$8,$4,$5,$6,$7,$8,$9}' temp1.bed > temp2.bed
intersectBed -wao -a temp2.bed -b CHH_col0_rep2_GSM980986.bed > temp3.bed
awk '{OFS="\t"}{if($14=="-1") print $1,$2,$3,$4,$5,$6,$7,$8,"NA"; else {print $1,$2,$3,$4,$5,$6,$7,$8,$14}}' temp3.bed > temp4.bed
awk '{OFS="\t"}{if($8==1) s+=1; print $0, s}' temp4.bed > WT_dataset.bed

intersectBed -wao -a temp2.bed -b CHH_nrpe1_GSM981040.bed > temp5.bed
awk '{OFS="\t"}{if($14=="-1") print $1,$2,$3,$4,$5,$6,$7,$8,"NA"; else {print $1,$2,$3,$4,$5,$6,$7,$8,$14}}' temp5.bed > temp6.bed
awk '{OFS="\t"}{if($8==1) s+=1; print $0, s}' temp6.bed > NRPE_dataset.bed

perl ../convert_2d_2.pl WT_dataset.bed 8 0 > input_WT.csv
perl ../convert_2d_2.pl NRPE_dataset.bed 8 0 > input_NRPE.csv
