#!/bin/sh
awk '{OFS="\t"}{if($6=="+") print $1, $2-300, $2+8995, $4, $5, $6, $7; else print $1, $3-8995, $3+300, $4, $5, $6, $7}' ../../SHORT_TRANSCRIPTS_12_1_15.bed > tss300.bed
sort -k7 -k1 -k2 -n tss300.bed > sort_tss300.bed
coverageBed -d -a sort_tss300.bed -b col0_AGO4.bed > Col_CHIP.bed
coverageBed -d -a sort_tss300.bed -b nrpe1_AGO4.bed > Nrpe_CHIP.bed
coverageBed -d -a sort_tss300.bed -b 23601.bed > Col_RIP.bed
coverageBed -d -a sort_tss300.bed -b 23602.bed > Nrpe_RIP.bed

awk '{OFS="\t"}{if($8==1) s+=1; print $0, s}' Col_CHIP.bed > Col_CHIP_num.bed 
awk '{OFS="\t"}{if($8==1) s+=1; print $0, s}' Nrpe_CHIP.bed > Nrpe_CHIP_num.bed 
awk '{OFS="\t"}{if($8==1) s+=1; print $0, s}' Col_RIP.bed > Col_RIP_num.bed 
awk '{OFS="\t"}{if($8==1) s+=1; print $0, s}' Nrpe_RIP.bed > Nrpe_RIP_num.bed 

mkdir temp
cd temp
awk '$6=="+"' ../Col_CHIP_num.bed > Col_CHIP_plus.bed
awk '$6=="-"' ../Col_CHIP_num.bed > Col_CHIP_minus.bed
sort -k10n -k8rn Col_CHIP_minus.bed > Col_CHIP_minus_rev_sort.bed
awk '{OFS="\t"}{s+=1; if($8!=1) print $1,$2,$3,$4,$5,$6,$7,s,$9,$10; else {print $1,$2,$3,$4,$5,$6,$7,s,$9,$10; s=0}}' Col_CHIP_minus_rev_sort.bed > Col_CHIP_minus_rev_order.bed
cat Col_CHIP_plus.bed Col_CHIP_minus_rev_order.bed > Col_CHIP_prep_nucleotide_coverage_perl1.bed
sort -k10 -k8 -n Col_CHIP_prep_nucleotide_coverage_perl1.bed > Col_CHIP_perl1.bed
perl ../../convert_2d_2.pl Col_CHIP_perl1.bed 8 0 > Col_CHIP_input.csv
perl ../../convert_2d_2.pl Col_CHIP_perl1.bed 8 1 > LENGTHS.csv

awk '$6=="+"' ../Nrpe_CHIP_num.bed > Nrpe_CHIP_plus.bed
awk '$6=="-"' ../Nrpe_CHIP_num.bed > Nrpe_CHIP_minus.bed
sort -k10n -k8rn Nrpe_CHIP_minus.bed > Nrpe_CHIP_minus_rev_sort.bed
awk '{OFS="\t"}{s+=1; if($8!=1) print $1,$2,$3,$4,$5,$6,$7,s,$9,$10; else {print $1,$2,$3,$4,$5,$6,$7,s,$9,$10; s=0}}' Nrpe_CHIP_minus_rev_sort.bed > Nrpe_CHIP_minus_rev_order.bed
cat Nrpe_CHIP_plus.bed Nrpe_CHIP_minus_rev_order.bed > Nrpe_CHIP_prep_nucleotide_coverage_perl1.bed
sort -k10 -k8 -n Nrpe_CHIP_prep_nucleotide_coverage_perl1.bed > Nrpe_CHIP_perl1.bed
perl ../../convert_2d_2.pl Nrpe_CHIP_perl1.bed 8 0 > Nrpe_CHIP_input.csv

awk '$6=="+"' ../Col_RIP_num.bed > Col_RIP_plus.bed
awk '$6=="-"' ../Col_RIP_num.bed > Col_RIP_minus.bed
sort -k10n -k8rn Col_RIP_minus.bed > Col_RIP_minus_rev_sort.bed
awk '{OFS="\t"}{s+=1; if($8!=1) print $1,$2,$3,$4,$5,$6,$7,s,$9,$10; else {print $1,$2,$3,$4,$5,$6,$7,s,$9,$10; s=0}}' Col_RIP_minus_rev_sort.bed > Col_RIP_minus_rev_order.bed
cat Col_RIP_plus.bed Col_RIP_minus_rev_order.bed > Col_RIP_prep_nucleotide_coverage_perl1.bed
sort -k10 -k8 -n Col_RIP_prep_nucleotide_coverage_perl1.bed > Col_RIP_perl1.bed
perl ../../convert_2d_2.pl Col_RIP_perl1.bed 8 0 > Col_RIP_input.csv
perl ../../convert_2d_2.pl Col_RIP_perl1.bed 8 1 > LENGTHS.csv

awk '$6=="+"' ../Nrpe_RIP_num.bed > Nrpe_RIP_plus.bed
awk '$6=="-"' ../Nrpe_RIP_num.bed > Nrpe_RIP_minus.bed
sort -k10n -k8rn Nrpe_RIP_minus.bed > Nrpe_RIP_minus_rev_sort.bed
awk '{OFS="\t"}{s+=1; if($8!=1) print $1,$2,$3,$4,$5,$6,$7,s,$9,$10; else {print $1,$2,$3,$4,$5,$6,$7,s,$9,$10; s=0}}' Nrpe_RIP_minus_rev_sort.bed > Nrpe_RIP_minus_rev_order.bed
cat Nrpe_RIP_plus.bed Nrpe_RIP_minus_rev_order.bed > Nrpe_RIP_prep_nucleotide_coverage_perl1.bed
sort -k10 -k8 -n Nrpe_RIP_prep_nucleotide_coverage_perl1.bed > Nrpe_RIP_perl1.bed
perl ../../convert_2d_2.pl Nrpe_RIP_perl1.bed 8 0 > Nrpe_RIP_input.csv