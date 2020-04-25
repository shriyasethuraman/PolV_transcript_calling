#!/bin/sh
awk '{OFS="\t"}{if($6=="+") print $1, $2-300, $2+8995, $4, $5, $6, $7; else print $1, $3-8995, $3+300, $4, $5, $6, $7}' ../../SHORT_TRANSCRIPTS_12_1_15.bed > tss300.bed
sort -k7 -k1 -k2 -n tss300.bed > sort_tss300.bed
coverageBed -d -a sort_tss300.bed -b ../24siRNA_GSM893118_wtedt.bed > Col_sense.bed
coverageBed -d -a sort_tss300.bed -b ../24siRNA_GSM893123_nrpd1.bed > Nrpd_sense.bed
coverageBed -d -a sort_tss300.bed -b ../24siRNA_GSM893115_nrpe_1.bed > Nrpe_sense.bed
coverageBed -d -a sort_tss300.bed -b ../24siRNA_GSM893112_Col_1.bed > Col_nrpe_control_sense.bed

awk '{OFS="\t"}{if($8==1) s+=1; print $0, s}' Col_sense.bed > Col_sense_num.bed 
awk '{OFS="\t"}{if($8==1) s+=1; print $0, s}' Nrpd_sense.bed > Nrpd_sense_num.bed 
awk '{OFS="\t"}{if($8==1) s+=1; print $0, s}' Nrpe_sense.bed > Nrpe_sense_num.bed 
awk '{OFS="\t"}{if($8==1) s+=1; print $0, s}' Col_nrpe_control_sense.bed > Col_nrpe_control_sense_num.bed

mkdir temp
cd temp
awk '$6=="+"' ../Col_sense_num.bed > Col_sense_plus.bed
awk '$6=="-"' ../Col_sense_num.bed > Col_sense_minus.bed
sort -k10n -k8rn Col_sense_minus.bed > Col_sense_minus_rev_sort.bed
awk '{OFS="\t"}{s+=1; if($8!=1) print $1,$2,$3,$4,$5,$6,$7,s,$9,$10; else {print $1,$2,$3,$4,$5,$6,$7,s,$9,$10; s=0}}' Col_sense_minus_rev_sort.bed > Col_sense_minus_rev_order.bed
cat Col_sense_plus.bed Col_sense_minus_rev_order.bed > Col_sense_prep_nucleotide_coverage_perl1.bed
sort -k10 -k8 -n Col_sense_prep_nucleotide_coverage_perl1.bed > Col_sense_perl1.bed
perl ../../convert_2d_2.pl Col_sense_perl1.bed 8 0 > Col_sense_input.csv
perl ../../convert_2d_2.pl Col_sense_perl1.bed 8 1 > LENGTHS.csv


awk '$6=="+"' ../Nrpd_sense_num.bed > Nrpd_sense_plus.bed
awk '$6=="-"' ../Nrpd_sense_num.bed > Nrpd_sense_minus.bed
sort -k10n -k8rn Nrpd_sense_minus.bed > Nrpd_sense_minus_rev_sort.bed
awk '{OFS="\t"}{s+=1; if($8!=1) print $1,$2,$3,$4,$5,$6,$7,s,$9,$10; else {print $1,$2,$3,$4,$5,$6,$7,s,$9,$10; s=0}}' Nrpd_sense_minus_rev_sort.bed > Nrpd_sense_minus_rev_order.bed
cat Nrpd_sense_plus.bed Nrpd_sense_minus_rev_order.bed > Nrpd_sense_prep_nucleotide_coverage_perl1.bed
sort -k10 -k8 -n Nrpd_sense_prep_nucleotide_coverage_perl1.bed > Nrpd_sense_perl1.bed
perl ../../convert_2d_2.pl Nrpd_sense_perl1.bed 8 0 > Nrpd_sense_input.csv

awk '$6=="+"' ../Nrpe_sense_num.bed > Nrpe_sense_plus.bed
awk '$6=="-"' ../Nrpe_sense_num.bed > Nrpe_sense_minus.bed
sort -k10n -k8rn Nrpe_sense_minus.bed > Nrpe_sense_minus_rev_sort.bed
awk '{OFS="\t"}{s+=1; if($8!=1) print $1,$2,$3,$4,$5,$6,$7,s,$9,$10; else {print $1,$2,$3,$4,$5,$6,$7,s,$9,$10; s=0}}' Nrpe_sense_minus_rev_sort.bed > Nrpe_sense_minus_rev_order.bed
cat Nrpe_sense_plus.bed Nrpe_sense_minus_rev_order.bed > Nrpe_sense_prep_nucleotide_coverage_perl1.bed
sort -k10 -k8 -n Nrpe_sense_prep_nucleotide_coverage_perl1.bed > Nrpe_sense_perl1.bed
perl ../../convert_2d_2.pl Nrpe_sense_perl1.bed 8 0 > Nrpe_sense_input.csv


awk '$6=="+"' ../Col_nrpe_control_sense_num.bed > Col_nrpe_control_plus.bed
awk '$6=="-"' ../Col_nrpe_control_sense_num.bed > Col_nrpe_control_minus.bed
sort -k10n -k8rn Col_nrpe_control_minus.bed > Col_nrpe_control_minus_rev_sort.bed
awk '{OFS="\t"}{s+=1; if($8!=1) print $1,$2,$3,$4,$5,$6,$7,s,$9,$10; else {print $1,$2,$3,$4,$5,$6,$7,s,$9,$10; s=0}}' Col_nrpe_control_minus_rev_sort.bed > Col_nrpe_control_minus_rev_order.bed
cat Col_nrpe_control_plus.bed Col_nrpe_control_minus_rev_order.bed > Col_nrpe_control_prep_nucleotide_coverage_perl1.bed
sort -k10 -k8 -n Col_nrpe_control_prep_nucleotide_coverage_perl1.bed > Col_nrpe_control_perl1.bed
perl ../../convert_2d_2.pl Col_nrpe_control_perl1.bed 8 0 > Col_nrpe_control_input.csv

