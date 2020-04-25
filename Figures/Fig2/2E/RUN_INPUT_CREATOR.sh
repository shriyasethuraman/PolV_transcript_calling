#!/bin/sh
awk '{OFS="\t"}{if($6=="+") print $1, $2-300, $2+8995, $4, $5, $6, $7; else print $1, $3-8995, $3+300, $4, $5, $6, $7}' ../../SHORT_TRANSCRIPTS_12_1_15.bed > tss300.bed
sort -k7 -k1 -k2 -n tss300.bed > sort_tss300.bed
/home/shriyas/Downloads/bedtools2/bin/coverageBed -d -a sort_tss300.bed -b histoneModFiles/H3_GSM701932.bed > H3_count.bed
/home/shriyas/Downloads/bedtools2/bin/coverageBed -d -a sort_tss300.bed -b histoneModFiles/H3K9me2_GSM701926.bed > H3K9me2_count.bed
/home/shriyas/Downloads/bedtools2/bin/coverageBed -d -a sort_tss300.bed -b histoneModFiles/H3K4me2_GSM701923.bed > H3K4me2_count.bed
awk '{OFS="\t"}{if($8==1) s+=1; print $0, s}' H3_count.bed > H3_count_num.bed
awk '{OFS="\t"}{if($8==1) s+=1; print $0, s}' H3K4me2_count.bed > H3K4me2_count_num.bed
awk '{OFS="\t"}{if($8==1) s+=1; print $0, s}' H3K9me2_count.bed > H3K9me2_count_num.bed

perl convert_2d_2.pl H3_count_num.bed 8 0 > H3_input.csv
perl convert_2d_2.pl H3K4me2_count_num.bed 8 0 > H3K4me2_input.csv
perl convert_2d_2.pl H3K9me2_count_num.bed 8 0 > H3K9me2_input.csv
perl convert_2d_2.pl H3_count_num.bed 8 1 > H3_LENGTHS.csv
