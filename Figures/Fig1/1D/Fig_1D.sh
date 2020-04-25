#!/bin/sh
cd temp_files
awk '{OFS="\t"}{print $1,$2-600,$3+600,$5,"255",$6}' ../summits_AW_PolVpeaks > peak_regions.bed
awk '{OFS="\t"}{print $0, $3-$2}' peak_regions.bed > peaks.bed

sed -e '/^[a-z]/d' ../nrpe1_POLV.bed > Chip_nrpe1.bed
sed -e '/^[a-z]/d' ../col0_POLV.bed > Chip_col0.bed

coverageBed -d -a peaks.bed -b Chip_col0.bed > CHIP_COL.bed
coverageBed -d -a peaks.bed -b Chip_nrpe1.bed > CHIP_NRPE.bed
coverageBed -d -a peaks.bed -b ../24993.bed > RIP_COL.bed
coverageBed -d -a peaks.bed -b ../24994.bed > RIP_NRPE.bed

awk '{OFS="\t"}{if($8==1) s+=1; print $0, s}' CHIP_COL.bed > CHIP_1.bed
awk '{OFS="\t"}{if($8==1) s+=1; print $0, s}' CHIP_NRPE.bed > CHIP_2.bed
awk '{OFS="\t"}{if($8==1) s+=1; print $0, s}' RIP_NRPE.bed > RIP_2.bed
awk '{OFS="\t"}{if($8==1) s+=1; print $0, s}' RIP_COL.bed > RIP_1.bed

perl ../convert_2d_2.pl CHIP_1.bed 8 0 > CHIP_input1.csv
perl ../convert_2d_2.pl CHIP_2.bed 8 0 > CHIP_input2.csv
perl ../convert_2d_2.pl RIP_1.bed 8 0 > RIP_input1.csv
perl ../convert_2d_2.pl RIP_2.bed 8 0 > RIP_input2.csv
perl ../convert_2d_2.pl CHIP_1.bed 8 1 > CHIP_LENGTHS.csv
perl ../convert_2d_2.pl RIP_1.bed 8 1 > RIP_LENGTHS.csv

