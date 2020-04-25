#!/bin/sh
#nucleotide-wise read count calc along the sense strands of the PolV transcripts
awk '{OFS="\t"}{if($6=="+") print $1, $2-300, $2+8995, $4, $5, $6, $7; else print $1, $3-8995, $3+300, $4, $5, $6, $7}' trial1.bed > tss300.bed #Extending the length of transcripts to that of the max. length = 8766+600
sort -k7 -k1 -k2 -n tss300.bed > sort_tss300.bed #sort them according to 1: Lt of transcript; 2:chr#; 3:start site

/home/shriyas/Downloads/bedtools2/bin/coverageBed -d -s -a sort_tss300.bed -b ../AllCol.bed > nucleotide_coverage_col0.bed #merged Col0 coverage
/home/shriyas/Downloads/bedtools2/bin/coverageBed -d -s -a sort_tss300.bed -b ../AllNrpe.bed > nucleotide_coverage_nrpe1.bed #merged Nrpe1 coverage
awk '{OFS="\t"}{if($8==1) s+=1; print $0, s}' nucleotide_coverage_col0.bed > nucleotide_coverage_col0_num.bed #Assign the transcript to which each nucleotide belogs
awk '{OFS="\t"}{if($8==1) s+=1; print $0, s}' nucleotide_coverage_nrpe1.bed > nucleotide_coverage_nrpe1_num.bed

mkdir sense_temp_files
cd sense_temp_files
awk '$6=="+"' ../nucleotide_coverage_col0_num.bed > plus.bed #Extract transcripts on the +strand
awk '$6=="-"' ../nucleotide_coverage_col0_num.bed > minus.bed #Extract transcripts on the -strand
sort -k10n -k8rn minus.bed > minus_rev_sort.bed #Reverse sort transcript based on nucleotide position
awk '{OFS="\t"}{s+=1; if($8!=1) print $1,$2,$3,$4,$5,$6,$7,s,$9,$10; else {print $1,$2,$3,$4,$5,$6,$7,s,$9,$10; s=0}}' minus_rev_sort.bed > minus_rev_order.bed #Fix the reversed order as the actual order
cat plus.bed minus_rev_order.bed > prep_nucleotide_coverage_perl1.bed #Concatenate the plus and corrected-minus transcripts again
sort -k10 -k8 -n prep_nucleotide_coverage_perl1.bed > nucleotide_coverage_perl1.bed #Sort back based on lt of Transcripts
perl ../../convert_2d_2.pl nucleotide_coverage_perl1.bed 8 0 > input1.csv #Obtain the right format of the files for heatmap plotting
perl ../../convert_2d_2.pl nucleotide_coverage_perl1.bed 8 1 > LENGTHS.csv

awk '$6=="+"' ../nucleotide_coverage_nrpe1_num.bed > plus_n.bed
awk '$6=="-"' ../nucleotide_coverage_nrpe1_num.bed > minus_n.bed
sort -k10n -k8rn minus_n.bed > minus_rev_sort_n.bed
awk '{OFS="\t"}{s+=1; if($8!=1) print $1,$2,$3,$4,$5,$6,$7,s,$9,$10; else {print $1,$2,$3,$4,$5,$6,$7,s,$9,$10; s=0}}' minus_rev_sort_n.bed > minus_rev_order_n.bed 
cat plus_n.bed minus_rev_order_n.bed > prep_nucleotide_coverage_perl2.bed
sort -k10 -k8 -n prep_nucleotide_coverage_perl2.bed > nucleotide_coverage_perl2.bed
perl ../../convert_2d_2.pl nucleotide_coverage_perl2.bed 8 0 > input2.csv

