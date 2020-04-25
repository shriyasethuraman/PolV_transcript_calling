#!/bin/sh
sed 1d Very_Significant_q0.05_transcripts.txt > sign1 #Next three lines alter the R-output to a linux-readable form
cut -f 2- sign1 > sign2
sed 's/"//g' sign2 > Sign.bed
awk '{OFS="\t"}{if(((($11+0.01)/($12+0.01))>4)&&((($13+0.01)/($14+0.01))>4)) print $0}' Sign.bed > double_enrich_sign.bed #Filtering for 4-fold enrichment in both repeats individually
awk '$11>2 && $13>2 && $2>=300' double_enrich_sign.bed > pre_final.bed #Filtering for more than 2 reads in Col0 of both reads and start site greater than 300
awk '$3<=9066 && $6=="-"' pre_final.bed > trim1.bed #Finding the transcripts that have start sites on the -strand lesser than the length of longest transcript
/home/shriyas/Downloads/bedtools2/bin/intersectBed -s -v -a pre_final.bed -b trim1.bed > FINAL.bed #Excluding above set of transcripts
awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$3-$2}' FINAL.bed > trial1.bed
awk -v max=0 '{if(($3-$2)>max){want=($3-$2); max=($3-$2)}}END{print want}' trial1.bed > max_len #maximum length
