#!/bin/sh

/home/shriyas/Downloads/bedtools2/bin/intersectBed -s -v -a 23601.bed -b genes_list.bed > Ago_col_nongenic.bed
sort -k1,1 -k2,2n -k3,3n -k6,6 -u Ago_col_nongenic.bed > x_ago.bed
/home/shriyas/Downloads/bedtools2/bin/intersectBed -s -v -a 23603.bed -b genes_list.bed > Ago_ago_nongenic.bed
sort -k1,1 -k2,2n -k3,3n -k6,6 -u Ago_ago_nongenic.bed > y_ago.bed

mkdir Ago_Transcripts
cd Ago_Transcripts

/home/shriyas/Downloads/bedtools2/bin/mergeBed -s -d 200 -i ../x_ago.bed > merged_col.txt #creating the merged transcripts with distance <=200
awk '{OFS="\t"}{s+=1; print $1,$2,$3,"trans"s,"255",$4}' merged_col.txt > merged_col.bed
/home/shriyas/Downloads/bedtools2/bin/intersectBed -s -c -F 1 -a merged_col.bed -b ../x_ago.bed > reads_col0.bed #Finding the deduplicated merged non-genic Col0 reads overlapping the transcripts
awk '$7>4 {print $0}' reads_col0.bed > greater_than_6.bed #Filtering for more than 6 Col0 reads
awk '{if((($3-$2)/100)<=$7) print $0}' greater_than_6.bed > v2_greater_than_1bp.bed #Filtering for transcripts with atleast 1read/100bp of transcript
/home/shriyas/Downloads/bedtools2/bin/intersectBed -s -F 1 -c -a v2_greater_than_1bp.bed -b ../y_ago.bed> col_ago4_read_counts_v2.bed #Finding the deduplicated merged non-genic ago4 reads overlapping the transcripts
awk '{OFS="\t"}{if((($7+0.01)/($8+0.01))>6) print $0}' col_ago4_read_counts_v2.bed > final_transcripts_v2.bed #Filtering for 4-fold enrichment in Col0 over ago4
sort -k1 -k2 -n final_transcripts_v2.bed > sort_transcripts_v2.bed
/home/shriyas/Downloads/bedtools2/bin/intersectBed -s -F 1 -c -a sort_transcripts_v2.bed -b ../23601.bed > col_transcripts.bed #Finding the number of duplicated merged Col0 reads overlapping the transcripts
/home/shriyas/Downloads/bedtools2/bin/intersectBed -s -F 1 -c -a col_transcripts.bed -b ../23603.bed > col_ago4_transcripts.bed #Finding the number of duplicated merged ago4 reads overlapping the transcripts
awk '((($9+0.01)/($10+0.01))>6) {print $0}' col_ago4_transcripts.bed > col_ago4_transcripts_final.bed #Filtering for 4-fold enrichment in Col0 over ago4 for all reads, incl. duplicated reads
/home/shriyas/Downloads/bedtools2/bin/intersectBed -s -v -a col_ago4_transcripts_final.bed -b ../genes_list.bed > NonGenic_PV_transcripts.bed #Filtering out the genic transcripts

awk '$2>=300' NonGenic_PV_transcripts.bed > pre_final.bed #Filtering for more than 2 reads in Col0 of both reads and start site greater than 300
awk '$3<=11175 && $6=="-"' pre_final.bed > trim1.bed #Finding the transcripts that have start sites on the -strand lesser than the length of longest transcript
/home/shriyas/Downloads/bedtools2/bin/intersectBed -s -v -a pre_final.bed -b trim1.bed > FINAL.bed #Excluding above set of transcripts
awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$3-$2}' FINAL.bed > trial1.bed
awk -v max=0 '{if(($3-$2)>max){want=($3-$2); max=($3-$2)}}END{print want}' trial1.bed > max_len #maximum length
