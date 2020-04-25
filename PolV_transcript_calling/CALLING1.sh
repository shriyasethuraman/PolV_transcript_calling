#!/bin/sh
cat 24993.bed rep2_col0.bed > AllCol.bed
cat 24994.bed rep2_nrpe1.bed > AllNrpe.bed
/home/shriyas/Downloads/bedtools2/bin/intersectBed -s -v -a AllCol.bed -b genes_list.bed > Col_nongenic.bed
sort -k1,1 -k2,2n -k3,3n -k6,6 -u Col_nongenic.bed > x.bed
/home/shriyas/Downloads/bedtools2/bin/intersectBed -s -v -a AllNrpe.bed -b genes_list.bed > Nrpe_nongenic.bed
sort -k1,1 -k2,2n -k3,3n -k6,6 -u Nrpe_nongenic.bed > y.bed

mkdir Correct_Transcripts
cd Correct_Transcripts

/home/shriyas/Downloads/bedtools2/bin/mergeBed -s -d 200 -i ../x.bed > merged_col.txt #creating the merged transcripts with distance <=200
awk '{OFS="\t"}{s+=1; print $1,$2,$3,"trans"s,"255",$4}' merged_col.txt > merged_col.bed
/home/shriyas/Downloads/bedtools2/bin/intersectBed -s -c -F 1 -a merged_col.bed -b ../x.bed > reads_col0.bed #Finding the deduplicated merged non-genic Col0 reads overlapping the transcripts
awk '$7>8 {print $0}' reads_col0.bed > greater_than_6.bed #Filtering for more than 6 Col0 reads
awk '{if((($3-$2)/100)<=$7) print $0}' greater_than_6.bed > v2_greater_than_1bp.bed #Filtering for transcripts with atleast 1read/100bp of transcript
/home/shriyas/Downloads/bedtools2/bin/intersectBed -s -F 1 -c -a v2_greater_than_1bp.bed -b ../y.bed> col_nrpe1_read_counts_v2.bed #Finding the deduplicated merged non-genic Nrpe1 reads overlapping the transcripts
awk '{OFS="\t"}{if((($7+0.01)/($8+0.01))>4) print $0}' col_nrpe1_read_counts_v2.bed > final_transcripts_v2.bed #Filtering for 4-fold enrichment in Col0 over Nrpe1
sort -k1 -k2 -n final_transcripts_v2.bed > sort_transcripts_v2.bed
/home/shriyas/Downloads/bedtools2/bin/intersectBed -s -F 1 -c -a sort_transcripts_v2.bed -b ../AllCol.bed > col_transcripts.bed #Finding the number of duplicated merged Col0 reads overlapping the transcripts
/home/shriyas/Downloads/bedtools2/bin/intersectBed -s -F 1 -c -a col_transcripts.bed -b ../AllNrpe.bed > col_nrpe_transcripts.bed #Finding the number of duplicated merged Nrpe1 reads overlapping the transcripts
awk '((($9+0.01)/($10+0.01))>4) {print $0}' col_nrpe_transcripts.bed > col_nrpe_transcripts_final.bed #Filtering for 4-fold enrichment in Col0 over Nrpe1 for all reads, incl. duplicated reads
/home/shriyas/Downloads/bedtools2/bin/intersectBed -s -v -a col_nrpe_transcripts_final.bed -b ../genes_list.bed > NonGenic_PV_transcripts.bed #Filtering out the genic transcripts
/home/shriyas/Downloads/bedtools2/bin/intersectBed -s -F 1 -c -a NonGenic_PV_transcripts.bed -b ../24993.bed > rep1Col.bed #Adding the individual read counts
/home/shriyas/Downloads/bedtools2/bin/intersectBed -s -F 1 -c -a rep1Col.bed -b ../24994.bed > rep1All.bed
/home/shriyas/Downloads/bedtools2/bin/intersectBed -s -F 1 -c -a rep1All.bed -b ../rep2_col0.bed > rep2Col.bed
/home/shriyas/Downloads/bedtools2/bin/intersectBed -s -F 1 -c -a rep2Col.bed -b ../rep2_nrpe1.bed > rep2All.bed
{ printf 'Chr\tStart\tStop\tDescription\tScore\tStrand\tMerged_dedup_Col\tMerged_dedup_Nrpe\tMerged_All_Col\tMerged_All_Nrpe\tRep1_Col\tRep1_Nrpe\tRep2_Col\tRep2_Nrpe\n'; cat rep2All.bed; } > head_rep2All.txt #Adding the header

