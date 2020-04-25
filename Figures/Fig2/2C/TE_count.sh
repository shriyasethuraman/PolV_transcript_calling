#!/bin/sh
awk '{OFS="\t"}{s+=1; if ($2=="TRUE") print $1,$3,$4,$1,"255","+"; else print $1,$3,$4,$1,"255","-"}' Transposons > Trans1
sed 1d Trans1 > Trans2
cut -c3 Trans2 > Trans_chr
cut -f2-6 Trans2 > Trans3
paste Trans_chr Trans3 > Trans.bed

