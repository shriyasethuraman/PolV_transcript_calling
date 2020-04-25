#!/bin/sh
intersectBed -c -s -a ../../SHORT_TRANSCRIPTS_12_1_15.bed -b rep2_col0.bed > temp1.bed
intersectBed -c -s -a temp1.bed -b rep2_nrpe1.bed > temp2.bed
intersectBed -c -s -a temp2.bed -b rep2_ago4.bed > temp3.bed
intersectBed -c -s -a temp3.bed -b rep2_idn2.bed > temp4.bed
