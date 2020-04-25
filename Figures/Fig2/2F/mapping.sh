cd ~/prog_run/RIP_manu_met1/

~/Downloads/sratoolkit.2.5.5-ubuntu64/bin/fastq-dump --outdir col0/ --split-files SRR1023819.sra
~/Downloads/sratoolkit.2.5.5-ubuntu64/bin/fastq-dump --outdir met1-3/ --split-files SRR1023820.sra

cd col0/
~/cifs-lab/Shriya/FastQC/fastqc SRR1023819_1.fastq 
#~/Downloads/trim_galore_zip/trim_galore --fastqc --length 50 SRR1023819_1.fastq

cutadapt --discard-trimmed -f fastq -a AGATCGGAAGAGC SRR1023819_1.fastq -o SRR1023819_1_no_adapter.fq 

#bowtie2 -p 4 -q -N 1 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -U SRR1023819_1.fastq -S mapped/untrimmed_col0.sam
#bowtie2 -p 4 -q -N 1 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -U SRR1023819_1.fastq -S mapped/trimmed_col0.sam
#bowtie2 -p 4 -q -N 1 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -U SRR1023819_1_no_adapter.fq -S mapped/NO_adapt_trimmed_col0.sam


bowtie -q -n 1 -m 1 --al --best -p 6 --sam ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie_index/bowtie_indx SRR1023819_1_no_adapter.fq mapped/NO_adapt_bowtie_trimmed_col0.sam 
samtools view -S -h -F 4 mapped/NO_adapt_bowtie_trimmed_col0.sam > mapped/aligned.sam
samtools view -bS mapped/aligned.sam > mapped/aligned.bam
bamToBed -i mapped/aligned.bam > mapped/aligned.bed
sort -k1 -k2 -n mapped/aligned.bed > mapped/sort_aligned.bed


cd ../met1-3/
~/cifs-lab/Shriya/FastQC/fastqc SRR1023820_1.fastq 

#~/Downloads/trim_galore_zip/trim_galore --fastqc --length 18 SRR1023820_1.fastq

cutadapt --discard-trimmed -f fastq -a AGATCGGAAGAGC SRR1023820_1.fastq -o SRR1023820_1_no_adapter.fq 

#bowtie2 -p 4 -q -N 1 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -U SRR1023820_1_trimmed.fq -S mapped/trimmed_met1-3.sam
#bowtie2 -p 4 -q -N 1 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -U SRR1023820_1_no_adapter.fq -S mapped/NO_adapt_trimmed_met1-3.sam

bowtie -n 1 -m 1 --al --best --sam -p 4 -q ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie_index/bowtie_indx SRR1023820_1_no_adapter.fq mapped/NO_adapt_bowtie_trimmed_met1-3.sam 
samtools view -S -h -F 4 mapped/NO_adapt_bowtie_trimmed_met1-3.sam > mapped/aligned.sam
samtools view -bS mapped/aligned.sam > mapped/aligned.bam
bamToBed -i mapped/aligned.bam > mapped/aligned.bed
sort -k1 -k2 -n mapped/aligned.bed > mapped/sort_aligned.bed

