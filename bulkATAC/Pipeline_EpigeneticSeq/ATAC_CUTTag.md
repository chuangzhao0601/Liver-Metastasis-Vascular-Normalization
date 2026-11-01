## Pre-processing Sample 

Use trimmomatic to remove Nextera Transposase and adapter sequence:

```
trimmomatic PE -threads 4 -phred33 input_forward.fq.gz input_reverse.fq.gz output_forward_paired.fq.gz output_forward_unpaired.fq.gz output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz ILLUMINACLIP:adapter.fa:2:30:10 SLIDINGWINDOW:5:20 MINLEN:36
```

## Install Bowtie2

We will focus on Bowtie2 and samtools here.

```
conda install bowtie2
conda install samtools
```

* Build reference

In order to align reads to a genome, the reference sequence must be indexed. It is a time- and memory-intense procedure, but it needs to be done only once for a given genome.

```
## mouse
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz
rm chromFa.tar.gz *random* *chrUn*
bowtie2-build chr1.fa,chr10.fa,chr11.fa,chr12.fa,chr13.fa,chr14.fa,chr15.fa,chr16.fa,chr17.fa,chr18.fa,chr19.fa,chr2.fa,chr3.fa,chr4.fa,chr5.fa,chr6.fa,chr7.fa,chr8.fa,chr9.fa,chrM.fa,chrX.fa,chrY.fa ../bowtie2_Ms/mm10
```

* Alignment

```
## Alignment
bowtie2 --local --very-sensitive --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 -x /root/vol_1T/Reference/bowtie2_Ms/mm10 -1 name_1.fastq.gz  -2 name_2.fastq.gz > sample_name.sam

## sam to bam and sort bam 
samtools view -@ 4 -Sb -q 10 sample.sam | samtools sort -@ 4 - > sample_sort.bam

## Index
samtools index sample_name_sort.bam

## visualize reads
samtools view sample_name_sort.bam | head
bamToBed -i sample_name_sort.bam | head
```

## Remove Duplicated Reads and Low Quality Mapped Reads

* install picard (one-time)

Download `picard.jar` from the [PICARD](https://broadinstitute.github.io/picard/) main page. Run `java -jar picard.jar -h` to test if picard has been installed successfully.

```
java -jar ~/picard.jar MarkDuplicates -I sample_sort.bam -O sample_sort_rmDup.bam -REMOVE_DUPLICATES true -M sample_RmDup_metrics.txt
bamToBed -i sample_sort_rmDup.bam | wc -l
```

## Remove Reads in blacklist

```
bedtools intersect -v -abam ${PWD##*/}_rmDup.bam -b /root/vol_1T/Reference/BlackList/mm10-blacklist.v2.bed > ${PWD##*/}_clean.bam
samtools index ${PWD##*/}_clean.bam
bamToBed -i ${PWD##*/}_clean.bam | wc -l
```