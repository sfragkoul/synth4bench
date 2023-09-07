#!/bin/sh
#
#This bach script calls NEAT in order to generate 10 individual synthetic data datasets.
#
#Input: fasta reference file
#Output: fastq files with pair end reads, "golden" bam file and bai index file, "golden" vcf file
#
#Individual files
echo  "Starting Run 1"
mkdir 1
python gen_reads.py -r TP53.fasta -R 150 -c 500 -M 0.1 -o 1 --pe 300 30 --bam --vcf


echo  "Starting Run 2"
mkdir 2
python gen_reads.py -r TP53.fasta -R 150 -c 500 -M 0.1 -o 2 --pe 300 30 --bam --vcf


echo  "Starting Run 3"
mkdir 3
python gen_reads.py -r TP53.fasta -R 150 -c 500 -M 0.1 -o 3 --pe 300 30 --bam --vcf


echo  "Starting Run 4"
mkdir 4
python gen_reads.py -r TP53.fasta -R 150 -c 500 -M 0.1 -o 4 --pe 300 30 --bam --vcf

echo  "Starting Run 5"
mkdir 5
python gen_reads.py -r TP53.fasta -R 150 -c 500 -M 0.1 -o 5 --pe 300 30 --bam --vcf


echo  "Starting Run 6"
mkdir 6
python gen_reads.py -r TP53.fasta -R 150 -c 500 -M 0.1 -o 6 --pe 300 30 --bam --vcf


echo  "Starting Run 7"
mkdir 7
python gen_reads.py -r TP53.fasta -R 150 -c 500 -M 0.1 -o 7 --pe 300 30 --bam --vcf


echo  "Starting Run 8"
mkdir 8
python gen_reads.py -r TP53.fasta -R 150 -c 500 -M 0.1 -o 8 --pe 300 30 --bam --vcf


echo  "Starting Run 9"
mkdir 9
python gen_reads.py -r TP53.fasta -R 150 -c 500 -M 0.1 -o 9 --pe 300 30 --bam --vcf


echo  "Starting Run 10"
mkdir 10
python gen_reads.py -r TP53.fasta -R 150 -c 500 -M 0.1 -o 10 --pe 300 30 --bam --vcf

#Merged file
echo  "Merging bam files"
samtools merge Merged2.bam 1/1_golden.bam 2/2_golden.bam 3/3_golden.bam 4/4_golden.bam 5/5_golden.bam 6/6_golden.bam 7/7_golden.bam 8/8_golden.bam 9/9_golden.bam 10/10_golden.bam
