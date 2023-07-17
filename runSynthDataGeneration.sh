#!/bin/sh
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