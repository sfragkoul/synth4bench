#!/bin/sh
echo  "Starting Run 1"
mkdir 1_TP53
python gen_reads.py -r testing/TP53/TP53.fasta -R 150 -c 500 -M 0.1 -o testing/TP53/1_TP53 --pe 300 30 --bam --vcf


echo  "Starting Run 2"
mkdir 2_TP53
python gen_reads.py -r testing/TP53/TP53.fasta -R 150 -c 500 -M 0.1 -o testing/TP53/2_TP53 --pe 300 30 --bam --vcf


echo  "Starting Run 3"
mkdir 3_TP53
python gen_reads.py -r testing/TP53/TP53.fasta -R 150 -c 500 -M 0.1 -o testing/TP53/3_TP53 --pe 300 30 --bam --vcf


echo  "Starting Run 4"
mkdir 4_TP53
python gen_reads.py -r testing/TP53/TP53.fasta -R 150 -c 500 -M 0.1 -o testing/TP53/4_TP53 --pe 300 30 --bam --vcf

echo  "Starting Run 5"
mkdir 5_TP53
python gen_reads.py -r testing/TP53/TP53.fasta -R 150 -c 500 -M 0.1 -o testing/TP53/5_TP53 --pe 300 30 --bam --vcf


echo  "Starting Run 6"
mkdir 6_TP53
python gen_reads.py -r testing/TP53/TP53.fasta -R 150 -c 500 -M 0.1 -o testing/TP53/6_TP53 --pe 300 30 --bam --vcf


echo  "Starting Run 7"
mkdir 7_TP53
python gen_reads.py -r testing/TP53/TP53.fasta -R 150 -c 500 -M 0.1 -o testing/TP53/7_TP53 --pe 300 30 --bam --vcf


echo  "Starting Run 8"
mkdir 8_TP53
python gen_reads.py -r testing/TP53/TP53.fasta -R 150 -c 500 -M 0.1 -o testing/TP53/8_TP53 --pe 300 30 --bam --vcf


echo  "Starting Run 9"
mkdir 9_TP53
python gen_reads.py -r testing/TP53/TP53.fasta -R 150 -c 500 -M 0.1 -o testing/TP53/9_TP53 --pe 300 30 --bam --vcf


echo  "Starting Run 10"
mkdir 10_TP53
python gen_reads.py -r testing/TP53/TP53.fasta -R 150 -c 500 -M 0.1 -o testing/TP53/10_TP53 --pe 300 30 --bam --vcf