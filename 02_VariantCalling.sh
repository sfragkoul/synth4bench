#!/bin/sh
#
#This bash script takes the Merged bam file and performs some preprocess steps before using Mutect2 to implement somatic variant calling.
#
#Input: Merged bam file, fasta reference file
#Output: processed bam files and vcf file with all variants that were detected
#
samtools sort -o Merged2.sorted.bam Merged2.bam 
samtools view -h -F 0x904 -b Merged2.sorted.bam > Merged2.sorted.uniq.bam 
samtools flagstat Merged2.sorted.uniq.bam > Merged2.align.stats.txt  
samtools addreplacerg -r '@RG\tID:Merged2\tSM:Merged2' Merged2.sorted.uniq.bam -o Merged2.sorted.uniq.rg.bam
samtools depth Merged2.sorted.uniq.rg.bam -o Merged2.depth.txt
samtools index Merged2.sorted.uniq.rg.bam
java -jar gatk-4.3.0.0.jar Mutect2 --reference TP53.fasta --input Merged2.sorted.uniq.rg.bam --tumor-sample Merged2  --output Merged2.vcf > Merged2.Mutect2.out 2>&1
