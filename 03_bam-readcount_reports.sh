#!/bin/sh
#
#This bash script calls bam-readcount to produce the bam reports with the genomic content at certain chromosomal positions.
#
#Input: fasta reference file, "golden" bam file
#Output: tsv file with the genomic content
#
echo "Index bam files"
samtools index 1/1_golden.bam
samtools index 2/2_golden.bam
samtools index 3/3_golden.bam
samtools index 4/4_golden.bam
samtools index 5/5_golden.bam
samtools index 6/6_golden.bam
samtools index 7/7_golden.bam
samtools index 8/8_golden.bam
samtools index 9/9_golden.bam
samtools index 10/10_golden.bam

echo  "Run bam-readcount"
bam-readcount/build/bin/bam-readcount -w 0 -f TP53.fasta 1/1_golden.bam > 1/1.tsv
bam-readcount/build/bin/bam-readcount -w 0 -f TP53.fasta 2/2_golden.bam > 2/2.tsv
bam-readcount/build/bin/bam-readcount -w 0 -f TP53.fasta 3/3_golden.bam > 3/3.tsv
bam-readcount/build/bin/bam-readcount -w 0 -f TP53.fasta 4/4_golden.bam > 4/4.tsv
bam-readcount/build/bin/bam-readcount -w 0 -f TP53.fasta 5/5_golden.bam > 5/5.tsv
bam-readcount/build/bin/bam-readcount -w 0 -f TP53.fasta 6/6_golden.bam > 6/6.tsv
bam-readcount/build/bin/bam-readcount -w 0 -f TP53.fasta 7/7_golden.bam > 7/7.tsv
bam-readcount/build/bin/bam-readcount -w 0 -f TP53.fasta 8/8_golden.bam > 8/8.tsv
bam-readcount/build/bin/bam-readcount -w 0 -f TP53.fasta 9/9_golden.bam > 9/9.tsv
bam-readcount/build/bin/bam-readcount -w 0 -f TP53.fasta 10/10_golden.bam > 10/10.tsv
bam-readcount/build/bin/bam-readcount -w 0 -f TP53.fasta Merged2.bam > Merged2_report.tsv
