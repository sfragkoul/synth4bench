#!/bin/sh
#
#This bash script is the basis of synth4bench workflow. It calls NEAT in order to generate 10 individual synthetic data datasets, 
#create one Merged bam file, performs some preprocess steps before implementing somatic variant calling and produces bam report files 
#with the genomic content at certain chromosomal positionsusing bam-readcount.
#
#Input:  fasta reference file
#Output: fastq files with pair end reads, "golden" bam file and bai index file, "golden" vcf file, 
#		 Merged bam file, processed bam files and vcf file with all variants that were detected, tsv file with the genomic content
#
#please replace all "path/to/generated/files/", "path/to/gatk/", "path/to/reference/" and  "path/to/bam-readcount/"  with desired folders

mkdir path/to/generated/files
mkdir path/to/generated/files/Plots

echo  "\n Starting Run 1"
mkdir path/to/generated/files/1
python gen_reads.py --rng 154 -r path/to/reference/TP53.fasta -M 0.1 -R 50 -c 100 -o path/to/generated/files/1/1 --pe 300 30 --bam --vcf

echo  "\n Starting Run 2"
mkdir path/to/generated/files/2
python gen_reads.py --rng 156 -r path/to/reference/TP53.fasta -M 0.1 -R 50 -c 100 -o path/to/generated/files/2/2 --pe 300 30 --bam --vcf


echo  "\n Starting Run 3"
mkdir path/to/generated/files/3
python gen_reads.py --rng 157 -r path/to/reference/TP53.fasta -M 0.1 -R 50 -c 100 -o path/to/generated/files/3/3 --pe 300 30 --bam --vcf


echo  "\n Starting Run 4"
mkdir path/to/generated/files/4
python gen_reads.py --rng 158 -r path/to/reference/TP53.fasta -M 0.1 -R 50 -c 100 -o path/to/generated/files/4/4 --pe 300 30 --bam --vcf

echo  "\n Starting Run 5"
mkdir path/to/generated/files/5
python gen_reads.py --rng 159 -r path/to/reference/TP53.fasta -M 0.1 -R 50 -c 100 -o path/to/generated/files/5/5 --pe 300 30 --bam --vcf


echo  "\n Starting Run 6"
mkdir path/to/generated/files/6
python gen_reads.py --rng 162 -r path/to/reference/TP53.fasta -M 0.1 -R 50 -c 100 -o path/to/generated/files/6/6 --pe 300 30 --bam --vcf


echo  "\n Starting Run 7"
mkdir path/to/generated/files/7
python gen_reads.py --rng 163 -r path/to/reference/TP53.fasta -M 0.1 -R 50 -c 100 -o path/to/generated/files/7/7 --pe 300 30 --bam --vcf


echo  "\n Starting Run 8"
mkdir path/to/generated/files/8
python gen_reads.py --rng 164 -r path/to/reference/TP53.fasta -M 0.1 -R 50 -c 100 -o path/to/generated/files/8/8 --pe 300 30 --bam --vcf


echo  "\n Starting Run 9"
mkdir path/to/generated/files/9
python gen_reads.py --rng 166 -r path/to/reference/TP53.fasta -M 0.1 -R 50 -c 100 -o path/to/generated/files/9/9 --pe 300 30 --bam --vcf


echo  "\n Starting Run 10"
mkdir path/to/generated/files/10
python gen_reads.py --rng 169 -r path/to/reference/TP53.fasta -M 0.1 -R 50 -c 100 -o path/to/generated/files/10/10 --pe 300 30 --bam --vcf

#Merged file
echo  "\n Merging bam files"
 samtools merge path/to/generated/files/Merged.bam path/to/generated/files/1/1_golden.bam path/to/generated/files/2/2_golden.bam path/to/generated/files/3/3_golden.bam path/to/generated/files/4/4_golden.bam path/to/generated/files/5/5_golden.bam path/to/generated/files/6/6_golden.bam path/to/generated/files/7/7_golden.bam path/to/generated/files/8/8_golden.bam path/to/generated/files/9/9_golden.bam path/to/generated/files/10/10_golden.bam

echo  "\n Merging vcf files"
bcftools index path/to/generated/files/1/1_golden.vcf.gz
bcftools index path/to/generated/files/2/2_golden.vcf.gz
bcftools index path/to/generated/files/3/3_golden.vcf.gz
bcftools index path/to/generated/files/4/4_golden.vcf.gz
bcftools index path/to/generated/files/5/5_golden.vcf.gz
bcftools index path/to/generated/files/6/6_golden.vcf.gz
bcftools index path/to/generated/files/7/7_golden.vcf.gz
bcftools index path/to/generated/files/8/8_golden.vcf.gz
bcftools index path/to/generated/files/9/9_golden.vcf.gz
bcftools index path/to/generated/files/10/10_golden.vcf.gz

bcftools merge path/to/generated/files/1/1_golden.vcf.gz path/to/generated/files/2/2_golden.vcf.gz path/to/generated/files/3/3_golden.vcf.gz path/to/generated/files/4/4_golden.vcf.gz  path/to/generated/files/5/5_golden.vcf.gz path/to/generated/files/6/6_golden.vcf.gz path/to/generated/files/7/7_golden.vcf.gz path/to/generated/files/8/8_golden.vcf.gz path/to/generated/files/9/9_golden.vcf.gz path/to/generated/files/10/10_golden.vcf.gz > path/to/generated/files/Merged_ground_truth.vcf


echo  "\n Variant Calling" 
samtools sort -o path/to/generated/files/Merged.sorted.bam path/to/generated/files/Merged.bam 
samtools view -h -F 0x904 -b path/to/generated/files/Merged.sorted.bam > path/to/generated/files/Merged.sorted.uniq.bam 
samtools flagstat path/to/generated/files/Merged.sorted.uniq.bam > path/to/generated/files/Merged.align.stats.txt  
samtools addreplacerg -r '@RG\tID:Merged\tSM:Merged' path/to/generated/files/Merged.sorted.uniq.bam -o path/to/generated/files/Merged.sorted.uniq.rg.bam
samtools depth path/to/generated/files/Merged.sorted.uniq.rg.bam -o path/to/generated/files/Merged.depth.txt
samtools index path/to/generated/files/Merged.sorted.uniq.rg.bam
java -jar path/to/gatk/gatk-4.3.0.0.jar Mutect2 --reference path/to/reference/TP53.fasta --input path/to/generated/files/Merged.sorted.uniq.rg.bam --tumor-sample Merged  --output path/to/generated/files/Merged_GATK.vcf > path/to/generated/files/Merged.Mutect.out 2>&1
 

echo  "\n bam-readcount reporting" 
samtools index path/to/generated/files/1/1_golden.bam
samtools index path/to/generated/files/2/2_golden.bam
samtools index path/to/generated/files/3/3_golden.bam
samtools index path/to/generated/files/4/4_golden.bam
samtools index path/to/generated/files/5/5_golden.bam
samtools index path/to/generated/files/6/6_golden.bam
samtools index path/to/generated/files/7/7_golden.bam
samtools index path/to/generated/files/8/8_golden.bam
samtools index path/to/generated/files/9/9_golden.bam
samtools index path/to/generated/files/10/10_golden.bam

  
#conda deactivate 
#conda activate bam-readcount
##cd path/to/bam-readcount/build/bin/bam-readcount
path/to/bam-readcount/build/bin/bam-readcount -w 0 -f path/to/reference/TP53.fasta path/to/generated/files/Merged.sorted.uniq.rg.bam > path/to/generated/files/Merged_report.tsv

path/to/bam-readcount/build/bin/bam-readcount -w 0 -f path/to/reference/TP53.fasta path/to/generated/files/1/1_golden.bam > path/to/generated/files/1/1_report.tsv
path/to/bam-readcount/build/bin/bam-readcount -w 0 -f path/to/reference/TP53.fasta path/to/generated/files/2/2_golden.bam > path/to/generated/files/2/2_report.tsv
path/to/bam-readcount/build/bin/bam-readcount -w 0 -f path/to/reference/TP53.fasta path/to/generated/files/3/3_golden.bam > path/to/generated/files/3/3_report.tsv
path/to/bam-readcount/build/bin/bam-readcount -w 0 -f path/to/reference/TP53.fasta path/to/generated/files/4/4_golden.bam > path/to/generated/files/4/4_report.tsv
path/to/bam-readcount/build/bin/bam-readcount -w 0 -f path/to/reference/TP53.fasta path/to/generated/files/5/5_golden.bam > path/to/generated/files/5/5_report.tsv
path/to/bam-readcount/build/bin/bam-readcount -w 0 -f path/to/reference/TP53.fasta path/to/generated/files/6/6_golden.bam > path/to/generated/files/6/6_report.tsv
path/to/bam-readcount/build/bin/bam-readcount -w 0 -f path/to/reference/TP53.fasta path/to/generated/files/7/7_golden.bam > path/to/generated/files/7/7_report.tsv
path/to/bam-readcount/build/bin/bam-readcount -w 0 -f path/to/reference/TP53.fasta path/to/generated/files/8/8_golden.bam > path/to/generated/files/8/8_report.tsv
path/to/bam-readcount/build/bin/bam-readcount -w 0 -f path/to/reference/TP53.fasta path/to/generated/files/9/9_golden.bam > path/to/generated/files/9/9_report.tsv
path/to/bam-readcount/build/bin/bam-readcount -w 0 -f path/to/reference/TP53.fasta path/to/generated/files/10/10_golden.bam > path/to/generated/files/10/10_report.tsv





 #gatk Mutect2 --R "./reference/TP53.fasta" --input "./Merged2.sorted.uniq.rg.bam" --tumor-sample Merged2  --output Merged2.vcf > Merged2.Mutect2.out 2>&1