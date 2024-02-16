#!/bin/sh
#
#This bash script is the basis of synth4bench workflow. It calls NEAT in order to generate 10 individual synthetic data datasets, 
#create one Merged bam file, performs some preprocess steps before implementing somatic variant calling and produces bam report files 
#with the genomic content at certain chromosomal positions using bam-readcount.
#
#Input:  fasta reference file
#Output: fastq files with pair end reads, "golden" bam file and bai index file, "golden" vcf file, 
#		 Merged bam file, processed bam files and vcf file with all variants that were detected, tsv file with the genomic content
#
#please replace all "path/to/files/"  with desired folders


mkdir path/to/files/1000_100
mkdir path/to/files/Plots

echo  "\n Starting Run 1"
mkdir path/to/files/1
python gen_reads.py --rng 213 -r path/to/files/TP53.fasta -M 0.1 -R 100 -c 100 -o path/to/files/1/1 --pe 300 30 --bam --vcf

echo  "\n Starting Run 2"
mkdir path/to/files/2
python gen_reads.py --rng 214 -r path/to/files/TP53.fasta -M 0.1 -R 100 -c 100 -o path/to/files/2/2 --pe 300 30 --bam --vcf


echo  "\n Starting Run 3"
mkdir path/to/files/3
python gen_reads.py --rng 215 -r path/to/files/TP53.fasta -M 0.1 -R 100 -c 100 -o path/to/files/3/3 --pe 300 30 --bam --vcf


echo  "\n Starting Run 4"
mkdir path/to/files/4
python gen_reads.py --rng 217 -r path/to/files/TP53.fasta -M 0.1 -R 100 -c 100 -o path/to/files/4/4 --pe 300 30 --bam --vcf

echo  "\n Starting Run 5"
mkdir path/to/files/5
python gen_reads.py --rng 218 -r path/to/files/TP53.fasta -M 0.1 -R 100 -c 100 -o path/to/files/5/5 --pe 300 30 --bam --vcf


echo  "\n Starting Run 6"
mkdir path/to/files/6
python gen_reads.py --rng 225 -r path/to/files/TP53.fasta -M 0.1 -R 100 -c 100 -o path/to/files/6/6 --pe 300 30 --bam --vcf


echo  "\n Starting Run 7"
mkdir path/to/files/7
python gen_reads.py --rng 226 -r path/to/files/TP53.fasta -M 0.1 -R 100 -c 100 -o path/to/files/7/7 --pe 300 30 --bam --vcf


echo  "\n Starting Run 8"
mkdir path/to/files/8
python gen_reads.py --rng 228 -r path/to/files/TP53.fasta -M 0.1 -R 100 -c 100 -o path/to/files/8/8 --pe 300 30 --bam --vcf


echo  "\n Starting Run 9"
mkdir path/to/files/9
python gen_reads.py --rng 229 -r path/to/files/TP53.fasta -M 0.1 -R 100 -c 100 -o path/to/files/9/9 --pe 300 30 --bam --vcf


echo  "\n Starting Run 10"
mkdir path/to/files/10
python gen_reads.py --rng 230 -r path/to/files/TP53.fasta -M 0.1 -R 100 -c 100 -o path/to/files/10/10 --pe 300 30 --bam --vcf

#Merged file
echo  "\n Merging bam files"
samtools merge path/to/files/Merged.bam path/to/files/1/1_golden.bam path/to/files/2/2_golden.bam path/to/files/3/3_golden.bam path/to/files/4/4_golden.bam path/to/files/5/5_golden.bam path/to/files/6/6_golden.bam path/to/files/7/7_golden.bam path/to/files/8/8_golden.bam path/to/files/9/9_golden.bam path/to/files/10/10_golden.bam

echo  "\n Merging vcf files"
bcftools index path/to/files/1/1_golden.vcf.gz
bcftools index path/to/files/2/2_golden.vcf.gz
bcftools index path/to/files/3/3_golden.vcf.gz
bcftools index path/to/files/4/4_golden.vcf.gz
bcftools index path/to/files/5/5_golden.vcf.gz
bcftools index path/to/files/6/6_golden.vcf.gz
bcftools index path/to/files/7/7_golden.vcf.gz
bcftools index path/to/files/8/8_golden.vcf.gz
bcftools index path/to/files/9/9_golden.vcf.gz
bcftools index path/to/files/10/10_golden.vcf.gz
bcftools merge path/to/files/1/1_golden.vcf.gz path/to/files/2/2_golden.vcf.gz path/to/files/3/3_golden.vcf.gz path/to/files/4/4_golden.vcf.gz  path/to/files/5/5_golden.vcf.gz path/to/files/6/6_golden.vcf.gz path/to/files/7/7_golden.vcf.gz path/to/files/8/8_golden.vcf.gz path/to/files/9/9_golden.vcf.gz path/to/files/10/10_golden.vcf.gz > path/to/files/Merged_ground_truth.vcf

bcftools norm path/to/files/Merged_ground_truth.vcf --output path/to/files/Merged_ground_truth_norm.vcf --output-type v -m "-"

echo  "\n Variant Calling" 
samtools sort -o path/to/files/Merged.sorted.bam path/to/files/Merged.bam 
samtools view -h -F 0x904 -b path/to/files/Merged.sorted.bam > path/to/files/Merged.sorted.uniq.bam 
samtools flagstat path/to/files/Merged.sorted.uniq.bam > path/to/files/Merged.align.stats.txt  
samtools addreplacerg -r '@RG\tID:Merged\tSM:Merged' path/to/files/Merged.sorted.uniq.bam -o path/to/files/Merged.sorted.uniq.rg.bam
samtools depth path/to/files/Merged.sorted.uniq.rg.bam -o path/to/files/Merged.depth.txt
samtools index path/to/files/Merged.sorted.uniq.rg.bam
#GATK
java -jar testing/gatk-4.3.0.0.jar Mutect2 --reference path/to/files/TP53.fasta --input path/to/files/Merged.sorted.uniq.rg.bam --tumor-sample Merged  --output path/to/files/Merged_GATK.vcf > path/to/files/Merged.Mutect.out 2>&1
bcftools norm path/to/files/Merged_GATK.vcf --output path/to/files/Merged_GATK_norm.vcf --output-type v -m "-"

#freebayes
freebayes --fasta-reference path/to/files/TP53.fasta --bam path/to/files/Merged.sorted.uniq.rg.bam > path/to/files/freebayes.vcf 
bcftools reheader --fai path/to/files/TP53.fasta.fai -o path/to/files/Merged_freebayes.vcf  path/to/files/freebayes.vcf
rm path/to/files/freebayes.vcf 
bcftools norm path/to/files/Merged_freebayes.vcf --output path/to/files/Merged_freebayes_norm.vcf --output-type v -m "-"


#vardict
vardict-java -G path/to/files/TP53.fasta -f 0.0001 -N Merged -b path/to/files/Merged.sorted.uniq.rg.bam -R hg38_knownGene_ENST00000610292.4:0-19080 > path/to/files/Merged_VarDict.txt
var2vcf_valid.pl -A path/to/files/Merged_VarDict.txt >  path/to/files/VarDict.vcf
rm path/to/files/Merged_VarDict.txt
bcftools reheader --fai path/to/files/TP53.fasta.fai -o path/to/files/Merged_VarDict.vcf  path/to/files/VarDict.vcf
rm path/to/files/VarDict.vcf
bcftools norm path/to/files/Merged_VarDict.vcf --output path/to/files/Merged_VarDict_norm.vcf --output-type v -m "-"

#varcan
samtools mpileup -f path/to/files/TP53.fasta path/to/files/Merged.sorted.uniq.rg.bam -a -o path/to/files/Merged.sorted.uniq.rg.mpileup
varscan pileup2cns path/to/files/Merged.sorted.uniq.rg.mpileup > path/to/files/VarScan.tsv
python Varscan2VCF/vscan_pileup2cns2vcf.py path/to/files/VarScan.tsv > path/to/files/Merged_VarScan.vcf
rm path/to/files/VarScan.tsv 
bcftools norm path/to/files/Merged_VarScan.vcf --output path/to/files/Merged_VarScan_norm.vcf --output-type v -m "-"

#lofreq
lofreq indelqual --dindel -f testing/TP53/TP53.fasta -o testing/TP53/read_length/1000_100/Merged_indels.sorted.uniq.rg.bam testing/TP53/read_length/1000_100/Merged.sorted.uniq.rg.bam
lofreq call -f testing/TP53/TP53.fasta --call-indels -o testing/TP53/read_length/1000_100/Lofreq.vcf testing/TP53/read_length/1000_100/Merged_indels.sorted.uniq.rg.bam
bcftools reheader --fai testing/TP53/TP53.fasta.fai -o testing/TP53/read_length/1000_100/Merged_Lofreq_norm.vcf  testing/TP53/read_length/1000_100/Lofreq.vcf
rm testing/TP53/read_length/1000_100/Lofreq.vcf

echo  "\n bam-readcount reporting" 
samtools index path/to/files/1/1_golden.bam
samtools index path/to/files/2/2_golden.bam
samtools index path/to/files/3/3_golden.bam
samtools index path/to/files/4/4_golden.bam
samtools index path/to/files/5/5_golden.bam
samtools index path/to/files/6/6_golden.bam
samtools index path/to/files/7/7_golden.bam
samtools index path/to/files/8/8_golden.bam
samtools index path/to/files/9/9_golden.bam
samtools index path/to/files/10/10_golden.bam

/home/sfragkoul/bam-readcount/build/bin/bam-readcount -w 0 -f path/to/files/TP53.fasta path/to/files/Merged.sorted.uniq.rg.bam > path/to/files/Merged_report.tsv
/home/sfragkoul/bam-readcount/build/bin/bam-readcount -w 0 -f path/to/files/TP53.fasta path/to/files/1/1_golden.bam > path/to/files/1/1_report.tsv
/home/sfragkoul/bam-readcount/build/bin/bam-readcount -w 0 -f path/to/files/TP53.fasta path/to/files/2/2_golden.bam > path/to/files/2/2_report.tsv
/home/sfragkoul/bam-readcount/build/bin/bam-readcount -w 0 -f path/to/files/TP53.fasta path/to/files/3/3_golden.bam > path/to/files/3/3_report.tsv
/home/sfragkoul/bam-readcount/build/bin/bam-readcount -w 0 -f path/to/files/TP53.fasta path/to/files/4/4_golden.bam > path/to/files/4/4_report.tsv
/home/sfragkoul/bam-readcount/build/bin/bam-readcount -w 0 -f path/to/files/TP53.fasta path/to/files/5/5_golden.bam > path/to/files/5/5_report.tsv
/home/sfragkoul/bam-readcount/build/bin/bam-readcount -w 0 -f path/to/files/TP53.fasta path/to/files/6/6_golden.bam > path/to/files/6/6_report.tsv
/home/sfragkoul/bam-readcount/build/bin/bam-readcount -w 0 -f path/to/files/TP53.fasta path/to/files/7/7_golden.bam > path/to/files/7/7_report.tsv
/home/sfragkoul/bam-readcount/build/bin/bam-readcount -w 0 -f path/to/files/TP53.fasta path/to/files/8/8_golden.bam > path/to/files/8/8_report.tsv
/home/sfragkoul/bam-readcount/build/bin/bam-readcount -w 0 -f path/to/files/TP53.fasta path/to/files/9/9_golden.bam > path/to/files/9/9_report.tsv
/home/sfragkoul/bam-readcount/build/bin/bam-readcount -w 0 -f path/to/files/TP53.fasta path/to/files/10/10_golden.bam > path/to/files/10/10_report.tsv
