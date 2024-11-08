mkdir ./results
mkdir ./Plots

echo Generating Synthetic files
echo Starting run 1
mkdir ./results/1
python /NEAT/gen_reads.py -r /reference/TP53.fasta --rng 213 -M 0.1 -R 100 -c 100 -o ./results/1/1 --pe 300 30 --bam --vcf

echo Starting run 2
mkdir ./results/2
python /NEAT/gen_reads.py -r /reference/TP53.fasta --rng 214 -M 0.1 -R 100 -c 100 -o ./results/2/2 --pe 300 30 --bam --vcf

echo Starting run 3
mkdir ./results/3
python /NEAT/gen_reads.py -r /reference/TP53.fasta --rng 215 -M 0.1 -R 100 -c 100 -o ./results/3/3 --pe 300 30 --bam --vcf

echo Starting run 4
mkdir ./results/4
python /NEAT/gen_reads.py -r /reference/TP53.fasta --rng 217 -M 0.1 -R 100 -c 100 -o ./results/4/4 --pe 300 30 --bam --vcf

echo Starting run 5
mkdir ./results/5
python /NEAT/gen_reads.py -r /reference/TP53.fasta --rng 218 -M 0.1 -R 100 -c 100 -o ./results/5/5 --pe 300 30 --bam --vcf

echo Merging bam files
samtools merge ./results/Merged_auto.bam ./results/1/1_golden.bam ./results/2/2_golden.bam ./results/3/3_golden.bam ./results/4/4_golden.bam ./results/5/5_golden.bam
echo Merging vcf files
bcftools index ./results/1/1_golden.vcf.gz 
bcftools index ./results/2/2_golden.vcf.gz 
bcftools index ./results/3/3_golden.vcf.gz 
bcftools index ./results/4/4_golden.vcf.gz 
bcftools index ./results/5/5_golden.vcf.gz 
bcftools merge ./results/1/1_golden.vcf.gz ./results/2/2_golden.vcf.gz ./results/3/3_golden.vcf.gz ./results/4/4_golden.vcf.gz ./results/5/5_golden.vcf.gz > ./results/Merged_auto_ground_truth.vcf
bcftools norm ./results/Merged_auto_ground_truth.vcf --output ./results/Merged_auto_ground_truth_norm.vcf --output-type v -m "-"
rm ./results/Merged_auto_ground_truth.vcf
