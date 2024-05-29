mkdir ./results
mkdir ./Plots

starting run 1
mkdir ./results/1
python gen_reads.py -r TP53.fasta --rng 213 -m 0.1 -r 100 -c 100 -o ./results/1/1 --pe 300 30 --bam --vcf

starting run 2
mkdir ./results/2
python gen_reads.py -r TP53.fasta --rng 214 -m 0.1 -r 100 -c 100 -o ./results/2/2 --pe 300 30 --bam --vcf

starting run 3
mkdir ./results/3
python gen_reads.py -r TP53.fasta --rng 215 -m 0.1 -r 100 -c 100 -o ./results/3/3 --pe 300 30 --bam --vcf

starting run 4
mkdir ./results/4
python gen_reads.py -r TP53.fasta --rng 217 -m 0.1 -r 100 -c 100 -o ./results/4/4 --pe 300 30 --bam --vcf

starting run 5
mkdir ./results/5
python gen_reads.py -r TP53.fasta --rng 218 -m 0.1 -r 100 -c 100 -o ./results/5/5 --pe 300 30 --bam --vcf

