
echo Common Pre-process Steps 
samtools sort -o ./results/results/Merged_auto.sorted.bam ./results/results/Merged_auto.bam
samtools view -h -F 0x904 -b ./results/results/Merged_auto.sorted.bam > ./results/results/Merged_auto.sorted.uniq.bam
samtools flagstat ./results/results/Merged_auto.sorted.uniq.bam > ./results/results/Merged_auto.align.stats.txt
samtools addreplacerg -r '@RG\tID:results/Merged_auto\tSM:results/Merged_auto' ./results/results/Merged_auto.sorted.uniq.bam -o ./results/results/Merged_auto.sorted.uniq.rg.bam 
samtools depth ./results/results/Merged_auto.sorted.uniq.rg.bam -o ./results/results/Merged_auto.depth.txt 
samtools index ./results/results/Merged_auto.sorted.uniq.rg.bam
samtools index ./results/1/1_golden.bamsamtools index ./results/2/2_golden.bamsamtools index ./results/3/3_golden.bamsamtools index ./results/4/4_golden.bamsamtools index ./results/5/5_golden.bam
 
echo bam-readcount reporting 
bam_readcount_path/bam-readcount/build/bin/bam-readcount -w 0 -f /mnt/d/sfragkoul/Synth_Data/Synthesizers/NEAT/testing/TP53/TP53.fasta ./results/results/Merged_auto.sorted.uniq.rg.bam > ./results/results/Merged_auto_report.tsv
bam_readcount_path/bam-readcount/build/bin/bam-readcount -w 0 -f /mnt/d/sfragkoul/Synth_Data/Synthesizers/NEAT/testing/TP53/TP53.fasta ./results/1/1_golden.bam > ./results/1/1_report.tsv
bam_readcount_path/bam-readcount/build/bin/bam-readcount -w 0 -f /mnt/d/sfragkoul/Synth_Data/Synthesizers/NEAT/testing/TP53/TP53.fasta ./results/2/2_golden.bam > ./results/2/2_report.tsv
bam_readcount_path/bam-readcount/build/bin/bam-readcount -w 0 -f /mnt/d/sfragkoul/Synth_Data/Synthesizers/NEAT/testing/TP53/TP53.fasta ./results/3/3_golden.bam > ./results/3/3_report.tsv
bam_readcount_path/bam-readcount/build/bin/bam-readcount -w 0 -f /mnt/d/sfragkoul/Synth_Data/Synthesizers/NEAT/testing/TP53/TP53.fasta ./results/4/4_golden.bam > ./results/4/4_report.tsv
bam_readcount_path/bam-readcount/build/bin/bam-readcount -w 0 -f /mnt/d/sfragkoul/Synth_Data/Synthesizers/NEAT/testing/TP53/TP53.fasta ./results/5/5_golden.bam > ./results/5/5_report.tsv

echo Variant Calling 
Printing lofreq commands
lofreq indelqual --dindel -f /mnt/d/sfragkoul/Synth_Data/Synthesizers/NEAT/testing/TP53/TP53.fasta -o ./results/results/Merged_auto_indels.sorted.uniq.rg.bam ./results/results/Merged_auto.sorted.uniq.rg.bam 
lofreq call -f /mnt/d/sfragkoul/Synth_Data/Synthesizers/NEAT/testing/TP53/TP53.fasta --call-indels -o ./results/Lofreq.vcf ./results/results/Merged_auto_indels.sorted.uniq.rg.bam 
bcftools reheader --fai /mnt/d/sfragkoul/Synth_Data/Synthesizers/NEAT/testing/TP53/TP53.fasta.fai -o ./results/results/Merged_auto_Lofreq_norm.vcf ./results/Lofreq.vcf
rm ./results/Lofreq.vcf
 
