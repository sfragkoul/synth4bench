
echo Common Pre-process Steps 
samtools sort -o ./results/Merged_auto.sorted.bam ./results/Merged_auto.bam
samtools view -h -F 0x904 -b ./results/Merged_auto.sorted.bam > ./results/Merged_auto.sorted.uniq.bam
samtools flagstat ./results/Merged_auto.sorted.uniq.bam > ./results/Merged_auto.align.stats.txt
samtools addreplacerg -r '@RG\tID:Merged_auto\tSM:Merged_auto' ./results/Merged_auto.sorted.uniq.bam -o ./results/Merged_auto.sorted.uniq.rg.bam 
samtools depth ./results/Merged_auto.sorted.uniq.rg.bam -o ./results/Merged_auto.depth.txt 
samtools index ./results/Merged_auto.sorted.uniq.rg.bam
samtools index ./results/1/1_golden.bam 
samtools index ./results/2/2_golden.bam 
samtools index ./results/3/3_golden.bam 
samtools index ./results/4/4_golden.bam 
samtools index ./results/5/5_golden.bam 

 
echo bam-readcount reporting 
/home/sfragkoul/bam-readcount/build/bin/bam-readcount -w 0 -f /mnt/d/sfragkoul/Synth_Data/Synthesizers/NEAT/testing/TP53/TP53.fasta ./results/Merged_auto.sorted.uniq.rg.bam > ./results/Merged_auto_report.tsv
/home/sfragkoul/bam-readcount/build/bin/bam-readcount -w 0 -f /mnt/d/sfragkoul/Synth_Data/Synthesizers/NEAT/testing/TP53/TP53.fasta ./results/1/1_golden.bam > ./results/1/1_report.tsv
/home/sfragkoul/bam-readcount/build/bin/bam-readcount -w 0 -f /mnt/d/sfragkoul/Synth_Data/Synthesizers/NEAT/testing/TP53/TP53.fasta ./results/2/2_golden.bam > ./results/2/2_report.tsv
/home/sfragkoul/bam-readcount/build/bin/bam-readcount -w 0 -f /mnt/d/sfragkoul/Synth_Data/Synthesizers/NEAT/testing/TP53/TP53.fasta ./results/3/3_golden.bam > ./results/3/3_report.tsv
/home/sfragkoul/bam-readcount/build/bin/bam-readcount -w 0 -f /mnt/d/sfragkoul/Synth_Data/Synthesizers/NEAT/testing/TP53/TP53.fasta ./results/4/4_golden.bam > ./results/4/4_report.tsv
/home/sfragkoul/bam-readcount/build/bin/bam-readcount -w 0 -f /mnt/d/sfragkoul/Synth_Data/Synthesizers/NEAT/testing/TP53/TP53.fasta ./results/5/5_golden.bam > ./results/5/5_report.tsv

echo Variant Calling 
#Printing vardict commands
vardict-java -G /mnt/d/sfragkoul/Synth_Data/Synthesizers/NEAT/testing/TP53/TP53.fasta -f 0.0001 -N Merged_auto -b ./results/Merged_auto.sorted.uniq.rg.bam -R hg38_knownGene_ENST00000610292.4:0-19080 > ./results/Merged_auto_VarDict.txt
var2vcf_valid.pl -A ./results/Merged_auto_VarDict.txt >  ./results/Merged_auto_VarDict_tmp.vcf
rm ./results/Merged_auto_VarDict.txt
bcftools reheader --fai /mnt/d/sfragkoul/Synth_Data/Synthesizers/NEAT/testing/TP53/TP53.fasta.fai -o ./results/Merged_auto_VarDict.vcf ./results/Merged_auto_VarDict_tmp.vcf
rm ./results/Merged_auto_VarDict_tmp.vcf
bcftools norm ./results/Merged_auto_VarDict.vcf --output ./results/Merged_auto_VarDict_norm.vcf --output-type v -m "-"
rm ./results/Merged_auto_VarDict.vcf
 
