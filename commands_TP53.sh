bwa index TP53.fasta
samtools faidx TP53.fasta
java -jar ../picard.jar CreateSequenceDictionary R=TP53.fasta O=TP53.dict
python gen_reads.py -r testing/TP53/TP53.fasta -R 150 -c 5000 -M 0.1 -o testing/TP53/TP53 --pe 300 30 --bam --vcf
samtools sort -o TP53_golden_edit.sorted.bam TP53_golden_edit.bam 
samtools view -h -F 0x904 -b TP53_golden_edit.sorted.bam > TP53_golden_edit.sorted.uniq.bam  
samtools flagstat TP53_golden_edit.sorted.uniq.bam > TP53_golden_edit.align.stats.txt  
samtools addreplacerg -r '@RG\tID:1\tSM:1' TP53_golden_edit.sorted.uniq.bam -o TP53_golden_edit.sorted.uniq.rg.bam
samtools depth TP53_golden_edit.sorted.uniq.rg.bam -o TP53_golden_edit.depth.txt
samtools index TP53_golden_edit.sorted.uniq.rg.bam
java -jar ../gatk-4.3.0.0.jar Mutect2 --reference TP53.fasta --input TP53_golden_edit.sorted.uniq.rg.bam --tumor-sample TP53_golden_edit  --output TP53_golden_edit_GATK_variants.vcf > TP53_golden_edit.Mutect2.out 2>&1
