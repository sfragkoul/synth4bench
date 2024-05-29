
function parse_yaml {
   local prefix=$2
   local s='[[:space:]]*' w='[a-zA-Z0-9_]*' fs=$(echo @|tr @ '\034')
   sed -ne "s|^\($s\):|\1|" \
        -e "s|^\($s\)\($w\)$s:$s[\"']\(.*\)[\"']$s\$|\1$fs\2$fs\3|p" \
        -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|\1$fs\2$fs\3|p"  $1 |
   awk -F$fs '{
      indent = length($1)/2;
      vname[indent] = $2;
      for (i in vname) {if (i > indent) {delete vname[i]}}
      if (length($3) > 0) {
         vn=""; for (i=0; i<indent; i++) {vn=(vn)(vname[i])("_")}
         printf("%s%s%s=\"%s\"\n", "'$prefix'",vn, $2, $3);
      }
   }'
}

eval $(parse_yaml parameters.yaml)

printf "\n"

# part 1 -------------------------------
printf "echo Common Pre-process Steps \n"
printf "samtools sort -o "
printf "${working_directory}/${folder}/${output_bam_merged}.sorted.bam "
printf "${working_directory}/${folder}/${output_bam_merged}.bam"
printf "\n"

printf "samtools view -h -F 0x904 -b "
printf "${working_directory}/${folder}/${output_bam_merged}.sorted.bam > "
printf "${working_directory}/${folder}/${output_bam_merged}.sorted.uniq.bam"
printf "\n"

printf "samtools flagstat "
printf "${working_directory}/${folder}/${output_bam_merged}.sorted.uniq.bam > "
printf "${working_directory}/${folder}/${output_bam_merged}.align.stats.txt"
printf "\n"

printf "samtools addreplacerg -r '@RG\x$(printf %x 92)tID:${output_bam_merged}\x$(printf %x 92)tSM:${output_bam_merged}' "
printf "${working_directory}/${folder}/${output_bam_merged}.sorted.uniq.bam -o "
printf "${working_directory}/${folder}/${output_bam_merged}.sorted.uniq.rg.bam "
printf "\n"

printf "samtools depth "
printf "${working_directory}/${folder}/${output_bam_merged}.sorted.uniq.rg.bam -o "
printf "${working_directory}/${folder}/${output_bam_merged}.depth.txt "
printf "\n"

printf "samtools index "
printf "${working_directory}/${folder}/${output_bam_merged}.sorted.uniq.rg.bam"
printf "\n"

for i in $( seq 1 1 $runs ); do
	printf "samtools index "
	printf "${working_directory}/${folder}/${i}/${i}_golden.bam"
done
printf "\n \n"

# part 2 -------------------------------
#bam-readcount reporting
printf "echo bam-readcount reporting \n"
printf "${path_to_bam_readcount}/bam-readcount/build/bin/bam-readcount -w 0 -f "
printf "${path_to_reference} "
printf "${working_directory}/${folder}/${output_bam_merged}.sorted.uniq.rg.bam > "
printf "${working_directory}/${folder}/${output_bam_merged}_report.tsv"
printf "\n"

for i in $( seq 1 1 $runs ); do
	printf "${path_to_bam_readcount}/bam-readcount/build/bin/bam-readcount -w 0 -f "
	printf "${path_to_reference} "
	printf "${working_directory}/${folder}/${i}/${i}_golden.bam > "
	printf "${working_directory}/${folder}/${i}/${i}_report.tsv"
	printf "\n"
done
printf "\n"

# part 3 -------------------------------
#call the variant callers

printf "echo Variant Calling \n"

#if ${caller}== GATK-Mutect2

printf "gatk Mutect2 --reference "
printf "${path_to_reference} --input "
printf "${working_directory}/${folder}/${output_bam_merged}.sorted.uniq.rg.bam --tumor-sample "
printf "${output_bam_merged} --output "
printf "${working_directory}/${folder}/${output_bam_merged}_GATK.vcf > "
printf "${working_directory}/${folder}/${output_bam_merged}.Mutect.out 2>&1"
printf "\n"
printf "bcftools norm ${working_directory}/${folder}/${output_bam_merged}_GATK.vcf --output "
printf "${working_directory}/${folder}/${output_bam_merged}_GATK_norm.vcf --output-type "
printf "v -m \"-\""
printf "rm ${working_directory}/${folder}/${output_bam_merged}_GATK.vcf"
printf "\n \n"

#if ${caller}== freebayes

printf "freebayes --fasta-reference "
printf "${path_to_reference} --bam "
printf "${working_directory}/${folder}/${output_bam_merged}.sorted.uniq.rg.bam > "
printf "${working_directory}/${folder}/freebayes.vcf"
printf "\n"
printf "bcftools reheader --fai "
printf "${path_to_reference}.fai -o "
printf "${working_directory}/${folder}/${output_bam_merged}_freebayes.vcf "
printf "${working_directory}/${folder}/freebayes.vcf"
printf "\n"
printf "rm ${working_directory}/${folder}/freebayes.vcf"
printf "\n"
printf "bcftools norm ${working_directory}/${folder}/${output_bam_merged}_freebayes.vcf --output "
printf "${working_directory}/${folder}/${output_bam_merged}_freebayes_norm.vcf --output-type "
printf "v -m \"-\""
printf "\n"
printf "rm ${working_directory}/${folder}/${output_bam_merged}_freebayes.vcf"
printf "\n \n"


#if ${caller}== lofreq

printf "lofreq indelqual --dindel -f "
printf "${path_to_reference} -o "
printf "${working_directory}/${folder}/${output_bam_merged}_indels.sorted.uniq.rg.bam "
printf "${working_directory}/${folder}/${output_bam_merged}.sorted.uniq.rg.bam "
printf "\n"
printf "lofreq call -f "
printf "${path_to_reference} --call-indels -o "
printf "${working_directory}/${folder}/Lofreq.vcf "
printf "${working_directory}/${folder}/${output_bam_merged}_indels.sorted.uniq.rg.bam "
printf "\n"
printf "bcftools reheader --fai "
printf "${path_to_reference}.fai -o "
printf "${working_directory}/${folder}/${output_bam_merged}_Lofreq_norm.vcf "
printf "${working_directory}/${folder}/Lofreq.vcf"
printf "\n"
printf "rm ${working_directory}/${folder}/Lofreq.vcf"
printf "\n \n"


#if ${caller}== varcan
#we need Varscan2VCF path

printf "samtools mpileup -f "
printf "${path_to_reference} "
printf "${working_directory}/${folder}/${output_bam_merged}.sorted.uniq.rg.bam -a -o "
printf "${working_directory}/${folder}/${output_bam_merged}.sorted.uniq.rg.mpileup"
printf "\n"
printf "varscan pileup2cns "
printf "${working_directory}/${folder}/${output_bam_merged}.sorted.uniq.rg.mpileup > "
printf "${working_directory}/${folder}/VarScan.tsv"
printf "\n"
printf "python path/to/Varscan2VCF/vscan_pileup2cns2vcf.py "
printf "${working_directory}/${folder}/VarScan.tsv > ${working_directory}/${folder}/${output_bam_merged}_VarScan.vcf"
printf "\n"
printf "rm ${working_directory}/${folder}/VarScan.tsv"
printf "\n"
printf "bcftools norm ${working_directory}/${folder}/${output_bam_merged}_VarScan.vcf --output "
printf "${working_directory}/${folder}/${output_bam_merged}_VarScan_norm.vcf --output-type "
printf "v -m \"-\""
printf "\n"
printf "rm ${working_directory}/${folder}/${output_bam_merged}_VarScan.vcf"
printf "\n \n"

#if ${caller}== vardict
#need to read this -R hg38_knownGene_ENST00000610292.4:0-19080

printf "vardict-java -G "
printf "${path_to_reference} -f 0.0001 -N ${output_bam_merged} -b "
printf "${working_directory}/${folder}/${output_bam_merged}.sorted.uniq.rg.bam -R "
printf "hg38_knownGene_ENST00000610292.4:0-19080 > ${working_directory}/${folder}/${output_bam_merged}_VarDict.txt"
printf "\n"
printf "var2vcf_valid.pl -A "
printf "${working_directory}/${folder}/${output_bam_merged}_VarDict.txt >  "
printf "${working_directory}/${folder}/${output_bam_merged}_VarDict.vcf"
printf "\n"
printf "rm ${working_directory}/${folder}/${output_bam_merged}_VarDict.txt"
printf "\n"
printf "bcftools reheader --fai "
printf "${path_to_reference}.fai -o "
printf "${working_directory}/${folder}/${output_bam_merged}_VarDict.vcf "
printf "${working_directory}/${folder}/VarDict.vcf"
printf "\n"
printf "rm ${working_directory}/${folder}/VarDict.vcf"
printf "\n"
printf "bcftools norm ${working_directory}/${folder}/${output_bam_merged}_VarDict.vcf --output "
printf "${working_directory}/${folder}/${output_bam_merged}_VarDict_norm.vcf --output-type "
printf "v -m \"-\""
printf "\n"
printf "rm ${working_directory}/${folder}/${output_bam_merged}_VarDict.vcf"
printf "\n \n"

