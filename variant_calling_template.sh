#!/bin/bash

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
	printf "${working_directory}/${folder}/${i}/${i}_golden.bam \n"
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

if [ "${caller}" == "Mutect2" ]; then

	bash mutect2_template.sh
	
elif [ "${caller}" == "Freebayes" ]; then

	bash freebayes_template.sh
	
elif [ "${caller}" == "LoFreq" ]; then

	bash lofreq_template.sh
	
elif [ "${caller}" == "VarScan" ]; then

	bash varscan_template.sh
	
elif [ "${caller}" == "VarDict" ]; then

	bash vardict_template.sh
	
else

	echo "Invalid caller option. Please provide one of the following: Mutect2, Freebayes, LoFreq, VarScan, VarDict"
	
fi
