
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

IFS=$IFS,
rngs=()
read -a rngs <<< "${rng}"

printf "mkdir ${working_directory}/${folder}\n"
printf "mkdir ${working_directory}/Plots\n"

printf "\n"

# part 1 ----------------------------------
printf "echo Generating Synthetic files\n"

for i in $( seq 1 1 $runs ); do
	printf "echo Starting run ${i}\n"
	printf "mkdir ${working_directory}/${folder}/${i}\n"
	printf "python ${path_to_neat}/gen_reads.py"
	printf " -r ${path_to_reference}"
	printf " --rng ${rngs[${i} - 1]}"
	printf " -m ${mutation_rate}"
	printf " -r ${read_length}"
	printf " -c ${coverage}"
	printf " -o ${working_directory}/${folder}/${i}/${i}"
	printf " --pe ${fragment_stats_1} ${fragment_stats_2}"
	printf " --bam --vcf\n"
	printf "\n"
done

# part 2 -------------------------------
#Merged file
printf "echo Merging bam files\n"

printf "samtools merge"
printf " ${working_directory}/${output_bam_merged}.bam"

for i in $( seq 1 1 $runs ); do
	
	printf " ${working_directory}/${folder}/${i}/${i}_golden.bam"
	
done

printf "\n"

for i in $( seq 1 1 $runs ); do
	printf "bcftools index"
	printf " ${working_directory}/${folder}/${i}/${i}_golden.vcf.gz \n"
done

printf "bcftools merge"


for i in $( seq 1 1 $runs ); do
	printf " ${working_directory}/${folder}/${i}/${i}_golden.vcf.gz"
done

printf " > ${working_directory}/${output_bam_merged}_ground_truth.vcf\n"	
	
	
printf "bcftools norm"
printf " ${working_directory}/${output_bam_merged}_ground_truth.vcf"
printf " --output ${working_directory}/${output_bam_merged}_ground_truth_norm.vcf"	
printf " --output-type v -m \"-\""
printf "\n"
