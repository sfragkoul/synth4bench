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

printf "#Printing Mutect2 commands"
printf "\n"

printf "gatk Mutect2"
printf " --reference ${path_to_reference}"													# what is this parameter
printf " --input ${working_directory}/${folder}/${output_bam_merged}.sorted.uniq.rg.bam"	# what is this parameter
printf " --tumor-sample ${output_bam_merged}"
printf " --output ${working_directory}/${folder}/${output_bam_merged}_Mutect2.vcf"
printf " > ${working_directory}/${folder}/${output_bam_merged}.Mutect2.out 2>&1"
printf "\n"

printf "bcftools norm ${working_directory}/${folder}/${output_bam_merged}_Mutect2.vcf"
printf " --output ${working_directory}/${folder}/${output_bam_merged}_Mutect2_norm.vcf"
printf " --output-type v"
printf " -m \"-\""
printf "\n"

printf "rm ${working_directory}/${folder}/${output_bam_merged}_Mutect2.vcf"
printf "\n \n"