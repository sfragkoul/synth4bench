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

printf "Printing vardict commands"
printf "\n"

printf "vardict-java"
printf " -G ${path_to_reference}"
printf " -f 0.0001"
printf " -N ${output_bam_merged}"
printf " -b ${working_directory}/${folder}/${output_bam_merged}.sorted.uniq.rg.bam"
printf " -R ${vardict_region_of_interest}"
printf " > ${working_directory}/${folder}/${output_bam_merged}_VarDict.txt"
printf "\n"

printf "var2vcf_valid.pl"
printf " -A ${working_directory}/${folder}/${output_bam_merged}_VarDict.txt"
printf " >  ${working_directory}/${folder}/${output_bam_merged}_VarDict.vcf"
printf "\n"

printf "rm ${working_directory}/${folder}/${output_bam_merged}_VarDict.txt"
printf "\n"

printf "bcftools reheader"
printf " --fai ${path_to_reference}.fai"
printf " -o ${working_directory}/${folder}/${output_bam_merged}_VarDict.vcf "
printf "${working_directory}/${folder}/VarDict.vcf"
printf "\n"

printf "rm ${working_directory}/${folder}/VarDict.vcf"
printf "\n"

printf "bcftools norm ${working_directory}/${folder}/${output_bam_merged}_VarDict.vcf"
printf " --output ${working_directory}/${folder}/${output_bam_merged}_VarDict_norm.vcf"
printf " --output-type v"
printf " -m \"-\""
printf "\n"

printf "rm ${working_directory}/${folder}/${output_bam_merged}_VarDict.vcf"
printf "\n \n"
