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


set -euo pipefail

#
# little helper to strip any stray CR (\r) from variables
#
trim_cr() {
  echo "${1//$'\r'/}"
}

#
# pull all keys in parameters.yaml into bash variables
#
function parse_yaml {
   local prefix=$2
   local s='[[:space:]]*' w='[a-zA-Z0-9_]*' \
         fs=$(echo @| tr @ '\034')
   sed -ne "s|^\($s\):|\1|" \
        -e "s|^\($s\)\($w\)$s:$s[\"']\(.*\)[\"']$s\$|\
\1$fs\2$fs\3|p" \
        -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|\
\1$fs\2$fs\3|p"  $1 |
   awk -F$fs '{
      indent = length($1)/2;
      vname[indent] = $2;
      for (i in vname) if (i > indent) delete vname[i];
      if (length($3) > 0) {
         vn=""; for (i=0; i<indent; i++) vn=vn vname[i] "_";
         printf("%s%s%s=\"%s\"\n",
                "'$prefix'", vn, $2, $3);
      }
   }'
}


# make sure folder has no stray CR
folder=$(trim_cr "$folder")

#
# synth-reads generation
#
gen_script="gen_${folder}.sh"
echo ">>> Generating $gen_script from synth_generation_template.sh"
bash bash/synth_generation_template.sh > "$gen_script"
chmod +x "$gen_script"

echo ">>> Running $gen_script"
bash "./$gen_script"

#
# variant-calling
#
vc_script="vc_${folder}.sh"
echo ">>> Generating $vc_script from variant_calling_template.sh"
bash bash/variant_calling_template.sh > "$vc_script"
chmod +x "$vc_script"

echo ">>> Running $vc_script"
bash "./$vc_script"


echo ">>> Running R analysis for each caller"

IFS=, 
read -ra callers_arr <<< "$callers"

for caller in "${callers_arr[@]}"; do
  echo ">>> S4BR analysis for $caller"
  Rscript R/S4BR.R \
    -c "$caller" \
    -r "$runs" \
    -w "${working_directory}/${folder}" \
    -m "$output_bam_merged"

  echo ">>> S4BR plotting for $caller"
  Rscript R/S4BR_plot.R \
    -c "$caller" \
    -w "${working_directory}/${folder}" \
    -m "$output_bam_merged"
done

echo "All done!"
