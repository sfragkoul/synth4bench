
# Synth4bench Documentation
This is the documentation of the scripts.

## Contents
- [Main scripts](https://github.com/sfragkoul/synth4bench/blob/main/docs/Documentation.md#main-scripts)
- [Extra scripts](https://github.com/sfragkoul/synth4bench/blob/main/docs/Documentation.md#extra-scripts)


Here follows the list of all scripts and their description.

## Main scripts

<p align="center"> 
<img src="https://github.com/sfragkoul/synth4bench/blob/main/images/merge.png" alt="Scehmatic for merging process of bam files" style="height: 200px; width: 650px;"/>
</p>

### Set parameters

`parameters.yaml` - File in which the user sets all the parameters of the workflow.


### Synthetic Data Generation

`synth_generation_template.sh` - This bash script prints the commands of the workflow for the synthetic data generation with parameters spesified by the user.

- Input files: fasta reference file
							 
- Output files: fastq files with pair end reads, ground truth "golden" bam file and bai index file, ground truth "golden" vcf file, Merged bam and vcf files for the ground truth.

### Variant Calling 

`variant_calling_template.sh` - This bash script prints the commands of the workflow for the variant calling process with parameters spesified by the user. 

- Input files:  ground truth Merged bam file
  			 
- Output files: vcf files from variant callers, bam-readcount tsv reports, files with stats for the ground truth Merged bam file


### Downstream Analysis

`S4BR.R` - A script, written in R, that calls the appropriate functions to perform the comparison between the ground truth and the caller.

- Input files:  ground truth files, caller vcf file
							 
- Output files: a tsv file with the comparison between the ground truth and the caller.

### Visualization

`S4BR_plot.R` - A script, written in R, that calls the appropriate functions to make visualizations to illustrate the comparison between the ground truth and the caller.

- Input files:  a tsv file with the comparison between the ground truth and the caller, ground truth vcf file, caller vcf file
							 
- Output files: multi planel figure and a Venn plot
 

`paper_plots.R` -  A script, written in R, that outputs the figures for the manuscript.

## Extra scripts

#### VarScan scripts

For the case of VarScan an extra step was required to convert its output to the standard VCF format. The script `vscan_pileup2cns2vcf.py` can be found [here](https://github.com/sfragkoul/Varscan2VCF).

#### Statistical Analysis scripts

`S4BR_read_pos.R` - A script to report all ground truth variants in each chromosomal position.

- Input files:  individual "golden" bam files produced by NEAT
							 
- Output files: reports in tsv format with the variants in each chromosomal position for each individual "golden" bam file

`S4BR_read_pos.R` - A script to report variants that were either  detected or not detected by each caller. The read positions were divided in bins to study possible correlated trends.

- Input files:  tsv file with reported variants from caller, reports in tsv format with the variants in each chromosomal position for each individual "golden" bam file 
							 
- Output files: reports with variants that were either detected or not detected by each caller. The read positions were divided in bins.

`stats_1_Prediction_Coverage_Length.R` - This script calculates the Kandall's tau coefficient of the Coverage and Read length with the modified Accuracy (mAC), False negative (FN) and False positive (FP) rates.

- Input files:  an xlsx file with a column containing the callers (named Caller), a column specifying the AC, FN, FP rate (named Variants), and a column for each coverage (in a Sheet named Coverage) or read length (in a Sheet named Read_length) with the respective rates. 
							 
- Output files: A xlsx file with the results in two sheets: Coverage and Read_length

`stats_2_Statistical_Analysis.Rmd` - This script runs the statistical analysis for a given dataset.

- Input files:  A csv file with variants that were either detected or not detected by each caller.
							 
- Output files: A word document with the results of the fatalistically analyses.

`stats_3_ROC_Curves.R` - This script produces a ROC curve for a given dataset.

- Input files:  A csv file with variants that were either detected or not detected by each caller.
							 
- Output files: A ROC curve.

`stats_4_ROC_Curves_Merged.R` - This script produces a merged figure with the ROC curves for a given caller, each curve corresponding to a specific Coverage and Read length.

- Input files:  A folder with the files containing all the datasets and the 5.3_ROC_Curves.R" file.
							 
- Output files: A figure with the ROC curves for each caller.

  </div>