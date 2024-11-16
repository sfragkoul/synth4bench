
# Synth4bench User Guide
This is the guide to help users execute our pipeline.

## Contents

- [Overview](https://github.com/sfragkoul/synth4bench/blob/main/UserGuide.md#overview)
- [Data Download](https://github.com/sfragkoul/synth4bench/blob/main/UserGuide.md#data-download)
- [Installation](https://github.com/sfragkoul/synth4bench/blob/main/UserGuide.md#installation)
- [Key Features](https://github.com/sfragkoul/synth4bench/blob/main/UserGuide.md#key-features)
- [Executing synth4bench](https://github.com/sfragkoul/synth4bench/blob/main/UserGuide.md#executing-synth4bench)
- [Requirements](https://github.com/sfragkoul/synth4bench/blob/main/UserGuide.md#requirements)


## Overview

<div align='justify'> Somatic variant calling algorithms are widely used to detect genomic alterations associated with cancer. Evaluating the performance of these algorithms can be challenging due to the lack of high-quality ground truth datasets. To address this issue, we developed a synthetic genomics data generation framework for benchmarking tumor-only somatic variant calling algorithms. We generated synthetic datasets based on <i>TP53</i> gene, using the NEATv3.3 (NExt-generation sequencing Analysis Toolkit version 3) simulator. Subsequently, we thoroughly evaluated the performance of variant calling algorithms using <strong> GATK-Mutect2, Freebayes, VarDict, VarScan and LoFreq </strong> on these datasets, and compared the results to the “golden” files produced by NEATv3.3 containing the actual variants. Synthetic datasets provide an excellent ground truth for studying the performance and behaviour of somatic variant calling algorithms, thus enabling researchers to evaluate and improve the accuracy of these algorithms for cancer genomics applications.</div>


## Data Download

The datasets required for this analysis are available for open access on [Zenodo](https://zenodo.org/records/10683211). 

This repository contains 10 synthetic genomics datasets, generated specifically for **benchmarking somatic variant callers**. Each dataset was produced with [NEAT v3](https://github.com/ncsa/NEAT/releases/tag/3.3), based on the *TP53* gene of Homo sapiens, and provides valuable resources for analyzing the effects of various NGS parameters on **tumor-only somatic variant calling** algorithms.

### Data Overview

These datasets explore two primary variables to observe their impact on somatic variant calling:

   - **Coverage:** Five datasets vary in coverage while keeping read length constant (150 bp). Coverage levels include 300x, 700x, 1000x, 3000x, and 5000x.

   - **Read Length:** Another set of five datasets vary in read length while maintaining a coverage level of 1000x. Read lengths include 50 bp, 100 bp, 150 bp, 170 bp, 200 bp, and 300 bp.

Each dataset contains paired-end reads, providing compatibility with most standard analysis pipelines for synthetic genomics data.

### Available Files

There are two compressed folders for download:

   - **Reference folders**: reference.rar (25.7 kB)

This folders includes reference sequences necessary for aligning or comparing generated synthetic reads.

   - **Synthetic folders:** synth_datasets.rar (55.9 MB)

This folders includes all 10 synthetic datasets, with filenames indicating coverage and read length details. Each dataset is structured and named according to the parameters it explores.

### Required Files for Analysis

To successfully **run the analysis**, users need the reference files located in the reference.rar folder.

## Installation

**Create the Conda Environment**:
To create the conda environment that was used for the analysis run the following command in your terminal **(Bash)**:

    conda env create -f environment.yml

This will install all the required dependencies specified in the environment.yml file.

**Activate the Conda Environment:**
After creating the environment, activate it with the following command in your terminal **(Bash)**:

    conda activate synth4bench

**Install NEATv3.3:**
Download the version [v3.3](https://github.com/ncsa/NEAT/releases/tag/3.3).
Once downloaded, in order to call the main script and view available options, run in your terminal **(Bash)**:

    python gen_reads.py --help
  
For any further details, please see the README.md file included in the downloaded NEATv3.3 version.

**Install bam-readcount:**
Follow the installation [instructions](https://github.com/genome/bam-readcount/tree/master?tab=readme-ov-file#build) from their official documentation.
After installation, verify that it has being installed properly with the following command in your terminal **(Bash)**:

    build/bin/bam-readcount --help

**Install R Package Dependencies:**
Install the required R packages by running the following command in your **R console**:

    install.packages(c("stringr", "data.table", "vcfR", "ggplot2", "ggvenn", "ggforce", "ggsci", "patchwork"))

**Download VarScan Extra Script:**
The extra script vscan_pileup2cns2vcf.py for VarScan can be found [here](https://github.com/sfragkoul/Varscan2VCF).


## Key Features

Once the environment and dependencies have been successfully installed, follow these steps to execute the analysis:

**Prepare the Parameters:**
Before running the analysis, you need to configure the parameters. Start by editing the `parameters.yaml` file with your specific analysis settings. This file contains various options for input files, output directories, and other customizable parameters for the execution of the scripts. Make sure to review and adjust each section based on your project requirements.

**Generate the Execution Scripts:**
To create customized execution scripts for data synthesis and variant calling, run the following commands in your **terminal (Bash)**:

   - **Generate the Synthesis Script**: This will create the `synth_generation_template.sh` script with the parameters specified in the `parameters.yaml` file:

    bash synth_generation_template.sh > desired_name.sh

Replace `desired_name.sh` with the desired name for your generated script.

   - **Generate the Variant Calling Script**: Similarly, this will create the `variant_calling_template.sh` script with the parameters from the `parameters.yaml` file:

    bash variant_calling_template.sh > desired_name.sh

Again, replace `desired_name.sh` with the name you would like to assign to this generated script.

Once you have generated the scripts, you should review them to ensure the parameters were correctly applied.

**Run the Scripts**: After generating the customized scripts, run them one by one by executing the following commands in your **terminal (Bash)**:

   - **Run the Synthesis Script**:

    bash synth_generation_run.sh

This will display the available options and details for the `S4BR.R` script.
   
   - **Run the Variant Calling Script**:

    bash variant_calling_run.sh

This command will execute the variant calling step.

**Check R Script Parameters**:
To review and check the parameters and usage for the R scripts involved in your analysis, run the following commands in your **terminal (Bash)**:

   - **Check Parameters for S4BR R Script**:

    Rscript R/S4BR.R --help

This will display the available options and details for the `S4BR.R` script.
   
   - **Check Parameters for S4BR Plotting R Script**:

    Rscript R/S4BR_plot.R --help

This will show you the available options and help documentation for the `S4BR_plot.R` script.

**Review Output**:
After executing the scripts, check the output files located in the directories specified in the `parameters.yaml` file. These will include any generated results from the **synthesis and variant calling** steps. You can now proceed to analyze or visualize the data as needed.


## Executing synth4bench

This guide will walk you through using Synth4bench with a practical example, from downloading data to executing the analysis and reviewing results. We assume you have the necessary environment and dependencies set up as outlined above.

**Set Up Your Analysis Folder**

   - **Create a New Analysis Folder:** Start by creating an analysis folder and navigating into it in your **terminal (Bash)**:

    mkdir synth4bench_analysis
    cd synth4bench_analysis

   - **Download the Data:** Download the reference files from Zenodo and extract the files into your analysis folder:
    
    wget https://zenodo.org/record/10683211/files/reference.rar
    unrar x reference.rar

**Configure Parameters:**

Edit `parameters.yaml` to set your working directory and specify paths to the output directory and reference files. Below is an example `parameters.yaml`:

    ### synth_generation_template.sh parameters ###

    # path to working directory
    working_directory: "/path/to/synth4bench_analysis"
    
    # folder for results
    folder: results 
    
    # path to NEAT scripts
    path_to_neat:  /path/to/NEAT-3.3
    
    # number of individual synth batches 
    runs: 5
    
    # seeds for each individual synth batch
    rng: 213,214,215,217,218
    
    # mutation rate parameter as taken from NEAT
    mutation_rate: 0.1
    
    # read length parameter as taken from NEAT
    read_length: 150
    
    # coverage parameter as taken from NEAT
    coverage: 100
    
    # fragment parameters as taken from NEAT
    fragment_stats:
     1: 300
     2: 30
     
    # path to reference file
    path_to_reference: /path/to/reference/TP53.fasta
    
    # name of final merged bam file
    output_bam_merged: Merged
   
   
    ### variant_calling_template.sh parameters ###
   
    # path to bam readcount scripts
    path_to_bam_readcount: /path/to/build/bin/bam_readcount
    
    # somatic variant caller name
    caller: VarScan
    
    # path to varscan_scripts_path extra scripts
    varscan_scripts_path: /path/to/varscan_extra_scripts
    
    # genomic region of interest as taken from VarScan
    vardict_region_of_interest: hg38_knownGene_ENST00000610292.4:0-19080

If you want to run the analysis at the synth4bench_analysis folder and you have downloaded all the files at this directory, you can change /path/to with a dot.

### **Run the Execution Scripts**

   - **Generate the Synthesis Script**: Run the following to create a customized synthesis script based on the parameters in `parameters.yaml`:

    bash synth_generation_template.sh > synth_generation_run.sh

   - **Generate the Variant Calling Script**: Create the variant calling script similarly:

    bash variant_calling_template.sh > variant_calling_run.sh

   - **Run the Synthesis Script**: Execute the synthesis script to generate synthetic reads:

    bash synth_generation_run.sh

This step will output synthetic FASTQ  and bam files ready for the variant calling analysis.
     
   - **Run the Variant Calling Script**: Execute the variant calling script to analyze the generated synthetic data:

    bash variant_calling_run.sh

This script will output a VCF file in the specified output_directory, containing the variant calls for the synthetic data.


### **Analyze the Results Using R Scripts**

   - Check **Parameters** for S4BR R Script: 

    Rscript R/S4BR.R --help
    
   
| **Options**  | **Description**   |
|:-----|:-|
|-v, --vcf_path | Directory path where VCF files are located.|
|-c, --caller | Choose caller name (Freebayes, Mutect2, LoFreq, VarDict, VarScan)| 
|-r, --runs |  Number of individual runs to produce synthetic data which will then be combined to form the final Merged ground truth file.|
|-w, --working_directory| Path of working directory.|
|-m, --merged_file | Indicate the name given to the final merged ground truth file.|
|-h, --help | Show this help message and exit.|


   - Run **S4BR.R** to perform the analysis:

    Rscript R/S4BR.R -v ./results -c VarScan -r 5 -w ./results -m test
    
### **Visualization Using R Scripts**
    
   - Check **Parameters** for S4BR_plot R Script: 

    Rscript R/S4BR_plot.R --help
    
| **Options**  | **Description**   |
|:-----|:-|
|-t, --gt_comparison | Directory path where Ground Truth vs Caller file tsv file is located.|
|-v, --vcf_path | Directory path where VCF files are located.|
|-g, --gt_path| Directory path where ground truth vcf file is located.|
|-c, --caller | Choose caller name (Freebayes, Mutect2, LoFreq, VarDict, VarScan)| 
|-w, --working_directory| Path of working directory.|
|-m, --merged_file | Indicate the name given to the final merged ground truth file.|
|-h, --help | Show this help message and exit.|

                
   - Run **S4BR_plot.R** to generate visualizations:
   
    Rscript R/S4BR_plot.R -t ./results -v ./results -g ./results -c VarScan -w . -m test

**Review the Output**
   
**Check Analysis Results:** Open the `analysis_output.tsv` file in the results directory. This file summarizes the variants found in the synthetic data.

**View Plots:** Navigate to the plots directory to review the generated visualizations. These plots provide insights into the variants called and their characteristics under different conditions.


## Requirements
