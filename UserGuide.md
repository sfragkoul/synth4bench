
#Synth4bench User Guide


## Contents



## Overview

<div align='justify'> Somatic variant calling algorithms are widely used to detect genomic alterations associated with cancer. Evaluating the performance of these algorithms can be challenging due to the lack of high-quality ground truth datasets. To address this issue, we developed a synthetic genomics data generation framework for benchmarking tumor-only somatic variant calling algorithms. We generated synthetic datasets based on <i>TP53</i> gene, using the NEATv3.3 (NExt-generation sequencing Analysis Toolkit version 3) simulator. Subsequently, we thoroughly evaluated the performance of variant calling algorithms using <strong> GATK-Mutect2, Freebayes, VarDict, VarScan and LoFreq </strong> on these datasets, and compared the results to the “golden” files produced by NEATv3.3 containing the actual variants. Synthetic datasets provide an excellent ground truth for studying the performance and behaviour of somatic variant calling algorithms, thus enabling researchers to evaluate and improve the accuracy of these algorithms for cancer genomics applications.</div>



### Key Features

## Synth4bench

### Synthetic Data Generation

### Variant Calling
### Visualization


## Data Download

The datasets required for this analysis are available for open access on [Zenodo](https://zenodo.org/records/10683211). 

This repository contains 10 synthetic genomics datasets, generated specifically for benchmarking somatic variant callers. Each dataset was produced with NEAT v3, based on the TP53 gene of Homo sapiens, and provides valuable resources for analyzing the effects of various NGS parameters on **tumor-only somatic variant calling** algorithms.

## Installation

1. **Create the Conda Environment**:
To create the conda environment that was used for the analysis run the following command in your terminal **(Bash)**:

       conda env create -f environment.yml

This will install all the required dependencies specified in the environment.yml file.

2. **Activate the Conda Environment:**
After creating the environment, activate it with the following command in your terminal **(Bash)**:

       conda activate synth4bench

3. **Install NEATv3.3:**
Download the version [v3.3](https://github.com/ncsa/NEAT/releases/tag/3.3).
Once downloaded, in order to call the main script and view available options, run in your terminal **(Bash)**:

       python gen_reads.py --help
  
For any further details, please see the README.md file included in the downloaded NEATv3.3 version.

4. **Install bam-readcount:**
Follow the installation [instructions](https://github.com/genome/bam-readcount/tree/master?tab=readme-ov-file#build) from their official documentation.
After installation, verify that it has being installed properly with the following command in your terminal **(Bash)**:

       build/bin/bam-readcount --help

5. **Install R Package Dependencies:**
Install the required R packages by running the following command in your **R console**:

       install.packages(c("stringr", "data.table", "vcfR", "ggplot2", "ggvenn", "ggforce", "ggsci", "patchwork"))

6. **Download VarScan Extra Script:**
The extra script vscan_pileup2cns2vcf.py for VarScan can be found [here](https://github.com/sfragkoul/Varscan2VCF).


## Execution

Once the environment and dependencies have been successfully installed, follow these steps to execute the analysis:

1. **Prepare the Parameters:**
Before running the analysis, you need to configure the parameters. Start by editing the `parameters.yaml` file with your specific analysis settings. This file contains various options for input files, output directories, and other customizable parameters for the execution of the scripts. Make sure to review and adjust each section based on your project requirements.

2. **Generate the Execution Scripts:**
To create customized execution scripts for data synthesis and variant calling, run the following commands in your **terminal (Bash)**:

   - **Generate the Synthesis Script**: This will create the `synth_generation_template.sh` script with the parameters specified in the `parameters.yaml` file:

         bash synth_generation_template.sh > desired_name.sh

     Replace `desired_name.sh` with the desired name for your generated script.

    - **Generate the Variant Calling Script**: Similarly, this will create the `variant_calling_template.sh` script with the parameters from the `parameters.yaml` file:

          bash variant_calling_template.sh > desired_name.sh

       Again, replace `desired_name.sh` with the name you would like to assign to this generated script.

Once you have generated the scripts, you should review them to ensure the parameters were correctly applied.

3. **Run the Scripts**: After generating the customized scripts, run them one by one by executing the following commands in your **terminal (Bash)**:

   - **Run the Synthesis Script**:

         bash desired_name.sh

    This will display the available options and details for the `S4BR.R` script.
   
   - **Check Parameters for S4BR Plotting R Script**:

         bash desired_name.sh

    This command will execute the variant calling step.

4. **Check R Script Parameters**:
To review and check the parameters and usage for the R scripts involved in your analysis, run the following commands in your **terminal (Bash)**:

   - **Check Parameters for S4BR R Script**:

         Rscript R/S4BR.R --help

    This will display the available options and details for the `S4BR.R` script.
   
   - **Check Parameters for S4BR Plotting R Script**:

         Rscript R/S4BR_plot.R --help

    This will show you the available options and help documentation for the `S4BR_plot.R` script.

5. **Review Output**:
After executing the scripts, check the output files located in the directories specified in the `parameters.yaml` file. These will include any generated results from the **synthesis and variant calling** steps. You can now proceed to analyze or visualize the data as needed.






## Using synth4bench


## Requirements
