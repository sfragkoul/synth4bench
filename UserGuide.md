
#Synth4bench User Guide


## Contents



## Overview

<div align='justify'> Somatic variant calling algorithms are widely used to detect genomic alterations associated with cancer. Evaluating the performance of these algorithms can be challenging due to the lack of high-quality ground truth datasets. To address this issue, we developed a synthetic genomics data generation framework for benchmarking tumor-only somatic variant calling algorithms. We generated synthetic datasets based on <i>TP53</i> gene, using the NEATv3.3 (NExt-generation sequencing Analysis Toolkit version 3) simulator. Subsequently, we thoroughly evaluated the performance of variant calling algorithms using <strong> GATK-Mutect2, Freebayes, VarDict, VarScan and LoFreq </strong> on these datasets, and compared the results to the “golden” files produced by NEATv3.3 containing the actual variants. Synthetic datasets provide an excellent ground truth for studying the performance and behaviour of somatic variant calling algorithms, thus enabling researchers to evaluate and improve the accuracy of these algorithms for cancer genomics applications.</div>



### Key Features

## Synth4bench

### Synthetic Data Generation

### Variant Calling


### Visualization


## Installation
To set up the environment and dependencies for synth4bench, follow these steps:

1. Create the Conda Environment
To create the conda environment that was used for the analysis run the following command in your terminal **(Bash)**:

       conda env create -f environment.yml

This will install all the required dependencies specified in the environment.yml file.

2. Activate the Conda Environment
After creating the environment, activate it with the following command in your terminal **(Bash)**:

        conda activate synth4bench

3. Install NEATv3.3
Download the version [v3.3](https://github.com/ncsa/NEAT/releases/tag/3.3).
Once downloaded, in order to call the main script and view available options, run in your terminal **(Bash)**:

       python gen_reads.py --help
  
For any further details, please see the README.md file included in the downloaded NEATv3.3 version.

4. Install bam-readcount
Follow the installation [instructions](https://github.com/genome/bam-readcount/tree/master?tab=readme-ov-file#build) from their official documentation.
After installation, verify that it has being installed properly with the following command in your terminal **(Bash)**:

       build/bin/bam-readcount --help

5. Install R Package Dependencies
Install the required R packages by running the following command in your **R console**:

       install.packages(c("stringr", "data.table", "vcfR", "ggplot2", "ggvenn", "ggforce", "ggsci", "patchwork"))

6. Download VarScan Extra Script
The extra script vscan_pileup2cns2vcf.py for VarScan can be found [here](https://github.com/sfragkoul/Varscan2VCF).

## Using synth4bench


## Requirements
