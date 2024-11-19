<p align="center"> 
<img src="https://github.com/sfragkoul/synth4bench/blob/main/images/synth4bench_logo_no_bg.png" alt="synth4bench logo" style="center; height: 80px; width:700px;"/>
</p>

## Abstract
<div align='justify'> Somatic variant calling algorithms are widely used to detect genomic alterations associated with cancer. Evaluating the performance of these algorithms can be challenging due to the lack of high-quality ground truth datasets. To address this issue, we developed a synthetic genomics data generation framework for benchmarking tumor-only somatic variant calling algorithms. We generated synthetic datasets based on <i>TP53</i> gene, using the NEATv3.3 (NExt-generation sequencing Analysis Toolkit version 3) simulator. Subsequently, we thoroughly evaluated the performance of variant calling algorithms using <strong> GATK-Mutect2, Freebayes, VarDict, VarScan and LoFreq </strong> on these datasets, and compared the results to the ground truth files produced by NEATv3.3 containing the real variants. Synthetic datasets provide an excellent ground truth for studying the performance and behaviour of somatic variant calling algorithms, thus enabling researchers to evaluate and improve the accuracy of these algorithms for cancer genomics applications.</div>

## Table of Contents

- [Abstract](https://github.com/sfragkoul/synth4bench/tree/main#abstract)
- [Motivation](https://github.com/sfragkoul/synth4bench/tree/main#motivation)
- [Description of Framework](https://github.com/sfragkoul/synth4bench/tree/main#description-of-framework)
- [Installation](https://github.com/sfragkoul/synth4bench/tree/main#installation)
- [Data Download](https://github.com/sfragkoul/synth4bench/tree/main#data-download)
- [Execution](https://github.com/sfragkoul/synth4bench/tree/main#execution)
- [Contribute](https://github.com/sfragkoul/synth4bench/tree/main#contribute)
- [Citation](https://github.com/sfragkoul/synth4bench/tree/main#citation)


## Motivation

<div align='justify'> Variant calling plays an important role in identifying genetic lesions. In the case of variants at low frequency (≤10%) identification becomes more difficult and the challenge that rises is the absence of a Ground Truth for reliable and consistent identification and benchmarking. </div>

## Description of Framework

<p align="center"> 
<img src="https://github.com/sfragkoul/synth4bench/blob/main/images/schematic.png" alt="synth4bench schematic" style="height: 100px; width:700px;"/>
</p>

<div align='justify'> Our framework focuses on addressing the challenge of variant calling, particularly for variants at low frequencies (≤10%). The main goal is to develop a reliable and consistent method for identifying genetic lesions, specifically in the context of cancer-associated genomic alterations. The absence of a ground truth, which refers to a reliable reference dataset with known variants, makes benchmarking and evaluating variant calling algorithms difficult. To overcome this challenge, the following steps are outlined in the framework:

1. Synth Data Generation: Synthetic genomics data is generated based on the *TP53* gene using the NEATv3.3 simulator in order to create synthetic datasets that mimic real cancer genome data. The "Ground Truth" is established by creating 10 individual datasets (each one of the same characteristics) containing Single Nucleotide Polymorphisms (SNPs) and Insertions/Deletions (INDELs). The genomic regions where variants accure with 100% Allele Frequency are chosen. The reason behind this choice is to avoid variants that are related to errors and products of noise. Then all these datasets are merged into one single file and the allele frequency is again measured at these genomic regions of interest.

2. Benchmarking Variant Callers: Somatic variant callers are evaluated using this synthetic Ground Truth dataset. GATK-Mutect2, Freebayes, VarDict, VarScan2 and LoFreq variant callers are assessed for their performance on our synthetic ground truth dataset. Their impact at low frequencies (≤10%) is explored, as these are particularly challenging to detect accurately.

The framework's overall aim is to provide a robust framework for evaluating the performance of tumor-only somatic variant calling algorithms by using synthetic datasets. By having a reliable ground truth, we can thoroughly test and improve the accuracy of variant calling algorithms for cancer genomics applications. This framework represents an essential step towards more precise and effective identification of genetic lesions associated with cancer and other diseases. </div>

## Data Download
All data are open and available in [Zenodo](https://zenodo.org/records/10683211).

## Installation
1.  To create the conda environment that was used for the analysis run `conda env create -f environment.yml` and to activate it run `conda activate synth4bench`.
2. To install NEATv3.3, download version [v3.3](https://github.com/ncsa/NEAT/releases/tag/3.3). To call the main script run the command `python gen_reads.py --help`. For any further info please see the README.md file from the downloaded files of version 3.3.
3. To install bam-readcount follow their [instructions](https://github.com/genome/bam-readcount/tree/master?tab=readme-ov-file#build) and then run `build/bin/bam-readcount --help` to see that it has being installed properly.
4. To install R packages dependencies run this command `install.packages(c("stringr", "data.table", "vcfR", "ggplot2", "ggvenn", "ggforce", "ggsci", "patchwork"))`.
5. The extra script `vscan_pileup2cns2vcf.py` for VarScan can be found [here](https://github.com/sfragkoul/Varscan2VCF).


## Execution
<div align='justify'> 
	
To run the bash scripts, fill in the parameters in the `parameters.yaml` file and then run:

- `bash synth_generation_template.sh > desired_name.sh` and
- `bash variant_calling_template.sh > desired_name.sh`

Run the following commands to check the paramaters of each R script:

- `Rscript R/S4BR.R --help`
- `Rscript R/S4BR_plot.R --help`

## Documentation

Here follows the list of all scripts and their description.

### Main scripts

<p align="center"> 
<img src="https://github.com/sfragkoul/synth4bench/blob/main/images/merge.png" alt="Scehmatic for merging process of bam files" style="height: 200px; width: 650px;"/>
</p>


`parameters.yaml` - File in which the user sets all the parameters of the workflow.

`synth_generation_template.sh` - This bash script prints the commands of the workflow for the synthetic data generation with parameters spesified by the user.

- Input files: fasta reference file
							 
- Output files: fastq files with pair end reads, ground truth "golden" bam file and bai index file, ground truth "golden" vcf file, Merged bam and vcf files for the ground truth.

`variant_calling_template.sh` - This bash script prints the commands of the workflow for the variant calling process with parameters spesified by the user. 

- Input files:  ground truth Merged bam file
  			 
- Output files: vcf files from variant callers, bam-readcount tsv reports, files with stats for the ground truth Merged bam file

`S4BR.R` - A script, written in R, that calls the appropriate functions to perform the comparison between the ground truth and the caller.

- Input files:  ground truth files, caller vcf file
							 
- Output files: a tsv file with the comparison between the ground truth and the caller.

`S4BR_plot.R` - A script, written in R, that calls the appropriate functions to make visualizations to illustrate the comparison between the ground truth and the caller.

- Input files:  a tsv file with the comparison between the ground truth and the caller, ground truth vcf file, caller vcf file
							 
- Output files: multi planel figure and a Venn plot
 

`paper_plots.R` -  A script, written in R, that outputs the figures for the manuscript.

### Extra scripts

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

## Contribute

We welcome and greatly appreciate any sort of feedback and/or contribution!

If you have any questions, please either open an issue [here](https://github.com/sfragkoul/synth4bench/issues/new) or write to us at `sfragkoul@certh.gr` or `inab.bioinformatics@lists.certh.gr`.

## Citation
Our work has been submitted to the *bioRxiv* preprint repository. If you use our work please cite as follows:

S.-C. Fragkouli, N. Pechlivanis, A. Anastasiadou, G. Karakatsoulis, A. Orfanou, P. Kollia, A. Agathangelidis, and F. E. Psomopoulos, *“Synth4bench: a framework for generating synthetic genomics data for the evaluation of tumor-only somatic variant calling algorithms.”* 2024, doi:[10.1101/2024.03.07.582313](https://www.biorxiv.org/content/10.1101/2024.03.07.582313v1).

## Related Publications
- <div align='justify'>  S.-C. Fragkouli, N. Pechlivanis, A. Anastasiadou, G. Karakatsoulis, A. Orfanou, P. Kollia, A. Agathangelidis, and F. Psomopoulos, <em>synth4bench: Benchmarking Somatic Variant Callers A Tale Unfolding In The Synthetic Genomics Feature Space</em>, <b>23rd European Conference On Computational Biology (ECCB24)</b>, Sep 2024, Turku, Finland doi: TBA </div>

- <div align='justify'> S.-C. Fragkouli, N. Pechlivanis, A. Anastasiadou, G. Karakatsoulis, A. Orfanou, P. Kollia, A. Agathangelidis, and F. Psomopoulos, <em>“Exploring Somatic Variant Callers' Behavior:  A Synthetic Genomics Feature Space Approach”</em>, <b>ELIXIR AHM24</b>, Jun 2024, Uppsala, Sweden, doi: <a href="https://doi.org/10.7490/f1000research.1119793.1">10.7490/f1000research.1119793.1</a></div>


- <div align='justify'> S.-C. Fragkouli, N. Pechlivanis, A. Orfanou, A. Anastasiadou, A. Agathangelidis and F. Psomopoulos, <em>Synth4bench: a framework for generating synthetic genomics data for the evaluation of somatic variant calling algorithms</em>, <b>17th Conference of Hellenic Society for Computational Biology and Bioinformatics (HSCBB)</b>, Oct 2023, Thessaloniki, Greece, doi:<a href="https://doi.org/10.5281/zenodo.8432060">10.5281/zenodo.8432060</a> </div>

  
- <div align='justify'>  S.-C. Fragkouli, N. Pechlivanis, A. Agathangelidis and F. Psomopoulos, <em>Synthetic Genomics Data Generation and Evaluation for the Use Case of Benchmarking Somatic Variant Calling Algorithms</em>, <b>31st Conference in Intelligent Systems For Molecular Biology and the 22nd European Conference On Computational Biology (ISΜB-ECCB23)</b>, Jul 2023, Lyon, France doi:<a href="https://doi.org/10.7490/f1000research.1119575.1">10.7490/f1000research.1119575.1</a> </div>

<div align='center'>
<img src="https://github.com/sfragkoul/synth4bench/blob/main/images/ECCB24_Poster.png" alt="ECCB2024 Poster" style="height: 550px; width:400px;"/>
</div>
