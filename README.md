# Synthetic Genomics Data Generation and Evaluation for the Use Case of Benchmarking Somatic Variant Calling Algorithms


This is the repository for the analysis that was presented as poster in the **31st Conference in Intelligent Systems For Molecular Biology and the 22nd European Conference On Computational Biology (ISΜB-ECCB23)**, Lyon, France.

<img src="/poster%20files/ISMBECCB2023_1226_Fragkouli_poster.png" alt="ISΜB-ECCB23 Poster" width="600"/>

All data are open and available in [Zenodo](https://zenodo.org/record/8095898).

## Abstract
Somatic variant calling algorithms are widely used to detect genomic alterations associated with cancer. However, evaluating the performance of these algorithms can be challenging due to the lack of high-quality ground truth datasets. To address this issue, we developed a synthetic genomics data generation and evaluation framework for benchmarking somatic variant calling algorithms. We generated synthetic datasets based on data from the TP53 gene, using the NEAT simulator. We then thoroughly evaluated the performance of GATK-Mutect2 on these datasets, and compared the results to the “golden” files produced by NEAT that contain the true variations. Our results demonstrate that the synthetic datasets generated using our framework can accurately capture the complexity and diversity of real cancer genome data. Moreover, the synthetic datasets provide an excellent ground truth for evaluating the performance of somatic variant calling algorithms. Our framework provides a valuable resource for testing the performance of somatic variant calling algorithms, enabling researchers to evaluate and improve the accuracy of these algorithms for cancer genomics applications.

## Motivation

Variant calling plays an important role in identifying genetic lesions. In the case of variants at low frequency (≤10%) identification becomes more difficult and the challenge that rises is the absence of a Ground Truth for reliable and consistent identification and benchmarking.

## Discription of Pipeline

Our pipeline focuses on addressing the challenge of variant calling, particularly for variants at low frequencies (≤10%). The main goal is to develop a reliable and consistent method for identifying genetic lesions, specifically in the context of cancer-associated genomic alterations. The absence of a ground truth, which refers to a reliable reference dataset with known variants, makes benchmarking and evaluating variant calling algorithms difficult. To overcome this challenge, the following steps are outlined in the pipeline:

1. Data Generation: Synthetic genomics data is generated based on the TP53 gene using the NEAT simulator in order to create synthetic datasets that mimic real cancer genome data.

2. Defining Ground Truth: The "Ground Truth" is established by creating 10 individual datasets (each one of average 500 coverage)containing Single Nucleotide Polymorphisms (SNPs) and Insertions/Deletions (INDELs). The genomic regions where variants accure with 100% Allele Frequency are chosen. The reason behind this choice is to avoid variants that are related to errors and products of noise. Then all these datasets are merged into one single file of 5000 coverage and the allele frequency is again measured at these genomic regions of interest.

3. Benchmarking Variant Callers: Somatic variant callers are evaluated using this synthetic Ground Truth dataset. The GATK-Mutect2 variant caller is specifically assessed for its performance on our synthetic dataset. Their impact at low frequencies (≤10%) is explored, as these are particularly challenging to detect accurately.

The pipeline's overall aim is to provide a robust framework for evaluating the performance of somatic variant calling algorithms by using synthetic datasets that closely resemble real cancer genome data. By having a reliable ground truth, we can thoroughly test and improve the accuracy of variant calling algorithms for cancer genomics applications. This pipeline represents an essential step towards more precise and effective identification of genetic lesions associated with cancer and other diseases.


## Installation
1. To install NEAT follow the intructions [here](https://github.com/ncsa/NEAT/blob/master/README.md#installation).
2. To install GATK-Mutect2 follow the intructions [here](https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4).
3. To install bam-readcount follow the instructions [here](https://github.com/genome/bam-readcount/tree/master#installation).
4. To create the conda environment that was used for the analysis run `conda env create -f environment.yml`


## Execution
Here follows the list of all scripts and their description:

`01_SynthDataGeneration.sh` - This bach script calls NEAT in order to generate 10 individual synthetic data datasets.

- Input: fasta reference file
							 
- Output: fastq files with pair end reads, "golden" bam file and bai index file, "golden" vcf file

`02_bam-readcount_reports.sh` - This bash script calls bam-readcount to produce the bam reports with the genomic content at certain chromosomal positions.

- Input: fasta reference file, "golden" bam file
							 
- Output: tsv file with the genomic content

`03_VariantCalling.sh` - This bash script takes the Merged bam file and performs some preprocess steps before using Mutect2 to implement somatic variant calling.

- Input: Merged bam file, fasta reference file
							 
- Output: processed bam files and vcf file with all variants that were detected

`04_load_bam_reports.R` - This R script compares the variants that Mutect2 reported against the ground truth. Firsty it identifies the variants with 100% Allele Frequency(AF) in the individual bam files and then caclulates their AF in the final Merged bam file.

- Input:  bam-readcount tsv reports, vcf file from Mutect2
							 
- Output: tsv file containing information regarding the ground truth variants

`05_clean_and_annotate.R` - This R script takes the tsv file and adds annotation information based on the gene chromosomal positions.

- Input:  tsv file containing information regarding the ground truth variants
							 
- Output: annotated tsv file containing information regarding the ground truth variants and annotation information

`06_bar_plots.R` - This R script produces the Allele Frequency and Coverage Barplots of the ground truth and the detected variants.

- Input: annotated tsv file
							 
- Output: Allele Frequency and Coverage Barplots

`06_bubble_plots.R` - This R script produces the bubble plot of the SNIPs of the ground truth and the detected variants.

- Input:  annotated tsv file
							 
- Output: bubble plot of the SNIPs

`06_density_plot.R` - This R script produces the Allele Frequency Density plots of Ground Truth and detected Variants per DNA Base.

- Input:  annotated tsv file
							 
- Output: Allele Frequency Density plots

`06_reference_barplot.R` - This R script produces the barplots of the genomic content of the reference.

- Input:  fasta reference file
							 
- Output: barplots of the genomic content of the reference

`07_mutation_overlap.R` - This R script produces the Venn plot of the Overall Variants

- Input:  ground truth vcf and Mutect2 vcf
							 
- Output: Venn plot of the Overall Variants


`08_patchwork.R` - This R script produces the final Figure of the Benchmarking of the poster.

- Input:  all produced plots
							 
- Output: final Figure for the poster


(The first 3 steps of the pipeline are optional since all the files are provided either in this repository or in the zenodo link.)


## Contribute

We welcome any sort of contribution!
If you have any developer-related questions, please open an issue or write us at sfragkoul@certh.gr.


## Citation
Styliani-Christina Fragkouli and Nikolaos Pechlivanis and Andreas Agathangelidis and Fotis Psomopoulos, *Synthetic Genomics Data Generation and Evaluation for the Use Case of Benchmarking Somatic Variant Calling Algorithms*, 31st Conference in Intelligent Systems For Molecular Biology and the 22nd European Conference On Computational Biology (ISΜB-ECCB23), 2023, doi: [10.7490/f1000research.1119575.1](https://doi.org/10.7490/f1000research.1119575.1)

