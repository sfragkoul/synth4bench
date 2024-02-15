# synth4bench: a framework for generating synthetic genomics data for the evaluation of somatic variant calling algorithms

## Abstract
<div align='justify'> Somatic variant calling algorithms are widely used to detect genomic alterations associated with cancer. Evaluating the performance of these algorithms can be challenging due to the lack of high-quality ground truth datasets. To address this issue, we developed a synthetic genomics data generation and evaluation framework for benchmarking somatic variant calling algorithms. We generated synthetic datasets based on sequence data from the TP53 gene, using the NEAT(NExt-generation sequencing Analysis Toolkit) simulator. Subsequently, we thoroughly evaluated the performance of variant calling algorithms using <strong>Mutect2, Freebayes, VarDict, VarScan and LoFreq</strong> on these datasets, and compared the results to the “golden” files produced by NEAT containing the actual variations. Our results demonstrate that the synthetic datasets generated using our framework can accurately capture the complexity and diversity of real cancer genomic data. Moreover, the synthetic datasets provide an excellent ground truth for evaluating the performance of somatic variant calling algorithms. Altogether, our framework provides a valuable resource for testing the performance of somatic variant calling algorithms, enabling researchers to evaluate and improve the accuracy of these algorithms for cancer genomics applications.</div>

## Table of Contents

- [Abstract](https://github.com/sfragkoul/synth4bench/tree/main#abstract)
- [Motivation](https://github.com/sfragkoul/synth4bench/tree/main#motivation)
- [Discription of Pipeline](https://github.com/sfragkoul/synth4bench/tree/main#discription-of-pipeline)
- [Installation](https://github.com/sfragkoul/synth4bench/tree/main#installation)
- [Data Download](https://github.com/sfragkoul/synth4bench/tree/main#data-download)
- [Execution](https://github.com/sfragkoul/synth4bench/tree/main#execution)
- [Contribute](https://github.com/sfragkoul/synth4bench/tree/main#contribute)
- [Citation](https://github.com/sfragkoul/synth4bench/tree/main#citation)


## Motivation

<div align='justify'> Variant calling plays an important role in identifying genetic lesions. In the case of variants at low frequency (≤10%) identification becomes more difficult and the challenge that rises is the absence of a Ground Truth for reliable and consistent identification and benchmarking. </div>

## Discription of Pipeline

<div align='justify'> Our pipeline focuses on addressing the challenge of variant calling, particularly for variants at low frequencies (≤10%). The main goal is to develop a reliable and consistent method for identifying genetic lesions, specifically in the context of cancer-associated genomic alterations. The absence of a ground truth, which refers to a reliable reference dataset with known variants, makes benchmarking and evaluating variant calling algorithms difficult. To overcome this challenge, the following steps are outlined in the pipeline:

1. Data Generation: Synthetic genomics data is generated based on the TP53 gene using the NEAT simulator in order to create synthetic datasets that mimic real cancer genome data.

2. Defining Ground Truth: The "Ground Truth" is established by creating 10 individual datasets (each one of average 500 coverage)containing Single Nucleotide Polymorphisms (SNPs) and Insertions/Deletions (INDELs). The genomic regions where variants accure with 100% Allele Frequency are chosen. The reason behind this choice is to avoid variants that are related to errors and products of noise. Then all these datasets are merged into one single file of 5000 coverage and the allele frequency is again measured at these genomic regions of interest.

3. Benchmarking Variant Callers: Somatic variant callers are evaluated using this synthetic Ground Truth dataset. The GATK-Mutect2, Freebayes, VarDict, VarScan2 and LoFreq variant callesr are assessed for their performance on our synthetic dataset. Their impact at low frequencies (≤10%) is explored, as these are particularly challenging to detect accurately.

The pipeline's overall aim is to provide a robust framework for evaluating the performance of somatic variant calling algorithms by using synthetic datasets that closely resemble real cancer genome data. By having a reliable ground truth, we can thoroughly test and improve the accuracy of variant calling algorithms for cancer genomics applications. This pipeline represents an essential step towards more precise and effective identification of genetic lesions associated with cancer and other diseases. </div>

## Data Download
All data are open and available in [Zenodo](https://zenodo.org/record/8095898).

## Installation
1. To install NEAT follow the intructions [here](https://github.com/ncsa/NEAT/blob/master/README.md#installation).
2. To install GATK-Mutect2 follow the intructions [here](https://anaconda.org/bioconda/gatk).
3. To install FreeBayes follow the intructions [here](https://anaconda.org/bioconda/freebayes).
4. To install VarDict follow the intructions [here](https://anaconda.org/bioconda/vardict).
5. To install VarScan follow the intructions [here](https://anaconda.org/bioconda/varscan).
6. To install LoFreq follow the intructions [here](https://anaconda.org/bioconda/lofreq).
7. To install SAMtools follow the intructions [here](https://anaconda.org/bioconda/samtools).
8. To install BCFtools follow the intructions [here](https://anaconda.org/bioconda/bcftools).
9. To install bam-readcount follow the instructions [here](https://github.com/genome/bam-readcount/tree/master#installation).
10. To create the conda environment that was used for the analysis run `conda env create -f environment.yml`
11. To install R packages dependencies run this command `install.packages(c("stringr", "data.table", "dplyr", "vcfR", "ggplot2", "ggvenn", "ggforce", "ggsci", "GenomicRanges", "systemPipeR", "AnnotationHub", "seqinr", "patchwork"))`


## Execution
<div align='justify'> Here follows the list of all scripts and their description:

`01_synth4bench.sh` - This bash script is the basis of synth4bench workflow. It calls NEAT in order to generate 10 individual synthetic data datasets, create one Merged bam file, performs some preprocess steps before implementing somatic variant calling and produces bam report files with the genomic content at certain chromosomal positionsusing bam-readcount. Please replace all `path/to/files/` with desired paths.

- Input: fasta reference file
							 
- Output: fastq files with pair end reads, "golden" bam file and bai index file, "golden" vcf file, Merged bam file, processed bam files and vcf file with all variants that were detected, tsv file with the genomic content

`02_load_bam_reports.R` - This R script compares the variants that Mutect2 reported against the ground truth. Firsty it identifies the variants with 100% Allele Frequency(AF) in the individual bam files and then caclulates their AF in the final Merged bam file.

- Input:  bam-readcount tsv reports, vcf file from Mutect2
							 
- Output: tsv file containing information regarding the ground truth variants

`03_clean_and_annotate.R` - This R script takes the tsv file and adds annotation information based on the gene chromosomal positions.

- Input:  tsv file containing information regarding the ground truth variants
							 
- Output: annotated tsv file containing information regarding the ground truth variants and annotation information

`04_patchwork_gatk.R` - This R script produces the final Figure of the Benchmarking of GATK.

- Input:  annotated tsv file, ground truth vcf, Mutect2 vcf
							 
- Output: final Figure for the Benchmarking of GATK </div>

`plot_helpers_gatk.R` - This R script incudes all necessary funtions for `04_patchwork_gatk.R` script.

`plot_libraries.R` - This R script incudes all necessary libraries for `04_patchwork_gatk.R` script.

### Extra scripts
For the case of VarScan an extra step was required to convert its output to the standard VCF format. The script that was developed can be found [here](https://github.com/sfragkoul/Varscan2VCF).

## Contribute

We welcome any sort of contribution!
If you have any developer-related questions, please open an issue or write us at sfragkoul@certh.gr.


## Related Publications
- Styliani-Christina Fragkouli, Nikolaos Pechlivanis, Aspasia Orfanou, Anastasia Anastasiadou, Andreas Agathangelidis and Fotis Psomopoulos, *Synth4bench: a framework for generating synthetic genomics data for the evaluation of somatic variant calling algorithms*, 17th Conference of Hellenic Society for Computational Biology and Bioinformatics (HSCBB), Oct 2023, Thessaloniki, Greece, doi: [10.5281/zenodo.8432060](https://doi.org/10.5281/zenodo.8432060)

  
- Styliani-Christina Fragkouli, Nikolaos Pechlivanis, Andreas Agathangelidis and Fotis Psomopoulos, *Synthetic Genomics Data Generation and Evaluation for the Use Case of Benchmarking Somatic Variant Calling Algorithms*, 31st Conference in Intelligent Systems For Molecular Biology and the 22nd European Conference On Computational Biology (ISΜB-ECCB23), Jul 2023, Lyon, France doi: [10.7490/f1000research.1119575.1](https://doi.org/10.7490/f1000research.1119575.1)
