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


## Using synth4bench


## Requirements
