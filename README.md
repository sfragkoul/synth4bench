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
- [Documentation](https://github.com/sfragkoul/synth4bench/tree/main#documentation)
- [Contribute](https://github.com/sfragkoul/synth4bench/tree/main#contribute)
- [Citation](https://github.com/sfragkoul/synth4bench/tree/main#citation)


## Motivation

<div align='justify'> Variant calling plays an important role in identifying genetic lesions. In the case of variants at low frequency (≤10%) identification becomes more difficult and the challenge that rises is the absence of a Ground Truth for reliable and consistent identification and benchmarking. </div>

## Description of Framework

<p align="center"> 
<img src="https://github.com/sfragkoul/synth4bench/blob/main/images/schematic.png" alt="synth4bench schematic" style="height: 270px; width:900px;"/>
</p>

<div align='justify'> Our framework focuses on addressing the challenge of variant calling, particularly for variants at low frequencies (≤10%). The main goal is to develop a reliable and consistent method for identifying genetic lesions, specifically in the context of cancer-associated genomic alterations. The absence of a ground truth, which refers to a reliable reference dataset with known variants, makes benchmarking and evaluating variant calling algorithms difficult. To overcome this challenge, the following steps are outlined in the framework:

1. **Synth Data Generation**: Synthetic genomics data is generated using the NEATv3.3 simulator in order to create datasets that mimic real genome data that plays the role of Ground Truth.

2. **Benchmarking Variant Callers**: Five somatic variant callers (i.e. GATK-Mutect2, Freebayes, VarDict, VarScan2 and LoFreq) are evaluated using our synthetic Ground Truth dataset.

 </div>

## Data Download
All data are open and available in [Zenodo](https://zenodo.org/records/10683211). For specific instructions please see our [**UserGuide**](https://github.com/sfragkoul/synth4bench/blob/main/docs/UserGuide.md#data-download).

## Installation
1.  To create the conda environment that was used for the analysis run `conda env create -f environment.yml` and to activate it run `conda activate synth4bench`.
2. To install NEATv3.3, download version [v3.3](https://github.com/ncsa/NEAT/releases/tag/3.3). To call the main script run the command `python gen_reads.py --help`. For any further info please see the README.md file from the downloaded files of version 3.3.
3. To install bam-readcount follow their [instructions](https://github.com/genome/bam-readcount/tree/master?tab=readme-ov-file#build) and then run `build/bin/bam-readcount --help` to see that it has being installed properly. If you face any problems with the installation during the `make` command please add the executable that can be found [here](https://github.com/sfragkoul/synth4bench/tree/main/bam-readcount) in the `bam-readcount\build\bin` folder.
4. The extra script `vscan_pileup2cns2vcf.py` for VarScan can be found [here](https://github.com/sfragkoul/Varscan2VCF).


## Execution
	
Just fill in the parameters in the `parameters.yaml` file and then run `bash s4b_run.sh`.

For detailed instructions regarding the execution please read our [**UserGuide**](https://github.com/sfragkoul/synth4bench/blob/main/docs/UserGuide.md#executing-synth4bench).

<div align='justify'> 

## Documentation
For more info regarding the documentation please visit [here](https://github.com/sfragkoul/synth4bench/blob/main/docs/Documentation.md).

## Contribute

We welcome and greatly appreciate any sort of feedback and/or contribution!

If you have any questions, please either open an issue [here](https://github.com/sfragkoul/synth4bench/issues/new) or send an email to `sfragkoul@certh.gr`.

## Citation
Our work has been submitted to the *bioRxiv* preprint repository. If you use our work please cite as follows:

S.-C. Fragkouli, N. Pechlivanis, A. Anastasiadou, G. Karakatsoulis, A. Orfanou, P. Kollia, A. Agathangelidis, and F. E. Psomopoulos, *“Synth4bench: a framework for generating synthetic genomics data for the evaluation of tumor-only somatic variant calling algorithms.”* 2024, doi:[10.1101/2024.03.07.582313](https://www.biorxiv.org/content/10.1101/2024.03.07.582313v1).

## Related Publications
- <div align='justify'>  S.-C. Fragkouli, N. Pechlivanis, A. Anastasiadou, G. Karakatsoulis, A. Orfanou, P. Kollia, A. Agathangelidis, and F. Psomopoulos, <em>synth4bench: Benchmarking Somatic Variant Callers A Tale Unfolding In The Synthetic Genomics Feature Space</em>, <b>23rd European Conference On Computational Biology (ECCB24)</b>, Sep 2024, Turku, Finland doi: <a href="https://zenodo.org/records/14186510">10.5281/zenodo.14186509</a> </div>

- <div align='justify'> S.-C. Fragkouli, N. Pechlivanis, A. Anastasiadou, G. Karakatsoulis, A. Orfanou, P. Kollia, A. Agathangelidis, and F. Psomopoulos, <em>“Exploring Somatic Variant Callers' Behavior:  A Synthetic Genomics Feature Space Approach”</em>, <b>ELIXIR AHM24</b>, Jun 2024, Uppsala, Sweden, doi: <a href="https://doi.org/10.7490/f1000research.1119793.1">10.7490/f1000research.1119793.1</a></div>


- <div align='justify'> S.-C. Fragkouli, N. Pechlivanis, A. Orfanou, A. Anastasiadou, A. Agathangelidis and F. Psomopoulos, <em>Synth4bench: a framework for generating synthetic genomics data for the evaluation of somatic variant calling algorithms</em>, <b>17th Conference of Hellenic Society for Computational Biology and Bioinformatics (HSCBB)</b>, Oct 2023, Thessaloniki, Greece, doi:<a href="https://doi.org/10.5281/zenodo.8432060">10.5281/zenodo.8432060</a> </div>

  
- <div align='justify'>  S.-C. Fragkouli, N. Pechlivanis, A. Agathangelidis and F. Psomopoulos, <em>Synthetic Genomics Data Generation and Evaluation for the Use Case of Benchmarking Somatic Variant Calling Algorithms</em>, <b>31st Conference in Intelligent Systems For Molecular Biology and the 22nd European Conference On Computational Biology (ISΜB-ECCB23)</b>, Jul 2023, Lyon, France doi:<a href="https://doi.org/10.7490/f1000research.1119575.1">10.7490/f1000research.1119575.1</a> </div>
