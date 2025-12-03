<p align="center">
  <img src="https://github.com/sfragkoul/synth4bench/blob/main/images/synth4bench_logo_no_bg.png" alt="synth4bench logo" height="80" width="700"/>
</p>

---

## Abstract

<div align='justify'> 
Somatic variant calling algorithms are essential for detecting genomic alterations. However, evaluating their performance can be challenging due to the lack of high-quality ground truth datasets. To address this issue, we developed a synthetic genomics data generation framework for benchmarking tumor-only somatic variant calling algorithms. We generated synthetic datasets using the NEAT simulator. Subsequently, we thoroughly evaluated the performance of variant calling algorithms including <strong>GATK-Mutect2, FreeBayes, VarDict, VarScan2 and LoFreq</strong> on these datasets, comparing results against the ground truth. Synthetic datasets provide an excellent ground truth for studying the performance and behavior of somatic variant calling algorithms.
</div>

---

## Table of Contents

- [Abstract](#abstract)
- [Description](#description)
- [Installation](#installation)
- [Data Download](#data-download)
- [Execution](#execution)
- [Documentation](#documentation)
- [Contribute](#contribute)
- [Citation](#citation)

---


## Description

<p align="center">
  <img src="https://github.com/sfragkoul/synth4bench/blob/main/images/schematic.png" alt="synth4bench schematic" height="270" width="900"/>
</p>

<div align='justify'>
Our framekwork addresses the challenge of variant calling, particularly for low-alle frequency variants (≤10%). The goal is to develop a reliable and consistent method for identifying genetic lesions in cancer-associated genomic alterations. The lack of ground truth datasets complicates benchmarking and evaluation. To overcome this, our framework includes:

1. **Synthetic Data Generation:** Using the NEAT v3.3 simulator, we generate synthetic genomics data that mimics real genome sequences, serving as a ground truth.
2. **Identifying three independent analyses**: to adress three variant types independently: SNVs true variants, SNVs noise and indels.
3. **Benchmarking Variant Callers:** We evaluate five somatic variant callers, GATK-Mutect2, Freebayes, VarDict, VarScan2, and LoFreq — using these synthetic datasets.
</div>

---

## Data Download

All data are openly available on [Zenodo](https://zenodo.org/records/16524193). For specific instructions, refer to our [**User Guide**](https://github.com/sfragkoul/synth4bench/blob/main/docs/UserGuide.md#data-download).

---

## Installation

1. **Create the Conda environment:**

    ```bash
    conda env create -f environment.yml
    conda activate synth4bench
    ```

2. **Install NEAT v3.3:**

    Download version [v3.3](https://github.com/ncsa/NEAT/releases/tag/3.3).  
    To call the main script:

    ```bash
    python gen_reads.py --help
    ```

    For further details, see the NEAT README included in the download.

3. **Install bam-readcount:**

    Follow their [installation instructions](https://github.com/genome/bam-readcount#build).  
    After building, verify installation:

    ```bash
    build/bin/bam-readcount --help
    ```

    If you encounter issues during the `make` process, you can alternatively use the executable available [here](https://github.com/sfragkoul/synth4bench/tree/main/bam-readcount) and place it in the `bam-readcount/build/bin` folder.

4. **Download VarScan Extra Script:**

    The extra script `vscan_pileup2cns2vcf.py` for VarScan is available [here](https://github.com/sfragkoul/Varscan2VCF).

---

## Execution

Simply configure your parameters in the `parameters.yaml` file, then execute:

```bash
bash s4b_run.sh
```

This single command generates synthetic data, runs variant calling for all selected tools, and performs downstream analysis and plotting.

For full execution instructions, see our [**User Guide**](https://github.com/sfragkoul/synth4bench/blob/main/docs/UserGuide.md#executing-synth4bench).

---

## Documentation

For further documentation, visit the [documentation page](https://github.com/sfragkoul/synth4bench/blob/main/docs/Documentation.md).

---

## Contribute

We welcome and greatly appreciate any feedback or contributions!

If you have questions, please open an issue [here](https://github.com/sfragkoul/synth4bench/issues/new) or email `sfragkoul@certh.gr`.

---

## Citation

Our work has been submitted to the *bioRxiv* preprint repository. If you use synth4bench, or any of our scripts/code, please cite:

> S.-C. Fragkouli, N. Pechlivanis, A. Anastasiadou, G. Karakatsoulis, A. Orfanou, P. Kollia, A. Agathangelidis, and F. E. Psomopoulos, *“Synth4bench: Synthetic Data Generation for Benchmarking Tumor-Only Somatic Variant Calling Algorithms”* 2025, doi:[10.1101/2024.03.07.582313](https://www.biorxiv.org/content/10.1101/2024.03.07.582313v2).

---

## Related Publications

- <div align='justify'> S.-C. Fragkouli, N. Pechlivanis, A. Anastasiadou, G. Karakatsoulis, A. Orfanou, P. Kollia, A. Agathangelidis, and F. Psomopoulos, <em>synth4bench: Benchmarking Somatic Variant Callers – A Tale Unfolding In The Synthetic Genomics Feature Space</em>, <b>23rd European Conference On Computational Biology (ECCB24)</b>, Sep 2024, Turku, Finland, doi: <a href="https://zenodo.org/records/14186510">10.5281/zenodo.14186509</a> </div>

- <div align='justify'> S.-C. Fragkouli, N. Pechlivanis, A. Anastasiadou, G. Karakatsoulis, A. Orfanou, P. Kollia, A. Agathangelidis, and F. Psomopoulos, <em>“Exploring Somatic Variant Callers' Behavior: A Synthetic Genomics Feature Space Approach”</em>, <b>ELIXIR AHM24</b>, Jun 2024, Uppsala, Sweden, doi: <a href="https://doi.org/10.7490/f1000research.1119793.1">10.7490/f1000research.1119793.1</a></div>

- <div align='justify'> S.-C. Fragkouli, N. Pechlivanis, A. Orfanou, A. Anastasiadou, A. Agathangelidis and F. Psomopoulos, <em>Synth4bench: a framework for generating synthetic genomics data for the evaluation of somatic variant calling algorithms</em>, <b>17th Conference of Hellenic Society for Computational Biology and Bioinformatics (HSCBB)</b>, Oct 2023, Thessaloniki, Greece, doi:<a href="https://doi.org/10.5281/zenodo.8432060">10.5281/zenodo.8432060</a> </div>

- <div align='justify'> S.-C. Fragkouli, N. Pechlivanis, A. Agathangelidis and F. Psomopoulos, <em>Synthetic Genomics Data Generation and Evaluation for the Use Case of Benchmarking Somatic Variant Calling Algorithms</em>, <b>31st Conference in Intelligent Systems For Molecular Biology and the 22nd European Conference On Computational Biology (ISΜB-ECCB23)</b>, Jul 2023, Lyon, France, doi:<a href="https://doi.org/10.7490/f1000research.1119575.1">10.7490/f1000research.1119575.1</a> </div>









