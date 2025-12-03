# Synth4bench User Guide

This guide helps users execute the synth4bench pipeline from data download to analysis and visualization.

---

## Contents

- [Data Download](#data-download)
- [Installation](#installation)
- [Key Features](#key-features)
- [Executing synth4bench](#executing-synth4bench)
- [Requirements](#requirements)

---

## Data Download

The datasets required for this analysis are **available for open access** on [Zenodo](https://zenodo.org/records/16524193).

This repository contains 10 synthetic genomics datasets, generated specifically for **benchmarking somatic variant callers**. Each dataset was produced with [NEAT v3](https://github.com/ncsa/NEAT/releases/tag/3.3), based on the *TP53* gene of Homo sapiens, and provides valuable resources for analyzing the effects of various NGS parameters on **tumor-only somatic variant calling** algorithms.

---

### Data Overview

These datasets explore two primary variables to observe their impact on somatic variant calling:

- **Coverage:** Five datasets vary in coverage while keeping read length constant (150 bp). Coverage levels include 300x, 700x, 1000x, 3000x, and 5000x.

- **Read Length:** Another set of five datasets varies in read length while maintaining a coverage level of 1000x. Read lengths include 50 bp, 100 bp, 150 bp, 170 bp, 200 bp, and 300 bp.

Each dataset contains paired-end reads, providing compatibility with most standard analysis pipelines for synthetic genomics data.

---

### Available Files

There are two compressed folders for download:

- **Reference folders**: `reference.rar` (25.7 kB)

  Contains reference sequences necessary for aligning or comparing generated synthetic reads.

- **Synthetic folders**: `synth_datasets.rar` (55.9 MB)

  Includes all 10 synthetic datasets, with filenames indicating coverage and read length details. Each dataset is structured and named according to the parameters it explores.

---

## Installation

### Create the Conda Environment

To create the conda environment used for the analysis, run:

```bash
conda env create -f environment.yml
```

---

### Activate the Conda Environment

After creating the environment, activate it with:

```bash
conda activate synth4bench
```

---

### Install NEAT v3.3

Download [NEAT v3.3](https://github.com/ncsa/NEAT/releases/tag/3.3).  
To view NEAT options:

```bash
python gen_reads.py --help
```

---

### Install bam-readcount

Follow the installation [instructions](https://github.com/genome/bam-readcount#build).  
Verify the install with:

```bash
build/bin/bam-readcount --help
```

---

### Install R Package Dependencies

Install the required R packages in your **R console**:

```r
install.packages(c(
  "stringr", "data.table", "vcfR", "ggplot2", 
  "ggvenn", "ggforce", "ggsci", "patchwork", 
  "optparse", "GenomicAlignments", "Rsamtools"
))
```

---

### Download VarScan Extra Script

The VarScan extra script `vscan_pileup2cns2vcf.py` is available [here](https://github.com/sfragkoul/Varscan2VCF).

---

## Key Features

Once the environment and dependencies are installed, follow these steps to execute synth4bench.

---

### Prepare Parameters

Edit the `parameters.yaml` file with your analysis settings. Example:

```yaml
working_directory: "/path/to/synth4bench_analysis"
folder: results
runs: 10
callers: "Mutect2,Freebayes,LoFreq,VarScan,VarDict"
output_bam_merged: Merged
```

Adjust paths and values for your system.

---

## Executing synth4bench

This section guides you from setup to analysis.

---

### Step 1 — Create an Analysis Folder

```bash
mkdir synth4bench_analysis
cd synth4bench_analysis
```

Download and extract your reference files as previously described.

---

### Step 2 — Configure parameters.yaml

Edit `parameters.yaml` to:

- define your working directory
- set output folder names
- specify coverage, read length, runs
- choose variant callers
- set paths to reference files

For local runs, you can replace `/path/to/` with `.`.

---

### Step 3 — Run the Full Pipeline

Instead of manually running each Bash template, you now only need one command:

```bash
bash s4b_run.sh
```

This performs:

✅ **Synthetic read generation**  
✅ **BAM pre-processing with samtools**  
✅ **Variant calling** for all callers you specified  
✅ **R analysis and plotting** for each caller

---

## How the Pipeline Works

`s4b_run.sh`:

- automatically parses your YAML file
- generates and executes:
  - `gen_<folder>.sh` → synthetic reads
  - `vc_<folder>.sh` → variant calling
- runs **R/S4BR.R** and **R/S4BR_plot.R** for each caller

Example of what’s executed for each caller:

```bash
Rscript R/S4BR.R \
  -c Mutect2 \
  -r 10 \
  -v ./results \
  -m Merged

Rscript R/S4BR_plot.R \
  -c Mutect2 \
  -w ./results \
  -m Merged
```

This loop repeats for each caller in your YAML configuration.

---

## Generated Scripts

`s4b_run.sh` produces:

```
bash/gen_<folder>.sh
bash/vc_<folder>.sh
```

These are **auto-generated** and do not need to be run manually—they are called automatically from `s4b_run.sh`.

> **Note:** These files are ignored in `.gitignore` and should not be committed to Git.

---

## Check R Script Parameters

To check help for either R script:

```bash
Rscript R/S4BR.R --help
Rscript R/S4BR_plot.R --help
```

---

### Options for `S4BR.R`

| Option            | Description                                                               |
|-------------------|----------------------------------------------------------------------------|
| `-c, --caller`    | Choose caller name (Freebayes, Mutect2, LoFreq, VarDict, VarScan)          |
| `-r, --runs`      | Number of individual synthetic runs combined to form the final merged file |
| `-v, --working_directory` | Path where files and results will be stored                      |
| `-m, --merged_file` | Name of the merged BAM file                                             |
| `-h, --help`      | Show help message and exit                                                |

---

### Options for `S4BR_plot.R`

| Option            | Description                                                               |
|-------------------|----------------------------------------------------------------------------|
| `-c, --caller`    | Choose caller name (Freebayes, Mutect2, LoFreq, VarDict, VarScan)          |
| `-w, --working_directory` | Path where files and results will be stored                      |
| `-m, --merged_file` | Name of the merged BAM file                                             |
| `-h, --help`      | Show help message and exit                                                |

---

## Analyze the Results

After running `s4b_run.sh`:

- Find variant calls and summary files under:
  ```
  <working_directory>/<folder>/
  ```
- Open files like:
  - `analysis_output.tsv` → summarizes called variants
  - plots from `S4BR_plot.R`

These outputs help assess variant caller performance under your simulated conditions.

---

## Requirements

- Linux/Unix shell environment (tested on Ubuntu, WSL)
- Conda for environment management
- R version ≥ 4.0
- NEAT v3.3 installed
- bam-readcount

---



