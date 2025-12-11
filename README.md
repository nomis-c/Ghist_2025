
[![license](https://img.shields.io/badge/license-GPL%20v3-black.svg?style=flat-square)](LICENSE)
# Ghist 2025

## Introduction

This repository contains a Snakemake workflow designed to detect selective sweeps for the GHIST 2025 competition (https://www.synapse.org/Synapse:syn65877330/wiki/631478). 

In the sweep detection category, 4 challenges exist:
 - singlesweep 
 - multiple sweeps
 - singlesweep with background selection
 - multiple sweeps with background selection

The goal is to find the right windows of the chromosome that contains the positive selection signals. VCF files are provided that needs to be analyzed. Final results were submitted as BED files, containing the regions with selective sweep signals.


## Prerequisites
To access the data analyzed in the competition, one needs a Synapse account to download the required datasets (https://accounts.synapse.org/?appId=synapse.org).

After cloning the repository (see Usage) VCF files should be stored as .gz files in each category folder:
- `GHIST_2025_singlesweep.15.final.vcf.gz` in `resources/data/single_sweep/`
- `GHIST_2025_multisweep.15.final.vcf.gz` in `resources/data/multiple_sweeps/`
- `GHIST_2025_singlesweep.growth_bg.15.final.vcf.gz` in `resources/data/single_sweep_background/`
- `GHIST_2025_multisweep.growth_bg.15.final.vcf.gz` in `resources/data/multiple_sweep_background/`

## Usage

1. Install Mambaforge (if not already installed):
   [Mambaforge installation guide](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html)

2. Clone this repository:

```
git clone https://github.com/nomis-c/Ghist_2025.git
cd  Ghist_2025
```

3. Create the environment:

```
mamba env create -f workflow/envs/env.yaml
```

4. Activate the environment:

```
mamba activate ghist
```

5. Run the analysis to replicate results on HPC (recommended):


```
snakemake -s workflow/Snakefile -c 1 --configfile config/singlesweep_final.yaml --profile config/slurm
snakemake -s workflow/Snakefile -c 1 --configfile config/multiplesweep_final.yaml --profile config/slurm
snakemake -s workflow/Snakefile -c 1 --configfile config/singlesweep_bg_final.yaml --profile config/slurm
snakemake -s workflow/Snakefile -c 1 --configfile config/multiplesweep_bg_final.yaml --profile config/slurm
```

6. Run the analysis to replicate results locally:

```
snakemake -s workflow/Snakefile -c 1 --configfile config/singlesweep_final.yaml
snakemake -s workflow/Snakefile -c 1 --configfile config/multiplesweep_final.yaml
snakemake -s workflow/Snakefile -c 1 --configfile config/singlesweep_bg_final.yaml
snakemake -s workflow/Snakefile -c 1 --configfile config/multiplesweep_bg_final.yaml
```


Users should adjust the `config.yaml` file in `config/slurm` according to their job scheduler.
