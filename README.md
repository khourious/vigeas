# VIral GEnome ASsembly pipelines for WGS

This repository contains scripts and files to run the bioinformatic analysis of whole genome sequencing of viruses using Illumina or Oxford Nanopore Technologies platforms.

Until now, this workflow was developed and tested for working with primer schemes available at ``vigeas-illumina -l`` or ``vigeas-nanopore -l``. **Tests with other primer schemes should be performed.**

## Setting up the pipeline

Download and install the pipeline from the github repo:

```sh
git clone --recursive https://github.com/khourious/vigeas.git; cd vigeas
chmod 700 -R INSTALL
bash INSTALL
```
