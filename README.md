# VIral GEnome ASsembly pipelines for WGS

This repository contains scripts and files to run the bioinformatic analysis of whole genome sequencing of viruses using Illumina or Oxford Nanopore Technologies platforms.

```
Usage: vigeas <command> or <miscellaneous>

Commands:
  ill   For Illumina Sequencing [*.fastQ data]
  ont   For ONT Sequencing [*.pod5 data]

Miscellaneous:
  clr3     List supported Clair3 models
  makedb   Create a BLAST database in this workflow -- for <vigeas ill -x hyb>
  panels   List available enrichment panels in this workflow -- for <vigeas ill -x hyb>
  primers  List available primer schemes in this workflow -- for <vigeas ill -x amp> and <vigeas ont -x bda>
  update   Update conda [micromamba] dependencies
```

## Setting up the pipeline

Download and install the pipeline from the github repo:

```sh
git clone --recursive https://github.com/khourious/vigeas.git; cd vigeas
chmod 700 -R INSTALL
bash INSTALL
```
