# VIral GEnome ASsembly pipelines for WGS

This repository contains scripts and files to run the bioinformatic analysis of whole genome sequencing of viruses using Illumina or Oxford Nanopore Technologies platforms.

```
Usage:   vigeas <command> [options]

Commands:
  -- ILLUMINA ASSEMBLY
       dnap     For libraries using Illumina DNA Prep or similar [*.fastQ data]
       rnapenr  For libraries using Illumina RNA Prep with Enrichment or similar [*.fastQ data]

  -- OXFORD NANOPORE TECHNOLOGIES ASSEMBLY
       ont      For libraries using ONT SQK-LSK109 or similar [*.fast5/*.pod5 data]

  -- Misc
       update   Update conda [micromamba] dependencies
       panels   List available enrichment panels in this workflow -- for <vigeas rnapenr>
       primers  List available primer schemes in this workflow -- for <vigeas dnap> and <vigeas ont>
```

## Setting up the pipeline

Download and install the pipeline from the github repo:

```sh
git clone --recursive https://github.com/khourious/vigeas.git; cd vigeas
chmod 700 -R INSTALL
bash INSTALL
```
