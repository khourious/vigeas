# VIral GEnome ASsembly pipelines for WGS

This repository contains scripts and files to run the bioinformatic analysis of whole genome sequencing of viruses using Illumina or Oxford Nanopore Technologies platforms.

```
Usage: vigeas <command> [options]

Commands:
  -- ILLUMINA ASSEMBLY
       ill-amp  For libraries using Amplicon chemistry [*.fastQ data]
       ill-hyb  For libraries using Hybrid-capture chemistry [*.fastQ data]

  -- OXFORD NANOPORE TECHNOLOGIES ASSEMBLY
       ont      For libraries using ONT SQK-LSK109 or similar [*.fast5/*.pod5 data]

  -- Misc
       makedb   Create a BLAST database in this workflow -- for <vigeas ill-hyb>
       panels   List available enrichment panels in this workflow -- for <vigeas ill-hyb>
       primers  List available primer schemes in this workflow -- for <vigeas ill-amp> and <vigeas ont>
       update   Update conda [micromamba] dependencies
```

## Setting up the pipeline

Download and install the pipeline from the github repo:

```sh
git clone --recursive https://github.com/khourious/vigeas.git; cd vigeas
chmod 700 -R INSTALL
bash INSTALL
```