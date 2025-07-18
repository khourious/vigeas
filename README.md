# VIral GEnome ASsembly pipelines for WGS

This repository contains scripts and files to run the bioinformatic analysis of whole genome sequencing of viruses using Illumina or Oxford Nanopore Technologies platforms.

## INSTALLATION

```sh
git clone --recursive https://github.com/khourious/vigeas.git && cd vigeas
chmod 700 -R INSTALL scripts
bash INSTALL
```

## USAGE

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
  update   Update conda/mamba dependencies
  version  Show last update information
```

## CITATION

## CONTRIBUTIONS
Thanks to:
- Ricardo Khouri ([Rkhour0](https://github.com/orgs/khourious/people/Rkhour0)) for mentorship and discussions on improving in the `-x amp` and `-x ill` workflows.
- Verity Hill ([ViralVerity](https://github.com/ViralVerity)) for her collaborative troubleshooting regarding low-coverage sequencing issues, primer trimming, and bioinformatic pipeline optimizations for DENV, ZIKV, and CHIKV data.
