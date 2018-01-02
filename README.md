# UKHLS_mQTL

This repository contains scripts used to perform the analyses reported in the manuscript Hannon et al which generted a dataset of DNA methylation quantitative trait loci and performed a series of downstream analyses. This readme provides details on which analyses each script descrbes. Filepaths have been removed and therefore these serve as a guide to how the analyses were performed. 

As the quality control and data preprocessing has been described elsewhere these scripts are not provided in this repo.

1. calculateMQTLs.r
This is an R script to load and test for  asocitions between all pairs of DNA methylation sites and genetic variants using the [MatrixEQTL package](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/). This script was produced with the help of the tutorials on the MatrixEQTL website. It implements a linear additive model and includes relevant covariates. The analysis has been streamlined by running the genetic variants on each chromosome separately. The specific chromosome is specific on the command line at execution e.g. to runt the analsis for all variants on chromosome 1 you would execute. 

```bash
Rscript calculateMQTLs.r 1
```