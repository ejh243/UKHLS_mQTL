# UKHLS_mQTL

This repository contains scripts used to perform the analyses reported in the manuscript Hannon et al which generated a dataset of DNA methylation quantitative trait loci and performed a series of downstream analyses. This readme provides details on which analyses each script describes. Filepaths have been removed and therefore these serve as a guide to how the analyses were performed. 

As the quality control and data preprocessing has been described elsewhere these scripts are not provided in this repo.

1. calculateMQTLs.r
This is an R script to load and test for  associations between all pairs of DNA methylation sites and genetic variants using the [MatrixEQTL package](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/). This script was produced with the help of the tutorials on the MatrixEQTL website. It implements a linear additive model and includes relevant covariates. The analysis has been streamlined by running the genetic variants on each chromosome separately. The chromosome is specified on the command line at execution e.g. to run the analysis for all variants on chromosome 1 you would execute. 

```bash
Rscript calculateMQTLs.r 1
```

2. collateMQTLOutput.r
This collates the QTL analyses run for each chromosome and annotates the DNA methylation sites with annotation provided by Illumina. It creates list of all DNA methylation sites with at least one cis QTL for other downstream analyses. It creates a table of number of QTL associations and characterizes heir effects on DNA  methylation for multiple P value thresholds for defining QTLs (Supplementary Table 1). Finally it creates a table summarizing the number of QTLs identified due to novel content on the EPIC array and shared with the 450K array (Supplementary Table 2).

3. createFilesToClumpMQTLs.r
Takes set of significant mQTLs from collateMQTLOutput.r and for all DNAm site with > 1 association creates GWAS files for use with PILNK's clumping function. Also creates script to execute clumping commands for each DNAm site.

4. collateClumpedQTLs.r
Collates clumped results for each DNAm site produced by PLINK and combines with DNAm sites that had only 1 association. 

5. groupMQTLs.r
Takes set of significant mQTLs from collateMQTLOutput.r and formulates groups of DNAm sites which are all associated with the same set of genetic variants and then groups these groups so that all genetic variants associated with the same set of DNAm sites are grouped together.  

6. calculateMQTLsCisSitesAll.r
This is an adaptation of calculateMQTLs.r to rerun the QTL analysis for all DNA sites with at least one significant cis QTL and save all cis P values for all genetic variants. It is again split by chromosome and then chunks within each chromosome. The chromosome and chunk number is specified on the command line at execution e.g. to run the analysis for all variants on chromosome 1 and DNAm sites in chunk 1 you would execute. 
```bash
Rscript calculateMQTLsCisSitesAll.r 1 1
```

7. colocOfNeighbouringDNAmSites.r
Takes all DNAm sites associated with at least one QTL and tests all pairs of DNAm sites located on the same chromosome and within 250kb of each other for co-localisation of genetic signals. Requires p values of all genetic variants in region (rather than just those significant) so need to rerun script calculateMQTLs.r to save all p values for these DNAm sites first. 

8. collateColocNeighbouringDNAmSites.r
Collates results from colocOfNeighbouringDNAmSites.r across chromosomes and performs enrichment tests across genic, CpG island and regulatory features. Produces table of convincing colocalising pairs of DNAm sites (Supplementary Table 3).

9. reformatMQTLsForSMR.r 
Takes set of significant mQTLs from collateMQTLOutput.r and creates files for SMR analysis as described [here]()


mQTL results were formatted for SMR using the following command:

```{bash}
smr_Linux --eqtl-flist US_Blood.flist --make-besd --out US_Blood 
```
SMR analyses for each trait were executed with the following command:

```{bash}
smr_Linux --bfile data_filtered_ForSMR --gwas-summary <GWAS file> --beqtl-summary US_Blood --out SMR_US_Blood_<GWAS trait> --thread-num 10;
```

10. collateSMRResultsAll.r
Collates output of SMR analysis and apply significance filters (Supplementary Table 5). Performs enrichment of location of associated DNAm sites in gens, regulatory features. Collates results of SMR analysis with blood eQTLs (using Westera et al data) and identifies overlap (Supplementary Table 6). 