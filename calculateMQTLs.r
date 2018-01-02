# takes chromosome specified on command line 
args <- commandArgs(trailingOnly = TRUE)
chr<-args

setwd("")
library(MatrixEQTL)
base.dir = find.package("MatrixEQTL")
useModel = modelLINEAR;  
pvOutputThreshold = 1e-8; ## set p value threshold of results to save
errorCovariance = numeric();

## specify filenames and paths
covariates_file_name = "US.txt"; 
SNP_file_name = paste("US_Genotypes_Imputed_Chr", chr, ".txt", sep = "")
expression_file_name<-"US_Methylation.txt"
output_file_name = paste("Output/US_chr", chr, ".txt", sep = "")

##load covariate data
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      
cvrt$fileOmitCharacters = "NA"; 
cvrt$fileSkipRows = 1;          
cvrt$fileSkipColumns = 1;       
cvrt$fileSliceSize = 2000;      
cvrt$LoadFile(covariates_file_name);

## load SNP data
snps = SlicedData$new();
snps$fileDelimiter = "\t";      
snps$fileOmitCharacters = "NA"; 
snps$fileSkipRows = 1;          
snps$fileSkipColumns = 1;       
snps$fileSliceSize = 2000;      
snps$LoadFile(SNP_file_name );

## load methylation data
gene = SlicedData$new();
gene$fileDelimiter = "\t";      
gene$fileOmitCharacters = "NA"; 
gene$fileSkipRows = 1;         
gene$fileSkipColumns = 1;       
gene$fileSliceSize = 2000;      
gene$LoadFile(expression_file_name);

## Run the analysis

me = Matrix_eQTL_main(
	snps = snps,
	gene = gene,
	cvrt = cvrt,
	output_file_name = output_file_name,
	pvOutputThreshold = pvOutputThreshold,
	useModel = useModel, 
	errorCovariance = errorCovariance, 
	verbose = TRUE,
	pvalue.hist = "qqplot",
	min.pv.by.genesnp = FALSE,
	noFDRsaveMemory = FALSE)

