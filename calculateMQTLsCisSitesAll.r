
# takes chromosome specified on command line 
args <- commandArgs(trailingOnly = TRUE)
chr<-args[1]
nchunk<-args[2]

setwd("MatrixEQTL/")

library(MatrixEQTL)
base.dir = find.package("MatrixEQTL")
useModel = modelLINEAR; 

## threshold of results to save
pvOutputThreshold = 1;
cisDist<-500000
errorCovariance = numeric();

covariates_file_name = "Covariates/US.txt"; 
SNP_file_name = paste("Genotypes/US_Genotypes_Chr", chr, ".txt", sep = "")
snps_location_file_name = paste("Genotypes/US_Genotypes_Imputed_MapInfo_Chr",chr,".txt", sep = "")
expression_file_name<-paste("Methylation/US_Methylation_cisProbes_chr", chr, "_", nchunk, ".txt", sep = "")
output_file_name = paste("Output/US_chr", chr, "_", nchunk, "cisProbe_allcis.txt", sep = "")

##load covariate data
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
cvrt$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
cvrt$LoadFile(covariates_file_name);

## load SNP data
snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
snps$LoadFile( SNP_file_name );

## load methylation data
gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
gene$LoadFile( expression_file_name);

## load snp location data
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
geno.map<-read.table(paste("Genotypes/US_eur_maf_geno_mind_hwe_chr", chr, ".bim", sep = ""), stringsAsFactors = FALSE)
colnames(geno.map)<-c("Chr", "id", "cm", "bp", "A1", "A2")
geno.map<-geno.map[match(snpspos$V1, geno.map$id),]
geno.map$id<-paste(gsub(":", ".", geno.map$id), geno.map$A1, sep = "_")
snpspos$V1<-geno.map$id

## load dnam site location data
gene_loc = read.table(paste("Methylation/US_Methylation_cisProbes_chr", chr, "_", nchunk, "_Info.txt",sep = ""), header = TRUE, stringsAsFactors = FALSE)
gene_loc<-cbind(gene_loc, gene_loc[,3])


## Run the analysis
me = Matrix_eQTL_main(
	snps = snps,
	gene = gene,
	cvrt = cvrt,
	pvOutputThreshold = 0,
	output_file_name.cis = output_file_name,
	pvOutputThreshold.cis = pvOutputThreshold,
	snpspos = snpspos, 
	genepos = gene_loc,
	cisDist = cisDist,
	useModel = useModel, 
	errorCovariance = errorCovariance, 
	verbose = TRUE,
	pvalue.hist = "qqplot",
	min.pv.by.genesnp = FALSE,
	noFDRsaveMemory = FALSE)

