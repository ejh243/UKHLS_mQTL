## for neighbouring DNAm sites, test if evidence of colocalisation signals

library(coloc)
library(data.table)

setwd("")

epicManifest<-read.csv("MethylationEPIC_v-1-0_B2.csv", skip = 7)

cisProbes<-read.csv("MatrixEQTL/AllCisProbes_1e-10.csv", stringsAsFactors = FALSE, row.names = 1)
epicManifest<-epicManifest[match(cisProbes[,1], epicManifest$IlmnID),]
epicManifest<-epicManifest[order(epicManifest$CHR, epicManifest$MAPINFO),]
epicManifest$IlmnID<-as.character(epicManifest$IlmnID) 


output<-NULL
for( chr in c(1:22)){
	gwas<-NULL
	file_names<-list.files("MatrixEQTL/Output", pattern = paste("US_chr", chr, "_", sep = ""))
	file_names<-file_names[grep("cisProbe_allcis.txt", file_names)]
	for(each in file_names){
		gwas<-rbind(gwas, fread(paste("MatrixEQTL/Output/", each, sep = "")))
		}
	freq<-fread(paste("Genotypes/MatrixEQTL/US_eur_maf_geno_mind_hwe_chr", chr, "_freq.frq", sep = ""))
	freq$SNP<-gsub(":", ".", freq$SNP)
	freq$SNP<-paste(freq$SNP, freq$A1, sep = "_")
	for(i in which(epicManifest$CHR == chr)){
		## as going through in order only need to test probes downstream to avoid repetition.
		region.end<-epicManifest$MAPINFO[i]
		region.start<-region.end-250000
		pairs<-epicManifest$IlmnID[which(epicManifest$CHR == chr & epicManifest$MAPINFO < region.end & epicManifest$MAPINFO >= region.start)]
		if(length(pairs) > 0){
			trait1<-gwas[which(gwas$gene == epicManifest$IlmnID[i]),]
			dataset1<-list(beta= trait1$beta, varbeta=(trait1$beta/trait1$"t-stat")^2,type = "quant",snp = trait1$SNP, MAF = freq$MAF[match(trait1$SNP, freq$SNP)],N=1159)
			for(each in pairs){
				trait2<-gwas[which(gwas$gene == each),]
				dataset2<-list(beta= trait2$beta, varbeta=(trait2$beta/trait2$"t-stat")^2,type = "quant",snp = trait2$SNP, MAF = freq$MAF[match(trait2$SNP, freq$SNP)],N=1159)
				my.res<-coloc.abf(dataset1, dataset2)
				output<-rbind(output, c("trait1" = as.character(epicManifest$IlmnID[i]), "trait2"=each, my.res$summary))
			}
		}
	}
	write.csv(output, paste("BayesianColoc_PairwiseCpgs_Chr", chr, ".csv", sep = ""))
}