## collate results from US mQTL analysis
library(data.table)
setwd("MatrixEQTL/")

epicManifest<-read.csv("MethylationEPIC_v-1-0_B2.csv", skip = 7)

allChr<-NULL
nprobes.cis<-0
for(chr in 1:22){
	output_file_name = paste("Output/US_chr", chr, ".txt", sep = "")
	tmp<-fread(output_file_name)
	tmp<-cbind(tmp, epicManifest[match(unlist(tmp[,2]), epicManifest$IlmnID),c("CHR", "MAPINFO", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_UCSC_CpG_Island", "Regulatory_Feature_Name", "Regulatory_Feature_Group", "DNase_Hypersensitivity_NAME", "OpenChromatin_NAME", "TFBS_Evidence_Count", "Methyl450_Loci", "SNP_ID", "SNP_DISTANCE", "SNP_MinorAlleleFrequency")])
	
	locs.tmp<-unlist(lapply(strsplit(unlist(tmp[,1]), "\\_"), head, n = 1))
	locs.tmp<-gsub("\\.SNP", "", locs.tmp)
	locs.tmp<-gsub("\\.INDEL", "", locs.tmp)
	snplocs<-matrix(unlist(strsplit(locs.tmp, "\\.")), ncol = 2, byrow = TRUE)
	snplocs[,1]<-gsub("chr", "", snplocs[,1])
	tmp<-cbind(snplocs, tmp)
	tmp$CHR<-as.character(tmp$CHR)
	
	nprobes.cis<-nprobes.cis+length(unique(tmp$gene[which(tmp$V1 == tmp$CHR & abs(as.numeric(tmp$V2)-tmp$MAPINFO) < 500000)]))
	allChr<-rbind(allChr, tmp)	
}
	
rm(tmp)


cisProbes<-allChr$gene[which(allChr$V1 == allChr$CHR & abs(as.numeric(allChr$V2)-allChr$MAPINFO) < 500000)]
cisProbes<-unique(cisProbes)
write.csv(cisProbes, "AllCisProbes_1e-10.csv")

allChr$V2<-as.numeric(allChr$V2)
save(allChr, file = "AllmQTL.rdata")


pThres<-c(5e-8/783967, 1e-13, 1e-12, 1e-11, 1e-10)
tab<-matrix(data = NA, nrow = length(pThres), ncol = 26)
colnames(tab)<-c("Threshold", "nMQTL", "nSNPs", "nProbes", "mean effect", "SD effect", "nSameChr","nSNPs", "nProbes","mean effect", "SD effect", "nwithin1MB","nSNPs", "nProbes","mean effect", "SD effect", "nwithin500kb","nSNPs", "nProbes","mean effect", "SD effect", "ntrans","nSNPs", "nProbes", "mean effect", "SD effect")
tab[,1]<-pThres

## Filter to autosomes
allChr<-allChr[which(allChr$CHR != "X" & allChr$CHR != "Y"),]

pThres<-c(5e-8/766714, 1e-13, 1e-12, 1e-11, 1e-10)
tab<-matrix(data = NA, nrow = length(pThres), ncol = 26)
colnames(tab)<-c("Threshold", "nMQTL", "nSNPs", "nProbes", "mean effect", "SD effect", "nSameChr","nSNPs", "nProbes","mean effect", "SD effect", "nwithin1MB","nSNPs", "nProbes","mean effect", "SD effect", "nwithin500kb","nSNPs", "nProbes","mean effect", "SD effect", "ntrans","nSNPs", "nProbes", "mean effect", "SD effect")
tab[,1]<-pThres
for(i in 1:length(pThres)){
	allChr.sub<-allChr[which(allChr$"p-value" < pThres[i]),]
	tab[i,2]<-nrow(allChr.sub)
	tab[i,3]<-length(unique(allChr.sub$SNP))
	tab[i,4]<-length(unique(allChr.sub$gene))
	tab[i,5]<-mean(abs(allChr.sub$beta))
	tab[i,6]<-sd(abs(allChr.sub$beta))
	tab[i,7]<-length(which(allChr.sub$V1 == allChr.sub$CHR))
	tab[i,8]<-length(unique(allChr.sub$SNP[which(allChr.sub$V1 == allChr.sub$CHR)]))
	tab[i,9]<-length(unique(allChr.sub$gene[which(allChr.sub$V1 == allChr.sub$CHR)]))
	tab[i,10]<-mean(abs(allChr.sub$beta[which(allChr.sub$V1 == allChr.sub$CHR)]))
	tab[i,11]<-sd(abs(allChr.sub$beta[which(allChr.sub$V1 == allChr.sub$CHR)]))	

	tab[i,12]<-length(which(allChr.sub$V1 == allChr.sub$CHR & abs(allChr.sub$V2 - allChr.sub$MAPINFO)< 1000000))
	tab[i,13]<-length(unique(allChr.sub$SNP[which(allChr.sub$V1 == allChr.sub$CHR & abs(allChr.sub$V2 - allChr.sub$MAPINFO)< 1000000)]))
	tab[i,14]<-length(unique(allChr.sub$gene[which(allChr.sub$V1 == allChr.sub$CHR & abs(allChr.sub$V2 - allChr.sub$MAPINFO)< 1000000)]))
	tab[i,15]<-mean(abs(allChr.sub$beta[which(allChr.sub$V1 == allChr.sub$CHR & abs(allChr.sub$V2 - allChr.sub$MAPINFO)< 1000000)]))
	tab[i,16]<-sd(abs(allChr.sub$beta[which(allChr.sub$V1 == allChr.sub$CHR & abs(allChr.sub$V2 - allChr.sub$MAPINFO)< 1000000)]))	

	tab[i,17]<-length(which(allChr.sub$V1 == allChr.sub$CHR & abs(allChr.sub$V2 - allChr.sub$MAPINFO)< 500000))
	tab[i,18]<-length(unique(allChr.sub$SNP[which(allChr.sub$V1 == allChr.sub$CHR & abs(allChr.sub$V2 - allChr.sub$MAPINFO)< 500000)]))
	tab[i,19]<-length(unique(allChr.sub$gene[which(allChr.sub$V1 == allChr.sub$CHR & abs(allChr.sub$V2 - allChr.sub$MAPINFO)< 500000)]))
	tab[i,20]<-mean(abs(allChr.sub$beta[which(allChr.sub$V1 == allChr.sub$CHR & abs(allChr.sub$V2 - allChr.sub$MAPINFO)< 500000)]))
	tab[i,21]<-sd(abs(allChr.sub$beta[which(allChr.sub$V1 == allChr.sub$CHR & abs(allChr.sub$V2 - allChr.sub$MAPINFO)< 500000)]))
	tab[i,22]<-length(which(allChr.sub$V1 != allChr.sub$CHR | abs(allChr.sub$V2 - allChr.sub$MAPINFO)>= 500000))
	tab[i,23]<-length(unique(allChr.sub$SNP[which(allChr.sub$V1 != allChr.sub$CHR | abs(allChr.sub$V2 - allChr.sub$MAPINFO)>= 500000)]))
	tab[i,24]<-length(unique(allChr.sub$gene[which(allChr.sub$V1 != allChr.sub$CHR | abs(allChr.sub$V2 - allChr.sub$MAPINFO)>= 500000)]))
	tab[i,25]<-mean(abs(allChr.sub$beta[which(allChr.sub$V1 != allChr.sub$CHR | abs(allChr.sub$V2 - allChr.sub$MAPINFO)>= 500000)]))
	tab[i,26]<-sd(abs(allChr.sub$beta[which(allChr.sub$V1 != allChr.sub$CHR | abs(allChr.sub$V2 - allChr.sub$MAPINFO)>= 500000)]))
}

write.csv(tab, "SupplementryTable1.csv")

allChr$UCSC_RefGene_Name<-as.character(allChr$UCSC_RefGene_Name)

tab<-matrix(data = NA, nrow = length(pThres), ncol = 10)
colnames(tab)<-c("Threshold", "nProbes", "nProbes_EPICSpecific", "nGenes", "nGenes_EPICnovel", "nGenes_EPICreplicate", "nProbes_intergenic", "nProbes_EPICnovel_intergenic", "nEPICspecificwithin1Kbof450k", "nEPICspecificwithin5Kbof450k")
tab[,1]<-pThres
window<-1000
for(i in 1:length(pThres)){
	allChr.sub<-allChr[which(allChr$"p-value" < pThres[i]),]

	tab[i,2]<-length(unique(allChr.sub$gene))
	tab[i,3]<-length(unique(allChr.sub$gene[is.na(allChr.sub$Methyl450_Loci)]))	
	tab[i,4]<-length(unique(unlist(strsplit(allChr.sub$UCSC_RefGene_Name, ";"))))
	genes.epic<-unique(unlist(strsplit(allChr.sub$UCSC_RefGene_Name[is.na(allChr.sub$Methyl450_Loci)], ";")))
	genes.450k<-unique(unlist(strsplit(allChr.sub$UCSC_RefGene_Name[which(allChr.sub$Methyl450_Loci == TRUE)], ";")))
	
	tab[i,5]<-sum(!genes.epic %in% genes.450k)
	tab[i,6]<-length(intersect(genes.epic, genes.450k))
	tab[i,7]<-length(unique(allChr.sub$gene[which(allChr.sub$UCSC_RefGene_Name == "")]))
	tab[i,8]<-length(unique(allChr.sub$gene[which(allChr.sub$UCSC_RefGene_Name == "" & is.na(allChr.sub$Methyl450_Loci))]))
	
	pos.450k<-unique(allChr.sub[which(allChr.sub$Methyl450_Loci == TRUE), c("CHR", "MAPINFO")])
	pos.epic<-unique(allChr.sub[is.na(allChr.sub$Methyl450_Loci), c("CHR", "MAPINFO")])
	gr.450k <- GRanges(seqnames = paste("chr", unlist(pos.450k[,1]), sep = ""), ranges = IRanges(unlist(pos.450k[,2])-window, width = 1+window*2),strand = "+")
	gr.450k<-reduce(gr.450k)
	gr.epic <- GRanges(seqnames = paste("chr", unlist(pos.epic[,1]), sep = ""), ranges = IRanges(unlist(pos.epic[,2])-1, width = 2),strand = "+")
	length(gr.450k)
	length(gr.epic)
	tab[i,9]<-length(intersect(gr.450k, gr.epic))
	gr.450k <- GRanges(seqnames = paste("chr", unlist(pos.450k[,1]), sep = ""), ranges = IRanges(unlist(pos.450k[,2])-5000, width = 1+5000*2),strand = "+")
	gr.450k<-reduce(gr.450k)
	tab[i,10]<-length(intersect(gr.450k, gr.epic))
}
write.csv(tab, "SupplementaryTable2.csv")