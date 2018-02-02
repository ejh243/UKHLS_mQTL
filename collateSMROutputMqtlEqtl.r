### collate output from SMR of mQTL against eQTL


## function to summaraize location of multriple DNAm sites assoicated with the same gene
charProbeDist<-function(term, cpgAnno){
	probes<-unlist(strsplit(term, ";"))
	locs<-cpgAnno[(cpgAnno$each %in% probes),]
	if(length(unique(locs$CHR)) > 1){
		diffChr<-1
		return(c(1, NA, NA, NA))
	} else {
		if(nrow(locs) > 1){
			return(c(0, min(locs$MAPINFO), max(locs$MAPINFO), (max(locs$MAPINFO)-min(locs$MAPINFO))/(nrow(locs)-1)))
		} else {
			return(c(0, locs$MAPINFO, locs$MAPINFO, NA))
		
		}
	}
}

epicManifest<-read.csv("MethylationEPIC_v-1-0_B2.csv", skip = 7)

setwd("")
cisProbes<-read.csv("MatrixEQTL/AllCisProbes_1e-10.csv", stringsAsFactors = FALSE, row.names = 1)
epicManifest<-epicManifest[match(cisProbes[,1], epicManifest$IlmnID),]


setwd("MatrixEQTL/SMR/mQTL_eQTL")

allFiles<-list.files()
allFiles<-allFiles[grep("smr", allFiles)]
out<-NULL
for(each in allFiles){
	tmp<-read.table(each, stringsAsFactors = FALSE, header = TRUE)
	if(nrow(tmp) > 0){
		out<-rbind(out, cbind(each, tmp))
	}
}

out$each<-gsub("\\.smr", "", out$each)
out<-cbind(out, epicManifest[match(out$each, epicManifest$IlmnID),
  c("CHR", "MAPINFO", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_UCSC_CpG_Island", "Regulatory_Feature_Name", 
  "Regulatory_Feature_Group", "DNase_Hypersensitivity_NAME", "OpenChromatin_NAME", "TFBS_Evidence_Count", "Methyl450_Loci", 
  "SNP_ID", "SNP_DISTANCE", "SNP_MinorAlleleFrequency")])

save(out, file = "AllmQTL_eQTL_ColocResults.rda")

table(table(out$each))
table(table(out$ProbeID))
par(mfrow = c(1,3))
hist(table(out$each), xlab = "Number of times DNAm site tested")
hist(table(out$ProbeID), xlab = "Number of times transcript tested")
hist(table(out$Gene), xlab = "Number of times gene tested")

write.csv(out, "../AllTested_mQTL_eQTL_pairs.csv")

## how many genes each DNAm site tested against?
summary(as.numeric(table(unlist(lapply(strsplit((unique(paste(out$each, out$Gene))), " "), head, n = 1)))))
sd(as.numeric(table(unlist(lapply(strsplit((unique(paste(out$each, out$Gene))), " "), head, n = 1)))))
summary(as.numeric(table(out$each)))
sd(as.numeric(table(out$each)))
summary(as.numeric(table(out$ProbeID)))
sd(as.numeric(table(out$ProbeID)))

write.csv(out[which(as.numeric(out$p_SMR) < 1.02e-07),], "../BonfSignif_mQTL_eQTL_pairs.csv")
## Supp Table 7
write.csv(out[which(as.numeric(out$p_SMR) < 1.02e-07 & as.numeric(out$p_HET) > 0.05),],"../BonfSignif_HetP_mQTL_eQTL_pairs.csv")

length(unique(out$Gene[which(as.numeric(out$p_SMR) < 1.02e-07)]))/length(unique(out$Gene))
length(unique(out$ProbeID[which(as.numeric(out$p_SMR) < 1.02e-07)]))/length(unique(out$ProbeID))
length(unique(out$each[which(as.numeric(out$p_SMR) < 1.02e-07)]))/length(unique(out$each))


signifPairs<-out[which(as.numeric(out$p_SMR) < 1.02e-07 & as.numeric(out$p_HET) > 0.05),]

table(table(signifPairs$each))
table(table(signifPairs$ProbeID))
par(mfrow = c(1,3))
hist(table(signifPairs$each), xlab = "Number of times DNAm site tested")
hist(table(signifPairs$ProbeID), xlab = "Number of times transcript tested")
hist(table(signifPairs$Gene), xlab = "Number of times gene tested")

summary(as.numeric(table(unlist(lapply(strsplit((unique(paste(signifPairs$each, signifPairs$Gene))), " "), head, n = 1)))))
sd(as.numeric(table(unlist(lapply(strsplit((unique(paste(signifPairs$each, signifPairs$Gene))), " "), head, n = 1)))))
summary(as.numeric(table(signifPairs$each)))
sd(as.numeric(table(signifPairs$each)))
summary(as.numeric(table(signifPairs$ProbeID)))
sd(as.numeric(table(signifPairs$ProbeID)))

length(unique(signifPairs$Gene))/length(unique(out$Gene))
length(unique(signifPairs$ProbeID))/length(unique(out$ProbeID))
length(unique(signifPairs$each))/length(unique(out$each))

geneFeatures<-unique(unlist(strsplit(as.character(out$UCSC_RefGene_Group), ";")))

backgroundAnno<-unique(out[,c("each", "CHR", "MAPINFO", "UCSC_RefGene_Name", "UCSC_RefGene_Group", 
  "Relation_to_UCSC_CpG_Island", "Regulatory_Feature_Name", "Regulatory_Feature_Group", "DNase_Hypersensitivity_NAME", 
  "OpenChromatin_NAME", "TFBS_Evidence_Count", "Methyl450_Loci", "SNP_ID", "SNP_DISTANCE", "SNP_MinorAlleleFrequency")])

signifAnno<-unique(signifPairs[,c("each", "CHR", "MAPINFO", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_UCSC_CpG_Island", 
  "Regulatory_Feature_Name", "Regulatory_Feature_Group", "DNase_Hypersensitivity_NAME", "OpenChromatin_NAME", "TFBS_Evidence_Count", 
  "Methyl450_Loci", "SNP_ID", "SNP_DISTANCE", "SNP_MinorAlleleFrequency")])

tab<-matrix(data = NA, nrow = length(geneFeatures)+1, ncol = 4)
rownames(tab)<-c(geneFeatures, "Intergenic")
colnames(tab)<-c("SignifPairs", "Background", "%SignifPairs", "%Background")
for(each in geneFeatures){
	tab[each, 1]<-length(grep(each, signifAnno$UCSC_RefGene_Group))
	tab[each, 2]<-length(grep(each, backgroundAnno$UCSC_RefGene_Group))
	
}

tab["Intergenic",1]<-length(which(signifAnno$UCSC_RefGene_Group == ""))
tab["Intergenic",2]<-length(which(backgroundAnno$UCSC_RefGene_Group == ""))

tab[,3]<-tab[,1]/length(unique(signifAnno$each))*100
tab[,4]<-tab[,2]/length(unique(backgroundAnno$each))*100

par(mfrow = c(1,1))
barplot(t(tab[,3:4]), beside = TRUE, col = c("blue", "lightgray"), ylab = "% annotated to gene feature")
legend("topright", c("Significant DNAm sites", "Background"), pch = 15, col = c("blue", "lightgray"))

vioData<-NULL
for(each in geneFeatures){
	tab[each, 1]<-length(grep(each, signifAnno$UCSC_RefGene_Group))
	tab[each, 2]<-length(grep(each, backgroundAnno$UCSC_RefGene_Group))
	
}

sepByFeature<-function(term, signifPairs){
	return(signifPairs$b_SMR[grep(term, signifPairs$UCSC_RefGene_Group)])

}

vioData<-sapply(geneFeatures, sepByFeature, signifPairs)
vioData[[8]]<-signifPairs$b_SMR[which(signifPairs$UCSC_RefGene_Group == "")]
vioplot(vioData[[1]], vioData[[2]], vioData[[3]], vioData[[4]], vioData[[5]], vioData[[6]], vioData[[7]], vioData[[8]], names = c(geneFeatures, "Intergenic"))

vioData.abs<-lapply(vioData, abs)
vioplot(vioData.abs[[1]], vioData.abs[[2]], vioData.abs[[3]], vioData.abs[[4]], vioData.abs[[5]], vioData.abs[[6]], vioData.abs[[7]], vioData.abs[[8]], names = c(geneFeatures, "Intergenic"))

annoMatch<-rep(NA, nrow(signifPairs))
for(i in 1:nrow(signifPairs)){
	if(length(grep(signifPairs$Gene[i], signifPairs$UCSC_RefGene_Name[i])) > 0){
		annoMatch[i]<-1
	} else{
		annoMatch[i]<-0
	}
}

## how many are not annotated to the same gene but are tested against
annoMatch.all<-rep(NA, nrow(out))
for(i in 1:nrow(out)){
	if(length(grep(out$Gene[i], out$UCSC_RefGene_Name[i])) > 0){
		annoMatch.all[i]<-1
	} else{
		annoMatch.all[i]<-0
	}
}

## are cpgs annotated to the same gene more likely to be associated?
sigInd<-(as.numeric(out$p_SMR) < 1.02e-07 & as.numeric(out$p_HET) > 0.05)
tab<-table(sigInd, annoMatch.all)

## are significant associations between different genes?
sig.opp<-out[which(sigInd == 1 & annoMatch.all == 0),]
sig.same<-out[which(sigInd == 1 & annoMatch.all == 1),]
length(intersect(unique(sig.opp$each), unique(sig.same$each))) 
table(unique(sig.same$each) %in% unique(sig.opp$each))
table(unique(sig.opp$each) %in% unique(sig.same$each))

## for cpgs associated with a diff gene, were they tested against themselves?

table(unique(sig.opp$each) %in% unique(out$each[which(annoMatch.all == 1)]))
table(unique(sig.opp$each) %in% unique(sig.same$each))
table(unique(sig.opp$each) %in% unique(out$each[which(annoMatch.all == 1)]))

## where multiple sites associated with the same transcript characterise the location/distribution of CpGs
tmpAg<-aggregate(signifPairs$each, by = list(signifPairs$ProbeID), paste, collapse = ";")
cpgAnno<-unique(signifPairs[,c("each", "CHR", "MAPINFO")])
nCPGs<-unlist(lapply(strsplit(tmpAg[,2], ";"), length))

data.probes<-sapply(tmpAg[,2], charProbeDist, cpgAnno)
data.probes<-t(data.probes)

par(mfrow = c(1,2))
hist(data.probes[which(nCPGs > 1),3]-data.probes[which(nCPGs > 1),2], xlab = "Range", breaks = 100)
mtext("A", line = 0.5, side = 3, adj = 0)
hist(data.probes[which(nCPGs > 1),4], xlab = "Density", breaks = 100)
mtext("B", line = 0.5, side = 3, adj = 0)

summary(data.probes[which(nCPGs > 1),3]-data.probes[which(nCPGs > 1),2])
sd(data.probes[which(nCPGs > 1),3]-data.probes[which(nCPGs > 1),2])
summary(data.probes[which(nCPGs > 1),4])
sd(data.probes[which(nCPGs > 1),4])
