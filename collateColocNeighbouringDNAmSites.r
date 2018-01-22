geneOverlap<-function(term){
	if(term[1] == "" | term[2] == ""){
		return(NA)
	} else {
		genes1<-strsplit(term[1], ";")
		genes2<-strsplit(term[2], ";")
		overlap<-intersect(genes1, genes2)
		if(length(overlap) > 0){
			return(1)
		} else {
			return(0)
		}
	}
}

setwd("")

epicManifest<-read.csv("MethylationEPIC_v-1-0_B2.csv", skip = 7)

output.all<-NULL
for(chr in 1:22){
	output<-read.csv(paste("BayesianColoc_PairwiseCpgs_Chr", chr, ".csv", sep = ""), row.names = 1, stringsAsFactors = FALSE)
	output.all<-rbind(output, output.all)
}

sum<-output.all[,7]+output.all[,8]
ratio<-output.all[,8]/output.all[,7]

output.all<-cbind(output.all, sum ,ratio)

output<-cbind(output.all, epicManifest[match(output.all$trait1, epicManifest$IlmnID), c("CHR", "MAPINFO", "UCSC_RefGene_Name", "UCSC_RefGene_Group")], epicManifest[match(output.all$trait2, epicManifest$IlmnID), c("CHR", "MAPINFO", "UCSC_RefGene_Name", "UCSC_RefGene_Group")])

rm(output.all)

dist<-abs(output[,12]-output[,16])
output<-cbind(output, dist)


output$UCSC_RefGene_Name<-as.character(output$UCSC_RefGene_Name)
output$UCSC_RefGene_Name.1<-as.character(output$UCSC_RefGene_Name.1)

geneMatch<-apply(output[,c(13,17)], 1, geneOverlap)

fisher.test(table(geneMatch, sum > 0.99 & ratio > 1))
fisher.test(table(geneMatch, sum > 0.99 & ratio > 5))

cpgMatch<-rep(NA, nrow(output))
cpgMatch[which(epicManifest[match(output$trait1, epicManifest$IlmnID),"HMM_Island"] != "" & epicManifest[match(output$trait1, epicManifest$IlmnID),"HMM_Island"] != "")]<-0
cpgMatch[which(
epicManifest[match(output$trait1, epicManifest$IlmnID),"HMM_Island"] == epicManifest[match(output$trait2, epicManifest$IlmnID),"HMM_Island"])]<-1

fisher.test(table(cpgMatch, sum > 0.99 & ratio > 1))
fisher.test(table(cpgMatch, sum > 0.99 & ratio > 5))


regMatch<-rep(NA, nrow(output))
regMatch[which(epicManifest[match(output$trait1, epicManifest$IlmnID),"Regulatory_Feature_Name"] != "" & epicManifest[match(output$trait1, epicManifest$IlmnID),"Regulatory_Feature_Name"] != "")]<-0
regMatch[which(
epicManifest[match(output$trait1, epicManifest$IlmnID),"Regulatory_Feature_Name"] == epicManifest[match(output$trait2, epicManifest$IlmnID),"Regulatory_Feature_Name"])]<-1

fisher.test(table(regMatch, sum > 0.99 & ratio > 1))
fisher.test(table(regMatch, sum > 0.99 & ratio > 5))

write.csv(output[which(sum > 0.99 & ratio > 1),], "SupplementaryTable3.csv")