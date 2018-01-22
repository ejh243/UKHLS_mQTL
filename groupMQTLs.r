
charProbeDist<-function(term){
	probes<-unlist(strsplit(term, ";"))
	locs<-cpgAnno[(cpgAnno$gene %in% probes),]
	if(length(unique(locs$CHR)) > 1){
		diffChr<-1
		return(c(1, NA, NA, NA))
	} else {
		return(c(0, min(locs$MAPINFO), max(locs$MAPINFO), (max(locs$MAPINFO)-min(locs$MAPINFO))/nrow(locs)))
	}
}

library(dplyr)

setwd("")
load("AllmQTL.rdata")
allChr<-allChr[which(allChr$CHR != "X" & allChr$CHR != "Y"),]
allChr$V2<-as.numeric(allChr$V2)
cisInd<-(allChr$V1 == allChr$CHR & abs(allChr$V2 - allChr$MAPINFO)< 500000)
allChr<-cbind(allChr, cisInd)


cpgGroups<-NULL
for(chr in 1:22){
	allChr.sub<-allChr[which(allChr$V1 == chr),]
	tmpAg<-aggregate(allChr.sub$gene, by = list(allChr.sub$SNP), paste, collapse = ";")
	tmpAg2<-aggregate(tmpAg$Group.1, by = list(tmpAg$x), paste, collapse = ";")
	cpgGroups<-rbind(cpgGroups, tmpAg2)
}

nCPGs<-unlist(lapply(strsplit(cpgGroups[,1], ";"), length))
nSNPs<-unlist(lapply(strsplit(cpgGroups[,2], ";"), length))

cpgAnno<-unique(allChr[,c("gene", "CHR", "MAPINFO")])

data.probes<-sapply(cpgGroups[,1], charProbeDist)
data.probes<-t(data.probes)

save(data.probes, cpgGroups, nCPGs, nSNPs, file = "cpgGroups.rda")

