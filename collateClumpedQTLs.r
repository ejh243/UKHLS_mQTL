## collate clumping results

setwd("")
load("AllmQTL.rdata")

## only need to do for cpgs with > 1 SNP
nSNPs<-table(allChr$gene)
nSNPs<-nSNPs[which(nSNPs > 1)]

clump.out<-NULL
for(each in names(nSNPs)){
	tmp<-read.table(paste("ToClump/", each, ".clumped", sep = ""), header = TRUE)
	tmp<-tmp[,1:11]
	clump.out<-rbind(clump.out, cbind(each, as.matrix(tmp)))

}

## add back in DNAm with only 1 genetic variant
nSNPs<-table(allChr$gene)
toAdd<-allChr[allChr$gene %in% names(nSNPs[which(nSNPs == 1)]),c("gene", "V1", "SNP", "V2", "p-value")]
toAdd<-cbind(toAdd[,1:2], "F"=NA, toAdd[,3:5], "TOTAL" = 1, "NSIG" = NA,  "S05" = NA, "S01" = NA, "S001" = NA, "S0001" = NA)
colnames(toAdd)<-colnames(clump.out)
toAdd<-as.matrix(toAdd)
clump.out<-rbind(clump.out, toAdd)
