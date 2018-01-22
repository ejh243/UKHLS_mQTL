## create files for clumping

setwd("")
load("AllmQTL.rdata")

## only need to do for cpgs with > 1 SNP

nSNPs<-table(allChr$gene)
nSNPs<-nSNPs[which(nSNPs > 1)]

for(each in names(nSNPs)){
	allChr.sub<-allChr[which(allChr$gene == each),]
	allChr.sub<-allChr.sub[,c("SNP","p-value")]
	allChr.sub$SNP<-gsub("\\.", ":", unlist(lapply(strsplit(allChr.sub$SNP, "_"), head, n = 1)))
	colnames(allChr.sub)<-c("SNP", "P")
	write.table(allChr.sub, paste("ToClump/", each, ".assoc", sep = ""), quote = FALSE, row.names = FALSE)
}

fileRef<-paste("plink2 --bfile Genotypes/Imputed/data_filtered --clump ToClump/",names(nSNPs),".assoc --clump-p1 1e-8 --clump-p2 1e-8 --clump-r2 0.1 --clump-kb 250 --out ToClump/", names(nSNPs), sep = "")
write.table(fileRef,"ClumpingScript.sh", quote = FALSE, row.names = FALSE)