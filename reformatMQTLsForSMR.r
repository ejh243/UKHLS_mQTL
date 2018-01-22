### reformat mQTL allChr for smr software


setwd("")

tmp<-NULL
for(chr in 1:22){

	geno.map<-read.table(paste("Genotypes/US_eur_maf_geno_mind_hwe_chr", chr, ".bim", sep = ""), stringsAsFactors = FALSE)
	colnames(geno.map)<-c("Chr", "id", "cm", "bp", "A1", "A2")
	snps<-read.table(paste("Genotypes/US_Genotypes_Imputed_MapInfo_MinGenoCount5_Chr", chr, ".txt", sep = ""), stringsAsFactors = FALSE, header = TRUE)
	geno.map<-geno.map[match(snps$V1, geno.map$id),]
	tmp<-rbind(tmp, geno.map)

}

tmp$id<-gsub(":SNP", "", tmp$id)
tmp$id<-gsub(":INDEL", "", tmp$id)
tmp$id<-gsub("chr", "", tmp$id)
geno.map<-tmp

write.table(tmp, "SMR/US_Blood.esi", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


## create list of duplicate snp ids to exclude
dupIds<-tmp$id[duplicated(tmp$id)]
write.table(dupIds, "SMR/DuplicateVariantsToExclude.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)



probes<-read.table("Methylation/US_Methylation_Info.txt", stringsAsFactors = FALSE, header = TRUE)
epicManifest<-read.csv("MethylationEPIC_v-1-0_B2.csv", skip = 7)
epicManifest<-epicManifest[match(probes$IlmnID, epicManifest$IlmnID),]

probes.tmp<-cbind(epicManifest[,c("CHR", "IlmnID")],0,epicManifest$MAPINFO, epicManifest$UCSC_RefGene_Name, "+")

probes.tmp[,5]<-as.character(probes.tmp[,5])
probes.tmp[which(probes.tmp[,5] == ""),5]<-"NA"
write.table(probes.tmp, "SMR/US_Blood.epi", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

allChr<-NULL
for(chr in 1:22){
	allChr_file_name = paste("Output/US_chr", chr, ".txt", sep = "")
	tmp<-fread(allChr_file_name)
	tmp<-cbind(tmp, epicManifest[match(unlist(tmp[,2]), epicManifest$IlmnID),c("CHR", "MAPINFO", "UCSC_RefGene_Name", "UCSC_RefGene_Group")])
	
	locs.tmp<-unlist(lapply(strsplit(unlist(tmp[,1]), "\\_"), head, n = 1))
	locs.tmp<-gsub("\\.SNP", "", locs.tmp)
	locs.tmp<-gsub("\\.INDEL", "", locs.tmp)
	snplocs<-matrix(unlist(strsplit(locs.tmp, "\\.")), ncol = 2, byrow = TRUE)
	snplocs[,1]<-gsub("chr", "", snplocs[,1])
	tmp<-cbind(snplocs, tmp)
	tmp$CHR<-as.character(tmp$CHR)
	
	tmp<-tmp[which(tmp$V1 == tmp$CHR & abs(as.numeric(tmp$V2)-tmp$MAPINFO) < 500000),]
	allChr<-rbind(allChr, tmp)	
}
rm(tmp)
allChr$SNP<-gsub("chr", "", allChr$SNP)
allChr$SNP<-gsub("\\.", ":", allChr$SNP)
allChr$SNP<-gsub("_.", "", allChr$SNP)
### limit to SNPs on the same chromosome
allChr<-allChr[which(allChr$V1 == allChr$CHR),]

freq.blood<-NULL
for(chr in 1:22){
	freq.tmp<-read.table(paste("Genotypes/US_eur_maf_geno_mind_hwe_chr", chr, "_freq.frq", sep = ""), header = TRUE, stringsAsFactors = FALSE)
	freq.blood<-rbind(freq.blood, freq.tmp)
}
rm(freq.tmp)
freq.blood$SNP<-gsub(":SNP", "", freq.blood$SNP)
freq.blood$SNP<-gsub(":INDEL", "", freq.blood$SNP)	
freq.blood$SNP<-gsub("chr", "", freq.blood$SNP)
	
dbSNP<-read.table("Genotypes/US_eur_maf_geno_mind_hwe_withdbSNPIDs.txt", stringsAsFactors = FALSE)
dbSNP$V2<-gsub("chr", "", dbSNP$V2)
allChr<-cbind(dbSNP$V7[match(allChr$SNP, dbSNP$V2)],allChr)
colnames(allChr)[1]<-"rsID"

allChr<-cbind(allChr, freq.blood$MAF[match(paste(allChr$V1, allChr$V2, sep = ":"), freq.blood$SNP)])
colnames(allChr)[14]<-"MAF"
allChr<-cbind(allChr, geno.map[match(paste(allChr$V1, allChr$V2, sep = ":"), geno.map$id), c("A1", "A2")])
gz1 = gzfile("MatrixEQTL/Unrelated/US_EPIC_AnnotatedwithrsID.txt.gz","w");
write.table(allChr, sep = "\t", quote = FALSE, row.names = FALSE, gz1);
close(gz1) 
  
### create a file for every single probe and a reference file
for(each in unique(allChr$gene)){
	sub<-allChr[which(allChr$gene == each),]
	res<-cbind(sub$V1, paste(sub$V1, sub$V2, sep = ":"), sub$V2, freq.blood[match(paste(sub$V1, sub$V2, sep = ":"), freq.blood$SNP),c("A1", "A2", "MAF")], sub$"beta", sub$beta/sub$"t-stat", sub$"p-value")
	colnames(res)<-c("Chr","SNP","Bp","A1","A2","Freq","Beta","se","p")
	res<-res[!res[,2] %in% dupIds,]
	write.table(res, paste("SMR/mQTLResults/", each, ".esd", sep = ""), row.names = FALSE, quote = FALSE, sep = "\t")

}

probes<-unique(allChr$gene)
fileRef<-cbind(probes.tmp[match(probes, probes.tmp$IlmnID),], paste("SMR/mQTLResults/", probes, ".esd", sep = ""))
colnames(fileRef)<-c("Chr","ProbeID","GeneticDistance","ProbeBp","Gene","Orientation","PathOfEsd")
fileRef$Gene<-as.character(fileRef$Gene)
fileRef$Gene[which(fileRef$Gene == "")]<-"NA"
write.table(fileRef, "SMR/US_Blood.flist", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
