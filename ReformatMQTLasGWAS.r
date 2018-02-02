
epicManifest<-read.csv(MethylationEPIC_v-1-0_B2.csv", skip = 7)

setwd("")
cisProbes<-read.csv("MatrixEQTL/AllCisProbes_1e-10.csv", stringsAsFactors = FALSE, row.names = 1)
epicManifest<-epicManifest[match(cisProbes[,1], epicManifest$IlmnID),]
epicManifest<-epicManifest[order(epicManifest$CHR, epicManifest$MAPINFO),]
epicManifest$IlmnID<-as.character(epicManifest$IlmnID) 

for(chr in 1:22){
  gwas<-NULL
  file_names<-list.files("MatrixEQTL/Output", pattern = paste("US_chr", chr, "_", sep = ""))
  file_names<-file_names[grep("withP10PCs_cisProbe_allcis.txt", file_names)]
  for(each in file_names){
      gwas<-rbind(gwas, fread(paste("MatrixEQTL/Output/", each, sep = "")))
    }
  bim<-fread(paste("Genotypes/MatrixEQTL/US_eur_maf_geno_mind_hwe_chr", chr, ".bim", sep = ""))
  bim$V2<-gsub(":", ".", bim$V2)
  bim$V2<-paste(bim$V2, bim$V5, sep = "_")
  bim<-bim[match(gwas$SNP, bim$V2),]

  freq<-fread(paste("Genotypes/MatrixEQTL/US_eur_maf_geno_mind_hwe_chr", chr, "_freq.frq", sep = ""))
  freq$SNP<-gsub(":", ".", freq$SNP)
  freq$SNP<-paste(freq$SNP, freq$A1, sep = "_")
  freq<-freq[match(gwas$SNP, freq$SNP),]

  locs.tmp<-unlist(lapply(strsplit(unlist(gwas[,1]), "\\_"), head, n = 1))
  locs.tmp<-gsub("\\.SNP", "", locs.tmp)
  locs.tmp<-gsub("\\.INDEL", "", locs.tmp)
  snplocs<-matrix(unlist(strsplit(locs.tmp, "\\.")), ncol = 2, byrow = TRUE)
  snplocs[,1]<-gsub("chr", "", snplocs[,1])

  gwas<-cbind(gwas, snplocs, bim, freq)

  for(i in which(epicManifest$CHR == chr)){
    trait1<-gwas[which(gwas$gene == epicManifest$IlmnID[i]),]
    out<-cbind(paste(trait1$V1, trait1$V2, sep = ":"), 
      trait1[,c("V5", "V6", "MAF", "beta")],trait1$beta/trait1$"t-stat", 
      trait1$"p-value", "1159")
    colnames(out)<-c("SNP", "A1", "A2", "freq", "b", "se", "p", "n")


    write.table(out, paste("MatrixEQTL/SMR/mQTL_GWASSumStats/", epicManifest$IlmnID[i], "_SMRFormat.txt", sep = ""), 
      quote = FALSE, row.names = FALSE)
  }

}

for(chr in 1:22){
fileRef<-"cd XXXX"
	for(i in which(epicManifest$CHR == chr)){
		fileRef<-rbind(fileRef, paste("smr_Linux --bfile refDataForSMR/WesteraeQTLOverlap 
      --gwas-summary MatrixEQTL/SMR/mQTL_GWASSumStats/", epicManifest$IlmnID[i], "_SMRFormat.txt 
      --beqtl-summary Westera_1000GFormat 
      --out MatrixEQTL/SMR/mQTL_eQTL/",epicManifest$IlmnID[i]," --thread-num 10", 
      sep = ""))
	}
	write.table(fileRef, paste("MatrixEQTL/SMR/Run_mQTL_eQTL_chr", chr, ".sh", sep = ""), quote = FALSE, row.names = FALSE)
}
