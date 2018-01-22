uniqueAnno<-function(row){
if(row != ""){
	return(paste(unique(unlist(strsplit(row, "\\;"))), collapse = ";"))
	} else {
	return(row)
	}
}
	
signif<-function(table, pval){
	return(table[which(table[,19] < pval),])
}

hetP<-function(table, pval = 0.05){
	return(table[which(table[,20] > pval),])
}

library(data.table)
setwd("")

blood.probes<-read.table("SMR/US_Blood.flist", stringsAsFactors = FALSE, header = TRUE)

AIR<-cbind(AIR, trait = "Acute insulin response")
AIR<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_AIR.smr", stringsAsFactors = FALSE, header = TRUE)
AIRadjBMISI<-cbind(AIRadjBMISI, trait = "Acute insulin response adj SI, BMI")
AIRadjBMISI<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_AIR_adj_BMI+SI.smr", stringsAsFactors = FALSE, header = TRUE)
AIRAdjSII<-cbind(AIRAdjSII, trait = "Acute insulin response adj SI")
AIRAdjSII<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_AIR_adj_SII.smr", stringsAsFactors = FALSE, header = TRUE)
als<-cbind(als, trait = "Amyloid lateral sclerosis")
als<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_ALS_lmm_GWAS.smr", stringsAsFactors = FALSE, header = TRUE)
apop_derm<-cbind(apop_derm, trait = "Atopic dermatitis (eczema)")
apop_derm<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_EAGLE_AD.smr", stringsAsFactors = FALSE, header = TRUE)
bmi<-cbind(bmi, trait = "BMI")
bmi<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_BMI_GWAS.smr", stringsAsFactors = FALSE, header = TRUE)
bpd<-cbind(bpd, trait = "Bipolar disorder")
bpd<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_BPD_GWAS.smr", stringsAsFactors = FALSE, header = TRUE)
cardio2<-cbind(cardio, trait = "Coronary artery disease")
cardio2<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_Cardiogram_2_GWAS.smr", stringsAsFactors = FALSE, header = TRUE)
cd<-cbind(cd, trait = "Crohn's disease")
cd<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_IBD_CD_GWAS.smr", stringsAsFactors = FALSE, header = TRUE)
Chronotype<-cbind(Chronotype, trait = "Chronotype")
Chronotype<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_Chronotype.smr", stringsAsFactors = FALSE, header = TRUE)
CKD<-cbind(CKD, trait = "Chronic Kidney Disease")
CKD<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_CKD.smr", stringsAsFactors = FALSE, header = TRUE)
cpd<-cbind(cpd, trait = "Cigarettes per day")
cpd<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_CPD_GWAS.smr", stringsAsFactors = FALSE, header = TRUE)
DI<-cbind(DI, trait = "Disposition index")
DI<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_DI.smr", stringsAsFactors = FALSE, header = TRUE)
DIadjBMI<-cbind(DIadjBMI, trait = "Disposition index adjBMI")
DIadjBMI<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_DI_adj_BMI.smr", stringsAsFactors = FALSE, header = TRUE)
eGFRcrea<-cbind(eGFRcrea, trait = "estimated glomerular filtration rate based on serum creatinine")
eGFRcrea<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_eGFRcrea_overall.smr", stringsAsFactors = FALSE, header = TRUE)
eGFRcys<-cbind(eGFRcys, trait = "estimated glomerular filtration rate cystatin C")
eGFRcys<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_eGFRcys_overall.smr", stringsAsFactors = FALSE, header = TRUE)
egg.bl<-cbind(egg.bl, trait = "Birth length")
egg.bl<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_EGG_BL_GWAS.smr", stringsAsFactors = FALSE, header = TRUE)
ever<-cbind(ever, trait = "Ever smoked")
ever<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_EverSmoked_GWAS.smr", stringsAsFactors = FALSE, header = TRUE)
extremebmi<-cbind(extremebmi, trait = "Extreme BMI")
extremebmi<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_EXTREME_BMI_GWAS.smr", stringsAsFactors = FALSE, header = TRUE)
extremeheight<-cbind(extremeheight, trait = "Extreme height")
extremeheight<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_EXTREME_HEIGHT_GWAS.smr", stringsAsFactors = FALSE, header = TRUE)
extremewhr<-cbind(extremewhr, trait = "Extreme waits hip ratio")
extremewhr<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_EXTREME_WHR_GWAS.smr", stringsAsFactors = FALSE, header = TRUE)
FatherAAD<-cbind(FatherAAD, trait = "Father age at death")
FatherAAD<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_father_age_death.smr", stringsAsFactors = FALSE, header = TRUE)
hdl<-cbind(hdl, trait = "HDL")
hdl<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_HDL_GWAS.smr", stringsAsFactors = FALSE, header = TRUE)
height<-cbind(height, trait = "Height")
height<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_Height_GWAS.smr", stringsAsFactors = FALSE, header = TRUE)
hip<-cbind(hip, trait = "Hip circumference")
hip<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_HIP_GWAS.smr", stringsAsFactors = FALSE, header = TRUE)
HR<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_HR.smr", stringsAsFactors = FALSE, header = TRUE)
ibd<-cbind(ibd, trait = "Inflammatory bowel disease")
ibd<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_IBD_GWAS.smr", stringsAsFactors = FALSE, header = TRUE, fill = TRUE)
ISR<-cbind(ISR, trait = "Insulin secretion rate")
ISR<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_ISR.smr", stringsAsFactors = FALSE, header = TRUE)
ISRadjBMI<-cbind(ISRadjBMI, trait = "Insulin secretion rate adjBMI")
ISRadjBMI<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_ISR_adj_BMI.smr", stringsAsFactors = FALSE, header = TRUE)
ldl<-cbind(ldl, trait = "LDL")
ldl<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_LDL_GWAS.smr", stringsAsFactors = FALSE, header = TRUE)
mdd<-cbind(mdd, trait = "Major depressive disorder")
mdd<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_MDD_GWAS.smr", stringsAsFactors = FALSE, header = TRUE)
migraine<-cbind(migraine, trait = "Migraine")
migraine<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_Migraine_GWAS.smr", stringsAsFactors = FALSE, header = TRUE)
MotherAAD<-cbind(MotherAAD, trait = "Mother age at death")
MotherAAD<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_mother_age_death.smr", stringsAsFactors = FALSE, header = TRUE)
obesity1<-cbind(obesity1, trait = "Obesity class 1")
obesity1<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_ObesityClass1_GWAS.smr", stringsAsFactors = FALSE, header = TRUE)
obesity2<-cbind(obesity2, trait = "Obesity class 2")
obesity2<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_ObesityClass2_GWAS.smr", stringsAsFactors = FALSE, header = TRUE)
obesity3<-cbind(obesity3, trait = "Obesity class 3")
obesity3<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_ObesityClass3_GWAS.smr", stringsAsFactors = FALSE, header = TRUE)
overweight<-cbind(overweight, trait = "Overweight")
overweight<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_Overweight_GWAS.smr", stringsAsFactors = FALSE, header = TRUE)
ParentsAAD<-cbind(ParentsAAD, trait = "Parents age at death")
ParentsAAD<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_parents_age_death.smr", stringsAsFactors = FALSE, header = TRUE)
PEAK<-cbind(PEAK, trait = "Peak insulin response")
PEAK<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_PEAK.smr", stringsAsFactors = FALSE, header = TRUE)
PEAKadjSI<-cbind(PEAKadjSI, trait = "Peak insulin response adj SI")
PEAKadjSI<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_PEAK_adj_SI.smr", stringsAsFactors = FALSE, header = TRUE)
PEAKadjSIBMI<-cbind(PEAKadjSIBMI, trait = "Peak insulin response adj SI BMI")
PEAKadjSIBMI<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_PEAK_adj_BMI+SI.smr", stringsAsFactors = FALSE, header = TRUE)
pg10f<-cbind(pg10f, trait = "Pubertal Growth - Height at age 10 (F)")
pg10f<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_Pubertal_growth_10F_GWAS.smr", stringsAsFactors = FALSE, header = TRUE)
pg10f12m<-cbind(pg10f12m, trait = "Pubertal Growth - Height at age 10 (F) and 12 (M)")
pg10f12m<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_Pubertal_growth_10F_12M_combined_GWAS.smr", stringsAsFactors = FALSE, header = TRUE)
pg12m<-cbind(pg12m, trait = "Pubertal Growth - Height at age 12 (M)")
pg12m<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_Pubertal_growth_12M_GWAS.smr", stringsAsFactors = FALSE, header = TRUE)
pgf<-cbind(pgf, trait = "Pubertal Growth - growth across the pubertal growth period (F)")
pgf<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_Pubertal_growth_PGF_GWAS.smr", stringsAsFactors = FALSE, header = TRUE)
pgfpgm<-cbind(pgfpgm, trait = "Pubertal Growth - growth across the pubertal growth period")
pgfpgm<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_Pubertal_growth_PGF_PGM_combined_GWAS.smr", stringsAsFactors = FALSE, header = TRUE)
pgm<-cbind(pgm, trait = "Pubertal Growth - growth across the pubertal growth period (M)")
pgm<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_Pubertal_growth_PGM_GWAS.smr", stringsAsFactors = FALSE, header = TRUE)
ptf<-cbind(ptf, trait = "Pubertal Growth - growth in late adolescence (F)")
ptf<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_Pubertal_growth_PTF_GWAS.smr", stringsAsFactors = FALSE, header = TRUE)
ptfptm<-cbind(ptfptm, trait = "Pubertal Growth - growth in late adolescence")
ptfptm<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_Pubertal_growth_PTF_PTM_combined_GWAS.smr", stringsAsFactors = FALSE, header = TRUE)
ptm<-cbind(ptm, trait = "Pubertal Growth - growth in late adolescence (M)")
ptm<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_Pubertal_growth_PTM_GWAS.smr", stringsAsFactors = FALSE, header = TRUE)
ra<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_RA_GWAS.smr", stringsAsFactors = FALSE, header = TRUE, fill = TRUE)
scz<-cbind(scz, trait = "Schizophrenia")
scz<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_SCZGWAS.smr", stringsAsFactors = FALSE, header = TRUE)
SleepDuration<-cbind(SleepDuration, trait = "Sleep Duration")
SleepDuration<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_SleepDuration.smr", stringsAsFactors = FALSE, header = TRUE)
t2d<-cbind(t2d, trait = "T2 Diabetes")
t2d<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_T2Diabetes_GWAS.smr", stringsAsFactors = FALSE, header = TRUE)
tannerf<-cbind(tannerf, trait = "Tanner stage (F)")
tannerf<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_Tanner_females_GWAS.smr", stringsAsFactors = FALSE, header = TRUE)
tannerfm<-cbind(tannerfm, trait = "Tanner stage")
tannerfm<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_Tanner_male_females_combined_GWAS.smr", stringsAsFactors = FALSE, header = TRUE)
tannerm<-cbind(tannerm, trait = "Tanner stage (M)")
tannerm<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_Tanner_males_GWAS.smr", stringsAsFactors = FALSE, header = TRUE)
tc<-cbind(tc, trait = "Total cholestoral")
tc<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_TC_GWAS.smr", stringsAsFactors = FALSE, header = TRUE)
tg<-cbind(tg, trait = "Triglycerides")
tg<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_TG_GWAS.smr", stringsAsFactors = FALSE, header = TRUE)
UACR<-cbind(UACR, trait = "Urinary albumin-to-creatinine ratio")
UACR<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_UACR_overall.smr", stringsAsFactors = FALSE, header = TRUE)
uc<-cbind(uc, trait = "Ulceritive colitis")
uc<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_IBD_UC_GWAS.smr", stringsAsFactors = FALSE, header = TRUE)
wc<-cbind(wc, trait = "Waist circumference")
wc<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_WC_GWAS.smr", stringsAsFactors = FALSE, header = TRUE)
wcadjbmi<-cbind(wcadjbmi, trait = "Waist circumference adjusted for BMI")
wcadjbmi<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_WCadjBMI_GWAS.smr", stringsAsFactors = FALSE, header = TRUE)
whr<-cbind(whr, trait = "Waist hip ratio")
whr<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_WHR_GWAS.smr", stringsAsFactors = FALSE, header = TRUE)
whradjbmi<-cbind(whradjbmi, trait = "Waist hip ratio adjusted for BMI")
whradjbmi<-read.table("MatrixEQTL/Unrelated/SMR/SMR_BloodmQTL_US_EPIC_WHRadjBMI_GWAS.smr", stringsAsFactors = FALSE, header = TRUE)

blood<-list(AIR,AIRadjBMISI,AIRAdjSII,als,apop_derm,bmi,bpd,cardio2,cd,Chronotype,CKD,cpd,DI,DIadjBMI,eGFRcrea,eGFRcys,egg.bl,ever,extremebmi,extremeheight,extremewhr,FatherAAD,hdl,height,hip,ibd,ISR,ISRadjBMI,ldl,mdd,migraine,MotherAAD,obesity1,obesity2,obesity3,overweight,ParentsAAD,PEAK,PEAKadjSI,PEAKadjSIBMI,pg10f,pg10f12m,pg12m,pgf,pgfpgm,pgm,ptf,ptfptm,ptm,scz,SleepDuration,t2d,tannerf,tannerfm,tannerm,tc,tg,UACR,uc,wc,wcadjbmi,whr,whradjbmi)

blood.signif<-lapply(blood, signif, 3.95e-7)
blood.signif.hetP<-lapply(blood.signif, hetP)

blood.all<-rbindlist(blood)
blood.signif<-rbindlist(blood.signif)
blood.signif.hetP<-rbindlist(blood.signif.hetP)
epicManifest<-read.csv("MethylationEPIC_v-1-0_B2.csv", skip = 7)

blood.signif$Gene<-as.character(epicManifest$UCSC_RefGene_Name)[match(blood.signif$ProbeID, epicManifest$IlmnID)]
blood.signif.hetP$Gene<-as.character(epicManifest$UCSC_RefGene_Name)[match(blood.signif.hetP$ProbeID, epicManifest$IlmnID)]

blood.signif.hetP<-cbind(blood.signif.hetP, epicManifest[match(blood.signif.hetP$ProbeID, epicManifest$IlmnID),c("CHR", "MAPINFO", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_UCSC_CpG_Island", "Regulatory_Feature_Name", "Regulatory_Feature_Group", "DNase_Hypersensitivity_NAME", "OpenChromatin_NAME", "TFBS_Evidence_Count", "Methyl450_Loci", "SNP_ID", "SNP_DISTANCE", "SNP_MinorAlleleFrequency")])

blood.all<-cbind(blood.all, epicManifest[match(blood.all$ProbeID, epicManifest$IlmnID),c("CHR", "MAPINFO", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_UCSC_CpG_Island", "Regulatory_Feature_Name", "Regulatory_Feature_Group", "DNase_Hypersensitivity_NAME", "OpenChromatin_NAME", "TFBS_Evidence_Count", "Methyl450_Loci", "SNP_ID", "SNP_DISTANCE", "SNP_MinorAlleleFrequency")])

blood.probes<-cbind(blood.probes, epicManifest[match(blood.probes$ProbeID, epicManifest$IlmnID),c("CHR", "MAPINFO", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_UCSC_CpG_Island", "Regulatory_Feature_Name", "Regulatory_Feature_Group", "DNase_Hypersensitivity_NAME", "OpenChromatin_NAME", "TFBS_Evidence_Count", "Methyl450_Loci", "SNP_ID", "SNP_DISTANCE", "SNP_MinorAlleleFrequency")])
	
blood.signif.hetP$Gene<-unlist(lapply(blood.signif.hetP$Gene, uniqueAnno))
write.csv(blood.signif.hetP, "SupplementaryTable5.csv")

length(which(blood.signif.hetP$Regulatory_Feature_Name != "" | blood.signif.hetP$DNase_Hypersensitivity_NAME != "" | blood.signif.hetP$OpenChromatin_NAME != ""))

length(which(blood.probes$Regulatory_Feature_Name != "" | blood.probes$DNase_Hypersensitivity_NAME != "" | blood.probes$OpenChromatin_NAME != ""))

length(which(blood.signif.hetP$Regulatory_Feature_Name == "" & blood.signif.hetP$DNase_Hypersensitivity_NAME == "" & blood.signif.hetP$OpenChromatin_NAME == "" & blood.signif.hetP$UCSC_RefGene_Name == ""))

## enriched for overlap with regulatory features

tab<-matrix(data = NA, nrow = 2, ncol = 2)

blood.probes.hetP<-blood.probes[blood.probes$ProbeID %in% blood.signif.hetP$ProbeID,]
tab[1,]<-table(blood.probes.hetP$Regulatory_Feature_Name != "" | blood.probes.hetP$DNase_Hypersensitivity_NAME != "" | blood.probes.hetP$OpenChromatin_NAME != "")
tab[2,]<-table(blood.probes$Regulatory_Feature_Name != "" | blood.probes$DNase_Hypersensitivity_NAME != "" | blood.probes$OpenChromatin_NAME != "") - tab[1,]
tab<-tab[,c(2,1)]
fisher.test(tab)

tab<-matrix(data = NA, nrow = 2, ncol = 2)
blood.probes.hetP<-blood.probes[blood.probes$ProbeID %in% blood.signif.hetP$ProbeID,]
tab[1,]<-table(blood.probes.hetP$Regulatory_Feature_Name != "")
tab[2,]<-table(blood.probes$Regulatory_Feature_Name != "") - tab[1,]
tab<-tab[,c(2,1)]
fisher.test(tab)

tab<-matrix(data = NA, nrow = 2, ncol = 2)
blood.probes.hetP<-blood.probes[blood.probes$ProbeID %in% blood.signif.hetP$ProbeID,]
tab[1,]<-table(blood.probes.hetP$DNase_Hypersensitivity_NAME != "")
tab[2,]<-table(blood.probes$DNase_Hypersensitivity_NAME != "") - tab[1,]
tab<-tab[,c(2,1)]
fisher.test(tab)

tab<-matrix(data = NA, nrow = 2, ncol = 2)
blood.probes.hetP<-blood.probes[blood.probes$ProbeID %in% blood.signif.hetP$ProbeID,]
tab[1,]<-table(blood.probes.hetP$OpenChromatin_NAME != "")
tab[2,]<-table(blood.probes$OpenChromatin_NAME != "") - tab[1,]
tab<-tab[,c(2,1)]
fisher.test(tab)


## enriched for overlap with genic features

tab<-matrix(data = NA, nrow = 2, ncol = 2)

tab[1,]<-table(blood.probes.hetP$UCSC_RefGene_Name != "")
tab[2,]<-table(blood.probes$UCSC_RefGene_Name != "") - tab[1,]
tab<-tab[,c(2,1)]
fisher.test(tab)


mQTL.genetraitsALL<-NULL
for(each in unique(blood.signif.hetP$trait)){
	sub<-blood.signif.hetP[which(blood.signif.hetP$trait == each),]
	mQTL.genetraitsALL<-c(mQTL.genetraitsALL, paste(unlist(strsplit(sub$Gene, ";")), each))
}

blood.signif.hetP.450k<-blood.signif.hetP[which(blood.signif.hetP$Methyl450_Loci == TRUE),]
mQTL.genetraitsALL.450k<-NULL
for(each in unique(blood.signif.hetP.450k$trait)){
	sub<-blood.signif.hetP[which(blood.signif.hetP.450k$trait == each),]
	mQTL.genetraitsALL.450k<-c(mQTL.genetraitsALL.450k, paste(unlist(strsplit(sub$Gene, ";")), each))
}

blood.signif.hetP.epic<-blood.signif.hetP[which(is.na(blood.signif.hetP$Methyl450_Loci) == TRUE),]
mQTL.genetraitsALL.epic<-NULL
for(each in unique(blood.signif.hetP.epic$trait)){
	sub<-blood.signif.hetP[which(blood.signif.hetP.epic$trait == each),]
	mQTL.genetraitsALL.epic<-c(mQTL.genetraitsALL.epic, paste(unlist(strsplit(sub$Gene, ";")), each))
}

##calc number of unique genes
length(unique(unlist(lapply(strsplit(unique(mQTL.genetraitsALL), " "), head, n = 1))))


##### collate eQTL

setwd("")

als<-read.table("SMR_westra_eqtl_ALS_lmm_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
scz<-read.table("SMR_westra_eqtl_SCZGWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
ra<-read.table("SMR_westra_eqtl_RA_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE, fill = TRUE)
bmi<-read.table("SMR_westra_eqtl_BMI_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
height<-read.table("SMR_westra_eqtl_Height_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
whradjbmi<-read.table("SMR_westra_eqtl_WHRadjBMI_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
cardio2<-read.table("SMR_westra_eqtl_Cardiogram_2_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
cardio<-read.table("SMR_westra_eqtl_Cardiogram_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
t2d<-read.table("SMR_westra_eqtl_T2Diabetes_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
t2d.2<-read.table("SMR_westra_eqtl_T2Diabetes_Metabochip_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
t2d.3<-read.table("SMR_westra_eqtl_T2Diabetes__2012DEC17_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
egg.bl<-read.table("SMR_westra_eqtl_EGG_BL_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
pg10f<-read.table("SMR_westra_eqtl_Pubertal_growth_10F_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
pg12m<-read.table("SMR_westra_eqtl_Pubertal_growth_12M_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
pg10f12m<-read.table("SMR_westra_eqtl_Pubertal_growth_10F_12M_combined_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
pgf<-read.table("SMR_westra_eqtl_Pubertal_growth_PGF_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
pgm<-read.table("SMR_westra_eqtl_Pubertal_growth_PGM_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
pgfpgm<-read.table("SMR_westra_eqtl_Pubertal_growth_PGF_PGM_combined_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
ptf<-read.table("SMR_westra_eqtl_Pubertal_growth_PTF_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
ptm<-read.table("SMR_westra_eqtl_Pubertal_growth_PTM_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
ptfptm<-read.table("SMR_westra_eqtl_Pubertal_growth_PTF_PTM_combined_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
tannerf<-read.table("SMR_westra_eqtl_Tanner_females_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
tannerm<-read.table("SMR_westra_eqtl_Tanner_males_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
tannerfm<-read.table("SMR_westra_eqtl_Tanner_male_females_combined_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
overweight<-read.table("SMR_westra_eqtl_Overweight_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
obesity1<-read.table("SMR_westra_eqtl_ObesityClass1_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
obesity2<-read.table("SMR_westra_eqtl_ObesityClass2_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
obesity3<-read.table("SMR_westra_eqtl_ObesityClass3_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
extremebmi<-read.table("SMR_westra_eqtl_EXTREME_BMI_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
extremeheight<-read.table("SMR_westra_eqtl_EXTREME_HEIGHT_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
extremewhr<-read.table("SMR_westra_eqtl_EXTREME_WHR_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
hip<-read.table("SMR_westra_eqtl_HIP_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
wc<-read.table("SMR_westra_eqtl_WC_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
wcadjbmi<-read.table("SMR_westra_eqtl_WCadjBMI_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
whr<-read.table("SMR_westra_eqtl_WHR_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
hdl<-read.table("SMR_westra_eqtl_HDL_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
ldl<-read.table("SMR_westra_eqtl_LDL_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
tc<-read.table("SMR_westra_eqtl_TC_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
tg<-read.table("SMR_westra_eqtl_TG_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
cd<-read.table("SMR_westra_eqtl_IBD_CD_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
ibd<-read.table("SMR_westra_eqtl_IBD_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE, fill = TRUE)
uc<-read.table("SMR_westra_eqtl_IBD_UC_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
bpd<-read.table("SMR_westra_eqtl_BPD_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
mdd<-read.table("SMR_westra_eqtl_MDD_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
cpd<-read.table("SMR_westra_eqtl_CPD_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
ever<-read.table("SMR_westra_eqtl_EverSmoked_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
migraine<-read.table("SMR_westra_eqtl_Migraine_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
ad12<-read.table("SMR_westra_eqtl_AD_Stage1_2Combined_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
ad<-read.table("SMR_westra_eqtl_AD_Stage1_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
## update
ISR<-read.table("SMR_westra_eqtl_ISR.smr", stringsAsFactors = FALSE, header = TRUE)
ISRadjBMI<-read.table("SMR_westra_eqtl_ISR_adj_BMI.smr", stringsAsFactors = FALSE, header = TRUE)
DI<-read.table("SMR_westra_eqtl_DI.smr", stringsAsFactors = FALSE, header = TRUE)
DIadjBMI<-read.table("SMR_westra_eqtl_DI_adj_BMI.smr", stringsAsFactors = FALSE, header = TRUE)
AIR<-read.table("SMR_westra_eqtl_AIR.smr", stringsAsFactors = FALSE, header = TRUE)
AIRAdjSII<-read.table("SMR_westra_eqtl_AIR_adj_SII.smr", stringsAsFactors = FALSE, header = TRUE)
AIRadjBMISI<-read.table("SMR_westra_eqtl_AIR_adj_BMI+SI.smr", stringsAsFactors = FALSE, header = TRUE)
apop_derm<-read.table("SMR_westra_eqtl_EAGLE_AD.smr", stringsAsFactors = FALSE, header = TRUE)

PEAK<-read.table("SMR_westra_eqtl_PEAK.smr", stringsAsFactors = FALSE, header = TRUE)
PEAKadjSI<-read.table("SMR_westra_eqtl_PEAK_adj_SI.smr", stringsAsFactors = FALSE, header = TRUE)
PEAKadjSIBMI<-read.table("SMR_westra_eqtl_PEAK_adj_BMI+SI.smr", stringsAsFactors = FALSE, header = TRUE)
HR<-read.table("SMR_westra_eqtl_HR.smr", stringsAsFactors = FALSE, header = TRUE)
Leptin<-read.table("SMR_westra_eqtl_LeptinNotAdjBMI.smr", stringsAsFactors = FALSE, header = TRUE)
LeptinadjBMI<-read.table("SMR_westra_eqtl_LeptinAdjBMI.smr", stringsAsFactors = FALSE, header = TRUE)
MotherAAD<-read.table("SMR_westra_eqtl_mother_age_death.smr", stringsAsFactors = FALSE, header = TRUE)
FatherAAD<-read.table("SMR_westra_eqtl_father_age_death.smr", stringsAsFactors = FALSE, header = TRUE)
ParentsAAD<-read.table("SMR_westra_eqtl_parents_age_death.smr", stringsAsFactors = FALSE, header = TRUE)
SleepDuration<-read.table("SMR_westra_eqtl_SleepDuration.smr", stringsAsFactors = FALSE, header = TRUE)
Chronotype<-read.table("SMR_westra_eqtl_Chronotype.smr", stringsAsFactors = FALSE, header = TRUE)
UACR<-read.table("SMR_westra_eqtl_UACR_overall.smr", stringsAsFactors = FALSE, header = TRUE)
eGFRcys<-read.table("SMR_westra_eqtl_eGFRcys_overall.smr", stringsAsFactors = FALSE, header = TRUE)
eGFRcrea<-read.table("SMR_westra_eqtl_eGFRcrea_overall.smr", stringsAsFactors = FALSE, header = TRUE)
CKD<-read.table("SMR_westra_eqtl_CKD.smr", stringsAsFactors = FALSE, header = TRUE)




als<-cbind(als, trait = "Amyloid lateral sclerosis")
scz<-cbind(scz, trait = "Schizophrenia")
bmi<-cbind(bmi, trait = "BMI")
height<-cbind(height, trait = "Height")
whradjbmi<-cbind(whradjbmi, trait = "Waist hip ratio adjusted for BMI")
cardio2<-cbind(cardio, trait = "Coronary artery disease (2015)")
cardio<-cbind(cardio, trait = "Coronary artery disease (2011)")
t2d<-cbind(t2d, trait = "T2 Diabetes (2015)")
t2d.2<-cbind(t2d.2, trait = "T2 Diabetes (Metabochip)")
t2d.3<-cbind(t2d.3, trait = "T2 Diabetes (2012)")
egg.bl<-cbind(egg.bl, trait = "Birth length")
pg10f<-cbind(pg10f, trait = "Pubertal Growth - Height at age 10 (F)")
pg12m<-cbind(pg12m, trait = "Pubertal Growth - Height at age 12 (M)")
pg10f12m<-cbind(pg10f12m, trait = "Pubertal Growth - Height at age 10 (F) and 12 (M)")
pgf<-cbind(pgf, trait = "Pubertal Growth - growth across the pubertal growth period (F)")
pgm<-cbind(pgm, trait = "Pubertal Growth - growth across the pubertal growth period (M)")
pgfpgm<-cbind(pgfpgm, trait = "Pubertal Growth - growth across the pubertal growth period")
ptf<-cbind(ptf, trait = "Pubertal Growth - growth in late adolescence (F)")
ptm<-cbind(ptm, trait = "Pubertal Growth - growth in late adolescence (M)")
ptfptm<-cbind(ptfptm, trait = "Pubertal Growth - growth in late adolescence")
tannerf<-cbind(tannerf, trait = "Tanner stage (F)")
tannerm<-cbind(tannerm, trait = "Tanner stage (M)")
tannerfm<-cbind(tannerfm, trait = "Tanner stage")
overweight<-cbind(overweight, trait = "Overweight")
obesity1<-cbind(obesity1, trait = "Obesity class 1")
obesity2<-cbind(obesity2, trait = "Obesity class 2")
obesity3<-cbind(obesity3, trait = "Obesity class 3")
extremebmi<-cbind(extremebmi, trait = "Extreme BMI")
extremeheight<-cbind(extremeheight, trait = "Extreme height")
extremewhr<-cbind(extremewhr, trait = "Extreme waits hip ratio")
hip<-cbind(hip, trait = "Hip circumference")
wc<-cbind(wc, trait = "Waist circumference")
wcadjbmi<-cbind(wcadjbmi, trait = "Waist circumference adjusted for BMI")
whr<-cbind(whr, trait = "Waist hip ratio")
hdl<-cbind(hdl, trait = "HDL")
ldl<-cbind(ldl, trait = "LDL")
tc<-cbind(tc, trait = "Total cholestoral")
tg<-cbind(tg, trait = "Triglycerides")
cd<-cbind(cd, trait = "Crohn's disease")
ibd<-cbind(ibd, trait = "Inflammatory bowel disease")
uc<-cbind(uc, trait = "Ulceritive colitis")
bpd<-cbind(bpd, trait = "Bipolar disorder")
mdd<-cbind(mdd, trait = "Major depressive disorder")
cpd<-cbind(cpd, trait = "Cigarettes per day")
ever<-cbind(ever, trait = "Ever smoked")
migraine<-cbind(migraine, trait = "Migraine")
ISR<-cbind(ISR, trait = "Insulin secretion rate")
ISRadjBMI<-cbind(ISRadjBMI, trait = "Insulin secretion rate adjBMI")
DI<-cbind(DI, trait = "Disposition index")
DIadjBMI<-cbind(DIadjBMI, trait = "Disposition index adjBMI")
AIR<-cbind(AIR, trait = "Acute insulin response")
AIRAdjSII<-cbind(AIRAdjSII, trait = "Acute insulin response adj SI")
AIRadjBMISI<-cbind(AIRadjBMISI, trait = "Acute insulin response adj SI, BMI")
apop_derm<-cbind(apop_derm, trait = "Atopic dermatitis (eczema)")
PEAK<-cbind(PEAK, trait = "Peak insulin response")
PEAKadjSI<-cbind(PEAKadjSI, trait = "Peak insulin response adj SI")
PEAKadjSIBMI<-cbind(PEAKadjSIBMI, trait = "Peak insulin response adj SI BMI")
MotherAAD<-cbind(MotherAAD, trait = "Mother age at death")
FatherAAD<-cbind(FatherAAD, trait = "Father age at death")
ParentsAAD<-cbind(ParentsAAD, trait = "Parents age at death")
SleepDuration<-cbind(SleepDuration, trait = "Sleep Duration")
Chronotype<-cbind(Chronotype, trait = "Chronotype")AIR<-cbind(AIR, trait = "Acute insulin response")
AIR<-read.table("SMR_westra_eqtl_AIR.smr", stringsAsFactors = FALSE, header = TRUE)
AIRadjBMISI<-cbind(AIRadjBMISI, trait = "Acute insulin response adj SI, BMI")
AIRadjBMISI<-read.table("SMR_westra_eqtl_AIR_adj_BMI+SI.smr", stringsAsFactors = FALSE, header = TRUE)
AIRAdjSII<-cbind(AIRAdjSII, trait = "Acute insulin response adj SI")
AIRAdjSII<-read.table("SMR_westra_eqtl_AIR_adj_SII.smr", stringsAsFactors = FALSE, header = TRUE)
als<-cbind(als, trait = "Amyloid lateral sclerosis")
als<-read.table("SMR_westra_eqtl_ALS_lmm_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
apop_derm<-cbind(apop_derm, trait = "Atopic dermatitis (eczema)")
apop_derm<-read.table("SMR_westra_eqtl_EAGLE_AD.smr", stringsAsFactors = FALSE, header = TRUE)
bmi<-cbind(bmi, trait = "BMI")
bmi<-read.table("SMR_westra_eqtl_BMI_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
bpd<-cbind(bpd, trait = "Bipolar disorder")
bpd<-read.table("SMR_westra_eqtl_BPD_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
cardio2<-cbind(cardio, trait = "Coronary artery disease (2015)")
cardio2<-read.table("SMR_westra_eqtl_Cardiogram_2_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
cd<-cbind(cd, trait = "Crohn's disease")
cd<-read.table("SMR_westra_eqtl_IBD_CD_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
Chronotype<-cbind(Chronotype, trait = "Chronotype")
Chronotype<-read.table("SMR_westra_eqtl_Chronotype.smr", stringsAsFactors = FALSE, header = TRUE)
CKD<-cbind(CKD, trait = "Chronic Kidney Disease")
CKD<-read.table("SMR_westra_eqtl_CKD.smr", stringsAsFactors = FALSE, header = TRUE)
cpd<-cbind(cpd, trait = "Cigarettes per day")
cpd<-read.table("SMR_westra_eqtl_CPD_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
DI<-cbind(DI, trait = "Disposition index")
DI<-read.table("SMR_westra_eqtl_DI.smr", stringsAsFactors = FALSE, header = TRUE)
DIadjBMI<-cbind(DIadjBMI, trait = "Disposition index adjBMI")
DIadjBMI<-read.table("SMR_westra_eqtl_DI_adj_BMI.smr", stringsAsFactors = FALSE, header = TRUE)
eGFRcrea<-cbind(eGFRcrea, trait = "estimated glomerular filtration rate based on serum creatinine")
eGFRcrea<-read.table("SMR_westra_eqtl_eGFRcrea_overall.smr", stringsAsFactors = FALSE, header = TRUE)
eGFRcys<-cbind(eGFRcys, trait = "estimated glomerular filtration rate cystatin C")
eGFRcys<-read.table("SMR_westra_eqtl_eGFRcys_overall.smr", stringsAsFactors = FALSE, header = TRUE)
egg.bl<-cbind(egg.bl, trait = "Birth length")
egg.bl<-read.table("SMR_westra_eqtl_EGG_BL_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
ever<-cbind(ever, trait = "Ever smoked")
ever<-read.table("SMR_westra_eqtl_EverSmoked_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
extremebmi<-cbind(extremebmi, trait = "Extreme BMI")
extremebmi<-read.table("SMR_westra_eqtl_EXTREME_BMI_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
extremeheight<-cbind(extremeheight, trait = "Extreme height")
extremeheight<-read.table("SMR_westra_eqtl_EXTREME_HEIGHT_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
extremewhr<-cbind(extremewhr, trait = "Extreme waits hip ratio")
extremewhr<-read.table("SMR_westra_eqtl_EXTREME_WHR_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
FatherAAD<-cbind(FatherAAD, trait = "Father age at death")
FatherAAD<-read.table("SMR_westra_eqtl_father_age_death.smr", stringsAsFactors = FALSE, header = TRUE)
hdl<-cbind(hdl, trait = "HDL")
hdl<-read.table("SMR_westra_eqtl_HDL_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
height<-cbind(height, trait = "Height")
height<-read.table("SMR_westra_eqtl_Height_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
hip<-cbind(hip, trait = "Hip circumference")
hip<-read.table("SMR_westra_eqtl_HIP_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
HR<-read.table("SMR_westra_eqtl_HR.smr", stringsAsFactors = FALSE, header = TRUE)
ibd<-cbind(ibd, trait = "Inflammatory bowel disease")
ibd<-read.table("SMR_westra_eqtl_IBD_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE, fill = TRUE)
ISR<-cbind(ISR, trait = "Insulin secretion rate")
ISR<-read.table("SMR_westra_eqtl_ISR.smr", stringsAsFactors = FALSE, header = TRUE)
ISRadjBMI<-cbind(ISRadjBMI, trait = "Insulin secretion rate adjBMI")
ISRadjBMI<-read.table("SMR_westra_eqtl_ISR_adj_BMI.smr", stringsAsFactors = FALSE, header = TRUE)
ldl<-cbind(ldl, trait = "LDL")
ldl<-read.table("SMR_westra_eqtl_LDL_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
mdd<-cbind(mdd, trait = "Major depressive disorder")
mdd<-read.table("SMR_westra_eqtl_MDD_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
migraine<-cbind(migraine, trait = "Migraine")
migraine<-read.table("SMR_westra_eqtl_Migraine_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
MotherAAD<-cbind(MotherAAD, trait = "Mother age at death")
MotherAAD<-read.table("SMR_westra_eqtl_mother_age_death.smr", stringsAsFactors = FALSE, header = TRUE)
obesity1<-cbind(obesity1, trait = "Obesity class 1")
obesity1<-read.table("SMR_westra_eqtl_ObesityClass1_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
obesity2<-cbind(obesity2, trait = "Obesity class 2")
obesity2<-read.table("SMR_westra_eqtl_ObesityClass2_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
obesity3<-cbind(obesity3, trait = "Obesity class 3")
obesity3<-read.table("SMR_westra_eqtl_ObesityClass3_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
overweight<-cbind(overweight, trait = "Overweight")
overweight<-read.table("SMR_westra_eqtl_Overweight_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
ParentsAAD<-cbind(ParentsAAD, trait = "Parents age at death")
ParentsAAD<-read.table("SMR_westra_eqtl_parents_age_death.smr", stringsAsFactors = FALSE, header = TRUE)
PEAK<-cbind(PEAK, trait = "Peak insulin response")
PEAK<-read.table("SMR_westra_eqtl_PEAK.smr", stringsAsFactors = FALSE, header = TRUE)
PEAKadjSI<-cbind(PEAKadjSI, trait = "Peak insulin response adj SI")
PEAKadjSI<-read.table("SMR_westra_eqtl_PEAK_adj_SI.smr", stringsAsFactors = FALSE, header = TRUE)
PEAKadjSIBMI<-cbind(PEAKadjSIBMI, trait = "Peak insulin response adj SI BMI")
PEAKadjSIBMI<-read.table("SMR_westra_eqtl_PEAK_adj_BMI+SI.smr", stringsAsFactors = FALSE, header = TRUE)
pg10f<-cbind(pg10f, trait = "Pubertal Growth - Height at age 10 (F)")
pg10f<-read.table("SMR_westra_eqtl_Pubertal_growth_10F_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
pg10f12m<-cbind(pg10f12m, trait = "Pubertal Growth - Height at age 10 (F) and 12 (M)")
pg10f12m<-read.table("SMR_westra_eqtl_Pubertal_growth_10F_12M_combined_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
pg12m<-cbind(pg12m, trait = "Pubertal Growth - Height at age 12 (M)")
pg12m<-read.table("SMR_westra_eqtl_Pubertal_growth_12M_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
pgf<-cbind(pgf, trait = "Pubertal Growth - growth across the pubertal growth period (F)")
pgf<-read.table("SMR_westra_eqtl_Pubertal_growth_PGF_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
pgfpgm<-cbind(pgfpgm, trait = "Pubertal Growth - growth across the pubertal growth period")
pgfpgm<-read.table("SMR_westra_eqtl_Pubertal_growth_PGF_PGM_combined_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
pgm<-cbind(pgm, trait = "Pubertal Growth - growth across the pubertal growth period (M)")
pgm<-read.table("SMR_westra_eqtl_Pubertal_growth_PGM_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
ptf<-cbind(ptf, trait = "Pubertal Growth - growth in late adolescence (F)")
ptf<-read.table("SMR_westra_eqtl_Pubertal_growth_PTF_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
ptfptm<-cbind(ptfptm, trait = "Pubertal Growth - growth in late adolescence")
ptfptm<-read.table("SMR_westra_eqtl_Pubertal_growth_PTF_PTM_combined_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
ptm<-cbind(ptm, trait = "Pubertal Growth - growth in late adolescence (M)")
ptm<-read.table("SMR_westra_eqtl_Pubertal_growth_PTM_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
ra<-read.table("SMR_westra_eqtl_RA_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE, fill = TRUE)
scz<-cbind(scz, trait = "Schizophrenia")
scz<-read.table("SMR_westra_eqtl_SCZGWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
SleepDuration<-cbind(SleepDuration, trait = "Sleep Duration")
SleepDuration<-read.table("SMR_westra_eqtl_SleepDuration.smr", stringsAsFactors = FALSE, header = TRUE)
t2d<-cbind(t2d, trait = "T2 Diabetes (2015)")
t2d<-read.table("SMR_westra_eqtl_T2Diabetes_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
tannerf<-cbind(tannerf, trait = "Tanner stage (F)")
tannerf<-read.table("SMR_westra_eqtl_Tanner_females_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
tannerfm<-cbind(tannerfm, trait = "Tanner stage")
tannerfm<-read.table("SMR_westra_eqtl_Tanner_male_females_combined_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
tannerm<-cbind(tannerm, trait = "Tanner stage (M)")
tannerm<-read.table("SMR_westra_eqtl_Tanner_males_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
tc<-cbind(tc, trait = "Total cholestoral")
tc<-read.table("SMR_westra_eqtl_TC_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
tg<-cbind(tg, trait = "Triglycerides")
tg<-read.table("SMR_westra_eqtl_TG_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
UACR<-cbind(UACR, trait = "Urinary albumin-to-creatinine ratio")
UACR<-read.table("SMR_westra_eqtl_UACR_overall.smr", stringsAsFactors = FALSE, header = TRUE)
uc<-cbind(uc, trait = "Ulceritive colitis")
uc<-read.table("SMR_westra_eqtl_IBD_UC_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
wc<-cbind(wc, trait = "Waist circumference")
wc<-read.table("SMR_westra_eqtl_WC_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
wcadjbmi<-cbind(wcadjbmi, trait = "Waist circumference adjusted for BMI")
wcadjbmi<-read.table("SMR_westra_eqtl_WCadjBMI_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
whr<-cbind(whr, trait = "Waist hip ratio")
whr<-read.table("SMR_westra_eqtl_WHR_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)
whradjbmi<-cbind(whradjbmi, trait = "Waist hip ratio adjusted for BMI")
whradjbmi<-read.table("SMR_westra_eqtl_WHRadjBMI_GWAS.chr.smr", stringsAsFactors = FALSE, header = TRUE)

UACR<-cbind(UACR, trait = "Urinary albumin-to-creatinine ratio")
eGFRcys<-cbind(eGFRcys, trait = "estimated glomerular filtration rate cystatin C")
eGFRcrea<-cbind(eGFRcrea, trait = "estimated glomerular filtration rate based on serum creatinine")
CKD<-cbind(CKD, trait = "Chronic Kidney Disease")

eqtl<-list(AIR,AIRadjBMISI,AIRAdjSII,als,apop_derm,bmi,bpd,cardio2,cd,Chronotype,CKD,cpd,DI,DIadjBMI,eGFRcrea,eGFRcys,egg.bl,ever,extremebmi,extremeheight,extremewhr,FatherAAD,hdl,height,hip,ibd,ISR,ISRadjBMI,ldl,mdd,migraine,MotherAAD,obesity1,obesity2,obesity3,overweight,ParentsAAD,PEAK,PEAKadjSI,PEAKadjSIBMI,pg10f,pg10f12m,pg12m,pgf,pgfpgm,pgm,ptf,ptfptm,ptm,scz,SleepDuration,t2d,tannerf,tannerfm,tannerm,tc,tg,UACR,uc,wc,wcadjbmi,whr,whradjbmi)


eqtl.signif<-lapply(eqtl, signif, 8.38e-6)
eqtl.signif.hetP<-lapply(eqtl.signif, hetP)

eqtl.signif<-rbindlist(eqtl.signif)
eqtl.signif.hetP<-rbindlist(eqtl.signif.hetP)

eqtl.probes<-read.table("Westera_1000GFormat.epi")

blood.genes<-unique(unlist(strsplit(as.character(blood.probes$Gene), ";")))
eqtl.genes<-unique(eqtl.probes$V5)

length(intersect(blood.genes, eqtl.genes))
overlapGenes<-intersect(blood.genes, eqtl.genes) ## overlap of tested genes
table(eqtl.probes$V5  %in% overlapGenes) 

count<-0
for(each in overlapGenes){

	count<-count+length(grep(paste(";", each, ";", sep = ""), paste(";", as.character(blood.probes$Gene), ";", sep = "")))
	}
	
mQTL.genetraits<-NULL
for(each in unique(blood.signif.hetP$trait)){
	sub<-blood.signif.hetP[which(blood.signif.hetP$trait == each),]
	mQTL.genetraits<-c(mQTL.genetraits, paste(unlist(unique(strsplit(sub$Gene, ";"))), each))
}
mQTL.genetraitsALL<-NULL
for(each in unique(blood.signif.hetP$trait)){
	sub<-blood.signif.hetP[which(blood.signif.hetP$trait == each),]
	mQTL.genetraitsALL<-c(mQTL.genetraitsALL, paste(unlist(strsplit(sub$Gene, ";")), each))
}

table(mQTL.genetraitsALL)[which(table(mQTL.genetraitsALL) > 2)]->multipleProbes

keep<-NULL
countSame<-NULL
countDiff<-NULL
for(each in names(multipleProbes)){
	gene<-unlist(strsplit(each, " "))[1]
	trait<-paste(unlist(strsplit(each, " "))[-1])
	keep<-rbind(keep, as.data.frame(blood.signif.hetP[intersect(which(blood.signif.hetP$trait == trait),grep(gene, blood.signif.hetP$Gene)),]))
	if(length(table(sign(blood.signif.hetP$b_SMR[intersect(which(blood.signif.hetP$trait == trait),grep(gene, blood.signif.hetP$Gene))]))) == 1){
		countSame<-c(countSame, each)
	} else {
		countDiff<-c(countDiff, each)
		
	}

}

keep<-cbind(keep, epicManifest[match(keep$ProbeID, epicManifest$IlmnID), c("UCSC_RefGene_Name", "UCSC_RefGene_Group")])

setwd("")

table(unlist(lapply(strsplit(unique(mQTL.genetraits), " "), head, n=1)) %in% eqtl.genes)

eQTL.genetraits<-paste(eqtl.signif.hetP$Gene, eqtl.signif.hetP$trait)

intersect(eQTL.genetraits, mQTL.genetraits)

length(intersect(unlist(unique(strsplit(blood.signif.hetP$Gene, ";"))), eqtl.genes))

table(unlist(strsplit(blood.signif.hetP$Gene, ";")) %in% intersect(unlist(unique(strsplit(blood.signif.hetP$Gene, ";"))), eqtl.genes))

testwithEQTL<-unlist(lapply(lapply(strsplit(blood.signif.hetP$Gene, ";"), intersect, eqtl.genes), length))

lookUp<-intersect(eQTL.genetraits, mQTL.genetraits)

eqtl.sub<-eqtl.signif.hetP[match(lookUp, paste(eqtl.signif.hetP$Gene, eqtl.signif.hetP$trait)),]

unique(unlist(lapply(strsplit(unique(mQTL.genetraits), " "), head , n = 1))[unlist(lapply(strsplit(unique(mQTL.genetraits), " "), head , n = 1)) %in% eqtl.genes])->eqtl.genestotest

count<-0
for(each in eqtl.genestotest){

	count<-count+length(grep(paste(";", each, ";", sep = ""), paste(";", as.character(epicManifest$UCSC_RefGene_Name)[match(unique(blood.signif.hetP$ProbeID), epicManifest$IlmnID)], ";", sep = "")))
	}


	
### for mQTL might be multiple probes
blood.signif.hetP$trait<-as.character(blood.signif.hetP$trait)
eqtl.sub$trait<-as.character(eqtl.sub$trait)

lookUp.table<-NULL

for(i in 1:nrow(eqtl.sub)){
	lookUp.table<-rbind(lookUp.table,cbind(blood.signif.hetP[intersect(which(blood.signif.hetP$trait == eqtl.sub$trait[i]), grep(eqtl.sub$Gene[i], blood.signif.hetP$Gene)),], eqtl.sub[i,]))

}

write.csv(lookUp.table, "SupplementaryTable6.csv")
