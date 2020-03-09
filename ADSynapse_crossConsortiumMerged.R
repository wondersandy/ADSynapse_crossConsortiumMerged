library(BiocManager)
install()
library(limma)
library(edgeR)
library(Glimma)
#library(AnnotationDbi)
options(java.parameters = "-Xmx8000m")
library(XLConnect)
library(RColorBrewer)
library(gplots)

## Import DGE summary stat (combined MSBB, ROSMAP, MAYO) 
dge.allMerged <- read.table("~/NKI/Database/differentialExpressionSummary_syn14237651_RNA.tsv", sep = "\t", 
                            header = T, stringsAsFactors = F)

dim(dge.allMerged) # 1330018      18
head(dge.allMerged)

table(dge.allMerged$Model)
#          APOE4       Diagnosis   Diagnosis.AOD   Diagnosis.Sex SourceDiagnosis 
#         345006          329422          130586          276756          248248 

table(dge.allMerged$Tissue)
#   CBE  DLPFC     FP    IFG    PHG    STG    TCX 
#221169 233760 163480 163480 163480 163480 221169

table(dge.allMerged$Comparison)
#                       1-0                       2-0                       2-1                       4-1                       4-2                AD-CONTROL 
#                    115002                    115002                    130586                     15584                     15584                    509618 
# AD-CONTROL.IN.FEMALE-MALE                  AD-OTHER               AD-PATH_AGE                       AOD                       CDR               FEMALE-MALE 
#                     15584                     99418                     34026                     15584                     65392                     31168 
#             OTHER-CONTROL          PATH_AGE-CONTROL               PSP-CONTROL 
#                     99418                     34026                     34026

table(dge.allMerged$Sex)
#     ALL  FEMALE    MALE 
# 1100014  115002  115002

table(dge.allMerged$Study)
#   MAYO   MSSM ROSMAP 
# 442338 653920 233760 

## To see all of above with one code
lapply(dge.allMerged[,c(1,2,3,17,18)], table)
# $Model

# APOE4       Diagnosis   Diagnosis.AOD   Diagnosis.Sex SourceDiagnosis 
# 345006          329422          130586          276756          248248 

# $Tissue

# CBE  DLPFC     FP    IFG    PHG    STG    TCX 
# 221169 233760 163480 163480 163480 163480 221169 

# $Comparison

# 1-0                       2-0                       2-1                       4-1                       4-2                AD-CONTROL 
# 115002                    115002                    130586                     15584                     15584                    509618 
# AD-CONTROL.IN.FEMALE-MALE                  AD-OTHER               AD-PATH_AGE                       AOD                       CDR               FEMALE-MALE 
# 15584                     99418                     34026                     15584                     65392                     31168 
# OTHER-CONTROL          PATH_AGE-CONTROL               PSP-CONTROL 
# 99418                     34026                     34026 

# $Sex

# ALL  FEMALE    MALE 
# 1100014  115002  115002 

# $Study

# MAYO   MSSM ROSMAP 
# 442338 653920 233760 

## Let's see how some characteristics meta data analyzed for ROSMAP study
#ifelse(dge.allMerged$Study == "ROSMAP", lapply(dge.allMerged[,c(1,2,3,17,18)], table), lapply(dge.allMerged[,c(1,2,3,17,18)], length)) # Didn't work

lapply(dge.allMerged[dge.allMerged$Study == "ROSMAP",c(1,2,3,17,18)], table)
# $Model

# APOE4       Diagnosis   Diagnosis.AOD   Diagnosis.Sex SourceDiagnosis 
# 46752           31168           31168           77920           46752 

# $Tissue

# DLPFC 
# 233760 

# $Comparison

# 1-0                       2-0                       2-1                       4-1                       4-2                AD-CONTROL 
# 15584                     15584                     31168                     15584                     15584                     77920 
# AD-CONTROL.IN.FEMALE-MALE                       AOD               FEMALE-MALE 
# 15584                     15584                     31168 

# $Sex

# ALL FEMALE   MALE 
# 202592  15584  15584 

# $Study

# ROSMAP 
# 233760 

## Subset ROSMAP data
dge.rosmap <- dge.allMerged[dge.allMerged$Study == "ROSMAP", ]
dim(dge.rosmap) # 233760     18
head(dge.rosmap)

table(dge.rosmap$Model)

table(dge.rosmap$Direction)
#  DOWN   NONE     UP 
#   485 232787    488 

sum(dge.rosmap$adj.P.Val < 0.05) # 13460

dge.rosmap[dge.rosmap$hgnc_symbol %in% c(grep("CTS|ATG|NEF|RAB", value = T, dge.rosmap$hgnc_symbol), "SQSTM1", "UVRAG", "MAPT", "BECN1"), -c(4,6,7,8,12,15,16)][1:20,]
dge.rosmap[dge.rosmap$hgnc_symbol %in% c(grep("CTS|ATG|NEF|RAB", value = T, dge.rosmap$hgnc_symbol), "SQSTM1", "UVRAG", "MAPT", "BECN1"), -c(4,6,7,8,12,15,16)][21:70,]
dge.rosmap[dge.rosmap$hgnc_symbol %in% c(grep("CTS|ATG|NEF|RAB", value = T, dge.rosmap$hgnc_symbol), "SQSTM1", "UVRAG", "MAPT", "BECN1"), -c(4,6,7,8,12,15,16)][71:120,]
dge.rosmap[dge.rosmap$hgnc_symbol %in% c(grep("CTS|ATG|NEF|RAB", value = T, dge.rosmap$hgnc_symbol), "SQSTM1", "UVRAG", "MAPT", "BECN1"), -c(4,6,7,8,12,15,16)][121:150,]
dge.rosmap[dge.rosmap$hgnc_symbol %in% c(grep("CTS|ATG|NEF|RAB", value = T, dge.rosmap$hgnc_symbol), "SQSTM1", "UVRAG", "MAPT", "BECN1"), -c(4,6,7,8,12,15,16)][151:200,]
dge.rosmap[dge.rosmap$hgnc_symbol %in% c(grep("CTS|ATG|NEF|RAB", value = T, dge.rosmap$hgnc_symbol), "SQSTM1", "UVRAG", "MAPT", "BECN1"), -c(4,6,7,8,12,15,16)][201:250,]

dge.rosmap[dge.rosmap$hgnc_symbol %in% c(grep("CTS|ATG|NEF|RAB|ULK|LAMP", value = T, dge.rosmap$hgnc_symbol), "SQSTM1", "UVRAG", "MAPT", "BECN1") & 
             dge.rosmap$Comparison == "AD-CONTROL", -c(4,6,7,8,12,15,16)][1:90,]
dge.rosmap[dge.rosmap$hgnc_symbol %in% c(grep("CTS|ATG|NEF|RAB|ULK|LAMP", value = T, dge.rosmap$hgnc_symbol), "SQSTM1", "UVRAG", "MAPT", "BECN1") & 
             dge.rosmap$Comparison == "AD-CONTROL", -c(4,6,7,8,12,15,16)][91:180,]
dge.rosmap[dge.rosmap$hgnc_symbol %in% c(grep("CTS|ATG|NEF|RAB|ULK|LAMP", value = T, dge.rosmap$hgnc_symbol), "SQSTM1", "UVRAG", "MAPT", "BECN1") & 
             dge.rosmap$Comparison == "AD-CONTROL", -c(4,6,7,8,12,15,16)][181:270,]
dge.rosmap[dge.rosmap$hgnc_symbol %in% c(grep("CTS|ATG|NEF|RAB|ULK|LAMP", value = T, dge.rosmap$hgnc_symbol), "SQSTM1", "UVRAG", "MAPT", "BECN1") & 
             dge.rosmap$Comparison == "AD-CONTROL", -c(4,6,7,8,12,15,16)][271:360,]
dge.rosmap[dge.rosmap$hgnc_symbol %in% c(grep("CTS|ATG|NEF|RAB|ULK|LAMP", value = T, dge.rosmap$hgnc_symbol), "SQSTM1", "UVRAG", "MAPT", "BECN1") & 
             dge.rosmap$Comparison == "AD-CONTROL", -c(4,6,7,8,12,15,16)][361:450,]
dge.rosmap[dge.rosmap$hgnc_symbol %in% c(grep("CTS|ATG|NEF|RAB|ULK|LAMP", value = T, dge.rosmap$hgnc_symbol), "SQSTM1", "UVRAG", "MAPT", "BECN1") & 
             dge.rosmap$Comparison == "AD-CONTROL", -c(4,6,7,8,12,15,16)][451:540,]
dge.rosmap[dge.rosmap$hgnc_symbol %in% c(grep("CTS|ATG|NEF|RAB|ULK|LAMP", value = T, dge.rosmap$hgnc_symbol), "SQSTM1", "UVRAG", "MAPT", "BECN1") & 
             dge.rosmap$Comparison == "AD-CONTROL", -c(4,6,7,8,12,15,16)][541:630,]
dge.rosmap[dge.rosmap$hgnc_symbol %in% c(grep("CTS|ATG|NEF|RAB|ULK|LAMP", value = T, dge.rosmap$hgnc_symbol), "SQSTM1", "UVRAG", "MAPT", "BECN1") & 
             dge.rosmap$Comparison == "AD-CONTROL", -c(4,6,7,8,12,15,16)][631:720,]


dge.rosmap[dge.rosmap$hgnc_symbol %in% "BECN1" & 
             dge.rosmap$Comparison == "AD-CONTROL", -c(4,6,7,8,12,15,16)]
dge.rosmap[dge.rosmap$hgnc_symbol %in% "BECN1", -c(4,6,7,8,12,15,16)]


lapply(dge.rosmap[, c(1,2,3,17,18)], table)
# $Model

# APOE4       Diagnosis   Diagnosis.AOD   Diagnosis.Sex SourceDiagnosis 
# 46752           31168           31168           77920           46752 

# $Tissue

# DLPFC 
# 233760 

# $Comparison

# 1-0                       2-0                       2-1                       4-1                       4-2                AD-CONTROL 
# 15584                     15584                     31168                     15584                     15584                     77920 
# AD-CONTROL.IN.FEMALE-MALE                       AOD               FEMALE-MALE 
# 15584                     15584                     31168 

# $Sex

# ALL FEMALE   MALE 
# 202592  15584  15584 

# $Study

# ROSMAP 
# 233760 

## 5 different models and 9 diff. comparisons. Will perform GSEA on all 9 comparisons. 
## Will make 9 diff DGE results files (subsets from dge.ROSMAP) representing 9 diff comparisons. Some comparisons overlap with the Model testes for e.g. 
## AD-CONTROL were tested for models "Diagnosis", "Diagnosis.AOD" and "Diagnosis.Sex". See below.

table(dge.rosmap$Model[dge.rosmap$Comparison == "AD-CONTROL"])
#     Diagnosis Diagnosis.AOD Diagnosis.Sex 
#         15584         15584         46752

levels(as.factor(dge.rosmap$Model))
# "APOE4"           "Diagnosis"       "Diagnosis.AOD"   "Diagnosis.Sex"   "SourceDiagnosis"
table(dge.rosmap$Comparison[dge.rosmap$Model %in% levels(as.factor(dge.rosmap$Model))])
# 1-0                       2-0                       2-1                       4-1                       4-2                AD-CONTROL 
# 15584                     15584                     31168                     15584                     15584                     77920 
# AD-CONTROL.IN.FEMALE-MALE                       AOD               FEMALE-MALE 
# 15584                     15584                     31168 

table(dge.rosmap$Model[dge.rosmap$Comparison %in% levels(as.factor(dge.rosmap$Comparison))])
#           APOE4       Diagnosis   Diagnosis.AOD   Diagnosis.Sex SourceDiagnosis 
#           46752           31168           31168           77920           46752

table(dge.rosmap$Model[dge.rosmap$Comparison == "1-0"])
# APOE4 
# 15584 

table(dge.rosmap$Model[dge.rosmap$Comparison == "2-0"])
# APOE4 
# 15584 

table(dge.rosmap$Model[dge.rosmap$Comparison == "2-1"])
# APOE4 SourceDiagnosis
# 15584           15584

table(dge.rosmap$Model[dge.rosmap$Comparison == "4-1"])
# SourceDiagnosis 
#           15584

table(dge.rosmap$Model[dge.rosmap$Comparison == "4-2"])
# SourceDiagnosis 
#           15584

## 1. Model: Diagnosis; Comparison: AD vs Control
dge.rosmap_ADvsCtr_dx <- dge.rosmap[dge.rosmap$Model == "Diagnosis" & dge.rosmap$Comparison == "AD-CONTROL", -c(6,7,8,12,15,16)] # dx stands for diagnosis
head(dge.rosmap_ADvsCtr_dx)
tail(dge.rosmap_ADvsCtr_dx)
dim(dge.rosmap_ADvsCtr_dx) # 15584    11
length(unique(dge.rosmap_ADvsCtr_dx$hgnc_symbol)) # 14386
sum(is.na(dge.rosmap_ADvsCtr_dx$hgnc_symbol)) # 0
sum(dge.rosmap_ADvsCtr_dx$hgnc_symbol == "") # 1199 + 14386 = 15585 (14386 unique has "" as a unique gene symbol and hence total is 15585)

length(unique(dge.rosmap_ADvsCtr_dx$ensembl_gene_id)) # 15582

## Import APList110 prepared on 12/18/18 (Autophagy genesets: Curated)
apList110_Hs <- loadWorkbook("~/NKI/Gene Sets/R dir_GeneSets/geneSets110_Hs withNum_121818.xlsx")
apList110_Hs <- readWorksheet(apList110_Hs, sheet = 1, header = T)
dim(apList110_Hs) # 891 110
head(apList110_Hs)
tail(apList110_Hs)
class(apList110_Hs) # "data.frame"
sapply(apList110_Hs, class) # All character

t.stat_ADvsCtr_dx <- as.numeric(dge.rosmap_ADvsCtr_dx$t)
names(t.stat_ADvsCtr_dx) <- dge.rosmap_ADvsCtr_dx$hgnc_symbol
# One-liner: t.stat_ADvsCtr_dx <- setNames(dge.rosmap_ADvsCtr_dx$t, dge.rosmap_ADvsCtr_dx$hgnc_symbol)

t.stat_ADvsCtr_dx[1:10]

idx.gs110_ADvsCtr_dx <- ids2indices(apList110_Hs, id = names(t.stat_ADvsCtr_dx))
str(idx.gs110_ADvsCtr_dx)

# fit CAMAERA
cam_gs110_ADvsCtr_dx <- cameraPR(statistic = t.stat_ADvsCtr_dx, index = idx.gs110_ADvsCtr_dx, inter.gene.cor = 0.01)
cam_gs110_ADvsCtr_dx
dir.create("GSEA")
dir.create("GSEA/dx_model")
write.table(cam_gs110_ADvsCtr_dx, file = "GSEA/dx_model/cam_gs110_ADvsCtr_dx.txt", sep = "\t", 
            row.names = T, col.names = NA, quote = F)

# Barcode plot
barcodeplot(t.stat_ADvsCtr_dx, index = idx.gs110_ADvsCtr_dx$MTORC1.substrates.20., main = "mTORC1 substrates (n=20) \nAD vs Ctr [Dx model]-DLPFC")
title(main = "P = 7.16x10-05; FDR = 0.0078", line = 0.1, cex.main = 0.7)
dev.print(pdf, "GSEA/dx_model/Rplot_barcode_mTORC1.Subs_ADvsCtr_dx.pdf", height = 6, width = 8)

barcodeplot(t.stat_ADvsCtr_dx, index = idx.gs110_ADvsCtr_dx$MTOR.cascade.32., main = "mTOR cascade (n=32) \nAD vs Ctr [Dx model]-DLPFC")
title(main = "P = 4.63x10-03; FDR = 0.15", line = 0.1, cex.main = 0.7)
dev.print(pdf, "GSEA/dx_model/Rplot_barcode_mTOR.cascade_ADvsCtr_dx.pdf", height = 6, width = 8)

barcodeplot(t.stat_ADvsCtr_dx, index = idx.gs110_ADvsCtr_dx$mTOR.s.transcriptional.factor.11., main = "mTOR' cascade's Transcriptional Factors (n=11) \nAD vs Ctr [Dx model]-DLPFC")
title(main = "P = 0.001; FDR = 0.16", line = 0.1, cex.main = 0.7)
dev.print(pdf, "GSEA/dx_model/Rplot_barcode_mTOR.TFs_ADvsCtr_dx.pdf", height = 6, width = 8)

barcodeplot(t.stat_ADvsCtr_dx, index = idx.gs110_ADvsCtr_dx$Lysosomal.Catabolic.enzymes.38., main = "Lysosomal Catabolic Enzymes (n=38): \nAD vs Ctr [Dx model]-DLPFC")
title(main = "P = 1.03x10-02; FDR = 0.16", line = 0.1, cex.main = 0.7)
dev.print(pdf, "GSEA/dx_model/Rplot_barcode_lysoCatEnz_ADvsCtr_dx.pdf", height = 6, width = 8)

barcodeplot(t.stat_ADvsCtr_dx, index = idx.gs110_ADvsCtr_dx$vacuolar.proton.transporting.V.type.ATPase.complex.26., main = "vATPase complex (n=26): \nAD vs Ctr [Dx model]-DLPFC")
title(main = "P = 0.014; FDR = 0.17", line = 0.1, cex.main = 0.7)
dev.print(pdf, "GSEA/dx_model/Rplot_barcode_vATPaseCmplx_ADvsCtr_dx.pdf", height = 6, width = 8)



## 2 Model: Diagnosis.AOD; Comparison: AD vs Control
dge.rosmap_ADvsCtr_dx.AOD <- dge.rosmap[dge.rosmap$Model == "Diagnosis.AOD" & dge.rosmap$Comparison == "AD-CONTROL", -c(6,7,8,12,15,16)] # dx stands for diagnosis
head(dge.rosmap_ADvsCtr_dx.AOD)
tail(dge.rosmap_ADvsCtr_dx.AOD)
dim(dge.rosmap_ADvsCtr_dx.AOD) # 15584    12
length(unique(dge.rosmap_ADvsCtr_dx.AOD$hgnc_symbol)) # 14386
sum(is.na(dge.rosmap_ADvsCtr_dx.AOD$hgnc_symbol)) # 0
sum(dge.rosmap_ADvsCtr_dx.AOD$hgnc_symbol == "") # 1199 + 14386 = 15585 (14386 unique has "" as a unique gene symbol and hence total is 15585)

length(unique(dge.rosmap_ADvsCtr_dx.AOD$ensembl_gene_id)) # 15582

#
t.stat_ADvsCtr_dx.AOD <- setNames(dge.rosmap_ADvsCtr_dx.AOD$t, dge.rosmap_ADvsCtr_dx.AOD$hgnc_symbol)

t.stat_ADvsCtr_dx.AOD[1:10]

idx.gs110_ADvsCtr_dx.AOD <- ids2indices(apList110_Hs, id = names(t.stat_ADvsCtr_dx.AOD))
str(idx.gs110_ADvsCtr_dx.AOD)

# fit CAMAERA
cam_gs110_ADvsCtr_dx.AOD <- cameraPR(statistic = t.stat_ADvsCtr_dx.AOD, index = idx.gs110_ADvsCtr_dx.AOD, inter.gene.cor = 0.01)
cam_gs110_ADvsCtr_dx.AOD
dir.create("GSEA/dx.AOD_model")
write.table(cam_gs110_ADvsCtr_dx.AOD, file = "GSEA/dx.AOD_model/cam_gs110_ADvsCtr_dx.AOD.txt", sep = "\t", 
            row.names = T, col.names = NA, quote = F)


# Barcode plot
barcodeplot(t.stat_ADvsCtr_dx.AOD, index = idx.gs110_ADvsCtr_dx.AOD$MTORC1.substrates.20., main = "mTORC1 substrates (n=20) \nAD vs Ctr [Dx.AOD model]-DLPFC")
title(main = "P = 9.27x10-05; FDR = 0.01", cex.main = 0.8, line = 0.01)
dev.print(pdf, "GSEA/dx.AOD_model/Rplot_barcode_mTORC1.subs_ADvsCtr_dx.AOD.pdf", height = 6, width = 8)

barcodeplot(t.stat_ADvsCtr_dx.AOD, index = idx.gs110_ADvsCtr_dx.AOD$MTOR.cascade.32., main = "mTOR Cascade (n=32) \nAD vs Ctr [Dx.AOD model]-DLPFC")
title(main = "P = 4.89x10-03; FDR = 0.13", cex.main = 0.8, line = 0.01)
dev.print(pdf, "GSEA/dx.AOD_model/Rplot_barcode_mTOR.cascade_ADvsCtr_dx.AOD.pdf", height = 6, width = 8)

barcodeplot(t.stat_ADvsCtr_dx.AOD, index = idx.gs110_ADvsCtr_dx.AOD$Lysosomal.Catabolic.enzymes.38., main = "Lysosomal Catabolic Enzymes (n=38): \nAD vs Ctr [Dx.AOD model]-DLPFC")
title(main = "P = 0.01; FDR = 0.16", line = 0.01, cex.main = 0.8)
dev.print(pdf, "GSEA/dx.AOD_model/Rplot_barcode_lysoCatEnz_ADvsCtr_dx.pdf", height = 6, width = 8)

barcodeplot(t.stat_ADvsCtr_dx.AOD, index = idx.gs110_ADvsCtr_dx.AOD$vacuolar.proton.transporting.V.type.ATPase.complex.26., main = "vATPase complex (n=26): \nAD vs Ctr [Dx.AOD model]-DLPFC")
title(main = "P = 0.011; FDR = 0.16", line = 0.01, cex.main = 0.8)
dev.print(pdf, "GSEA/dx.AOD_model/Rplot_barcode_vATPaseCmplx_ADvsCtr_dx.pdf", height = 6, width = 8)



## 3 Model: Diagnosis.SEX; Comparison: AD vs Control
dge.rosmap_ADvsCtr_dx.sex <- dge.rosmap[dge.rosmap$Model == "Diagnosis.Sex" & dge.rosmap$Comparison == "AD-CONTROL", -c(6,7,8,12,15,16)] # dx stands for diagnosis
head(dge.rosmap_ADvsCtr_dx.sex)
tail(dge.rosmap_ADvsCtr_dx.sex)
dim(dge.rosmap_ADvsCtr_dx.sex) # 46752    12
15584*3 # 46752 # One each for "male", "female" and "all" sex
#length(unique(dge.rosmap_ADvsCtr_dx.sex$hgnc_symbol)) # 14386
#sum(is.na(dge.rosmap_ADvsCtr_dx.AOD$hgnc_symbol)) # 0
sum(dge.rosmap_ADvsCtr_dx.sex$hgnc_symbol == "") # 3597
3597/3 # 1199

#length(unique(dge.rosmap_ADvsCtr_dx.AOD$ensembl_gene_id)) # 15582

# Subsetting for each sex category
# 3.1) Male
dge.rosmap_ADvsCtr_dx.sexM <- dge.rosmap[dge.rosmap$Model == "Diagnosis.Sex" & dge.rosmap$Comparison == "AD-CONTROL" 
                                         & dge.rosmap$Sex == "MALE", -c(6,7,8,12,15,16)] # dx stands for diagnosis
head(dge.rosmap_ADvsCtr_dx.sexM)
tail(dge.rosmap_ADvsCtr_dx.sexM)
dim(dge.rosmap_ADvsCtr_dx.sexM) # 15584    12

t.stat_ADvsCtr_dx.sexM <- setNames(dge.rosmap_ADvsCtr_dx.sexM$t, dge.rosmap_ADvsCtr_dx.sexM$hgnc_symbol)

t.stat_ADvsCtr_dx.sexM[1:10]

idx.gs110_ADvsCtr_dx.sexM <- ids2indices(apList110_Hs, id = names(t.stat_ADvsCtr_dx.sexM))
str(idx.gs110_ADvsCtr_dx.sexM)

# fit CAMAERA
cam_gs110_ADvsCtr_dx.sexM <- cameraPR(statistic = t.stat_ADvsCtr_dx.sexM, index = idx.gs110_ADvsCtr_dx.sexM, inter.gene.cor = 0.01)
cam_gs110_ADvsCtr_dx.sexM
dir.create("GSEA/dx.sex_model")
write.table(cam_gs110_ADvsCtr_dx.sexM, file = "GSEA/dx.sex_model/cam_gs110_ADvsCtr_dx.sexM.txt", sep = "\t", 
            row.names = T, col.names = NA, quote = F)
# Barcode plot
#barcodeplot(t.stat_ADvsCtr_dx.AOD, index = idx.gs110_ADvsCtr_dx.AOD$MTORC1.substrates.20., main = "mTORC1 substrates (n=20): AD vs Ctr from Dx.AOD model-DLPFC")

# 3.2) Female
dge.rosmap_ADvsCtr_dx.sexF <- dge.rosmap[dge.rosmap$Model == "Diagnosis.Sex" & dge.rosmap$Comparison == "AD-CONTROL" 
                                         & dge.rosmap$Sex == "FEMALE", -c(6,7,8,12,15,16)] # dx stands for diagnosis
head(dge.rosmap_ADvsCtr_dx.sexF)
tail(dge.rosmap_ADvsCtr_dx.sexF)
dim(dge.rosmap_ADvsCtr_dx.sexF) # 15584    12

t.stat_ADvsCtr_dx.sexF <- setNames(dge.rosmap_ADvsCtr_dx.sexF$t, dge.rosmap_ADvsCtr_dx.sexF$hgnc_symbol)

t.stat_ADvsCtr_dx.sexF[1:10]

idx.gs110_ADvsCtr_dx.sexF <- ids2indices(apList110_Hs, id = names(t.stat_ADvsCtr_dx.sexF))
str(idx.gs110_ADvsCtr_dx.sexF)

# fit CAMAERA
cam_gs110_ADvsCtr_dx.sexF <- cameraPR(statistic = t.stat_ADvsCtr_dx.sexF, index = idx.gs110_ADvsCtr_dx.sexF, inter.gene.cor = 0.01)
cam_gs110_ADvsCtr_dx.sexF
write.table(cam_gs110_ADvsCtr_dx.sexF, file = "GSEA/dx.sex_model/cam_gs110_ADvsCtr_dx.sexF.txt", sep = "\t", 
            row.names = T, col.names = NA, quote = F)

# Barcode plot
barcodeplot(t.stat_ADvsCtr_dx.sexF, index = idx.gs110_ADvsCtr_dx.sexF$MTORC1.substrates.20., main = "mTORC1 substrates (n=20) \nAD vs Ctr from Dx.Sex(F) model-DLPFC")
title(main = "P = 0.00011; FDR = 0.0068", line = 0.01, cex.main = 0.8)
dev.print(pdf, "GSEA/dx.sex_model/Rplot_barcode_mTORC1.subs_ADvsCtr_dx.sexF.pdf", height = 6, width = 8)

barcodeplot(t.stat_ADvsCtr_dx.sexF, index = idx.gs110_ADvsCtr_dx.sexF$vacuolar.proton.transporting.V.type.ATPase.complex.26., main = "vATPase cmpx (n=26) \nAD vs Ctr [Dx.Sex(F) model]-DLPFC")
title(main = "P = 0.00012; FDR = 0.0068", line = 0.01, cex.main = 0.8)
dev.print(pdf, "GSEA/dx.sex_model/Rplot_barcode_vATPaseCmplx_ADvsCtr_dx.sexF.pdf", height = 6, width = 8)

barcodeplot(t.stat_ADvsCtr_dx.sexF, index = idx.gs110_ADvsCtr_dx.sexF$Lysosomal.acidification.14., main = "Lysosomal Acidification (n=14) \nAD vs Ctr [Dx.Sex(F) model]-DLPFC")
title(main = "P = 0.00029; FDR = 0.010", line = 0.01, cex.main = 0.8)
dev.print(pdf, "GSEA/dx.sex_model/Rplot_barcode_lysoAcidification_ADvsCtr_dx.sexF.pdf", height = 6, width = 8)

# 3.3) All
dge.rosmap_ADvsCtr_dx.sexAll <- dge.rosmap[dge.rosmap$Model == "Diagnosis.Sex" & dge.rosmap$Comparison == "AD-CONTROL" 
                                         & dge.rosmap$Sex == "ALL", -c(6,7,8,12,15,16)] # dx stands for diagnosis
head(dge.rosmap_ADvsCtr_dx.sexAll)
tail(dge.rosmap_ADvsCtr_dx.sexAll)
dim(dge.rosmap_ADvsCtr_dx.sexAll) # 15584    12

t.stat_ADvsCtr_dx.sexAll <- setNames(dge.rosmap_ADvsCtr_dx.sexAll$t, dge.rosmap_ADvsCtr_dx.sexAll$hgnc_symbol)

t.stat_ADvsCtr_dx.sexAll[1:10]

idx.gs110_ADvsCtr_dx.sexAll <- ids2indices(apList110_Hs, id = names(t.stat_ADvsCtr_dx.sexAll))
str(idx.gs110_ADvsCtr_dx.sexAll)

# fit CAMAERA
cam_gs110_ADvsCtr_dx.sexAll <- cameraPR(statistic = t.stat_ADvsCtr_dx.sexAll, index = idx.gs110_ADvsCtr_dx.sexAll, inter.gene.cor = 0.01)
cam_gs110_ADvsCtr_dx.sexAll
write.table(cam_gs110_ADvsCtr_dx.sexAll, file = "GSEA/dx.sex_model/cam_gs110_ADvsCtr_dx.sexALL.txt", sep = "\t", 
            row.names = T, col.names = NA, quote = F)

cam2_gs110_ADvsCtr_dx.sexAll <- cameraPR(statistic = t.stat_ADvsCtr_dx.sexAll, index = idx.gs110_ADvsCtr_dx.sexAll)
cam2_gs110_ADvsCtr_dx.sexAll

# Barcode plot
barcodeplot(t.stat_ADvsCtr_dx.sexAll, index = idx.gs110_ADvsCtr_dx.sexAll$MTORC1.substrates.20., main = "mTORC1 substrates (n=20) \nAD vs Ctr [Dx.Sex(ALL)] model-DLPFC")
title(main = "P = 0.00020; FDR = 0.023", line = 0.01, cex.main = 0.8)
dev.print(pdf, "GSEA/dx.sex_model/Rplot_barcode_mTORC1.subs_ADvsCtr_dx.sexALL.pdf", height = 6, width = 8)

barcodeplot(t.stat_ADvsCtr_dx.sexAll, index = idx.gs110_ADvsCtr_dx.sexAll$MTOR.cascade.32., main = "mTOR Cascade (n=32) \nAD vs Ctr [Dx.Sex(ALL)] model-DLPFC")
title(main = "P = 0.0029; FDR = 0.081", line = 0.01, cex.main = 0.8)
dev.print(pdf, "GSEA/dx.sex_model/Rplot_barcode_mTOR.cascade_ADvsCtr_dx.sexALL.pdf", height = 6, width = 8)

#barcodeplot(t.stat_ADvsCtr_dx.sexAll, index = idx.gs110_ADvsCtr_dx.sexAll$vacuolar.proton.transporting.V.type.ATPase.complex.26., main = "vATPase cmpx (n=26): AD vs Ctr [Dx.Sex(ALL) model]-DLPFC")
barcodeplot(t.stat_ADvsCtr_dx.sexAll, index = idx.gs110_ADvsCtr_dx.sexAll$Lysosomal.Catabolic.enzymes.38., main = "Lysosomal Catabolic Enzymes (n=38) \nAD vs Ctr [Dx.Sex(ALL) model]-DLPFC")
title(main = "P = 0.0029; FDR = 0.081", line = 0.01, cex.main = 0.8)
dev.print(pdf, "GSEA/dx.sex_model/Rplot_barcode_lysoCatEnz_ADvsCtr_dx.sexALL.pdf", height = 6, width = 8)


## 4 Model: SourceDiagnosis; Comparison: 4-1
dge.rosmap_4vs1_SrcDx <- dge.rosmap[dge.rosmap$Model == "SourceDiagnosis" & dge.rosmap$Comparison == "4-1" , -c(6,7,8,12,15,16)] # dx stands for diagnosis
head(dge.rosmap_4vs1_SrcDx)
tail(dge.rosmap_4vs1_SrcDx)
dim(dge.rosmap_4vs1_SrcDx) # 15584    12

t.stat_4vs1_SrcDx <- setNames(dge.rosmap_4vs1_SrcDx$t, dge.rosmap_4vs1_SrcDx$hgnc_symbol)

t.stat_4vs1_SrcDx[1:10]

idx.gs110_4vs1_SrcDx <- ids2indices(apList110_Hs, id = names(t.stat_4vs1_SrcDx))
str(idx.gs110_4vs1_SrcDx)

# fit CAMAERA
cam_gs110_4vs1_SrcDx <- cameraPR(statistic = t.stat_4vs1_SrcDx, index = idx.gs110_4vs1_SrcDx, inter.gene.cor = 0.01)
cam_gs110_4vs1_SrcDx
dir.create("GSEA/SrcDx_model")
write.table(cam_gs110_4vs1_SrcDx, file = "GSEA/SrcDx_model/cam_gs110_4vs1_SrcDx.txt", sep = "\t", 
            row.names = T, col.names = NA, quote = F)
cam2_gs110_4vs1_SrcDx <- cameraPR(statistic = t.stat_4vs1_SrcDx, index = idx.gs110_4vs1_SrcDx)
cam2_gs110_4vs1_SrcDx

cam3_gs110_4vs1_SrcDx <- cameraPR(statistic = t.stat_4vs1_SrcDx, index = idx.gs110_4vs1_SrcDx, inter.gene.cor = 0.001)
cam3_gs110_4vs1_SrcDx

# Barcode plot
barcodeplot(t.stat_4vs1_SrcDx, index = idx.gs110_4vs1_SrcDx$MTORC1.substrates.20., main = "mTORC1 substrates (n=20) \n4 vs 1 [SourceDiagnosis model]-DLPFC")
title(main = "P = 0.00020; FDR = 0.023", line = 0.01, cex.main = 0.8)
dev.print(pdf, "GSEA/SrcDx_model/Rplot_barcode_mTORC1.subs_4vs1_SrcDx.pdf", height = 6, width = 8)



## 4 Model: SourceDiagnosis; Comparison: 4-2
dge.rosmap_4vs2_SrcDx <- dge.rosmap[dge.rosmap$Model == "SourceDiagnosis" & dge.rosmap$Comparison == "4-2" , -c(6,7,8,12,15,16)] # dx stands for diagnosis
head(dge.rosmap_4vs2_SrcDx)
tail(dge.rosmap_4vs2_SrcDx)
dim(dge.rosmap_4vs2_SrcDx) # 15584    12

t.stat_4vs2_SrcDx <- setNames(dge.rosmap_4vs2_SrcDx$t, dge.rosmap_4vs2_SrcDx$hgnc_symbol)

t.stat_4vs2_SrcDx[1:10]

idx.gs110_4vs2_SrcDx <- ids2indices(apList110_Hs, id = names(t.stat_4vs2_SrcDx))
str(idx.gs110_4vs2_SrcDx)

# fit CAMAERA
cam_gs110_4vs2_SrcDx <- cameraPR(statistic = t.stat_4vs2_SrcDx, index = idx.gs110_4vs2_SrcDx, inter.gene.cor = 0.01)
cam_gs110_4vs2_SrcDx
#dir.create("GSEA/SrcDx_model")
write.table(cam_gs110_4vs2_SrcDx, file = "GSEA/SrcDx_model/cam_gs110_4vs2_SrcDx.txt", sep = "\t", 
            row.names = T, col.names = NA, quote = F)
#cam2_gs110_4vs1_SrcDx <- cameraPR(statistic = t.stat_4vs1_SrcDx, index = idx.gs110_4vs1_SrcDx)
#cam2_gs110_4vs1_SrcDx

#cam3_gs110_4vs1_SrcDx <- cameraPR(statistic = t.stat_4vs1_SrcDx, index = idx.gs110_4vs1_SrcDx, inter.gene.cor = 0.001)
#cam3_gs110_4vs1_SrcDx

# Barcode plot
barcodeplot(t.stat_4vs2_SrcDx, index = idx.gs110_4vs2_SrcDx$MTOR.cascade.32., main = "mTOR Cascade (n=32) \n4 vs 2 [SourceDiagnosis] model-DLPFC")
title(main = "P = 0.00029; FDR = 0.031", line = 0.01, cex.main = 0.8)
dev.print(pdf, "GSEA/SrcDx_model/Rplot_barcode_mTOR.cascade_4vs2_SrcDx.pdf", height = 6, width = 8)

barcodeplot(t.stat_4vs2_SrcDx, index = idx.gs110_4vs2_SrcDx$endosome.to.lysosome.transport.10., main = "Endosome to Lysosome Transport (n=10) \n4 vs 2 [SourceDiagnosis] model-DLPFC")
title(main = "P = 0.00092; FDR = 0.031", line = 0.01, cex.main = 0.8)
dev.print(pdf, "GSEA/SrcDx_model/Rplot_barcode_endoToLysoTrans_4vs2_SrcDx.pdf", height = 6, width = 8)

barcodeplot(t.stat_4vs2_SrcDx, index = idx.gs110_4vs2_SrcDx$All.PI3K.genes.44., main = "All PI3K genes (n=44) \n4 vs 2 [SourceDiagnosis] model-DLPFC")
title(main = "P = 0.00096; FDR = 0.031", line = 0.01, cex.main = 0.8)
dev.print(pdf, "GSEA/SrcDx_model/Rplot_barcode_allPI3Kgenes_4vs2_SrcDx.pdf", height = 6, width = 8)

barcodeplot(t.stat_4vs2_SrcDx, index = idx.gs110_4vs2_SrcDx$PI3K.complex.12., main = "PI3K Complex (n=12) \n4 vs 2 [SourceDiagnosis] model-DLPFC")
title(main = "P = 0.0011; FDR = 0.031", line = 0.01, cex.main = 0.8)
dev.print(pdf, "GSEA/SrcDx_model/Rplot_barcode_PI3KComplex_4vs2_SrcDx.pdf", height = 6, width = 8)

barcodeplot(t.stat_4vs2_SrcDx, index = idx.gs110_4vs2_SrcDx$MTORC1.substrates.20., main = "mTORC1 substrates (n=19) \n4 vs 2 [SourceDiagnosis] model-DLPFC")
title(main = "P = 0.0044; FDR = 0.082", line = 0.01, cex.main = 0.8)
dev.print(pdf, "GSEA/SrcDx_model/Rplot_barcode_mTORC1.subs_4vs2_SrcDx.pdf", height = 6, width = 8)



## 5 Model: SourceDiagnosis; Comparison: 2-1
dge.rosmap_2vs1_SrcDx <- dge.rosmap[dge.rosmap$Model == "SourceDiagnosis" & dge.rosmap$Comparison == "2-1" , -c(6,7,8,12,15,16)] # dx stands for diagnosis
head(dge.rosmap_2vs1_SrcDx)
tail(dge.rosmap_2vs1_SrcDx)
dim(dge.rosmap_2vs1_SrcDx) # 15584    12

t.stat_2vs1_SrcDx <- setNames(dge.rosmap_2vs1_SrcDx$t, dge.rosmap_2vs1_SrcDx$hgnc_symbol)

t.stat_2vs1_SrcDx[1:10]

idx.gs110_2vs1_SrcDx <- ids2indices(apList110_Hs, id = names(t.stat_2vs1_SrcDx))
str(idx.gs110_2vs1_SrcDx)

# fit CAMAERA
cam_gs110_2vs1_SrcDx <- cameraPR(statistic = t.stat_2vs1_SrcDx, index = idx.gs110_2vs1_SrcDx, inter.gene.cor = 0.01)
cam_gs110_2vs1_SrcDx
#dir.create("GSEA/SrcDx_model")
write.table(cam_gs110_2vs1_SrcDx, file = "GSEA/SrcDx_model/cam_gs110_2vs1_SrcDx.txt", sep = "\t", 
            row.names = T, col.names = NA, quote = F)
#cam2_gs110_4vs1_SrcDx <- cameraPR(statistic = t.stat_4vs1_SrcDx, index = idx.gs110_4vs1_SrcDx)
#cam2_gs110_4vs1_SrcDx

#cam3_gs110_4vs1_SrcDx <- cameraPR(statistic = t.stat_4vs1_SrcDx, index = idx.gs110_4vs1_SrcDx, inter.gene.cor = 0.001)
#cam3_gs110_4vs1_SrcDx

# Barcode plot
#barcodeplot(t.stat_4vs2_SrcDx, index = idx.gs110_4vs2_SrcDx$MTOR.cascade.32., main = "mTOR Cascade (n=32) \n4 vs 2 [SourceDiagnosis] model-DLPFC")
#title(main = "P = 0.00029; FDR = 0.031", line = 0.01, cex.main = 0.8)
#dev.print(pdf, "GSEA/SrcDx_model/Rplot_barcode_mTOR.cascade_4vs2_SrcDx.pdf", height = 6, width = 8)



## 6 Model: APOE4; Comparison: 2-0
dge.rosmap_2vs0_apoe4 <- dge.rosmap[dge.rosmap$Model == "APOE4" & dge.rosmap$Comparison == "2-0" , -c(6,7,8,12,15,16)] # dx stands for diagnosis
head(dge.rosmap_2vs0_apoe4)
tail(dge.rosmap_2vs0_apoe4)
dim(dge.rosmap_2vs0_apoe4) # 15584    12

t.stat_2vs0_apoe4 <- setNames(dge.rosmap_2vs0_apoe4$t, dge.rosmap_2vs0_apoe4$hgnc_symbol)

t.stat_2vs0_apoe4[1:10]

idx.gs110_2vs0_apoe4 <- ids2indices(apList110_Hs, id = names(t.stat_2vs0_apoe4))
str(idx.gs110_2vs0_apoe4)

# fit CAMAERA
cam_gs110_2vs0_apoe4 <- cameraPR(statistic = t.stat_2vs0_apoe4, index = idx.gs110_2vs0_apoe4, inter.gene.cor = 0.01)
cam_gs110_2vs0_apoe4
dir.create("GSEA/APOE4_model")
write.table(cam_gs110_2vs0_apoe4, file = "GSEA/APOE4_model/cam_gs110_2vs0_APOE4.txt", sep = "\t", 
            row.names = T, col.names = NA, quote = F)
#cam2_gs110_4vs1_SrcDx <- cameraPR(statistic = t.stat_4vs1_SrcDx, index = idx.gs110_4vs1_SrcDx)
#cam2_gs110_4vs1_SrcDx

#cam3_gs110_4vs1_SrcDx <- cameraPR(statistic = t.stat_4vs1_SrcDx, index = idx.gs110_4vs1_SrcDx, inter.gene.cor = 0.001)
#cam3_gs110_4vs1_SrcDx

# Barcode plot
barcodeplot(t.stat_2vs0_apoe4, index = idx.gs110_2vs0_apoe4$Transcriptional.factors.104., main = "Transcriptional Factors (n=104) \n2 vs 0 [APOE4 model]-DLPFC")
title(main = "P = 4.26x10-06; FDR = 0.00046", line = 0.01, cex.main = 0.8)
dev.print(pdf, "GSEA/APOE4_model/Rplot_barcode_TFs_2vs0_APOE4.pdf", height = 6, width = 8)

barcodeplot(t.stat_2vs0_apoe4, index = idx.gs110_2vs0_apoe4$NfKb.pathway.23., main = "NFKB pathway (n=23) \n2 vs 0 [APOE4 model]-DLPFC")
title(main = "P = 1.13x10-05; FDR = 0.00062", line = 0.01, cex.main = 0.8)
dev.print(pdf, "GSEA/APOE4_model/Rplot_barcode_NFKBpathway_2vs0_APOE4.pdf", height = 6, width = 8)


## 7 Model: APOE4; Comparison: 1-0
dge.rosmap_1vs0_apoe4 <- dge.rosmap[dge.rosmap$Model == "APOE4" & dge.rosmap$Comparison == "1-0" , -c(6,7,8,12,15,16)] # dx stands for diagnosis
head(dge.rosmap_1vs0_apoe4)
tail(dge.rosmap_1vs0_apoe4)
dim(dge.rosmap_1vs0_apoe4) # 15584    12

t.stat_1vs0_apoe4 <- setNames(dge.rosmap_1vs0_apoe4$t, dge.rosmap_1vs0_apoe4$hgnc_symbol)

t.stat_1vs0_apoe4[1:10]

idx.gs110_1vs0_apoe4 <- ids2indices(apList110_Hs, id = names(t.stat_1vs0_apoe4))
str(idx.gs110_1vs0_apoe4)

# fit CAMAERA
cam_gs110_1vs0_apoe4 <- cameraPR(statistic = t.stat_1vs0_apoe4, index = idx.gs110_1vs0_apoe4, inter.gene.cor = 0.01)
cam_gs110_1vs0_apoe4
#dir.create("GSEA/APOE4_model")
write.table(cam_gs110_1vs0_apoe4, file = "GSEA/APOE4_model/cam_gs110_1vs0_APOE4.txt", sep = "\t", 
            row.names = T, col.names = NA, quote = F)
#cam2_gs110_4vs1_SrcDx <- cameraPR(statistic = t.stat_4vs1_SrcDx, index = idx.gs110_4vs1_SrcDx)
#cam2_gs110_4vs1_SrcDx

#cam3_gs110_4vs1_SrcDx <- cameraPR(statistic = t.stat_4vs1_SrcDx, index = idx.gs110_4vs1_SrcDx, inter.gene.cor = 0.001)
#cam3_gs110_4vs1_SrcDx

# Barcode plot

## 8 Model: APOE4; Comparison: 2-1
dge.rosmap_2vs1_apoe4 <- dge.rosmap[dge.rosmap$Model == "APOE4" & dge.rosmap$Comparison == "2-1" , -c(6,7,8,12,15,16)] # dx stands for diagnosis
head(dge.rosmap_2vs1_apoe4)
tail(dge.rosmap_2vs1_apoe4)
dim(dge.rosmap_2vs1_apoe4) # 15584    12

t.stat_2vs1_apoe4 <- setNames(dge.rosmap_2vs1_apoe4$t, dge.rosmap_2vs1_apoe4$hgnc_symbol)

t.stat_2vs1_apoe4[1:10]

idx.gs110_2vs1_apoe4 <- ids2indices(apList110_Hs, id = names(t.stat_2vs1_apoe4))
str(idx.gs110_2vs1_apoe4)

# fit CAMAERA
cam_gs110_2vs1_apoe4 <- cameraPR(statistic = t.stat_2vs1_apoe4, index = idx.gs110_2vs1_apoe4, inter.gene.cor = 0.01)
cam_gs110_2vs1_apoe4
#dir.create("GSEA/APOE4_model")
write.table(cam_gs110_2vs1_apoe4, file = "GSEA/APOE4_model/cam_gs110_2vs1_APOE4.txt", sep = "\t", 
            row.names = T, col.names = NA, quote = F)
#cam2_gs110_4vs1_SrcDx <- cameraPR(statistic = t.stat_4vs1_SrcDx, index = idx.gs110_4vs1_SrcDx)
#cam2_gs110_4vs1_SrcDx

#cam3_gs110_4vs1_SrcDx <- cameraPR(statistic = t.stat_4vs1_SrcDx, index = idx.gs110_4vs1_SrcDx, inter.gene.cor = 0.001)
#cam3_gs110_4vs1_SrcDx

# Barcode plot






























































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































