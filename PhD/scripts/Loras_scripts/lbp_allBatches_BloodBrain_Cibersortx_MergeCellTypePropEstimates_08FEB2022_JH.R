#merge all of the cell type estimates from Cibersort into one and then do analyses with the CBC blood count measures that we have 


scp -r /Users/liharl02/Documents/Thesis/LivingBrainProjectDa-ListOfLBPSubjects_DATA_2022-02-01_1517.csv liharl02@chimera.hpc.mssm.edu:/sc/arion/projects/psychgen/lbp/data/emr/

# scp -r /Users/liharl02/Downloads/Postmortem_Brain_Samples_Chart_10JAN2017_LEL.xlsx liharl02@chimera.hpc.mssm.edu:/sc/arion/work/liharl02/
# scp -r /Users/liharl02/Downloads/Postmortem_Brain_Samples_Chart_LEL_8APR2021.csv liharl02@chimera.hpc.mssm.edu:/sc/arion/work/liharl02/


module load R/4.0.3
R

rm(list=ls())
options(stringsAsFactors=F)
library(data.table)
library(openxlsx)
###############################################################################################################

setwd("/sc/arion/projects/psychgen/lbp/files/lbp_cibersortx/alldata")

covs <- readRDS("/sc/arion/projects/psychgen/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_BLOODandBRAIN_CovariatesTable_PlusBankPMI_776Samples_NoBLOOD582_NoBRAIN734_NoBRAIN722_18JAN2022.RDS")

#clinical data
clinical <- read.xlsx("/sc/arion/projects/psychgen/lbp/data/emr/msd494/MSD494LBPMRN112018.xlsx", sheet = 11, colNames = T, rowNames = F)
head(clinical)
#   MEDICAL_RECORD_NUMBER RESULT_DATE_TIME           ORDER_CODE
# 1               1734475        31-DEC-10 3002-ER VENOUS PANEL
# 2               1734475        31-DEC-10 3002-ER VENOUS PANEL
# 3               1734475        31-DEC-10 3002-ER VENOUS PANEL
# 4               1734475        31-DEC-10 3002-ER VENOUS PANEL
# 5               1734475        31-DEC-10 3002-ER VENOUS PANEL
# 6               1734475        31-DEC-10 3002-ER VENOUS PANEL
#                RESULT_CODE VALUE UNIT_OF_MEASURE
# 1 106-WB UREA NITROGEN-VEN    18           Mg/Dl
# 2     116-WB GLUCOSE - VEN   108           Mg/Dl
# 3  110-WB CREATININE - VEN   1.1           Mg/Dl
# 4    104-WB CHLORIDE - VEN   103           Meq/L
# 5     111-HEMATOCRIT - VEN    51               %
# 6   113-WB POTASSIUM - VEN   6.1           Meq/L

mrnmap <- fread("/sc/arion/projects/psychgen/lbp/data/emr/LivingBrainProjectDa-ListOfLBPSubjects_DATA_2022-02-01_1517.csv", header=T)
#this file has personal info of patients

common <- intersect(clinical$MEDICAL_RECORD_NUMBER, mrnmap$mrn)
length(common)
#[1] 167
length(unique(clinical$MEDICAL_RECORD_NUMBER))
#[1] 168
length(unique(mrnmap$mrn))
#[1] 307

labs <- merge(clinical, mrnmap[, c(1,8)], by.x="MEDICAL_RECORD_NUMBER", by.y="mrn", all.x=TRUE)

length(unique(labs$subject_id))
#[1] 168
length(unique(labs$MEDICAL_RECORD_NUMBER))
#[1] 168
length(common)
#[1] 167

labs[!labs$MEDICAL_RECORD_NUMBER %in% common, c("MEDICAL_RECORD_NUMBER", "subject_id")]      
#MEDICAL_RECORD_NUMBER subject_id
#7344419       			<NA>

#the NA counts as a unique value, if you look at unique(labs$subject_id), it's in there

#remove the participant who doesn't have a pt_id

labs2 <- labs[labs$MEDICAL_RECORD_NUMBER %in% common, ]

grep("CBC", unique(labs2$ORDER_CODE), value=T)
#[1] "203-CBC & PLT & DIFF"   "202-CBC & PLT"          "20226-CBC & PLT & DIFF"

cbc <- grep("CBC", unique(labs2$ORDER_CODE), value=T)
labs3 <- labs2[labs2$ORDER_CODE %in% cbc, ]

unique(labs3$RESULT_CODE)

# [1] "2610-MEAN PLT VOLUME"       "2651-NEUTROPHIL %"         
#  [3] "2605-MEAN CORP. VOLUME"     "2693-MONOCYTE #"           
#  [5] "2601-WHITE BLOOD CELL"      "2606-MEAN CORP. HGB"       
#  [7] "2607-MEAN CORP. HGB CONC."  "2653-MONOCYTE %"           
#  [9] "2603-HEMOGLOBIN"            "2692-LYMPHOCYTE #"         
# [11] "5520-NUCLEATE RBC%"         "2691-NEUTROPHIL #"         
# [13] "2694-EOSINOPHIL #"          "2654-EOSINOPHIL %"         
# [15] "2604-HEMATOCRIT"            "2695-BASOPHIL #"           
# [17] "2655-BASOPHIL %"            "2638-PLATELET"             
# [19] "2602-RED BLOOD CELL"        "2608-RED DISTRIB. WIDTH"   
# [21] "5519-NRBC#"                 "2652-LYMPHOCYTE %"         
# [23] "20245-NEUTROPHIL"           "20247-BASOPHIL"            
# [25] "20246-EOSINOPHIL"           "20243-LYMPHOCYTE"          
# [27] "20244-MONOCYTE"             "2640-RBC MORPHOLOGY"       
# [29] "2639-RBC MORPHOLOGY"        "2641-RBC MORPHOLOGY"       
# [31] "2642-RBC MORPHOLOGY"        "20915-ABS NEUTROPHIL COUNT"
# [33] "2611-PLATELET ESTIMATE"     "2644-WBC MORPHOLOGY"       
# [35] "2645-WBC MORPHOLOGY"        "20237-HEMOGLOBIN"          
# [37] "20241-MEAN CORP. HGB CONC." "20243-LYMPHOCYTE %"        
# [39] "20240-MEAN CORP. HGB"       "20248-LYMPHOCYTE #"        
# [41] "20236-RED BLOOD CELL"       "20252-BASOPHIL #"          
# [43] "20254-MEAN PLT VOLUME"      "20244-MONOCYTE %"          
# [45] "20235-WHITE BLOOD CELL"     "20239-MEAN CORP. VOLUME"   
# [47] "20249-MONOCYTE #"           "20246-EOSINOPHIL %"        
# [49] "20250-NEUTROPHIL #"         "20245-NEUTROPHIL %"        
# [51] "20247-BASOPHIL %"           "20251-EOSINOPHIL #"        
# [53] "20238-HEMATOCRIT"           "20242-PLATELET"            
# [55] "20253-RED DISTRIB. WIDTH"   "2616-ATYPICAL LYMPHOCYTES" 
# [57] "2618-BAND CELL"             "2643-RBC MORPHOLOGY"       
# [59] "2620-MYELOCYTE"             "2619-METAMYELOCYTE"     


# lymphocytes -- T cell, B cell, NK cell, Plasma Cell
# monocytes -- differentiate into dendritic cell, macrophage but don't include DCs and macrophages in estimate

dim(labs3)
#[1] 11770     7

labs4 <- labs3[labs3$RESULT_CODE %in% c("2692-LYMPHOCYTE #", "2652-LYMPHOCYTE %", "2651-NEUTROPHIL %", "2691-NEUTROPHIL #", "2653-MONOCYTE %", "2693-MONOCYTE #", "2601-WHITE BLOOD CELL",  "2638-PLATELET", "2694-EOSINOPHIL #", "2654-EOSINOPHIL %", "2695-BASOPHIL #", "2655-BASOPHIL %", "20245-NEUTROPHIL", "20247-BASOPHIL", "20246-EOSINOPHIL", "20243-LYMPHOCYTE", "20244-MONOCYTE", "20915-ABS NEUTROPHIL COUNT", "2611-PLATELET ESTIMATE", "20243-LYMPHOCYTE %", "20248-LYMPHOCYTE #", "20252-BASOPHIL #", "20244-MONOCYTE %", "20235-WHITE BLOOD CELL", "20249-MONOCYTE #", "20246-EOSINOPHIL %", "20250-NEUTROPHIL #", "20245-NEUTROPHIL %", "20247-BASOPHIL %", "20251-EOSINOPHIL #", "20242-PLATELET", "2616-ATYPICAL LYMPHOCYTES"), ]

dim(labs4)
#[1] 6432    7

# lbpBatch1.fqc[, V3:=gsub(".Aligned.out.bam", "", fixed=T, V3)]
# lbpBatch1.ftc[, V1:=tstrsplit(gsub("Sample_", "", V1), split="-", fixed=T, keep=1L)]
# lbpBatch1.ftc = dcast( lbpBatch1.ftc, V1 ~ V2, value.var="V3", fill=NA) #fill empty values with NA when using fun.aggregate

#remove all spaces and change % to percent
labs4 <- data.table(labs4)
labs4[, RESULT_CODE:=gsub("%", "_PCT", fixed=T, RESULT_CODE)]
labs4[, RESULT_CODE:=gsub("#", "_COUNT", fixed=T, RESULT_CODE)]
labs4[, RESULT_CODE:=gsub(" ", "", fixed=T, RESULT_CODE)]


###############################################################################################################
#Cibersort outputs

#1) Lake -- brain single cell
lakeout <- fread("/sc/arion/projects/psychgen/lbp/files/lbp_cibersortx/output_lake/CIBERSORTx_Results.txt")

dim(lakeout)
#[1] 533   9

colnames(lakeout)
# [1] "Mixture"     "GLU"         "GABA"        "AST"         "ODC"        
# [6] "MG"          "P-value"     "Correlation" "RMSE"    

#2) LM22 -- blood facs
lm22out <- fread("/sc/arion/projects/psychgen/lbp/files/lbp_cibersortx/output_LM22/output2/CIBERSORTx_Results.txt")

dim(lm22out)
#[1] 243  26

colnames(lm22out)
#  [1] "Mixture"                      "B cells naive"               
#  [3] "B cells memory"               "Plasma cells"                
#  [5] "T cells CD8"                  "T cells CD4 naive"           
#  [7] "T cells CD4 memory resting"   "T cells CD4 memory activated"
#  [9] "T cells follicular helper"    "T cells gamma delta"         
# [11] "T cells regulatory (Tregs)"   "NK cells resting"            
# [13] "NK cells activated"           "Monocytes"                   
# [15] "Macrophages M0"               "Macrophages M1"              
# [17] "Macrophages M2"               "Dendritic cells resting"     
# [19] "Dendritic cells activated"    "Mast cells resting"          
# [21] "Mast cells activated"         "Eosinophils"                 
# [23] "Neutrophils"                  "P-value"                     
# [25] "Correlation"                  "RMSE"        

#From Noam:
# mer_LM22$Lymphocyte=mer_LM22$`B cells naive` + mer_LM22$`B cells memory` + mer_LM22$`Plasma cells` + mer_LM22$`T cells CD8` + mer_LM22$`T cells CD4 naive` + mer_LM22$`T cells CD4 memory resting` + mer_LM22$`T cells CD4 memory activated` + mer_LM22$`T cells follicular helper` + mer_LM22$`T cells gamma delta` + mer_LM22$`T cells regulatory (Tregs)`


#lymphocytes
lm22out$Lymphocyte_total=lm22out$`B cells naive` + lm22out$`B cells memory` + lm22out$`Plasma cells` + lm22out$`T cells CD8` + lm22out$`T cells CD4 naive` + lm22out$`T cells CD4 memory resting` + lm22out$`T cells CD4 memory activated` + lm22out$`T cells follicular helper` + lm22out$`T cells gamma delta` + lm22out$`T cells regulatory (Tregs)` + lm22out$`NK cells resting` + lm22out$`NK cells activated`

#Also monocytes is their own category not inclusive of macrophages and dendritic cells 


#3) Wilk -- single cell
wilkout <- fread("/sc/arion/projects/psychgen/lbp/files/lbp_cibersortx/output_wilk/output5/CIBERSORTx_Results.txt")

dim(wilkout)
#[1] 243  24

#PB = Peripheral Blood B Cell

colnames(wilkout)
#  [1] "Mixture"                   "RBC"                      
#  [3] "IgM PB"                    "IgG PB"                   
#  [5] "IgA PB"                    "CD14 Monocyte"            
#  [7] "CD8m T"                    "CD4m T"                   
#  [9] "CD4n T"                    "B"                        
# [11] "Platelet"                  "IFN-stim CD4 T"           
# [13] "NK"                        "Neutrophil"               
# [15] "CD16 Monocyte"             "Proliferative Lymphocytes"
# [17] "gd T"                      "pDC"                      
# [19] "Developing Neutrophil"     "SC & Eosinophil"          
# [21] "DC"                        "P-value"                  
# [23] "Correlation"               "RMSE"  

#From Noam:
# mer_wilk$Lymphocyte=mer_wilk$`CD8m T` + mer_wilk$`CD4m T` + mer_wilk$`B` + mer_wilk$`IFN-stim CD4 T` + mer_wilk$`Proliferative Lymphocytes` + mer_wilk$`gd T` + mer_wilk$`IgM PB` + mer_wilk$`IgG PB` + mer_wilk$`IgA PB`
# mer_wilk$Monocyte=mer_wilk$`CD14 Monocyte` + mer_wilk$`CD16 Monocyte`

#CD4nT +  NK missing from lymphocytes above


wilkout$Lymphocyte_total=wilkout$`CD8m T` + wilkout$`CD4m T` + wilkout$`B` + wilkout$`IFN-stim CD4 T` + wilkout$`Proliferative Lymphocytes` + wilkout$`gd T` + wilkout$`IgM PB` + wilkout$`IgG PB` + wilkout$`IgA PB` + wilkout$`NK` + wilkout$`CD4n T`

#gM PB + IgG PB + IgA PB + CD8mT + CD4nT + CD4nT + B + NK + gd T + IFN-stim CD4 T+ Proliferative Lymphocytes = lymphocytes

wilkout$Monocyte_total=wilkout$`CD14 Monocyte` + wilkout$`CD16 Monocyte`


#4) SCP424
scpout <- fread("/sc/arion/projects/psychgen/lbp/files/lbp_cibersortx/output_SCP424/output4/CIBERSORTx_Results.txt")

dim(scpout)
#[1] 243  14

colnames(scpout)
#  [1] "Mixture"                     "CD4+_T_cell"                
#  [3] "Cytotoxic_T_cell"            "Natural_killer_cell"        
#  [5] "CD16+_monocyte"              "CD14+_monocyte"             
#  [7] "Megakaryocyte"               "B_cell"                     
#  [9] "Dendritic_cell"              "Plasmacytoid_dendritic_cell"
# [11] "Unassigned"                  "P-value"                    
# [13] "Correlation"                 "RMSE"      

#From Noam:
# mer_SCP424$Lymphocyte=mer_SCP424$`CD4+_T_cell` + mer_SCP424$`Cytotoxic_T_cell` + mer_SCP424$`B_cell`
# mer_SCP424$Monocyte=mer_SCP424$`CD16+_monocyte` + mer_SCP424$`CD14+_monocyte` 

scpout$Lymphocyte_total=scpout$`CD4+_T_cell` + scpout$`Cytotoxic_T_cell` + scpout$`B_cell` + scpout$`Natural_killer_cell`
scpout$Monocyte_total=scpout$`CD16+_monocyte` + scpout$`CD14+_monocyte` 

#5) Cibersort2
figout <- fread("/sc/arion/projects/psychgen/lbp/files/lbp_cibersortx/output_Fig2b-WholeBlood_RNAseq/output3/CIBERSORTx_Results.txt")

colnames(figout)
# [1] "Mixture"     "T cells CD8" "Monocytes"   "T cells CD4" "NKT cells"  
# [6] "B cells"     "NK cells"    "P-value"     "Correlation" "RMSE"  


dim(figout)
#[1] 243  10

#From Noam:
#mer_NSCLC_PBMC$Lymphocyte=mer_NSCLC_PBMC$`T cells CD8` + mer_NSCLC_PBMC$`T cells CD4` + mer_NSCLC_PBMC$`B cells`

figout$lymphocyte_total=figout$`T cells CD8` + figout$`T cells CD4` + figout$`B cells` + figout$`NKT cells` + figout$`NK cells`

###############################################################################################################

labs4$RESULT_CODE <- as.factor(labs4$RESULT_CODE)
labs4$ORDER_CODE <- as.factor(labs4$ORDER_CODE)
labs4$subject_id <- as.factor(labs4$subject_id)

###############################################################################################################
#ok so, remove the dates such that it's only one date per patient for everything -- need to find the dates closest to the surgeries


master=fread("/sc/arion/projects/psychgen/lbp/files/sema4_bulk_rna_sample_sheet/Bulk_RNA_Isolation_Mastertable_BRAINANDBLOOD_forSEMA4_awcFormatted.tsv",header=T)

#labmaster <- merge(master[, c(6,16)], labs4, by.x="iid", by.y="subject_id", all.y=T)
#can't do this either because of dups


master[grep("PT-0034", master$iid), c("iid", "collection_date")]
#        iid collection_date
# 1: PT-0034      2014-10-13
# 2: PT-0034      2014-09-15
# 3: PT-0034      2014-10-13
# 4: PT-0034      2014-09-15

#unique(labs4[grep("PT-0034", labs4$subject_id),])

unique(labs4[grep("PT-0034", labs4$subject_id), c("subject_id", "RESULT_DATE_TIME")])
#   subject_id RESULT_DATE_TIME
# 1:    PT-0034        19-AUG-14
# 2:    PT-0034        24-SEP-14
# 3:    PT-0034        17-DEC-18

labs4$newDate <- format(strptime(labs4$RESULT_DATE_TIME, format = "%d-%b-%y"), "%Y-%m-%d")

unique(labs4[grep("PT-0034", labs4$subject_id), c("subject_id", "newDate")])
#    subject_id    newDate
# 1:    PT-0034 2014-08-19
# 2:    PT-0034 2014-09-24
# 3:    PT-0034 2018-12-17

master[grep("PT-0034", master$iid), c("iid", "collection_date")]
#        iid collection_date
# 1: PT-0034      2014-10-13
# 2: PT-0034      2014-09-15
# 3: PT-0034      2014-10-13
# 4: PT-0034      2014-09-15

################################################################################
#HERE: next is to choose the date in labs4 closest to the collection date; remove everything that is not in the range of years and then not in the same month go from there -- if the date is not within a few days; after this is done, go back to dcasting, and then merging with the cell type counts

range(labs4$newDate)
#[1] "2005-04-04" "2018-12-19"

labs4$newDate2 <- as.Date(labs4$newDate)

class(master$collection_date)
#[1] "POSIXct" "POSIXt" 

master$newDate <- as.Date(as.character(master$collection_date))

#df[df$date >= "some date" & df$date <= "some date", ]

#labs5 <- labs4[labs4$newDate2 >= "2017-04-01" & labs4$newDate2 <= "2019-03-30", ] 
#4/13/17 - 3/12/19 

#length(unique(master[master$newDate >= "2017-04-01" & master$newDate <= "2019-03-30", ]$iid))
#[1] 56


length(unique(master[master$living==1, ]$iid))
#[1] 172

mliving <- master[master$living==1, ]

mliv2 <- mliving[, colnames(master) %in% c("LBPSEMA4_ID", "iid", "living", "collection_date", "newDate"), with=FALSE]

mliv2[, Year:=tstrsplit(newDate, split="-", fixed=T, keep=1L)]

# lbpBatch1.ftc[, V1:=tstrsplit(gsub("Sample_", "", V1), split="-", fixed=T, keep=1L)]

mliv2$Year <- as.numeric(mliv2$Year)

summary(mliv2$Year)
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
   # 2013    2015    2016    2016    2017    2019       1


labs4[, Year:=tstrsplit(newDate2, split="-", fixed=T, keep=1L)]

labs5 <- labs4[labs4$Year >= "2013" & labs4$Year <= "2019", ] # only results from 2013 to 2019, the range of the collection dates in the master table

labs5[, Month:=tstrsplit(newDate2, split="-", fixed=T, keep=2L)]
labs5[, Day:=tstrsplit(newDate2, split="-", fixed=T, keep=3L)]


#length(intersect(unique(master[master$newDate >= "2017-04-01" & master$newDate <= "2019-03-30", ]$iid), labs5$subject_id))

length(intersect(unique(mliv2$iid), labs5$subject_id))
#[1] 128

ids <- intersect(unique(mliv2$iid), labs5$subject_id)

livlabs <- labs5[subject_id %in% ids, ] #only iids in common with master
livlabs <- livlabs[, 2:ncol(livlabs)]

livlabs <- livlabs[order(subject_id),] 

#need to take blood only
bloodids <- grep("BLOOD", mliv2$LBPSEMA4_ID, value=T)
mliv3 <- mliv2[LBPSEMA4_ID %in% bloodids, ] #only collection dates for blood samples 

test <- table(mliv3$iid)>=2 
#test2 <- table(mliv3$iid)
#t3 <- test2[test2==2]
#dupiid <- rownames(test2[test2==2]) #subjects with 2 blood samples 

dupiid <- rownames(table(mliv3$iid)[table(mliv3$iid)==2]) #subjects with 2 blood samples 

#####################################################
#for subjects with just one blood sample, determine closest date of lab result to collection date of the blood sample 
#then, do that for subjects with two blood samples

oneiid <- rownames(table(mliv3$iid)[table(mliv3$iid)==1])

length(oneiid)
#[1] 66
length(dupiid)
#[1] 89
89 + 66
#[1] 155
length(unique(mliv3$iid))
#[1] 155

mliv3$newDate <- as.Date(as.character(mliv3$newDate))

length(intersect(livlabs$subject_id, oneiid))
#[1] 48



library(birk)

#blood1 <- data.table("iid"=oneiid, "Date"=as.Date(character()))

#here, creating a table with the blood iids that have only one sample each and populating the table with the date closest to the date of surgery; for iids that have two blood samples, the output is "NOT THERE" 
blood1 <- data.table("iid"=oneiid, "SurgeryDate1" = "SurgeryDate", "LabDate1"="DATE")

for (i in oneiid) {
	if (!i %in% livlabs$subject_id) {
	blood1[iid==i]$LabDate1 <- "NOT THERE"
	}
	else {
	x <- livlabs[subject_id==i]$newDate2
	y <- mliv3[iid==i]$newDate
	#blood1[iid==i]$Date<-as.Date(x[which.closest(x, y)])
	#blood1[iid==i]$SurgeryDate1<-mliv3[iid==i]$newDate #some issue where it converts the date to that weird numbers format; as.character stops it from doing that
	blood1[iid==i]$SurgeryDate1<-as.character(y)
	blood1[iid==i]$LabDate1<-as.character(x[which.closest(x, y)])
	}
}

blood1$SurgeryDate1<-as.Date(blood1$SurgeryDate1)

length(which(blood1$Date=="NOT THERE"))
#[1] 18
66-18
#[1] 48

# for (i in oneiid) {
# 	print(oneiid)
# 	x <- livlabs[subject_id==i]$newDate2
# 	y <- mliv3[iid==i]$newDate
# 	#blood1[iid==i]$Date<-as.Date(x[which.closest(x, y)])
# 	blood1[iid==i]$Date<-as.character(x[which.closest(x, y)])
# }

#####################################################
#for subjects with two blood samples, (i.e., two hemispheres of brain taken on two different days) determine closest date of lab to each date of collection

length(dupiid)
#[1] 89 -- should match number of rows of blood2

length(intersect(livlabs$subject_id, dupiid))
#[1] 67

# blood2 <- data.table("iid"=dupiid, "LBPSEMA4_ID" = mliv3[iid %in% dupiid]$LBPSEMA4_ID, "Date1"="DATE1", "Date2"="DATE2")

#do the same thing for the iids that have two samples each as was done for the one-sample iids except with both dates per iid
blood2 <- data.table("iid"=dupiid, "SurgeryDate1" = "SurgeryDate1", "LabDate1"="DATE1", "SurgeryDate2" = "SurgeryDate2", "LabDate2"="DATE2")

blood2 <- blood2[order(iid),]

#blood2$Dup <- rep(c(1,2), len=dim(blood2)[1])

for (i in dupiid) {
	if (!i %in% livlabs$subject_id) {
	blood2[iid==i]$LabDate1 <- "NOT THERE"
	blood2[iid==i]$LabDate2 <- "NOT THERE"
	}
	else {
	x <- livlabs[subject_id==i]$newDate2
	y <- mliv3[iid==i]$newDate
	blood2[iid==i]$SurgeryDate1<-as.character(y[1])
	blood2[iid==i]$SurgeryDate2<-as.character(y[2])
	blood2[iid==i]$LabDate1<-as.character(x[which.closest(x, y[1])])
	blood2[iid==i]$LabDate2<-as.character(x[which.closest(x, y[2])])
	}
}

blood2$SurgeryDate1<-as.Date(blood2$SurgeryDate1)
blood2$SurgeryDate2<-as.Date(blood2$SurgeryDate2)

length(which(blood2$LabDate1=="NOT THERE"))
#[1] 22

89 - 22
#[1] 67 -- checks out 

# if (!blood2[iid==i]$Date1 %in% c("DATE1", "NOT THERE") {
# 	blood2[iid==i]$Date2 <-
# }

#colnames(blood1)[2] <- "Date1"


all.blood <- merge(blood1, blood2, all=TRUE) #now all of the one-sample iids have <NA> for the SurgeryDate2 and Date2


length(oneiid) + length(dupiid)
#[1] 155

length(intersect(livlabs$subject_id, dupiid)) + length(intersect(livlabs$subject_id,oneiid))
#[1] 115
#(67 + 48 [1] 115)
length(which(all.blood$LabDate1=="NOT THERE"))
#[1] 40
155-40
#[1] 115 -- checks out

all.blood2 <- all.blood[!LabDate1=="NOT THERE", ]

dim(all.blood2)
#[1] 115   5


##########################################################
#Figure out how the closest dates match the sample IDs

bids <- all.blood2$iid #all of the blood iids that have labs

bliv <- mliv3[iid %in% bids, c("iid", "LBPSEMA4_ID", "newDate")] #subset mliv3 for the blood iids that have labs
#now do closest date to bliv$newDate and all.blood2

bliv$labDate <- "DATE"
bliv <- bliv[order(iid), ]

##########################################################################
#DON'T RUN
##########################################################################
for (i in bliv$iid) {
	if (is.na(all.blood2[iid==i]$Date2)) {
		bliv[iid==i]$labDate <- all.blood2[iid==i]$Date1
	} #stop here, the top part is for the one-sample ids
# 	else{
# 		x <- as.Date(c(all.blood2[iid==i]$Date1, all.blood2[iid==i]$Date2))
# 		bliv[iid==i]$labDate <- as.character(x[which.closest(x, bliv[iid==i]$newDate)])
# 	}
# }

#this works, but needs to be split in two so there are two dataframes with unique sample ids per iid

##########################################################################
bliv.d1 <- bliv[duplicated(bliv$iid),] #the blood iids with two samples each, this is just taking one of the samples from the two samples for each of these iids
length(unique(bliv.d1$iid))
#[1] 67

bliv.d2 <- bliv[!LBPSEMA4_ID %in% bliv.d1$LBPSEMA4_ID, ] #the bloods iids with one sample each and the second sample from the blood iids with two samples, the idea being that bliv.d1 and bliv.d2 both have only one sample for each iid

dim(bliv.d1)
#[1] 67  4
dim(bliv.d2)
#[1] 115   4
length(unique(bliv.d1$iid))
#[1] 67
dim(bliv)
#[1] 182   4
115+67
#[1] 182

bliv.d1 <- bliv.d1[order(iid), ]
bliv.d2 <- bliv.d2[order(iid), ]

for (i in bliv.d1$iid) {
	x <- as.Date(c(all.blood2[iid==i]$LabDate1, all.blood2[iid==i]$LabDate2))
	bliv.d1[iid==i]$labDate <- as.character(x[which.closest(x, bliv.d1[iid==i]$newDate)])
	}

bliv.d1$dateDif <- bliv.d1$newDate - as.Date(bliv.d1$labDate)


for (i in bliv.d2$iid) {
	if (is.na(all.blood2[iid==i]$LabDate2)) {
		bliv.d2[iid==i]$labDate <- all.blood2[iid==i]$LabDate1
	}
	else{
		x <- as.Date(c(all.blood2[iid==i]$LabDate1, all.blood2[iid==i]$LabDate2))
		bliv.d2[iid==i]$labDate <- as.character(x[which.closest(x, bliv.d2[iid==i]$newDate)])
	}
}

bliv.d2$dateDif <- bliv.d2$newDate - as.Date(bliv.d2$labDate)



labs4[labs4$subject_id=="PT-0181",]

bliv.d2[iid=="PT-0181",]


head(bliv.d1)
#        iid      LBPSEMA4_ID    newDate    labDate   dateDif
# 1: PT-0018 LBPSEMA4BLOOD548 2014-07-25 2017-02-21 -942 days
# 2: PT-0021 LBPSEMA4BLOOD556 2014-06-13 2014-06-12    1 days
# 3: PT-0022 LBPSEMA4BLOOD395 2014-06-19 2014-06-04   15 days
# 4: PT-0023 LBPSEMA4BLOOD305 2014-07-11 2014-06-25   16 days
# 5: PT-0024 LBPSEMA4BLOOD728 2014-07-09 2014-06-25   14 days
# 6: PT-0026 LBPSEMA4BLOOD773 2014-06-30 2014-06-18   12 days
head(bliv.d2)
#        iid      LBPSEMA4_ID    newDate    labDate   dateDif
# 1: PT-0018 LBPSEMA4BLOOD484 2014-06-25 2017-02-21 -972 days
# 2: PT-0019 LBPSEMA4BLOOD326 2014-04-28 2014-04-16   12 days
# 3: PT-0020 LBPSEMA4BLOOD294 2014-05-19 2014-06-04  -16 days
# 4: PT-0021 LBPSEMA4BLOOD216 2014-05-12 2014-04-18   24 days
# 5: PT-0022 LBPSEMA4BLOOD370 2014-05-19 2014-06-04  -16 days
# 6: PT-0023 LBPSEMA4BLOOD300 2014-06-11 2014-06-25  -14 days


mblood <- rbind(bliv.d1, bliv.d2)

dim(bliv.d1)
#[1] 67  5
dim(bliv.d2)
#[1] 115   5
dim(mblood)
#[1] 182   5 -- checks out, correct number of blood samples, same as in bliv

colnames(mblood)[3] <- "surgeryDate"
mblood <- mblood[order(iid),]

mblood$labDate <- as.Date(mblood$labDate)

saveRDS(mblood, file="/sc/arion/projects/psychgen/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_BLOODOnly_SurgeryDates_ClinicalLabsDates_Mapping_182BloodSamples_08FEB2022.RDS")
###############################################################################################################
#Map pt-id to sample id -- DONE
#Merge labs4 with the lymphocyte and monocyte counts -- go back to the dcasting

livlabs2 <- livlabs[livlabs$newDate2 %in% mblood$labDate, ]

dim(livlabs)
#[1] 5174    9

dim(livlabs2)
#[1] 2034    9

bloodlabs <- merge(mblood, livlabs2[, c(1:5, 8)], by.x="labDate", by.y="newDate2", all=TRUE) #can't do this because of the date dups in livlabs; need to dcast to make each result code its own column


#dcast the labs4 table so that each RESULT_CODE is a column
###############################################################################################################
library(reshape2)


#need to dcast this to make the two levels of value into separate columns 
# variable would be the column that has the values that we want to make wide, so that would be RESULT_CODE
# value.var = "VALUE"
# also need to do something with the UNIT_OF_MEASURE although i don't actually need this 


labcast <- dcast(setDT(livlabs2), subject_id + rowid(subject_id) + newDate2 + ORDER_CODE + UNIT_OF_MEASURE ~ RESULT_CODE, value.var = c("VALUE"))

labcast <- data.table(labcast)

dim(labcast)
#[1] 2034   25


saveRDS(livlabs2, file="/sc/arion/projects/psychgen/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_BLOODOnly_ClinicalLabsDates_TMP2_182BloodSamples_08FEB2022.RDS")


test2 <- livlabs2[subject_id=="PT-0010",]
test2[grep("20243-LYMPHOCYTE", test2$RESULT_CODE),]
#    MEDICAL_RECORD_NUMBER RESULT_DATE_TIME           ORDER_CODE      RESULT_CODE VALUE UNIT_OF_MEASURE subject_id
# 1:               7287778        23-FEB-15 203-CBC & PLT & DIFF 20243-LYMPHOCYTE  10.5               %    PT-0010
# 2:               7287778        18-APR-17 203-CBC & PLT & DIFF 20243-LYMPHOCYTE  16.6               %    PT-0010

#the same patient has two different results for this, but that is because they are on different dates, so i am confused abut what the problem is

test2[grep("2694-EOSINOPHIL_COUNT", test2$RESULT_CODE),]
#    MEDICAL_RECORD_NUMBER RESULT_DATE_TIME           ORDER_CODE           RESULT_CODE VALUE UNIT_OF_MEASURE subject_id
# 1:               7287778        05-MAR-14 203-CBC & PLT & DIFF 2694-EOSINOPHIL_COUNT   0.0        X10 3/Ul    PT-0010
# 2:               7287778        23-FEB-15 203-CBC & PLT & DIFF 2694-EOSINOPHIL_COUNT   0.0             X10    PT-0010
# 3:               7287778        18-APR-17 203-CBC & PLT & DIFF 2694-EOSINOPHIL_COUNT   0.1             X10    PT-0010


# data1.melt<-plyr::rename(x=data1.melt,
#                        replace=c(Var1="ProteinNumber",
#                                  Var2="SampleID",
#                                  value="Expression"))

# colnames(data4.cast)[colnames(data4.cast)=="BodyMassIndex"] <- "bmi" 
# table(duplicated(data4.cast$SampleID))
# colSums(is.na(data4.cast))
# which(is.na(data4.cast$bmi)) #155


#Internet answer for my problem:
# the ID/TIME combination does not indicate a unique row. In fact, there are two rows with each ID/TIME combinations. reshape2 assumes a single value for each possible combination of the variables and will apply a summary function to create a single variable is there are multiple entries. That is why there is the warning

# Aggregation function missing: defaulting to length

# You can get something that works if you add another variable which breaks that redundancy.

# my.df$cycle <- rep(1:2, each=num.id*num.time)
# dcast(melt(my.df, id.vars=c("cycle", "ID", "TIME")), cycle+ID~variable+TIME)

# #This works because cycle/ID/time now uniquely defines a row in my.df.

# Update: Do NOT need to do the above! just do rowid()
#rowid() function generates unique ids directly in the formula for the function dcast! so you don't need to do the random unique number generation thing that doesn't seem to work


#livlabs3$test <- rep(c(1:5), len=dim(livlabs2)[1])
#livlabs3$test <- sample(100, size = nrow(livlabs3), replace = TRUE)

#blood2$Dup <- rep(c(1,2), len=dim(blood2)[1])

#######################################
#this problem is occurring because for some, there are the same patient and same test result on the same date, but two different values or more than two
#also, the table has duplicate rows, need to remove them
#for loop, create a column that is result index, for any instance where you have these trios of subject id, new date 2 and result code, and when that trio 

dim(livlabs2)
#[1] 2034    9
dim(unique(livlabs2))
#[1] 1989    9

livlabs2[,.N,list(subject_id,newDate2,RESULT_CODE)][N>1]
livlabs2[,.N,list(subject_id,newDate2,ORDER_CODE, RESULT_CODE, UNIT_OF_MEASURE, VALUE)][N>1]

#setkeyv(livlabs2, c("subject_id", "newDate2")) #sorts the data based on the key variables, confused about what this means but ok
#The code above first sets two keys for the data.table. The key acts as an identifier and the data are automatically sorted based on the key variables. This is one of the reasons why the data.table package can be so fast at doing many of its tasks. UPDATE AS OF 2017: this solution no longer works, as data.table no longer considers unique() in keys. The option unique(, by = c(keys)) 

livlabs3 <- unique(livlabs2)
dim(livlabs3)
#[1] 1989    9

livlabs3[,rInd:=1]
trios = livlabs3[,.N,list(subject_id,newDate2,RESULT_CODE)][N>1] #N=2 for all; 21 total
trios[,id:=paste(subject_id,newDate2,RESULT_CODE)] #e.g., "PT-0089 2016-03-09 2638-PLATELET"
xt=livlabs3[paste(subject_id,newDate2,RESULT_CODE) %in% trios$id]
dim(xt)
#[1] 42 10 -- checks out, 21x2=42
xx=livlabs3[!paste(subject_id,newDate2,RESULT_CODE) %in% trios$id]

for (i in trios$id){
  cur <- xt[paste(subject_id,newDate2,RESULT_CODE)==i]
  cur[,rInd:=.I]
  livlabs4 <- rbind(xx, cur) #rbind the rest of livlabs3 that isn't in trios with the trios rows that have 1 and 2 as the ind to amke them unique while the others just have 1 as the rInd
}

range(livlabs4$rInd)
#[1] 1 2

labcast2 <- data.table(dcast(livlabs4, subject_id + newDate2 + rInd + ORDER_CODE + UNIT_OF_MEASURE ~ RESULT_CODE, value.var="VALUE"))

saveRDS(livlabs4, file="/sc/arion/projects/psychgen/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_BLOODOnly_CBC_ClinicalLabResults_09FEB2022.RDS")

saveRDS(labcast2, file="/sc/arion/projects/psychgen/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_BLOODOnly_CBC_ClinicalLabResults_LongFormat_dcast_09FEB2022.RDS")


##############################################
xt[paste(subject_id,newDate2,RESULT_CODE)==i] #takes all rows and columns of xt for which when those 3 are pasted they equal that id in trios

#     RESULT_DATE_TIME           ORDER_CODE   RESULT_CODE VALUE UNIT_OF_MEASURE subject_id    newDate   newDate2 Year rInd
# 1:        09-MAR-16 203-CBC & PLT & DIFF 2638-PLATELET   303             X10    PT-0089 2016-03-09 2016-03-09 2016    1
# 2:        09-MAR-16 203-CBC & PLT & DIFF 2638-PLATELET   279             X10    PT-0089 2016-03-09 2016-03-09 2016    1


#Orig from Alex:
# x = readRDS("/sc/arion/projects/psychgen/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_BLOODOnly_ClinicalLabsDates_TMP2_182BloodSamples_08FEB2022.RDS")
# x=unique(x)
# x[,rIndex:=1]
# trios = x[,.N,list(subject_id,newDate2,RESULT_CODE)][N>1]
# trios[,id:=paste(subject_id,newDate2,RESULT_CODE)]
# xt=x[paste(subject_id,newDate2,RESULT_CODE) %in% trios$id]
# xx=x[!paste(subject_id,newDate2,RESULT_CODE) %in% trios$id]
# for (i in trios$id){
#   cur <- xt[paste(subject_id,newDate2,RESULT_CODE)==i]
#   cur[,rIndex:=.I]
#   xx <- rbind(xx, cur)
# }
# dcast(xx, subject_id+newDate2+rIndex ~ RESULT_CODE, value.var="VALUE")

#     ID TIME X  Y
# 1   A    1  1 16
# 2   B    1  2 17
# 3   C    1  3 18
# 4   A    2  4 19
# 5   B    2  5 20
# 6   C    2  6 21
# 7   A    3  7 22
# 8   B    3  8 23
# 9   C    3  9 24
# 10  A    4 10 25
# 11  B    4 11 26
# 12  C    4 12 27
# 13  A    5 13 28
# 14  B    5 14 29
# 15  C    5 15 30

#   ID X_1 X_2 X_3 X_4 X_5 Y_1 Y_2 Y_3 Y_4 Y_5
# 1  A   1   4   7  10  13  16  19  22  25  28
# 2  B   2   5   8  11  14  17  20  23  26  29
# 3  C   3   6   9  12  15  18  21  24  27  30
##############################################







# de.compare <- merge(hde, pm.de, by="gene", suffixes = c(".harvardpmide", ".livpmde"))

# lbpBatch1.all = merge(merge(merge(merge(merge(merge(lbpBatch1.fqc,lbpBatch1.str),lbpBatch1.qc1),lbpBatch1.ftc),lbpBatch1.qc3),lbpBatch1.qc4),lbpBatch1.qc5)


setfacl -R -m u:beckmn01:rwx /sc/arion/projects/psychgen/lbp/files/lbp_cibersortx

