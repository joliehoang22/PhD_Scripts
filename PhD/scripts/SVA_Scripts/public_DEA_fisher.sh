library(biomaRt)
library(data.table)
library(ggplot2)
library(readr)
library(readxl)
library(dplyr)

######################## BP ############################################################
#################### public DEGs datasets ##############################################
public_bp<- fread("/sc/arion/projects/mscic1/results/jolie/public_DEA/public_BP_signatures.csv", data.table = FALSE)
# Filter out specific groups
public_bp <- public_bp[!(public_bp$group %in% c("Iwamoto_2004", "Wang_2013", "Nurnberger_2014")), ]
dim(public_bp) #4342 x 3, after filtering, it's 1465 x 3

colnames(public_bp)[colnames(public_bp) == "group"] <- "dataset"

ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
#ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast", host = "http://useast.ensembl.org")
result_public_bp <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                filters = "hgnc_symbol", #hugo gene symbol
                values = public_bp$symbol,
                mart = ensembl)
dim(result_public_bp) #2495 x 2 -- why is there 1000 less? @Noam, new: 3129 x 2; 1326 x 2
head(result_public_bp)

public_bp=merge(unique(result_public_bp[,1:2]),public_bp,by.x="hgnc_symbol",by.y="symbol",all.y=TRUE)

# table(public_bp$dataset)
#       Choi_2011        Chu_2009    Iwamoto_2004 Nurnberger_2014     Padmos_2008 
#             462              81              58              10              33 
#       Park_2022       Ryan_2006       Wang_2013 
#             942             101            3078 

######################## hbcc ##########################################################
bp_hbcc_no_sva <-fread("/sc/arion/projects/mscic1/results/jolie/CMC/BP/hbcc_no_sva_bp_11302024.txt",data.table=FALSE)# 19017 genes
#bp_hbcc_sva_be <-fread("/sc/arion/projects/mscic1/results/jolie/CMC/BP/hbcc_sva_be_bp.txt",data.table=FALSE)
bp_hbcc_sva_be <-fread("/sc/arion/projects/mscic1/results/jolie/CMC/BP/hbcc_sva_be_bp_12012024.txt",data.table=FALSE)

bp_hbcc_sva_leek <- bp_hbcc_no_sva ##leek is the same as no sva because leek chooses 0 SVs

table(bp_hbcc_no_sva$adj.P.Val<0.05,bp_hbcc_no_sva$logFC<0)

##EXTRACTING ONLY THE SIGNIFICANT DEGS defined as adj.P.Val<=0.25 AND RENAMING COLUMNS
mel_bp_hbcc_no_sva <-data.frame(geneID=bp_hbcc_no_sva$X, DE=as.numeric(bp_hbcc_no_sva$adj.P.Val<=0.05))
colnames(mel_bp_hbcc_no_sva)[colnames(mel_bp_hbcc_no_sva) == "geneID"]<-"my_ensembl_id"
colnames(mel_bp_hbcc_no_sva)[colnames(mel_bp_hbcc_no_sva) == "DE"]<-"my_DE"
#remove numbers after decimal place
mel_bp_hbcc_no_sva$my_ensembl_id <- sub("\\..*", "", mel_bp_hbcc_no_sva$my_ensembl_id)
head(mel_bp_hbcc_no_sva)

mel_bp_hbcc_sva_leek <- mel_bp_hbcc_no_sva

mel_bp_hbcc_sva_be <-data.frame(geneID=bp_hbcc_sva_be$X, DE=as.numeric(bp_hbcc_sva_be$adj.P.Val<=0.25))
colnames(mel_bp_hbcc_sva_be)[colnames(mel_bp_hbcc_sva_be) == "geneID"]<-"my_ensembl_id"
colnames(mel_bp_hbcc_sva_be)[colnames(mel_bp_hbcc_sva_be) == "DE"]<-"my_DE"
mel_bp_hbcc_sva_be$my_ensembl_id <- sub("\\..*", "", mel_bp_hbcc_sva_be$my_ensembl_id)
head(mel_bp_hbcc_sva_be)

res=list()
for(dataset in unique(public_bp$dataset)){
  subset=public_bp[public_bp$dataset==dataset,] #subset = new df where all rows from public_bp where the value in the dataset column matches the current dataset value being iterated over and all columns are included
  merged_df <- merge(subset, mel_bp_hbcc_no_sva, by.x = "ensembl_gene_id", by.y = "my_ensembl_id", all.y = TRUE)
  print(dim(merged_df))
  merged_df$dataset[!is.na(merged_df$dataset)] <- TRUE #dataset = TRUE while it has ensembl_gene_id
  merged_df$dataset[is.na(merged_df$dataset)] <- FALSE #dataset = FALSE if that dataset is NA for the ensembl_gene_id in which my_DE = 0 or 1 
  table(merged_df$my_DE, merged_df$dataset) 
  res[[dataset]]=fisher.test(table(merged_df$my_DE, merged_df$dataset), alternative="greater")
} #this fisher test is only done for one dataset, in this case, Nurnberger_2014 
results=data.frame(DElist=names(res), pval=sapply(res,function(x){x$p.value}), OR=sapply(res,function(x){x$estimate})) #create a df with 3 columns: DElist, pval and  OR 
results$adjustment="noSVA"
allResults=results
#park and chu are significant - which other papers are park using?
#                          DElist         pval        OR adjustment
# Park_2022             Park_2022 2.277808e-42 4.0625669      noSVA ##significant
# Choi_2011             Choi_2011 1.119936e-01 1.3061825      noSVA ##1
# Iwamoto_2004       Iwamoto_2004 8.894314e-01 0.4527529      noSVA 
# Wang_2013             Wang_2013 9.989326e-01 0.7598885      noSVA 
# Chu_2009               Chu_2009 8.056932e-03 2.6790045      noSVA ##significant
# Padmos_2008         Padmos_2008 6.789354e-01 0.9062489      noSVA ##2
# Ryan_2006             Ryan_2006 3.241962e-01 1.4324737      noSVA ##3
# Nurnberger_2014 Nurnberger_2014 1.000000e+00 0.0000000      noSVA ##dont count

#                  DElist         pval        OR adjustment
# Park_2022     Park_2022 2.277808e-42 4.0625669      noSVA ##
# Choi_2011     Choi_2011 1.119936e-01 1.3061825      noSVA ##
# Chu_2009       Chu_2009 8.056932e-03 2.6790045      noSVA ##
# Padmos_2008 Padmos_2008 6.789354e-01 0.9062489      noSVA ##
# Ryan_2006     Ryan_2006 3.241962e-01 1.4324737      noSVA ##

res=list()
for(dataset in unique(public_bp$dataset)){
  subset=public_bp[public_bp$dataset==dataset,] 
  merged_df <- merge(subset, mel_bp_hbcc_sva_be, by.x = "ensembl_gene_id", by.y = "my_ensembl_id", all.y = TRUE)
  print(dim(merged_df))
  merged_df$dataset[!is.na(merged_df$dataset)] <- TRUE #dataset = TRUE while it has ensembl_gene_id
  merged_df$dataset[is.na(merged_df$dataset)] <- FALSE #dataset = FALSE if that dataset is NA for the ensembl_gene_id in which my_DE = 0 or 1 
  table(merged_df$my_DE, merged_df$dataset) 
  res[[dataset]]=fisher.test(table(merged_df$my_DE, merged_df$dataset), alternative="greater")
} #this fisher test is only done for one dataset, in this case, Nurnberger_2014 
results=data.frame(DElist=names(res), pval=sapply(res,function(x){x$p.value}), OR=sapply(res,function(x){x$estimate})) #create a df with 3 columns: DElist, pval and  OR 
results$adjustment="SVA_be"
allResults=rbind(allResults, results)
#                           DElist         pval        OR adjustment
# Park_2022              Park_2022 2.277808e-42 4.0625669      noSVA
# Choi_2011              Choi_2011 1.119936e-01 1.3061825      noSVA
# Iwamoto_2004        Iwamoto_2004 8.894314e-01 0.4527529      noSVA
# Wang_2013              Wang_2013 9.989326e-01 0.7598885      noSVA
# Chu_2009                Chu_2009 8.056932e-03 2.6790045      noSVA
# Padmos_2008          Padmos_2008 6.789354e-01 0.9062489      noSVA
# Ryan_2006              Ryan_2006 3.241962e-01 1.4324737      noSVA
# Nurnberger_2014  Nurnberger_2014 1.000000e+00 0.0000000      noSVA ##doesn't count
# Park_20221             Park_2022 8.307918e-14 2.6285780     SVA_be
# Choi_20111             Choi_2011 6.466397e-02 1.4565598     SVA_be
# Iwamoto_20041       Iwamoto_2004 4.469056e-01 1.3557843     SVA_be
# Wang_20131             Wang_2013 1.367322e-13 1.8902050     SVA_be
# Chu_20091               Chu_2009 7.430444e-02 2.1489372     SVA_be
# Padmos_20081         Padmos_2008 5.482380e-01 1.3101945     SVA_be
# Ryan_20061             Ryan_2006 3.326589e-01 1.5131768     SVA_be
# Nurnberger_20141 Nurnberger_2014 1.000000e+00 0.0000000     SVA_be

#                   DElist         pval        OR adjustment
# Park_2022      Park_2022 2.277808e-42 4.0625669      noSVA
# Choi_2011      Choi_2011 1.119936e-01 1.3061825      noSVA
# Chu_2009        Chu_2009 8.056932e-03 2.6790045      noSVA
# Padmos_2008  Padmos_2008 6.789354e-01 0.9062489      noSVA
# Ryan_2006      Ryan_2006 3.241962e-01 1.4324737      noSVA
# Park_20221     Park_2022 8.307918e-14 2.6285780     SVA_be
# Choi_20111     Choi_2011 6.466397e-02 1.4565598     SVA_be
# Chu_20091       Chu_2009 7.430444e-02 2.1489372     SVA_be
# Padmos_20081 Padmos_2008 5.482380e-01 1.3101945     SVA_be
# Ryan_20061     Ryan_2006 3.326589e-01 1.5131768     SVA_be

##for AD: p-values without SVA are smaller than p-values with SVA
# for BP, this is true for Park and Chu 2009

#pval w sva decreases for:
#Choi_2011 (diff noSVA-SVA_be: 0.2195465)
#Wang_2013 (diff noSVA-SVA_be: 0.9999575)

#pval w sva increases for: 
#Chu_2009 (diff noSVA-SVA_be: -0.05861409) 
#Padmos_2008 (diff noSVA-SVA_be: -0.1561181)
#Ryan_2006 (diff noSVA-SVA_be: -0.049793)

res=list()
for(dataset in unique(public_bp$dataset)){
  subset=public_bp[public_bp$dataset==dataset,] 
  merged_df <- merge(subset, mel_bp_hbcc_sva_leek, by.x = "ensembl_gene_id", by.y = "my_ensembl_id", all.y = TRUE)
  print(dim(merged_df))
  merged_df$dataset[!is.na(merged_df$dataset)] <- TRUE #dataset = TRUE while it has ensembl_gene_id
  merged_df$dataset[is.na(merged_df$dataset)] <- FALSE #dataset = FALSE if that dataset is NA for the ensembl_gene_id in which my_DE = 0 or 1 
  table(merged_df$my_DE, merged_df$dataset) 
  res[[dataset]]=fisher.test(table(merged_df$my_DE, merged_df$dataset), alternative="greater")
} #this fisher test is only done for one dataset, in this case, Nurnberger_2014 
results=data.frame(DElist=names(res), pval=sapply(res,function(x){x$p.value}), OR=sapply(res,function(x){x$estimate})) #create a df with 3 columns: DElist, pval and  OR 
results$adjustment="SVA_leek"
allResults=rbind(allResults, results)
allResults_bp_hbcc=allResults


######################## BP-mssm #######################################################
bp_mssm_no_sva <-fread("/sc/arion/projects/mscic1/results/jolie/CMC/BP/cmc_no_sva_bp.txt",data.table=FALSE)
bp_mssm_sva_be <-fread("/sc/arion/projects/mscic1/results/jolie/CMC/BP/cmc_sva_be_bp.txt",data.table=FALSE)
bp_mssm_sva_leek <- bp_mssm_no_sva ##leek is the same as no sva because leek chooses 0 SVs

mel_bp_mssm_no_sva <-data.frame(geneID=bp_mssm_no_sva$X, DE=as.numeric(bp_mssm_no_sva$adj.P.Val<=0.25))
colnames(mel_bp_mssm_no_sva)[colnames(mel_bp_mssm_no_sva) == "geneID"]<-"my_ensembl_id"
colnames(mel_bp_mssm_no_sva)[colnames(mel_bp_mssm_no_sva) == "DE"]<-"my_DE"
mel_bp_mssm_no_sva$my_ensembl_id <- sub("\\..*", "", mel_bp_mssm_no_sva$my_ensembl_id)

mel_bp_mssm_sva_leek <- mel_bp_mssm_no_sva

mel_bp_mssm_sva_be <-data.frame(geneID=bp_mssm_sva_be$X, DE=as.numeric(bp_mssm_sva_be$adj.P.Val<=0.25))
colnames(mel_bp_mssm_sva_be)[colnames(mel_bp_mssm_sva_be) == "geneID"]<-"my_ensembl_id"
colnames(mel_bp_mssm_sva_be)[colnames(mel_bp_mssm_sva_be) == "DE"]<-"my_DE"
mel_bp_mssm_sva_be$my_ensembl_id <- sub("\\..*", "", mel_bp_mssm_sva_be$my_ensembl_id)
head(mel_bp_mssm_sva_be)

res=list()
for(dataset in unique(public_bp$dataset)){
  subset=public_bp[public_bp$dataset==dataset,] #subset = new df where all rows from public_bp where the value in the dataset column matches the current dataset value being iterated over and all columns are included
  merged_df <- merge(subset, mel_bp_mssm_no_sva, by.x = "ensembl_gene_id", by.y = "my_ensembl_id", all.y = TRUE)
  print(dim(merged_df))
  merged_df$dataset[!is.na(merged_df$dataset)] <- TRUE #dataset = TRUE while it has ensembl_gene_id
  merged_df$dataset[is.na(merged_df$dataset)] <- FALSE #dataset = FALSE if that dataset is NA for the ensembl_gene_id in which my_DE = 0 or 1 
  table(merged_df$my_DE, merged_df$dataset) 
  res[[dataset]]=fisher.test(table(merged_df$my_DE, merged_df$dataset), alternative="greater")
} #this fisher test is only done for one dataset, in this case, Nurnberger_2014 
results=data.frame(DElist=names(res), pval=sapply(res,function(x){x$p.value}), OR=sapply(res,function(x){x$estimate})) #create a df with 3 columns: DElist, pval and  OR 
results$adjustment="noSVA"
allResults=results

##Choi and Chu are significant!
#                          DElist       pval        OR adjustment
# Choi_2011             Choi_2011 0.02475225 2.5531681      noSVA
# Wang_2013             Wang_2013 0.97393893 0.6173575      noSVA
# Chu_2009               Chu_2009 0.01623135 5.9635318      noSVA
# Padmos_2008         Padmos_2008 1.00000000 0.0000000      noSVA
# Ryan_2006             Ryan_2006 1.00000000 0.0000000      noSVA
# Nurnberger_2014 Nurnberger_2014 1.00000000 0.0000000      noSVA

##                 DElist       pval       OR adjustment
# Park_2022     Park_2022 0.92321315 0.559167      noSVA
# Choi_2011     Choi_2011 0.02475225 2.553168      noSVA
# Chu_2009       Chu_2009 0.01623135 5.963532      noSVA
# Padmos_2008 Padmos_2008 1.00000000 0.000000      noSVA
# Ryan_2006     Ryan_2006 1.00000000 0.000000      noSVA

res=list()
for(dataset in unique(public_bp$dataset)){
  subset=public_bp[public_bp$dataset==dataset,] 
  merged_df <- merge(subset, mel_bp_mssm_sva_be, by.x = "ensembl_gene_id", by.y = "my_ensembl_id", all.y = TRUE)
  print(dim(merged_df))
  merged_df$dataset[!is.na(merged_df$dataset)] <- TRUE #dataset = TRUE while it has ensembl_gene_id
  merged_df$dataset[is.na(merged_df$dataset)] <- FALSE #dataset = FALSE if that dataset is NA for the ensembl_gene_id in which my_DE = 0 or 1 
  table(merged_df$my_DE, merged_df$dataset) 
  res[[dataset]]=fisher.test(table(merged_df$my_DE, merged_df$dataset), alternative="greater")
} #this fisher test is only done for one dataset, in this case, Nurnberger_2014 
results=data.frame(DElist=names(res), pval=sapply(res,function(x){x$p.value}), OR=sapply(res,function(x){x$estimate})) #create a df with 3 columns: DElist, pval and  OR 
results$adjustment="SVA_be"
allResults=rbind(allResults, results)
#                           DElist        pval        OR adjustment
# Choi_2011              Choi_2011 0.024752254 2.5531681      noSVA
# Wang_2013              Wang_2013 0.973938932 0.6173575      noSVA
# Chu_2009                Chu_2009 0.016231348 5.9635318      noSVA
# Padmos_2008          Padmos_2008 1.000000000 0.0000000      noSVA
# Ryan_2006              Ryan_2006 1.000000000 0.0000000      noSVA
# Nurnberger_2014  Nurnberger_2014 1.000000000 0.0000000      noSVA
# Choi_20111             Choi_2011 0.020801695 2.6521268     SVA_be
# Wang_20131             Wang_2013 0.192171458 1.2322426     SVA_be
# Chu_20091               Chu_2009 0.001813727 8.3465031     SVA_be
# Padmos_20081         Padmos_2008 0.126493285 7.8975960     SVA_be
# Ryan_20061             Ryan_2006 1.000000000 0.0000000     SVA_be
# Nurnberger_20141 Nurnberger_2014 1.000000000 0.0000000     SVA_be

#                   DElist        pval       OR adjustment
# Park_2022      Park_2022 0.923213146 0.559167      noSVA
# Choi_2011      Choi_2011 0.024752254 2.553168      noSVA
# Chu_2009        Chu_2009 0.016231348 5.963532      noSVA
# Padmos_2008  Padmos_2008 1.000000000 0.000000      noSVA
# Ryan_2006      Ryan_2006 1.000000000 0.000000      noSVA
# Park_20221     Park_2022 0.142537334 1.512584     SVA_be
# Choi_20111     Choi_2011 0.020801695 2.652127     SVA_be ##more sig w SVA?!
# Chu_20091       Chu_2009 0.001813727 8.346503     SVA_be #more sig w SVA?!
# Padmos_20081 Padmos_2008 0.126493285 7.897596     SVA_be
# Ryan_20061     Ryan_2006 1.000000000 0.000000     SVA_be

res=list()
for(dataset in unique(public_bp$dataset)){
  subset=public_bp[public_bp$dataset==dataset,] 
  merged_df <- merge(subset, mel_bp_mssm_sva_leek, by.x = "ensembl_gene_id", by.y = "my_ensembl_id", all.y = TRUE)
  print(dim(merged_df))
  merged_df$dataset[!is.na(merged_df$dataset)] <- TRUE #dataset = TRUE while it has ensembl_gene_id
  merged_df$dataset[is.na(merged_df$dataset)] <- FALSE #dataset = FALSE if that dataset is NA for the ensembl_gene_id in which my_DE = 0 or 1 
  table(merged_df$my_DE, merged_df$dataset) 
  res[[dataset]]=fisher.test(table(merged_df$my_DE, merged_df$dataset), alternative="greater")
} #this fisher test is only done for one dataset, in this case, Nurnberger_2014 
results=data.frame(DElist=names(res), pval=sapply(res,function(x){x$p.value}), OR=sapply(res,function(x){x$estimate})) #create a df with 3 columns: DElist, pval and  OR 
results$adjustment="SVA_leek"
allResults_bp_mssm=rbind(allResults, results)

#allResults_bp_mssm and allResults_bp_hbcc
allResults_bp_mssm$data="BP_MSSM"
allResults_bp_hbcc$data="BP_HBCC"

allResults=rbind(allResults_bp_mssm,allResults_bp_hbcc)
allResults$SVA="none"
allResults$SVA[allResults$adjustment=="SVA_be"]="be"
allResults$SVA[allResults$adjustment=="SVA_leek"]="leek"
#ggplot(allResults,aes(x= data, y=-log10(pval), fill= SVA)) + geom_boxplot()
# ggplot(allResults,aes(x= data, y=OR, fill= adjustment)) + geom_boxplot(notch=TRUE)

# wilcox.test(allResults$pval[allResults$adjustment=="noSVA" & allResults$data=="BP_MSSM"],allResults$pval[allResults$adjustment=="SVA_be" & allResults$data=="BP_MSSM"],alternative="less") #p-value = 0.8012 NOT SIGNIFICANT

# wilcox.test(allResults$pval[allResults$adjustment=="noSVA" & allResults$data=="BP_HBCC"],allResults$pval[allResults$adjustment=="SVA_be" & allResults$data=="BP_HBCC"],alternative="less") #p-value = 0.5 NOT SIGNIFICANT

# library(tidyverse)
# library(ggpubr)
# library(rstatix)

# allResults$log10pval=-log10(allResults$pval)
# allResults$log10pval[!is.finite(allResults$log10pval)]=-log10(.Machine$double.xmin)

# stat.test1 <- allResults[allResults$data=="BP_MSSM",] %>% wilcox_test(formula = log10pval ~ SVA, alternative="less")
# stat.test1$data="BP_MSSM"
# attr(stat.test1, "args")$data=allResults
# stat.test1 = stat.test1 %>% add_xy_position(x = "data", dodge = 0.8)

# stat.test2 <- allResults[allResults$data=="BP_HBCC",] %>% wilcox_test(formula = log10pval ~ SVA, alternative="less")
# stat.test2$data="BP_HBCC"
# attr(stat.test2, "args")$data=allResults
# stat.test2 = stat.test2 %>% add_xy_position(x = "data", dodge = 0.8)

# stat.test2$xmin = stat.test2$xmin 
# stat.test2$xmax = stat.test2$xmax 

# stat.test=rbind(stat.test1,stat.test2)
# stat.test$data=NULL
# stat.test$p.adj=p.adjust(stat.test$p,method="fdr")
# stat.test$p.adj.signif[stat.test$p.adj<0.05]="*"
# stat.test$p.adj.signif[stat.test$p.adj<0.005]="**"
# stat.test$p.adj.signif[stat.test$p.adj<0.0005]="***"
# # stat.test = stat.test %>% add_xy_position(x = "data", dodge = 0.8)

# #p1<-ggboxplot(allResults, x = "data", y = "log10pval", fill="SVA") + stat_pvalue_manual(stat.test, step.increase = 0.0,label = "p.adj") + theme_bw()

# #ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/BP_DEA_fisher_wilcoxon.pdf", plot = p1, width = 8, height = 5)

# Filter out "leek"
allResults <- allResults %>% filter(SVA != "leek")

# Calculate -log10(p-value), handle infinite values
allResults$log10pval <- -log10(allResults$pval)
allResults$log10pval[!is.finite(allResults$log10pval)] <- -log10(.Machine$double.xmin)

# BP_MSSM pairwise test
stat.test1 <- allResults[allResults$data == "BP_MSSM",] %>%
  pairwise_wilcox_test(log10pval ~ SVA, p.adjust.method = "fdr", alternative = "less")
stat.test1$data <- "BP_MSSM"
attr(stat.test1, "args")$data <- allResults
stat.test1 <- stat.test1 %>% add_xy_position(x = "data", dodge = 0.8)

# BP_HBCC pairwise test
stat.test2 <- allResults[allResults$data == "BP_HBCC",] %>%
  pairwise_wilcox_test(log10pval ~ SVA, p.adjust.method = "fdr", alternative = "less")
stat.test2$data <- "BP_HBCC"
attr(stat.test2, "args")$data <- allResults
stat.test2 <- stat.test2 %>% add_xy_position(x = "data", dodge = 0.8)

# Combine stats
stat.test <- rbind(stat.test1, stat.test2)
#stat.test$data <- NULL
stat.test$p.adj.signif[stat.test$p.adj < 0.05] <- "*"
stat.test$p.adj.signif[stat.test$p.adj < 0.005] <- "**"
stat.test$p.adj.signif[stat.test$p.adj < 0.0005] <- "***"

allResults <- allResults %>%
  filter(SVA %in% c("be", "none")) %>%
  mutate(
    SVA = factor(SVA, levels = c("none", "be")),  
    facet_label = case_when(
      data == "BP_MSSM" ~ "BP: CMC",
      data == "BP_HBCC" ~ "BP: CMC-HBCC"
    )
  )
allResults$pval <- ifelse(allResults$pval == 0, 1e-300, allResults$pval)

# Strip whitespace and make sure it's character
# allResults$data <- trimws(as.character(allResults$data))

# # Now filter again
# bp_cmc <- allResults[allResults$data == "BP_MSSM", ]

# Split datasets
bp_cmc <- subset(allResults, grepl("BP: CMC", facet_label))
bp_hbcc <- subset(allResults, grepl("BP: CMC-HBCC", facet_label))

bp_cmc$facet_label <- factor("BP: CMC", levels = c("BP: CMC", "BP: CMC-HBCC"))
bp_hbcc$facet_label <- factor("BP: CMC-HBCC", levels = c("BP: CMC", "BP: CMC-HBCC"))

library(scales)

# Y-axis breaks
breaks_bp <- c(1, 1e-10, 1e-20, 1e-30,1e-40,1e-50)

p_bp_cmc <- ggplot(bp_cmc, aes(x = SVA, y = pval)) +
  geom_boxplot(outlier.shape = NA, fill = NA, color = "black", size = 1, width = 0.4) +
  geom_line(aes(group = DElist), color = "gray40", alpha = 0.25) +
  geom_point(aes(group = DElist), color = "black", shape = 16, size = 4, alpha = 0.5) +
  scale_y_neglog10(
    limits = c(1, 1e-10), 
    breaks = breaks_bp,
    labels = scales::trans_format("log10", math_format(10^.x))
  ) +
  facet_wrap(~facet_label, scales = "fixed") +
  scale_x_discrete(labels = c("none" = "No SVA", "be" = "SVA 'BE'")) +
  theme_bw(base_size = 16) +
  theme(
    plot.margin = margin(5, 5, 5, 5),
    legend.position = "none",
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text = element_text(size = 16),
    strip.text = element_text(face = "bold", size = 18)
  ) +
  labs(y = "Fisher's Exact Test p-value", x = "SVA Methods")

p_bp_hbcc <- ggplot(bp_hbcc, aes(x = SVA, y = pval)) +
  geom_boxplot(outlier.shape = NA, fill = NA, color = "black", size = 1, width = 0.4) +
  geom_line(aes(group = DElist), color = "gray40", alpha = 0.25) +
  geom_point(aes(group = DElist), color = "black", shape = 16, size = 4, alpha = 0.5) +
  scale_y_neglog10(
    limits = c(1, 1e-45),
    breaks = breaks_bp,
    labels = scales::trans_format("log10", math_format(10^.x))
  ) +
  facet_wrap(~facet_label, scales = "fixed") +  # <- changed from "free_y" to "fixed"
  scale_x_discrete(labels = c("none" = "No SVA", "be" = "SVA 'BE'")) +
  theme_bw(base_size = 16) +
  theme(
    plot.margin = margin(5, 5, 5, 5),
    legend.position = "none",
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text = element_text(size = 16),
    strip.text = element_text(face = "bold", size = 18)
  ) +
  labs(y = "Fisher's Exact Test p-value", x = "SVA Methods")

combo <- p_bp_cmc + p_bp_hbcc
ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/bp_public_DEA.pdf", plot = combo, width = 9, height = 6)

# # Plot 2: BP: CMC-HBCC
# p_bp_hbcc <- ggplot(bp_hbcc, aes(x = SVA, y = pval)) +
#   geom_boxplot(outlier.shape = NA, fill = NA, color = "black", size = 1, width = 0.4) +
#   geom_line(aes(group = DElist), color = "gray40", alpha = 0.25, position = position_nudge(x = 0)) +
#   geom_point(aes(group = DElist), color = "black", shape = 16, size = 4, alpha = 0.5, position = position_nudge(x = 0)) +
#   scale_y_neglog10(breaks = breaks_bp, labels = scales::label_log()) +
#   facet_wrap(~facet_label, scales = "free_y") +
#   scale_x_discrete(labels = c("none" = "No SVA", "be" = "SVA 'BE'")) +
#   theme_bw(base_size = 16) +
#   theme(
#     plot.margin = margin(5, 5, 5, 5),
#     legend.position = "none",
#     axis.title.x = element_text(size = 18),
#     axis.title.y = element_text(size = 18),
#     axis.text = element_text(size = 16),
#     strip.text = element_text(face = "bold", size = 18)
#   ) +
#   labs(y = "Fisher's Exact Test p-value", x = "SVA Methods")

# combo <- p_bp_cmc + p_bp_hbcc


######################## SZ ############################################################
#################### public DEGs datasets ##############################################
## this dataset contains 3 datasets
public_sz_done<- fread("/sc/arion/projects/mscic1/results/jolie/public_DEA/public_SZ_signatures.csv", data.table = FALSE)
dim(public_sz_done) #1176 x 4
colnames(public_sz_done)[colnames(public_sz_done) == "group"] <- "dataset"
colnames(public_sz_done)[colnames(public_sz_done) == "ensemblID"] <- "ensembl_gene_id"
colnames(public_sz_done)[colnames(public_sz_done) == "symbol"] <- "hgnc_symbol"
# #reorder columns:
public_sz_done <- public_sz_done[, c("hgnc_symbol", "ensembl_gene_id", "DE", "dataset")]

##need to load in Hwang and Wu separately and convert gene symbols to ensembl IDs
public_sz<- fread("/sc/arion/projects/mscic1/results/jolie/public_DEA/public_SZ_signatures2.csv", data.table = FALSE)
dim(public_sz) #866 x 3
colnames(public_sz)[colnames(public_sz) == "group"] <- "dataset"

ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
#ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast", host = "http://useast.ensembl.org")
result_public_sz <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                          filters = "hgnc_symbol", #hugo gene symbol
                          values = public_sz$symbol,
                          mart = ensembl)
##in case emsembl doesn't work, load this in for result_public_sz
# write.table(
#   result_public_sz,
#   file = "/sc/arion/projects/mscic1/results/jolie/public_DEA/result_public_sz_backup.txt",
#   sep = "\t",
#   quote = FALSE,
#   row.names = FALSE
# )

dim(result_public_sz) #853 x 2 -- 866-853 = 13 less aka 13 gene symbol that doesnt have an ensembl id?
head(result_public_sz)

# colnames(result_public_sz)[colnames(result_public_sz) == "ensembl_gene_id"] <- "public_ensembl_id"
# result_public_sz$public_DE <- 1
# head(result_public_sz) 

# #check for duplicates: 
# any(duplicated(result_public_sz$public_ensembl_id))

##need to check this
public_sz=merge(unique(result_public_sz[,1:2]),public_sz,by.x="hgnc_symbol",by.y="symbol",all.y=TRUE)

public_sz<-rbind(public_sz,public_sz_done)
dim(public_sz) #2129 x 4

any(duplicated(public_sz$public_ensembl_id)) #FALSE

# public_sz2 <- public_sz[!duplicated(public_sz$public_ensembl_id), ]
# table(public_sz2$DE) #1889 unique ensebl IDs that are not duplicated

######################## sz-hbcc #######################################################
sz_hbcc_no_sva <-fread("/sc/arion/projects/mscic1/results/jolie/CMC/hbcc_no_sva.txt",data.table=FALSE)
sz_hbcc_sva_be <-fread("/sc/arion/projects/mscic1/results/jolie/CMC/hbcc_sva_be.txt",data.table=FALSE)
sz_hbcc_sva_leek <- sz_hbcc_no_sva ##leek is the same as no sva because leek chooses 0 SVs

##EXTRACTING ONLY THE SIGNIFICANT DEGS defined as adj.P.Val<=0.25 AND RENAMING COLUMNS
mel_sz_hbcc_no_sva <-data.frame(my_ensembl_id=sz_hbcc_no_sva$X, my_DE=as.numeric(sz_hbcc_no_sva$adj.P.Val<=0.25))
mel_sz_hbcc_no_sva$my_ensembl_id <- sub("\\..*", "", mel_sz_hbcc_no_sva$my_ensembl_id)
head(mel_sz_hbcc_no_sva)
table(mel_sz_hbcc_no_sva$my_DE) ##5652 DEGs out of 19086 genes

mel_sz_hbcc_sva_leek <- mel_sz_hbcc_no_sva

mel_sz_hbcc_sva_be <-data.frame(my_ensembl_id=sz_hbcc_sva_be$X, my_DE=as.numeric(sz_hbcc_sva_be$adj.P.Val<=0.25))
mel_sz_hbcc_sva_be$my_ensembl_id <- sub("\\..*", "", mel_sz_hbcc_sva_be$my_ensembl_id)
head(mel_sz_hbcc_sva_be)

table(public_sz$dataset)
#   Collado_Torres_2019            Hwang_2013 Manchia_2017_GSE12649 
#                   293                   159                   400 
# Manchia_2017_GSE25673               Wu_2012 
#                   483                   794 

res=list()
for(dataset in unique(public_sz$dataset)){
  subset=public_sz[public_sz$dataset==dataset,] 
  merged_df <- merge(subset, mel_sz_hbcc_no_sva, by.x = "ensembl_gene_id", by.y = "my_ensembl_id", all.y = TRUE)
  print(dim(merged_df))
  merged_df$dataset[!is.na(merged_df$dataset)] <- TRUE #dataset = TRUE while it has ensembl_gene_id
  merged_df$dataset[is.na(merged_df$dataset)] <- FALSE #dataset = FALSE if that dataset is NA for the ensembl_gene_id in which my_DE = 0 or 1 
  table(merged_df$my_DE, merged_df$dataset) 
  res[[dataset]]=fisher.test(table(merged_df$my_DE, merged_df$dataset), alternative="greater")
} #this fisher test is only done for one dataset, in this case, Nurnberger_2014 
results=data.frame(DElist=names(res), pval=sapply(res,function(x){x$p.value}), OR=sapply(res,function(x){x$estimate})) #create a df with 3 columns: DElist, pval and  OR 
results$adjustment="noSVA"
allResults=results

#very nice! 
#                                      DElist         pval       OR adjustment
# Wu_2012                             Wu_2012 4.822464e-05 1.432085      noSVA
# Hwang_2013                       Hwang_2013 1.976770e-03 1.841156      noSVA
# Collado_Torres_2019     Collado_Torres_2019 2.369661e-25 3.572490      noSVA
# Manchia_2017_GSE25673 Manchia_2017_GSE25673 1.716710e-06 1.594557      noSVA
# Manchia_2017_GSE12649 Manchia_2017_GSE12649 1.564360e-24 2.890049      noSVA

res=list()
for(dataset in unique(public_sz$dataset)){
  subset=public_sz[public_sz$dataset==dataset,] 
  merged_df <- merge(subset, mel_sz_hbcc_sva_be, by.x = "ensembl_gene_id", by.y = "my_ensembl_id", all.y = TRUE)
  print(dim(merged_df))
  merged_df$dataset[!is.na(merged_df$dataset)] <- TRUE #dataset = TRUE while it has ensembl_gene_id
  merged_df$dataset[is.na(merged_df$dataset)] <- FALSE #dataset = FALSE if that dataset is NA for the ensembl_gene_id in which my_DE = 0 or 1 
  table(merged_df$my_DE, merged_df$dataset) 
  res[[dataset]]=fisher.test(table(merged_df$my_DE, merged_df$dataset), alternative="greater")
} #this fisher test is only done for one dataset, in this case, Nurnberger_2014 
results=data.frame(DElist=names(res), pval=sapply(res,function(x){x$p.value}), OR=sapply(res,function(x){x$estimate})) #create a df with 3 columns: DElist, pval and  OR 
results$adjustment="SVA_be"
allResults=rbind(allResults, results)
#                                       DElist         pval        OR adjustment
# Wu_2012                              Wu_2012 4.822464e-05 1.4320855      noSVA
# Hwang_2013                        Hwang_2013 1.976770e-03 1.8411562      noSVA
# Collado_Torres_2019      Collado_Torres_2019 2.369661e-25 3.5724903      noSVA
# Manchia_2017_GSE25673  Manchia_2017_GSE25673 1.716710e-06 1.5945573      noSVA
# Manchia_2017_GSE12649  Manchia_2017_GSE12649 1.564360e-24 2.8900490      noSVA
# Wu_20121                             Wu_2012 2.478920e-02 1.6946425     SVA_be
# Hwang_20131                       Hwang_2013 6.325063e-01 0.9362483     SVA_be
# Collado_Torres_20191     Collado_Torres_2019 8.699338e-12 5.4991056     SVA_be
# Manchia_2017_GSE256731 Manchia_2017_GSE25673 2.977475e-02 1.7424763     SVA_be
# Manchia_2017_GSE126491 Manchia_2017_GSE12649 3.026183e-01 1.2355308     SVA_be

#compared to pvalues without SVA, pvalues WITHOUT SVA for 5 out of 5 datasets are smaller

res=list()
for(dataset in unique(public_sz$dataset)){
  subset=public_sz[public_sz$dataset==dataset,] 
  merged_df <- merge(subset, mel_sz_hbcc_sva_leek, by.x = "ensembl_gene_id", by.y = "my_ensembl_id", all.y = TRUE)
  print(dim(merged_df))
  merged_df$dataset[!is.na(merged_df$dataset)] <- TRUE #dataset = TRUE while it has ensembl_gene_id
  merged_df$dataset[is.na(merged_df$dataset)] <- FALSE #dataset = FALSE if that dataset is NA for the ensembl_gene_id in which my_DE = 0 or 1 
  table(merged_df$my_DE, merged_df$dataset) 
  res[[dataset]]=fisher.test(table(merged_df$my_DE, merged_df$dataset), alternative="greater")
} #this fisher test is only done for one dataset, in this case, Nurnberger_2014 
results=data.frame(DElist=names(res), pval=sapply(res,function(x){x$p.value}), OR=sapply(res,function(x){x$estimate})) #create a df with 3 columns: DElist, pval and  OR 
results$adjustment="SVA_leek"
allResults_sz_hbcc=rbind(allResults, results)

######################## sz-mssm #######################################################
sz_mssm_no_sva <-fread("/sc/arion/projects/mscic1/results/jolie/CMC/mssm_no_sva.txt",data.table=FALSE)
sz_mssm_sva_be <-fread("/sc/arion/projects/mscic1/results/jolie/CMC/mssm_sva_be.txt",data.table=FALSE)
sz_mssm_sva_leek <- sz_mssm_no_sva ##leek = no sva because leek chooses 0 SVs

mel_sz_mssm_no_sva <-data.frame(my_ensembl_id=sz_mssm_no_sva$X, my_DE=as.numeric(sz_mssm_no_sva$adj.P.Val<=0.25))
mel_sz_mssm_no_sva$my_ensembl_id <- sub("\\..*", "", mel_sz_mssm_no_sva$my_ensembl_id)
head(mel_sz_mssm_no_sva)
table(mel_sz_mssm_no_sva$my_DE) ##1,371 DEGs out of 19,086 genes

mel_sz_mssm_sva_leek <- mel_sz_mssm_no_sva

mel_sz_mssm_sva_be <-data.frame(my_ensembl_id=sz_mssm_sva_be$X, my_DE=as.numeric(sz_mssm_sva_be$adj.P.Val<=0.25))
mel_sz_mssm_sva_be$my_ensembl_id <- sub("\\..*", "", mel_sz_mssm_sva_be$my_ensembl_id)
head(mel_sz_mssm_sva_be)
table(mel_sz_mssm_sva_be$my_DE)  #2,053 DEGs out of 19086 DEGs

res=list()
for(dataset in unique(public_sz$dataset)){
  subset=public_sz[public_sz$dataset==dataset,] 
  merged_df <- merge(subset, mel_sz_mssm_no_sva, by.x = "ensembl_gene_id", by.y = "my_ensembl_id", all.y = TRUE)
  print(dim(merged_df))
  merged_df$dataset[!is.na(merged_df$dataset)] <- TRUE #dataset = TRUE while it has ensembl_gene_id
  merged_df$dataset[is.na(merged_df$dataset)] <- FALSE #dataset = FALSE if that dataset is NA for the ensembl_gene_id in which my_DE = 0 or 1 
  table(merged_df$my_DE, merged_df$dataset) 
  res[[dataset]]=fisher.test(table(merged_df$my_DE, merged_df$dataset), alternative="greater")
} #this fisher test is only done for one dataset, in this case, Nurnberger_2014 
results=data.frame(DElist=names(res), pval=sapply(res,function(x){x$p.value}), OR=sapply(res,function(x){x$estimate})) #create a df with 3 columns: DElist, pval and  OR 
results$adjustment="noSVA"
allResults=results

#                                      DElist         pval       OR adjustment
# Wu_2012                             Wu_2012 7.516362e-04 1.612160      noSVA ##
# Hwang_2013                       Hwang_2013 5.564854e-03 2.267370      noSVA ##
# Collado_Torres_2019     Collado_Torres_2019 2.560875e-14 3.610138      noSVA ##
# Manchia_2017_GSE25673 Manchia_2017_GSE25673 1.914321e-05 1.916940      noSVA
# Manchia_2017_GSE12649 Manchia_2017_GSE12649 1.738204e-01 1.209551      noSVA #not significant

res=list()
for(dataset in unique(public_sz$dataset)){
  subset=public_sz[public_sz$dataset==dataset,] 
  merged_df <- merge(subset, mel_sz_mssm_sva_be, by.x = "ensembl_gene_id", by.y = "my_ensembl_id", all.y = TRUE)
  print(dim(merged_df))
  merged_df$dataset[!is.na(merged_df$dataset)] <- TRUE #dataset = TRUE while it has ensembl_gene_id
  merged_df$dataset[is.na(merged_df$dataset)] <- FALSE #dataset = FALSE if that dataset is NA for the ensembl_gene_id in which my_DE = 0 or 1 
  table(merged_df$my_DE, merged_df$dataset) 
  res[[dataset]]=fisher.test(table(merged_df$my_DE, merged_df$dataset), alternative="greater")
} #this fisher test is only done for one dataset, in this case, Nurnberger_2014 
results=data.frame(DElist=names(res), pval=sapply(res,function(x){x$p.value}), OR=sapply(res,function(x){x$estimate})) #create a df with 3 columns: DElist, pval and  OR 
results$adjustment="SVA_be"
allResults=rbind(allResults, results)
##p-values are consistently higher with SVA (4 out of 5) than without SVA
#                                       DElist         pval       OR adjustment
# Wu_2012                              Wu_2012 7.516362e-04 1.612160      noSVA
# Hwang_2013                        Hwang_2013 5.564854e-03 2.267370      noSVA
# Collado_Torres_2019      Collado_Torres_2019 2.560875e-14 3.610138      noSVA
# Manchia_2017_GSE25673  Manchia_2017_GSE25673 1.914321e-05 1.916940      noSVA
# Manchia_2017_GSE12649  Manchia_2017_GSE12649 1.738204e-01 1.209551      noSVA
# Wu_20121                             Wu_2012 8.605996e-04 1.497812     SVA_be
# Hwang_20131                       Hwang_2013 1.230172e-01 1.450339     SVA_be
# Collado_Torres_20191     Collado_Torres_2019 4.428277e-05 1.941961     SVA_be
# Manchia_2017_GSE256731 Manchia_2017_GSE25673 5.435990e-07 1.918612     SVA_be
# Manchia_2017_GSE126491 Manchia_2017_GSE12649 1.312700e-02 1.414946     SVA_be #becomes significant

#for SZ, pvalues are larger WITH SVA for 2 out of 5 datasets compared to pvalues without SVA

res=list()
for(dataset in unique(public_sz$dataset)){
  subset=public_sz[public_sz$dataset==dataset,] 
  merged_df <- merge(subset, mel_sz_mssm_sva_leek, by.x = "ensembl_gene_id", by.y = "my_ensembl_id", all.y = TRUE)
  print(dim(merged_df))
  merged_df$dataset[!is.na(merged_df$dataset)] <- TRUE #dataset = TRUE while it has ensembl_gene_id
  merged_df$dataset[is.na(merged_df$dataset)] <- FALSE #dataset = FALSE if that dataset is NA for the ensembl_gene_id in which my_DE = 0 or 1 
  table(merged_df$my_DE, merged_df$dataset) 
  res[[dataset]]=fisher.test(table(merged_df$my_DE, merged_df$dataset), alternative="greater")
} #this fisher test is only done for one dataset, in this case, Nurnberger_2014 
results=data.frame(DElist=names(res), pval=sapply(res,function(x){x$p.value}), OR=sapply(res,function(x){x$estimate})) #create a df with 3 columns: DElist, pval and  OR 
results$adjustment="SVA_leek"
allResults_sz_mssm=rbind(allResults, results)

#allResults_bp_mssm and allResults_bp_hbcc
allResults_sz_mssm$data="SZ_MSSM"
allResults_sz_hbcc$data="SZ_HBCC"

allResults=rbind(allResults_sz_mssm,allResults_sz_hbcc)
allResults$SVA="none"
allResults$SVA[allResults$adjustment=="SVA_be"]="be"
allResults$SVA[allResults$adjustment=="SVA_leek"]="leek"
#ggplot(allResults,aes(x= data, y=-log10(pval), fill= SVA)) + geom_boxplot()
# ggplot(allResults,aes(x= data, y=OR, fill= adjustment)) + geom_boxplot(notch=TRUE)

# Filter out "leek"
allResults <- allResults %>% filter(SVA != "leek")

# Calculate -log10(p-value), handle infinite values
allResults$log10pval <- -log10(allResults$pval)
allResults$log10pval[!is.finite(allResults$log10pval)] <- -log10(.Machine$double.xmin)

# BP_MSSM pairwise test
stat.test1 <- allResults[allResults$data == "SZ_MSSM",] %>%
  pairwise_wilcox_test(log10pval ~ SVA, p.adjust.method = "fdr", alternative = "less")
stat.test1$data <- "SZ_MSSM"
attr(stat.test1, "args")$data <- allResults
stat.test1 <- stat.test1 %>% add_xy_position(x = "data", dodge = 0.8)

# BP_HBCC pairwise test
stat.test2 <- allResults[allResults$data == "SZ_HBCC",] %>%
  pairwise_wilcox_test(log10pval ~ SVA, p.adjust.method = "fdr", alternative = "less")
stat.test2$data <- "SZ_HBCC"
attr(stat.test2, "args")$data <- allResults
stat.test2 <- stat.test2 %>% add_xy_position(x = "data", dodge = 0.8)

# Combine stats
stat.test <- rbind(stat.test1, stat.test2)
#stat.test$data <- NULL
stat.test$p.adj.signif[stat.test$p.adj < 0.05] <- "*"
stat.test$p.adj.signif[stat.test$p.adj < 0.005] <- "**"
stat.test$p.adj.signif[stat.test$p.adj < 0.0005] <- "***"

allResults <- allResults %>%
  filter(SVA %in% c("be", "none")) %>%
  mutate(
    SVA = factor(SVA, levels = c("none", "be")),  
    facet_label = case_when(
      data == "SZ_MSSM" ~ "SZ: CMC",
      data == "SZ_HBCC" ~ "SZ: CMC-HBCC"
    )
  )
allResults$pval <- ifelse(allResults$pval == 0, 1e-300, allResults$pval)

# Split datasets
sz_cmc <- subset(allResults, grepl("SZ: CMC", facet_label))
sz_hbcc <- subset(allResults, grepl("SZ: CMC-HBCC", facet_label))

sz_cmc$facet_label <- factor("SZ: CMC", levels = c("SZ: CMC", "SZ: CMC-HBCC"))
sz_hbcc$facet_label <- factor("SZ: CMC-HBCC", levels = c("SZ: CMC", "SZ: CMC-HBCC"))

library(scales)

sz_cmc <- filter(sz_cmc, data == "SZ_MSSM")

# Y-axis breaks
breaks_sz <- c(1, 1e-10, 1e-20, 1e-30)

p_sz_cmc <- ggplot(sz_cmc, aes(x = SVA, y = pval)) +
  geom_boxplot(outlier.shape = NA, fill = NA, color = "black", size = 1, width = 0.4) +
  geom_line(aes(group = DElist), color = "gray40", alpha = 0.25) +
  geom_point(aes(group = DElist), color = "black", shape = 16, size = 4, alpha = 0.5) +
  scale_y_neglog10(
    limits = c(1, 1e-30), 
    breaks = breaks_sz,
    labels = scales::trans_format("log10", math_format(10^.x))
  ) +
  facet_wrap(~facet_label, scales = "fixed") +
  scale_x_discrete(labels = c("none" = "No SVA", "be" = "SVA 'BE'")) +
  theme_bw(base_size = 16) +
  theme(
    plot.margin = margin(5, 5, 5, 5),
    legend.position = "none",
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text = element_text(size = 16),
    strip.text = element_text(face = "bold", size = 18)
  ) +
  labs(y = "Fisher's Exact Test p-value", x = "SVA Methods")

p_sz_hbcc <- ggplot(sz_hbcc, aes(x = SVA, y = pval)) +
  geom_boxplot(outlier.shape = NA, fill = NA, color = "black", size = 1, width = 0.4) +
  geom_line(aes(group = DElist), color = "gray40", alpha = 0.25) +
  geom_point(aes(group = DElist), color = "black", shape = 16, size = 4, alpha = 0.5) +
  scale_y_neglog10(
    limits = c(1, 1e-30),
    breaks = breaks_sz,
    labels = scales::trans_format("log10", math_format(10^.x))
  ) +
  facet_wrap(~facet_label, scales = "fixed") +  # <- changed from "free_y" to "fixed"
  scale_x_discrete(labels = c("none" = "No SVA", "be" = "SVA 'BE'")) +
  theme_bw(base_size = 16) +
  theme(
    plot.margin = margin(5, 5, 5, 5),
    legend.position = "none",
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text = element_text(size = 16),
    strip.text = element_text(face = "bold", size = 18)
  ) +
  labs(y = "Fisher's Exact Test p-value", x = "SVA Methods")

combo <- p_sz_cmc + p_sz_hbcc
ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/sz_public_DEA.pdf", plot = combo, width = 9, height = 6)

combo_sz_and_bp <- (p_bp_cmc + p_bp_hbcc) / (p_sz_cmc + p_sz_hbcc)

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/sz_and_bp_public_DEA.pdf", plot = combo_sz_and_bp, width = 9, height = 9)

##############################################################################################
################################ ARCHIVED ####################################################
##############################################################################################

wilcox.test(allResults$pval[allResults$adjustment=="noSVA" & allResults$data=="SZ_MSSM"],allResults$pval[allResults$adjustment=="SVA_be" & allResults$data=="SZ_MSSM"],alternative="less") #p-value = 0.4206 NOT SIGNIFICANT
wilcox.test(allResults$pval[allResults$adjustment=="noSVA" & allResults$data=="SZ_HBCC"],allResults$pval[allResults$adjustment=="SVA_be" & allResults$data=="SZ_HBCC"],alternative="less") #p-value = 0.02778 - SIGNIFICANT

library(tidyverse)
library(ggpubr)
library(rstatix)

allResults$log10pval=-log10(allResults$pval)
allResults$log10pval[!is.finite(allResults$log10pval)]=-log10(.Machine$double.xmin)

stat.test1 <- allResults[allResults$data=="SZ_MSSM",] %>% wilcox_test(formula = log10pval ~ SVA, alternative="less")
stat.test1$data="SZ_MSSM"
attr(stat.test1, "args")$data=allResults
stat.test1 = stat.test1 %>% add_xy_position(x = "data", dodge = 0.8)

stat.test2 <- allResults[allResults$data=="SZ_HBCC",] %>% wilcox_test(formula = log10pval ~ SVA, alternative="less")
stat.test2$data="SZ_HBCC"
attr(stat.test2, "args")$data=allResults
stat.test2 = stat.test2 %>% add_xy_position(x = "data", dodge = 0.8)

stat.test2$xmin = stat.test2$xmin 
stat.test2$xmax = stat.test2$xmax 

stat.test=rbind(stat.test1,stat.test2)
stat.test$data=NULL
stat.test$p.adj=p.adjust(stat.test$p,method="fdr")
stat.test$p.adj.signif[stat.test$p.adj<0.05]="*"
stat.test$p.adj.signif[stat.test$p.adj<0.005]="**"
stat.test$p.adj.signif[stat.test$p.adj<0.0005]="***"
# stat.test = stat.test %>% add_xy_position(x = "data", dodge = 0.8)

p1<-ggboxplot(allResults, x = "data", y = "log10pval", fill="SVA") + stat_pvalue_manual(stat.test, step.increase = 0.0,label = "p.adj") + theme_bw()

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/SZ_DEA_fisher_wilcoxon.pdf", plot = p1, width = 8, height = 5)

############################## ARCHIVE ##########################################
############################## ARCHIVE ##########################################
############################## ARCHIVE ##########################################
############################## ARCHIVE ##########################################
############################## ARCHIVE ##########################################
############################## ARCHIVE ##########################################

######################## PD ############################################################
#################### public DEGs datasets ##############################################
public_pd<- fread("/sc/arion/projects/mscic1/results/jolie/public_DEA/public_PD_signatures.csv", data.table = FALSE)
dim(public_pd) #2127 x 3
colnames(public_pd)[colnames(public_pd) == "group"] <- "dataset"

ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
#ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast", host = "http://useast.ensembl.org")
result_public_pd <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                          filters = "hgnc_symbol", #hugo gene symbol
                          values = public_pd$symbol,
                          mart = ensembl)
dim(result_public_pd) #1681 x 2 -- 446 less? @Noam
head(result_public_pd)

##need to check this
public_pd=merge(unique(result_public_pd[,1:2]),public_pd,by.x="hgnc_symbol",by.y="symbol",all.y=TRUE)

colnames(result_public_pd)[colnames(result_public_pd) == "ensembl_gene_id"] <- "public_ensembl_id"
result_public_pd$public_DE <- 1
head(result_public_pd) 

#check for duplicates: 
any(duplicated(result_public_pd$public_ensembl_id)) #FALSE

##picking two PD datasets, ones with the biggest DEGs list - GSE8397 and GSE49036 (highest numbers of DEGs
GSE8397_none <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE8397_DEA_results_without_sva_final.txt", header = TRUE, sep = "\t")
mel_GSE8397_none <-data.frame(my_ensembl_id=GSE8397_none$X, my_DE=as.numeric(GSE8397_none$adj.P.Val<=0.25))
table(mel_GSE8397_none$my_DE)
#     0     1 
# 18826  3457
GSE8397_sva_be <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE8397_DEA_results_with_sva_be_final.txt", header = TRUE, sep = "\t")

GSE49036_none <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE49036_DEA_results_without_sva_final.txt", header = TRUE, sep = "\t")
mel_GSE49036_none <-data.frame(my_ensembl_id=GSE49036_none$X, my_DE=as.numeric(GSE49036_none$adj.P.Val<=0.25))
table(mel_GSE49036_none$my_DE)
#     0     1 
# 47458  7217 
#need to do leek for this one
GSE49036_sva_be <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE49036_DEA_results_with_sva_be_final.txt", header = TRUE, sep = "\t")


GSE7621_none <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE7621_DEA_results_without_sva_final.txt", header = TRUE, sep = "\t")
mel_GSE7621_none <-data.frame(my_ensembl_id=GSE7621_none$X, my_DE=as.numeric(GSE7621_none$adj.P.Val<=0.25))
table(mel_GSE7621_none$my_DE)
#     0     1 
# 53822   496 

GSE24378_none <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE24378_DEA_results_without_sva_final.txt", header = TRUE, sep = "\t")
mel_GSE24378_none <-data.frame(my_ensembl_id=GSE24378_none$X, my_DE=as.numeric(GSE24378_none$adj.P.Val<=0.25))
table(mel_GSE24378_none$my_DE)
#    0     1 
# 61354     5

GSE20292_none <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE20292_DEA_results_without_sva_final.txt", header = TRUE, sep = "\t")
mel_GSE20292_none <-data.frame(my_ensembl_id=GSE20292_none$X, my_DE=as.numeric(GSE20292_none$adj.P.Val<=0.25))
table(mel_GSE20292_none$my_DE)
#     0     1 
# 21564   719 

GSE20141_none <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE20141_DEA_results_without_sva_final.txt", header = TRUE, sep = "\t")
mel_GSE20141_none <-data.frame(my_ensembl_id=GSE20141_none$X, my_DE=as.numeric(GSE20141_none$adj.P.Val<=0.25))
table(mel_GSE20141_none$my_DE)
#     0     1 
# 54674     1 

GSE20163_none <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE20163_DEA_results_without_sva_final.txt", header = TRUE, sep = "\t")
mel_GSE20163_none <-data.frame(my_ensembl_id=GSE20163_none$X, my_DE=as.numeric(GSE20163_none$adj.P.Val<=0.25))
table(mel_GSE20163_none$my_DE)
#     0     1 
# 21860   423

GSE20164_none <- read.table("/sc/arion/projects/mscic1/results/jolie/GEO/GSE20164_DEA_results_without_sva_final.txt", header = TRUE, sep = "\t")
mel_GSE20164_none <-data.frame(my_ensembl_id=GSE20164_none$X, my_DE=as.numeric(GSE20164_none$adj.P.Val<=0.25))
table(mel_GSE20164_none$my_DE)
#     0 
# 22283


######################################################################################
mel_GSE8397_none
mel_GSE49036_none
#
##need to match probe 
head(mel_GSE49036_none)
  my_ensembl_id my_DE
1         22106     1
2         35468     1
3          2225     1
4         12929     1
5         54169     1
6         32992     1

res=list()
for(dataset in unique(public_pd$dataset)){
  subset=public_pd[public_pd$dataset==dataset,] 
  merged_df <- merge(subset, mel_GSE8397_none, by.x = "ensembl_gene_id", by.y = "my_ensembl_id", all.y = TRUE)
  print(dim(merged_df))
  merged_df$dataset[!is.na(merged_df$dataset)] <- TRUE #dataset = TRUE while it has ensembl_gene_id
  merged_df$dataset[is.na(merged_df$dataset)] <- FALSE #dataset = FALSE if that dataset is NA for the ensembl_gene_id in which my_DE = 0 or 1 
  table(merged_df$my_DE, merged_df$dataset) 
  res[[dataset]]=fisher.test(table(merged_df$my_DE, merged_df$dataset), alternative="greater")
} #this fisher test is only done for one dataset, in this case, Nurnberger_2014 
results=data.frame(DElist=names(res), pval=sapply(res,function(x){x$p.value}), OR=sapply(res,function(x){x$estimate})) #create a df with 3 columns: DElist, pval and  OR 
results$adjustment="noSVA"
allResults=results

##############################ORIGINAL ATTEMPT##########################################
public_SZ<- fread("/sc/arion/projects/mscic1/results/jolie/public_DEA/Clifton_2022_DEA_SZ.csv", data.table = FALSE)
dim(public_SZ) #16278 x 3

public_BP<- fread("/sc/arion/projects/mscic1/results/jolie/public_DEA/Ryan_2006_DEA_BP.csv", data.table = FALSE)
dim(public_BP) #79 x 2

#data <- read_excel("/sc/arion/projects/mscic1/results/jolie/public_DEA/Simunovic_2008_DEA_PD.xls")
#public_PD<- read.csv("/sc/arion/projects/mscic1/results/jolie/public_DEA/Simunovic_2008_DEA_PD.csv")


######SCHIZOPHRENIA#######################################################################
##########################################################################################
public_SZ<- fread("/sc/arion/projects/mscic1/results/jolie/public_DEA/Clifton_2022_DEA_SZ.csv", data.table = FALSE)
dim(public_SZ) #32003 x 4

##already have ensembl ID!!! YAY!!

colnames(public_SZ)[colnames(public_SZ) == "Ensembl_ID"] <- "public_ensembl_id"
public_SZ$public_DE <- 1

table(public_SZ$public_DE) #16278 1s ####list all DEGs of public datasets as 1

sz_hbcc_no_sva <-fread("/sc/arion/projects/mscic1/results/jolie/CMC/hbcc_no_sva.txt",data.table=FALSE)
sz_mssm_no_sva <-fread("/sc/arion/projects/mscic1/results/jolie/CMC/mssm_no_sva.txt",data.table=FALSE)

mel_sz_hbcc_no_sva <-data.frame(geneID=sz_hbcc_no_sva$X, DE=as.numeric(sz_hbcc_no_sva$adj.P.Val<=0.25))
mel_sz_hbcc_no_sva$geneID <- gsub("\\.[0-9]+$", "", mel_sz_hbcc_no_sva$geneID)
colnames(mel_sz_hbcc_no_sva)[colnames(mel_sz_hbcc_no_sva) == "geneID"]<-"my_ensembl_id"
colnames(mel_sz_hbcc_no_sva)[colnames(mel_sz_hbcc_no_sva) == "DE"]<-"my_DE"
#table(mel_msbb_sva_leek$my_DE)

#############################without SVA#############################
merged_df <- merge(public_SZ, mel_sz_hbcc_no_sva, by.x = "public_ensembl_id", by.y = "my_ensembl_id", all.y = TRUE) #by.x and by.y specifies the columns in result_ad and mel_msbb_no_sva to merge by.
dim(merged_df) #19086 x 5 with all.y=TRUE
dim(merged_df) # 20508 x 5 will all=TRUE

#fill missing values for hgnc_symbol, gwas_DE, and my_DE
merged_df$Gene_symbol[is.na(merged_df$Gene_symbol)] <- NA
merged_df$public_DE[is.na(merged_df$public_DE)] <- 0
merged_df$my_DE[is.na(merged_df$my_DE)] <- 0

table(merged_df$my_DE, merged_df$public_DE) 
  #   FALSE  TRUE
  # 0  3346 11510
  # 1   884  4768
merged_df$public_DE <- ifelse(merged_df$public_DE == 1, TRUE, FALSE)

fisher.test(table(merged_df$my_DE, merged_df$public_DE), alternative="greater")
p_value1 <- fisher.test(table(merged_df$my_DE, merged_df$public_DE), alternative="greater")$p.value

###result of all=TRUE in merge
# 	Fisher's Exact Test for Count Data
# data:  table(merged_df$my_DE, merged_df$public_DE)
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#  1.463418      Inf
# sample estimates:
# odds ratio 
#   1.567925 

###result of all.y=TRUE in merge
# 	Fisher's Exact Test for Count Data
# data:  table(merged_df$my_DE, merged_df$public_DE)
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#  1.669259      Inf
# sample estimates:
# odds ratio 
#   1.788923 

#############################with SVA#############################
sz_hbcc_sva_be <-fread("/sc/arion/projects/mscic1/results/jolie/CMC/hbcc_sva_be.txt",data.table=FALSE) 
sz_mssm_sva_be <-fread("/sc/arion/projects/mscic1/results/jolie/CMC/mssm_sva_be.txt",data.table=FALSE)
#did not do leek, need to do leek 

mel_sz_hbcc_sva_be <-data.frame(geneID=sz_hbcc_sva_be$X, DE=as.numeric(sz_hbcc_sva_be$adj.P.Val<=0.25))
mel_sz_hbcc_sva_be$geneID <- gsub("\\.[0-9]+$", "", mel_sz_hbcc_sva_be$geneID)
colnames(mel_sz_hbcc_sva_be)[colnames(mel_sz_hbcc_sva_be) == "geneID"]<-"my_ensembl_id"
colnames(mel_sz_hbcc_sva_be)[colnames(mel_sz_hbcc_sva_be) == "DE"]<-"my_DE"

merged_df <- merge(public_SZ, mel_sz_hbcc_sva_be, by.x = "public_ensembl_id", by.y = "my_ensembl_id", all.y = TRUE) #by.x and by.y specifies the columns in result_ad and mel_msbb_no_sva to merge by.
dim(merged_df) #19086 x 5 with all.y=TRUE
dim(merged_df) # 20508 x 5 will all=TRUE

#fill missing values for hgnc_symbol, gwas_DE, and my_DE
merged_df$Gene_symbol[is.na(merged_df$Gene_symbol)] <- NA
merged_df$public_DE[is.na(merged_df$public_DE)] <- 0
merged_df$my_DE[is.na(merged_df$my_DE)] <- 0

table(merged_df$my_DE, merged_df$public_DE) 
#   FALSE  TRUE
# 0  4159 14524
# 1    71   332
merged_df$public_DE <- ifelse(merged_df$public_DE == 1, TRUE, FALSE)

fisher.test(table(merged_df$my_DE, merged_df$public_DE), alternative="greater")
p_value2 <- fisher.test(table(merged_df$my_DE, merged_df$public_DE), alternative="greater")$p.value
# 	Fisher's Exact Test for Count Data
# data:  table(merged_df$my_DE, merged_df$public_DE)
# p-value = 0.01376
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#  1.072424      Inf
# sample estimates:
# odds ratio 
#   1.338975 

# Perform a Wilcoxon test on the two p-values
wilcox_test_result <- wilcox.test(c(p_value1, p_value2),alternative="less"); wilcox_test_result

# 	Wilcoxon signed rank exact test
# data:  c(p_value1, p_value2)
# V = 3, p-value = 1
# alternative hypothesis: true location is less than 0



######BIPOLAR DISORDER#######################################################################
##########################################################################################
public_BP<- fread("/sc/arion/projects/mscic1/results/jolie/public_DEA/Ryan_2006_DEA_BP.csv", data.table = FALSE)
dim(public_BP) #79 x 2


bp_hbcc_no_sva <-fread("/sc/arion/projects/mscic1/results/jolie/CMC/BP/hbcc_no_sva_bp.txt",data.table=FALSE)
bp_mssm_no_sva <-fread("/sc/arion/projects/mscic1/results/jolie/CMC/BP/cmc_no_sva_bp.txt",data.table=FALSE)

bp_hbcc_sva_be <-fread("/sc/arion/projects/mscic1/results/jolie/CMC/BP/hbcc_sva_be_bp.txt",data.table=FALSE) 
bp_mssm_sva_be <-fread("/sc/arion/projects/mscic1/results/jolie/CMC/BP/cmc_sva_be_bp.txt",data.table=FALSE)
#did not do leek, need to do leek 























