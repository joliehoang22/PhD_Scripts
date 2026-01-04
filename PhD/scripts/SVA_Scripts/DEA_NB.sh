##NOAM'S PUBLIC DEA##
###########AD
#Zhang B, Gaiteri C, Bodea LG, et al. Integrated systems approach identifies genetic
#219 nodes and networks in late-onset Alzheimer's disease. Cell. 2013;153(3):707-720.

#Mostafavi S, Gaiteri C, Sullivan SE, et al. A molecular network of the aging human brain provides insights into the pathology and cognitive decline of Alzheimer's disease. Nat Neurosci. 2018;21(6):811-819.

library(data.table)
library(ggplot2)
library(biomaRt)

public_AD <- fread("/sc/arion/projects/mscic1/results/jolie/public_DEA_comparison/public_AD_signatures.tsv", data.table = FALSE)
dim(public_AD) #32003 x 4
length(unique(public_AD$Symbol)) #11261 unique genes
public_AD$dataset=unlist(lapply(strsplit(public_AD$shortName,".",fixed=TRUE),function(x){x[1]}))
table(public_AD$dataset) #20 datasets? Zhang_Atrophy_CB, Zhang_Atrophy_PFC,    Zhang_Braak_CB

colnames(public_AD)[colnames(public_AD) == "logFC.or.Cor"] <- "logFC"

## new OCT 29 2025 - REMOVE MYERS
public_AD <- public_AD[!(public_AD$dataset %in% "Myers_CB"), ]
public_AD <- public_AD[!(public_AD$dataset %in% "Myers_TC"), ]

table(public_AD$dataset)
#find ensembl IDs for these gene symbols
# ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# #ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast", host = "http://useast.ensembl.org")
# result_public_ad <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
#                 filters = "hgnc_symbol", #hugo gene symbol
#                 values = public_AD$Symbol,
#                 mart = ensembl)

# dim(result_public_ad) #26802 x 2
# head(result_public_ad)

# write.table(
#   result_public_ad,
#   file = "/sc/arion/projects/mscic1/results/jolie/public_DEA_comparison/ensembl_to_hgnc_ad.txt",
#   sep = "\t",
#   quote = FALSE,
#   row.names = FALSE
# )

result_public_ad <- read.table(
  "/sc/arion/projects/mscic1/results/jolie/public_DEA_comparison/ensembl_to_hgnc_ad.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE)

public_AD=merge(unique(result_public_ad[,1:2]),public_AD,by.x="hgnc_symbol",by.y="Symbol",all.y=TRUE)

colnames(result_public_ad)[colnames(result_public_ad) == "ensembl_gene_id"] <- "public_ensembl_id"
result_public_ad$public_DE <- 1
head(result_public_ad) 

any(duplicated(result_public_ad$public_ensembl_id)) #TRUE

##there are duplicates in public_ensembl_id and hgnc_symbol - KEEP? i feel like we should remove duplicates for public_ensembl_id @Noam

table(result_public_ad$public_DE) #26,802 1s WITH duplicates 

result_public_ad <- result_public_ad[!duplicated(result_public_ad$public_ensembl_id), ]
table(result_public_ad$public_DE)# 11,197 1s WITHOUT duplicates!
#          Allen_CB          Allen_TC      Avramopoulos           Blalock 
#              1275                58              1026              2791 
#         Colangelo             Liang        Miller_CA1        Miller_CA3 
#                31              1060               406              1232 
#      Mostafavi_Ab  Mostafavi_ClinAD  Mostafavi_CogDec     Mostafavi_NFT 
#              4286              2971              4752              2485 
# Mostafavi_PathoAD          Myers_CB          Myers_TC             Satoh 
#              2056                31               711               520 
#  Zhang_Atrophy_CB Zhang_Atrophy_PFC    Zhang_Braak_CB   Zhang_Braak_PFC 
#               156              2363               979              2814 
              
######TEST ASSOCIATION BETWEEN PUBLIC DEA AND MY DE WITHOUT SVA IN MSBB FOR AD#############
############################################################################################
msbb_no_sva <-fread("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/MSBBv30/msbb_no_sva_1.11.txt", data.table = FALSE) #26618 x 7
msbb_be <-fread("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/MSBBv30/msbb_sva_be_1.11.txt", data.table = FALSE) #26618 x 7
msbb_leek <-fread("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/MSBBv30/msbb_sva_leek_1.11.txt", data.table = FALSE) #26618 x 7

table(msbb_no_sva$adj.P.Val<0.05,msbb_no_sva$logFC<0)

mel_msbb_no_sva <-data.frame(geneID=msbb_no_sva$X, DE=as.numeric(msbb_no_sva$adj.P.Val<=0.25))
colnames(mel_msbb_no_sva)[colnames(mel_msbb_no_sva) == "geneID"]<-"my_ensembl_id"
colnames(mel_msbb_no_sva)[colnames(mel_msbb_no_sva) == "DE"]<-"my_DE"

mel_msbb_sva_be <-data.frame(geneID=msbb_be$X, DE=as.numeric(msbb_be$adj.P.Val<=0.25))
colnames(mel_msbb_sva_be)[colnames(mel_msbb_sva_be) == "geneID"]<-"my_ensembl_id"
colnames(mel_msbb_sva_be)[colnames(mel_msbb_sva_be) == "DE"]<-"my_DE"

mel_msbb_sva_leek <-data.frame(geneID=msbb_leek$X, DE=as.numeric(msbb_leek$adj.P.Val<=0.25))
colnames(mel_msbb_sva_leek)[colnames(mel_msbb_sva_leek) == "geneID"]<-"my_ensembl_id"
colnames(mel_msbb_sva_leek)[colnames(mel_msbb_sva_leek) == "DE"]<-"my_DE"

table(mel_msbb_sva_be$my_DE) ##these numbers seem different ~
#     0     1 
# 26375   243 
# dataset="Mostafavi2018NN.B.amyloid_Negative_Cor.in.DLPFC"

res=list()
for(dataset in unique(public_AD$dataset)){
	subset=public_AD[public_AD$dataset==dataset,]
	merged_df <- merge(subset, mel_msbb_no_sva, by.x = "ensembl_gene_id", by.y = "my_ensembl_id", all.y = TRUE) #by.x and by.y specifies the columns in result_ad and mel_msbb_no_sva to merge by.
	print(dim(merged_df))
	merged_df$dataset[!is.na(merged_df$dataset)] <- TRUE
	merged_df$dataset[is.na(merged_df$dataset)] <- FALSE
	table(merged_df$my_DE, merged_df$dataset) 
	res[[dataset]]=fisher.test(table(merged_df$my_DE, merged_df$dataset), alternative="greater")
}
results=data.frame(DElist=names(res), pval=sapply(res,function(x){x$p.value}), OR=sapply(res,function(x){x$estimate}))
results$adjustment="noSVA"
allResults=results

##all pval in results are significant!!!

res=list()
for(dataset in unique(public_AD$dataset)){
	subset=public_AD[public_AD$dataset==dataset,]
	merged_df <- merge(subset, mel_msbb_sva_be, by.x = "ensembl_gene_id", by.y = "my_ensembl_id", all.y = TRUE) #by.x and by.y specifies the columns in result_ad and mel_msbb_sva_be to merge by.
	print(dim(merged_df))
	merged_df$dataset[!is.na(merged_df$dataset)] <- TRUE
	merged_df$dataset[is.na(merged_df$dataset)] <- FALSE
	table(merged_df$my_DE, merged_df$dataset) 
	res[[dataset]]=fisher.test(table(merged_df$my_DE, merged_df$dataset), alternative="greater")
}
results=data.frame(DElist=names(res), pval=sapply(res,function(x){x$p.value}), OR=sapply(res,function(x){x$estimate}))
results$adjustment="SVA_be"
allResults=rbind(allResults, results)
#                               DElist         pval        OR adjustment
# Mostafavi_NFT          Mostafavi_NFT 1.099575e-13 1.5470137      noSVA
# Zhang_Braak_CB        Zhang_Braak_CB 7.298167e-15 1.9157063      noSVA
# Miller_CA1                Miller_CA1 5.048644e-04 1.5480588      noSVA
# Mostafavi_CogDec    Mostafavi_CogDec 1.342860e-52 1.9401354      noSVA
# Mostafavi_PathoAD  Mostafavi_PathoAD 1.428855e-40 2.2292591      noSVA
# Mostafavi_Ab            Mostafavi_Ab 8.787707e-69 2.1806448      noSVA
# Miller_CA3                Miller_CA3 1.289184e-30 2.2876930      noSVA
# Mostafavi_ClinAD    Mostafavi_ClinAD 2.164611e-60 2.3088796      noSVA
# Zhang_Braak_PFC      Zhang_Braak_PFC 1.994296e-51 2.1572493      noSVA
# Zhang_Atrophy_PFC  Zhang_Atrophy_PFC 1.869062e-60 2.4265080      noSVA
# Satoh                          Satoh 7.102492e-03 1.3452368      noSVA
# Allen_CB                    Allen_CB 2.263803e-22 2.0250859      noSVA
# Liang                          Liang 4.993905e-03 1.2514368      noSVA
# Zhang_Atrophy_CB    Zhang_Atrophy_CB 1.063438e-01 1.3505682      noSVA
# Blalock                      Blalock 3.905167e-15 1.5158969      noSVA
# Avramopoulos            Avramopoulos 3.488707e-22 2.1385607      noSVA
# Allen_TC                    Allen_TC 3.118412e-01 1.2518860      noSVA
# Colangelo                  Colangelo 4.324716e-01 1.1993608      noSVA
# Mostafavi_NFT1         Mostafavi_NFT 7.425048e-05 2.1345493     SVA_be
# Zhang_Braak_CB1       Zhang_Braak_CB 1.972149e-01 1.3672424     SVA_be
# Miller_CA11               Miller_CA1 4.754935e-01 1.1302439     SVA_be
# Mostafavi_CogDec1   Mostafavi_CogDec 1.212477e-03 1.6430398     SVA_be
# Mostafavi_PathoAD1 Mostafavi_PathoAD 3.045146e-03 1.8516684     SVA_be
# Mostafavi_Ab1           Mostafavi_Ab 7.874436e-03 1.5209085     SVA_be
# Miller_CA31               Miller_CA3 1.184256e-01 1.4348379     SVA_be
# Mostafavi_ClinAD1   Mostafavi_ClinAD 1.724078e-02 1.5338684     SVA_be
# Zhang_Braak_PFC1     Zhang_Braak_PFC 1.169197e-04 1.9793040     SVA_be
# Zhang_Atrophy_PFC1 Zhang_Atrophy_PFC 8.821927e-06 2.2807458     SVA_be
# Satoh1                         Satoh 1.799914e-02 2.2483246     SVA_be
# Allen_CB1                   Allen_CB 3.286988e-03 2.0304870     SVA_be
# Liang1                         Liang 8.314909e-02 1.5491960     SVA_be
# Zhang_Atrophy_CB1   Zhang_Atrophy_CB 7.162401e-01 0.7972836     SVA_be
# Blalock1                     Blalock 1.761889e-04 1.9380010     SVA_be
# Avramopoulos1           Avramopoulos 9.498036e-01 0.5421661     SVA_be
# Allen_TC1                   Allen_TC 4.019702e-01 1.9773890     SVA_be
# Colangelo1                 Colangelo 1.000000e+00 0.0000000     SVA_be

res=list()
for(dataset in unique(public_AD$dataset)){
	subset=public_AD[public_AD$dataset==dataset,]
	merged_df <- merge(subset, mel_msbb_sva_leek, by.x = "ensembl_gene_id", by.y = "my_ensembl_id", all.y = TRUE) #by.x and by.y specifies the columns in result_ad and mel_msbb_sva_leek to merge by.
	print(dim(merged_df))
	merged_df$dataset[!is.na(merged_df$dataset)] <- TRUE
	merged_df$dataset[is.na(merged_df$dataset)] <- FALSE
	table(merged_df$my_DE, merged_df$dataset) 
	res[[dataset]]=fisher.test(table(merged_df$my_DE, merged_df$dataset), alternative="greater")
}
results=data.frame(DElist=names(res), pval=sapply(res,function(x){x$p.value}), OR=sapply(res,function(x){x$estimate}))
results$adjustment="SVA_leek"
allResults=rbind(allResults, results)

# ggplot(allResults,aes(x=OR, color=adjustment, fill=adjustment)) + geom_density(alpha=0.4)
# ggplot(allResults,aes(x=-log10(pval), color=adjustment, fill=adjustment)) + geom_density(alpha=0.4)

######TEST ASSOCIATION BETWEEN PUBLIC DEA AND MY DE WITHOUT SVA IN ROSMAP FOR AD#############
############################################################################################
rosmap_no_sva <-fread("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/rosmap_no_sva_1.11.txt", data.table = FALSE) #26618 x 7
rosmap_be <-fread("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/rosmap_sva_be_1.16.txt", data.table = FALSE) #26618 x 7
rosmap_leek <-fread("/sc/arion/projects/mscic1/results/jolie/ROSMAP-MSBB/Rosmapv30/rosmap_sva_leek_1.16.txt", data.table = FALSE) #26618 x 7

table(rosmap_no_sva$adj.P.Val<0.05,rosmap_no_sva$logFC<0)

mel_rosmap_no_sva <-data.frame(geneID=rosmap_no_sva$X, DE=as.numeric(rosmap_no_sva$adj.P.Val<=0.25))
colnames(mel_rosmap_no_sva)[colnames(mel_rosmap_no_sva) == "geneID"]<-"my_ensembl_id"
colnames(mel_rosmap_no_sva)[colnames(mel_rosmap_no_sva) == "DE"]<-"my_DE"
mel_rosmap_sva_be <-data.frame(geneID=rosmap_be$X, DE=as.numeric(rosmap_be$adj.P.Val<=0.25))
colnames(mel_rosmap_sva_be)[colnames(mel_rosmap_sva_be) == "geneID"]<-"my_ensembl_id"
colnames(mel_rosmap_sva_be)[colnames(mel_rosmap_sva_be) == "DE"]<-"my_DE"
mel_rosmap_sva_leek <-data.frame(geneID=rosmap_leek$X, DE=as.numeric(rosmap_leek$adj.P.Val<=0.25))
colnames(mel_rosmap_sva_leek)[colnames(mel_rosmap_sva_leek) == "geneID"]<-"my_ensembl_id"
colnames(mel_rosmap_sva_leek)[colnames(mel_rosmap_sva_leek) == "DE"]<-"my_DE"

table(mel_rosmap_sva_leek$my_DE)
#     0     1 
# 14337  5139 

# dataset="Mostafavi2018NN.B.amyloid_Negative_Cor.in.DLPFC"

res=list()
for(dataset in unique(public_AD$dataset)){
	subset=public_AD[public_AD$dataset==dataset,]
	merged_df <- merge(subset, mel_rosmap_no_sva, by.x = "ensembl_gene_id", by.y = "my_ensembl_id", all.y = TRUE) #by.x and by.y specifies the columns in result_ad and mel_rosmap_no_sva to merge by.
	print(dim(merged_df))
	merged_df$dataset[!is.na(merged_df$dataset)] <- TRUE
	merged_df$dataset[is.na(merged_df$dataset)] <- FALSE
	table(merged_df$my_DE, merged_df$dataset) 
	res[[dataset]]=fisher.test(table(merged_df$my_DE, merged_df$dataset), alternative="greater")
}
results=data.frame(DElist=names(res), pval=sapply(res,function(x){x$p.value}), OR=sapply(res,function(x){x$estimate}))
results$adjustment="noSVA"
allResults_rosmap=results

res=list()
for(dataset in unique(public_AD$dataset)){
	subset=public_AD[public_AD$dataset==dataset,]
	merged_df <- merge(subset, mel_rosmap_sva_be, by.x = "ensembl_gene_id", by.y = "my_ensembl_id", all.y = TRUE) #by.x and by.y specifies the columns in result_ad and mel_rosmap_sva_be to merge by.
	print(dim(merged_df))
	merged_df$dataset[!is.na(merged_df$dataset)] <- TRUE
	merged_df$dataset[is.na(merged_df$dataset)] <- FALSE
	table(merged_df$my_DE, merged_df$dataset) 
	res[[dataset]]=fisher.test(table(merged_df$my_DE, merged_df$dataset), alternative="greater")
}
results=data.frame(DElist=names(res), pval=sapply(res,function(x){x$p.value}), OR=sapply(res,function(x){x$estimate}))
results$adjustment="SVA_be"
allResults_rosmap=rbind(allResults_rosmap, results)

res=list()
for(dataset in unique(public_AD$dataset)){
	subset=public_AD[public_AD$dataset==dataset,]
	merged_df <- merge(subset, mel_rosmap_sva_leek, by.x = "ensembl_gene_id", by.y = "my_ensembl_id", all.y = TRUE) #by.x and by.y specifies the columns in result_ad and mel_rosmap_sva_leek to merge by.
	print(dim(merged_df))
	merged_df$dataset[!is.na(merged_df$dataset)] <- TRUE
	merged_df$dataset[is.na(merged_df$dataset)] <- FALSE
	table(merged_df$my_DE, merged_df$dataset) 
	res[[dataset]]=fisher.test(table(merged_df$my_DE, merged_df$dataset), alternative="greater")
}
results=data.frame(DElist=names(res), pval=sapply(res,function(x){x$p.value}), OR=sapply(res,function(x){x$estimate}))
results$adjustment="SVA_leek"
allResults_rosmap=rbind(allResults_rosmap, results)

allResults_rosmap$data="ROSMAP"
allResults$data="MSBB"

allResults=rbind(allResults,allResults_rosmap)
allResults$SVA="none"
allResults$SVA[allResults$adjustment=="SVA_be"]="be"
allResults$SVA[allResults$adjustment=="SVA_leek"]="leek"

# ggplot(allResults_rosmap,aes(x=OR, color=adjustment, fill=adjustment)) + geom_density(alpha=0.4)
# ggplot(allResults_rosmap,aes(x=-log10(pval), color=adjustment, fill=adjustment)) + geom_density(alpha=0.4)
# ggplot(allResults_rosmap,aes(x= adjustment, y=-log10(pval))) + geom_boxplot()

####################### MOVE TO Figure6_20250622.sh or Figure6.sh ################# 
###################################################################################
#ggplot(allResults,aes(x= data, y=-log10(pval), fill= SVA)) + geom_boxplot()
# ggplot(allResults,aes(x= data, y=OR, fill= adjustment)) + geom_boxplot(notch=TRUE)
wilcox.test(allResults$pval[allResults$adjustment=="noSVA" & allResults$data=="MSBB"],allResults$pval[allResults$adjustment=="SVA_be" & allResults$data=="MSBB"],alternative="less")
wilcox.test(allResults$pval[allResults$adjustment=="noSVA" & allResults$data=="ROSMAP"],allResults$pval[allResults$adjustment=="SVA_be" & allResults$data=="ROSMAP"],alternative="less")

library(tidyverse)
library(ggpubr)
library(rstatix)

allResults$log10pval=-log10(allResults$pval)
allResults$log10pval[!is.finite(allResults$log10pval)]=-log10(.Machine$double.xmin)

stat.test1 <- allResults[allResults$data=="MSBB",] %>% wilcox_test(formula = log10pval ~ SVA, alternative="less")
stat.test1$data="MSBB"
attr(stat.test1, "args")$data=allResults
stat.test1 = stat.test1 %>% add_xy_position(x = "data", dodge = 0.8)
stat.test2 <- allResults[allResults$data=="ROSMAP",] %>% wilcox_test(formula = log10pval ~ SVA, alternative="less")
stat.test2$data="ROSMAP"
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

ggboxplot(allResults, x = "data", y = "log10pval", fill="SVA") + stat_pvalue_manual(stat.test, step.increase = 0.0,label = "p.adj") + theme_bw() # Add pairwise comparisons p-value