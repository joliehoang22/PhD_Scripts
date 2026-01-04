#+NAME: DONE_livpm_de_by_anesthesia_dosing
#+BEGIN_SRC R

# setup
  rm(list=ls())
  options(stringsAsFactors=F)
  suppressMessages(library(readxl))
  suppressMessages(library(data.table))
  suppressMessages(library(variancePartition))
  suppressMessages(library(limma))
  suppressMessages(library(edgeR))
  suppressMessages(library(Glimma))
  suppressMessages(library(BiocParallel))
  suppressMessages(library(ggplot2))
  suppressMessages(library(ggthemes))
  suppressMessages(library(readxl))
  Sys.setenv(OMP_NUM_THREADS = 6)
  setwd("/sc/arion/projects/psychgen/lbp/")

# read in lbp data
  lbp <- readRDS("/sc/arion/projects/psychgen/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_Covs-featureCounts-vobjDream-Resids-LivPmDE_FINALModel_onlyBRAIN_518Samples_Excluding-Outliers-MislabeledSamples-BadSamples_19JUL2021.RDS")
  vob <- lbp$vobjDream
  met <- lbp$covariates

# read in sample sheet
  myfile <- "/sc/arion/projects/psychgen/lbp/files/sema4_bulk_rna_sample_sheet/Bulk_RNA_Isolation_Mastertable_BRAINANDBLOOD_forSEMA4_awcFormatted.tsv"
  dt <- fread(myfile, na=c("na","","NA"))[grep("brain", tissue, ignore.case=T)]
  dt <- dt[living==1,.(iid, sid=LBPSEMA4_ID, collection_date)]
  dt[,collection_date:=as.Date(collection_date)]

# read in anesthesia dosing data
  dop <- as.data.table(read_excel("./data/emr/lbp_lel2021_clinical/lel2021_anesthesia_data_pull_11APR2022.xlsx"), na=c(""))
  dop <- dop[,c("subject_id", "brain_surgery_date", "Dexmedetomidine_bolus", "Fentanyl_bolus", "Propofol_bolus", 
                "Dexmedetomidine_infusion", "Propofol_infusion"), with=F]
  dop[,brain_surgery_date:=as.Date(brain_surgery_date)]

# format
  dt <- merge(dt, dop, by.x=c("iid","collection_date"), by.y=c("subject_id","brain_surgery_date"), all.x=T)
  dt[is.na(Propofol_bolus),Propofol_bolus:=0]
  dt[is.na(Propofol_infusion),Propofol_infusion:=0]
  dt[is.na(Dexmedetomidine_bolus),Dexmedetomidine_bolus:=0]
  dt[is.na(Dexmedetomidine_infusion),Dexmedetomidine_infusion:=0]
  dt[is.na(Fentanyl_bolus),Fentanyl_bolus:=0]
  dt[,propofol:=Propofol_infusion+Propofol_bolus]
  dt[,dexmedetomidine:=Dexmedetomidine_infusion+Dexmedetomidine_bolus]
  dt[,fentanyl:=Fentanyl_bolus]
  dt <- dt[,.(iid, sid, propofol, dexmedetomidine, fentanyl)]
  dt[,x:=propofol+dexmedetomidine+fentanyl]
  dt[,MISS:=FALSE]
  dt[x==0,MISS:=TRUE]
  dt[,x:=NULL]
  dt2 <- dt[MISS==FALSE]


  table(dt2$MISS)
  #FALSE 
  #281

# define drug dosing groups
  ##
  ## view distributions
  ##
  cor(dt2$dexmedetomidine, dt2$propofol) #[1] 0.1184813 -- Alex; but Lora is [1] 0.06835994
  ggplot(dt2, aes(propofol)) + geom_histogram(fill="white", color="black", bins=50) + theme_base()
  dev.new()
  ggplot(dt2, aes(fentanyl)) + geom_histogram(fill="white", color="black", bins=50) + theme_base()
  dev.new()
  ggplot(dt2, aes(dexmedetomidine)) + geom_histogram(fill="white", color="black", bins=50) + theme_base()

# bin living samples by dose

  ## propofol
  dt2 <- dt2[order(propofol)]
  dt2[,propquartile := floor( 1 + 4 * (.I-1) / .N)]
  dt2[,propquartile := floor( 1 + 4 * (.I-1) / .N)]

  ## dexmedetomidine
  dt2 <- dt2[order(dexmedetomidine)]
  dt2[,dexquartile := floor( 1 + 4 * (.I-1) / .N)]
  dt2[,dexquartile := floor( 1 + 4 * (.I-1) / .N)]

  ## fentanyl
  dt2 <- dt2[order(fentanyl)]
  dt2[,fenquartile := floor( 1 + 4 * (.I-1) / .N)]
  dt2[,fenquartile := floor( 1 + 4 * (.I-1) / .N)]

  ## list
  lp1 <- dt2[propquartile==1]$sid
  lp2 <- dt2[propquartile==2]$sid
  lp3 <- dt2[propquartile==3]$sid
  lp4 <- dt2[propquartile==4]$sid
  ld1 <- dt2[dexquartile==1]$sid
  ld2 <- dt2[dexquartile==2]$sid
  ld3 <- dt2[dexquartile==3]$sid
  ld4 <- dt2[dexquartile==4]$sid
  lf1 <- dt2[fenquartile==1]$sid
  lf2 <- dt2[fenquartile==2]$sid
  lf3 <- dt2[fenquartile==3]$sid
  lf4 <- dt2[fenquartile==4]$sid

# postmortem bins

  ## all
  pmt <- met[mymet_postmortem==1]$SAMPLE_ISMMS

  ## split into 4s -- randomly sample 
  set.seed(676)
  hf1 <- sample(pmt, size=round(length(pmt)/2, 0), replace=FALSE)
  hf2 <- pmt[!pmt %in% hf1]
  pm1 <- sample(hf1, size=round(length(hf1)/2, 0), replace=FALSE)
  pm2 <- hf1[!hf1 %in% pm1]
  pm3 <- sample(hf2, size=round(length(hf2)/2, 0), replace=FALSE)
  pm4 <- hf2[!hf2 %in% pm3]

# counts
  length(pm1) #[1] 61
  length(pm2) #[1] 61
  length(pm3) #[1] 60
  length(pm4) #[1] 61
  length(lp1) #[1] 71
  length(lp2) #[1] 70
  length(lp3) #[1] 70
  length(lp4) #[1] 70
  length(ld1) #[1] 71
  length(ld2) #[1] 70
  length(ld3) #[1] 70
  length(ld4) #[1] 70
  length(lf1) #[1] 71
  length(lf2) #[1] 70
  length(lf3) #[1] 70
  length(lf4) #[1] 70


length(unique(met[met$SAMPLE_ISMMS %in% lp1, ]$SAMPLE_ISMMS))
# [1] 70
length(unique(met[met$SAMPLE_ISMMS %in% lp2, ]$SAMPLE_ISMMS))
# [1] 68
length(unique(met[met$SAMPLE_ISMMS %in% lp3, ]$SAMPLE_ISMMS))
# [1] 70
length(unique(met[met$SAMPLE_ISMMS %in% lp4, ]$SAMPLE_ISMMS))
# [1] 67
length(unique(met[met$SAMPLE_ISMMS %in% ld1, ]$SAMPLE_ISMMS))
# [1] 70
length(unique(met[met$SAMPLE_ISMMS %in% ld2, ]$SAMPLE_ISMMS))
# [1] 68
length(unique(met[met$SAMPLE_ISMMS %in% ld3, ]$SAMPLE_ISMMS))
# [1] 70
length(unique(met[met$SAMPLE_ISMMS %in% ld4, ]$SAMPLE_ISMMS))
# [1] 67
length(unique(met[met$SAMPLE_ISMMS %in% lf1, ]$SAMPLE_ISMMS))
# [1] 70
length(unique(met[met$SAMPLE_ISMMS %in% lf2, ]$SAMPLE_ISMMS))
# [1] 68
length(unique(met[met$SAMPLE_ISMMS %in% lf3, ]$SAMPLE_ISMMS))
# [1] 70
length(unique(met[met$SAMPLE_ISMMS %in% lf4, ]$SAMPLE_ISMMS))
# [1] 67


length(unique(met[met$SAMPLE_ISMMS %in% lp1, ]$IID_ISMMS))
#[1] 61
length(unique(met[met$SAMPLE_ISMMS %in% lp2, ]$IID_ISMMS))
# [1] 63
length(unique(met[met$SAMPLE_ISMMS %in% lp3, ]$IID_ISMMS))
# [1] 59
length(unique(met[met$SAMPLE_ISMMS %in% lp4, ]$IID_ISMMS))
# [1] 58

length(unique(met[met$SAMPLE_ISMMS %in% ld1, ]$IID_ISMMS))
# [1] 61
length(unique(met[met$SAMPLE_ISMMS %in% ld2, ]$IID_ISMMS))
# [1] 63
length(unique(met[met$SAMPLE_ISMMS %in% ld3, ]$IID_ISMMS))
# [1] 59
length(unique(met[met$SAMPLE_ISMMS %in% ld4, ]$IID_ISMMS))
# [1] 58

length(unique(met[met$SAMPLE_ISMMS %in% lf1, ]$IID_ISMMS))
# [1] 61
length(unique(met[met$SAMPLE_ISMMS %in% lf2, ]$IID_ISMMS))
# [1] 63
length(unique(met[met$SAMPLE_ISMMS %in% lf3, ]$IID_ISMMS))
# [1] 59
length(unique(met[met$SAMPLE_ISMMS %in% lf4, ]$IID_ISMMS))
# [1] 58


#dt[, .N, by=a]

# update metadata
  met <- lbp$covariates
  met[,proDE:="notassigned"]
  met[,dexDE:="notassigned"]
  met[,fenDE:="notassigned"]
  met[SAMPLE_ISMMS %in% lp1, proDE:="LIV_PQ1"]
  met[SAMPLE_ISMMS %in% lp2, proDE:="LIV_PQ2"]
  met[SAMPLE_ISMMS %in% lp3, proDE:="LIV_PQ3"]
  met[SAMPLE_ISMMS %in% lp4, proDE:="LIV_PQ4"]
  met[SAMPLE_ISMMS %in% ld1, dexDE:="LIV_DQ1"]
  met[SAMPLE_ISMMS %in% ld2, dexDE:="LIV_DQ2"]
  met[SAMPLE_ISMMS %in% ld3, dexDE:="LIV_DQ3"]
  met[SAMPLE_ISMMS %in% ld4, dexDE:="LIV_DQ4"]
  met[SAMPLE_ISMMS %in% lf1, fenDE:="LIV_FQ1"]
  met[SAMPLE_ISMMS %in% lf2, fenDE:="LIV_FQ2"]
  met[SAMPLE_ISMMS %in% lf3, fenDE:="LIV_FQ3"]
  met[SAMPLE_ISMMS %in% lf4, fenDE:="LIV_FQ4"]
  met[SAMPLE_ISMMS %in% pm1, proDE:="PM_PQ1"]
  met[SAMPLE_ISMMS %in% pm2, proDE:="PM_PQ2"]
  met[SAMPLE_ISMMS %in% pm3, proDE:="PM_PQ3"]
  met[SAMPLE_ISMMS %in% pm4, proDE:="PM_PQ4"]
  met[SAMPLE_ISMMS %in% pm1, dexDE:="PM_DQ1"]
  met[SAMPLE_ISMMS %in% pm2, dexDE:="PM_DQ2"]
  met[SAMPLE_ISMMS %in% pm3, dexDE:="PM_DQ3"]
  met[SAMPLE_ISMMS %in% pm4, dexDE:="PM_DQ4"]
  met[SAMPLE_ISMMS %in% pm1, fenDE:="PM_FQ1"]
  met[SAMPLE_ISMMS %in% pm2, fenDE:="PM_FQ2"]
  met[SAMPLE_ISMMS %in% pm3, fenDE:="PM_FQ3"]
  met[SAMPLE_ISMMS %in% pm4, fenDE:="PM_FQ4"]
  met[,proDE:=as.factor(proDE)]
  met[,dexDE:=as.factor(dexDE)]
  met[,fenDE:=as.factor(fenDE)]
  met <- as.data.frame(met)
  rownames(met) <- met$SAMPLE_ISMMS


length(unique(met[met$SAMPLE_ISMMS %in% lp1, ]$IID_ISMMS))

# formula
  form1 <- ~0 + proDE + (1|mymet_sex) + mymet_rin + neuronal + RNASeqMetrics_MEDIAN_3PRIME_BIAS + 
              RNASeqMetrics_PCT_MRNA_BASES + (1|IID_ISMMS) + (1|mymet_depletionbatch) + 
              InsertSizeMetrics_MEDIAN_INSERT_SIZE + AlignmentSummaryMetrics_STRAND_BALANCE_FIRST_OF_PAIR
  form2 <- ~0 + dexDE + (1|mymet_sex) + mymet_rin + neuronal + RNASeqMetrics_MEDIAN_3PRIME_BIAS + 
              RNASeqMetrics_PCT_MRNA_BASES + (1|IID_ISMMS) + (1|mymet_depletionbatch) + 
              InsertSizeMetrics_MEDIAN_INSERT_SIZE + AlignmentSummaryMetrics_STRAND_BALANCE_FIRST_OF_PAIR
  form3 <- ~0 + fenDE + (1|mymet_sex) + mymet_rin + neuronal + RNASeqMetrics_MEDIAN_3PRIME_BIAS + 
              RNASeqMetrics_PCT_MRNA_BASES + (1|IID_ISMMS) + (1|mymet_depletionbatch) + 
              InsertSizeMetrics_MEDIAN_INSERT_SIZE + AlignmentSummaryMetrics_STRAND_BALANCE_FIRST_OF_PAIR

# sanity check
  identical(rownames(met), colnames(vob$E)) #[1] TRUE

# contrasts
  p1Con <- getContrast(vob, form1, met, c(paste0("proDE","PM_PQ1"), paste0("proDE","LIV_PQ1")))
  p2Con <- getContrast(vob, form1, met, c(paste0("proDE","PM_PQ2"), paste0("proDE","LIV_PQ2")))
  p3Con <- getContrast(vob, form1, met, c(paste0("proDE","PM_PQ3"), paste0("proDE","LIV_PQ3")))
  p4Con <- getContrast(vob, form1, met, c(paste0("proDE","PM_PQ4"), paste0("proDE","LIV_PQ4")))
  d1Con <- getContrast(vob, form2, met, c(paste0("dexDE","PM_DQ1"), paste0("dexDE","LIV_DQ1")))
  d2Con <- getContrast(vob, form2, met, c(paste0("dexDE","PM_DQ2"), paste0("dexDE","LIV_DQ2")))
  d3Con <- getContrast(vob, form2, met, c(paste0("dexDE","PM_DQ3"), paste0("dexDE","LIV_DQ3")))
  d4Con <- getContrast(vob, form2, met, c(paste0("dexDE","PM_DQ4"), paste0("dexDE","LIV_DQ4")))
  f1Con <- getContrast(vob, form3, met, c(paste0("fenDE","PM_FQ1"), paste0("fenDE","LIV_FQ1")))
  f2Con <- getContrast(vob, form3, met, c(paste0("fenDE","PM_FQ2"), paste0("fenDE","LIV_FQ2")))
  f3Con <- getContrast(vob, form3, met, c(paste0("fenDE","PM_FQ3"), paste0("fenDE","LIV_FQ3")))
  f4Con <- getContrast(vob, form3, met, c(paste0("fenDE","PM_FQ4"), paste0("fenDE","LIV_FQ4")))
  L1 <- cbind(p1Con, p2Con, p3Con, p4Con)
  L2 <- cbind(d1Con, d2Con, d3Con, d4Con)
  L3 <- cbind(f1Con, f2Con, f3Con, f4Con)

# de
  fitmm1 <- dream( vob, form1, met, L1, BPPARAM = MulticoreParam(5))
  fitmm2 <- dream( vob, form2, met, L2, BPPARAM = MulticoreParam(5))
  fitmm3 <- dream( vob, form3, met, L3, BPPARAM = MulticoreParam(5))

# format de results
  p1De <- topTable(fitmm1, coef="p1Con", number=nrow(vob))
  p2De <- topTable(fitmm1, coef="p2Con", number=nrow(vob))
  p3De <- topTable(fitmm1, coef="p3Con", number=nrow(vob))
  p4De <- topTable(fitmm1, coef="p4Con", number=nrow(vob))
  d1De <- topTable(fitmm2, coef="d1Con", number=nrow(vob))
  d2De <- topTable(fitmm2, coef="d2Con", number=nrow(vob))
  d3De <- topTable(fitmm2, coef="d3Con", number=nrow(vob))
  d4De <- topTable(fitmm2, coef="d4Con", number=nrow(vob))
  f1De <- topTable(fitmm3, coef="f1Con", number=nrow(vob))
  f2De <- topTable(fitmm3, coef="f2Con", number=nrow(vob))
  f3De <- topTable(fitmm3, coef="f3Con", number=nrow(vob))
  f4De <- topTable(fitmm3, coef="f4Con", number=nrow(vob))
  p1De <- data.table(gene = rownames(p1De), p1De)[order(logFC)]
  p2De <- data.table(gene = rownames(p2De), p2De)[order(logFC)]
  p3De <- data.table(gene = rownames(p3De), p3De)[order(logFC)]
  p4De <- data.table(gene = rownames(p4De), p4De)[order(logFC)]
  d1De <- data.table(gene = rownames(d1De), d1De)[order(logFC)]
  d2De <- data.table(gene = rownames(d2De), d2De)[order(logFC)]
  d3De <- data.table(gene = rownames(d3De), d3De)[order(logFC)]
  d4De <- data.table(gene = rownames(d4De), d4De)[order(logFC)]
  f1De <- data.table(gene = rownames(f1De), f1De)[order(logFC)]
  f2De <- data.table(gene = rownames(f2De), f2De)[order(logFC)]
  f3De <- data.table(gene = rownames(f3De), f3De)[order(logFC)]
  f4De <- data.table(gene = rownames(f4De), f4De)[order(logFC)]
  p1De[, DEG:="NOTDEG"]
  p2De[, DEG:="NOTDEG"]
  p3De[, DEG:="NOTDEG"]
  p4De[, DEG:="NOTDEG"]
  d1De[, DEG:="NOTDEG"]
  d2De[, DEG:="NOTDEG"]
  d3De[, DEG:="NOTDEG"]
  d4De[, DEG:="NOTDEG"]
  f1De[, DEG:="NOTDEG"]
  f2De[, DEG:="NOTDEG"]
  f3De[, DEG:="NOTDEG"]
  f4De[, DEG:="NOTDEG"]
  p1De[adj.P.Val<0.05, DEG:="DEG"]
  p2De[adj.P.Val<0.05, DEG:="DEG"]
  p3De[adj.P.Val<0.05, DEG:="DEG"]
  p4De[adj.P.Val<0.05, DEG:="DEG"]
  d1De[adj.P.Val<0.05, DEG:="DEG"]
  d2De[adj.P.Val<0.05, DEG:="DEG"]
  d3De[adj.P.Val<0.05, DEG:="DEG"]
  d4De[adj.P.Val<0.05, DEG:="DEG"]
  f1De[adj.P.Val<0.05, DEG:="DEG"]
  f2De[adj.P.Val<0.05, DEG:="DEG"]
  f3De[adj.P.Val<0.05, DEG:="DEG"]
  f4De[adj.P.Val<0.05, DEG:="DEG"]
  p1De[logFC<0, LFC:="NEGLFC"]
  p2De[logFC<0, LFC:="NEGLFC"]
  p3De[logFC<0, LFC:="NEGLFC"]
  p4De[logFC<0, LFC:="NEGLFC"]
  d1De[logFC<0, LFC:="NEGLFC"]
  d2De[logFC<0, LFC:="NEGLFC"]
  d3De[logFC<0, LFC:="NEGLFC"]
  d4De[logFC<0, LFC:="NEGLFC"]
  f1De[logFC<0, LFC:="NEGLFC"]
  f2De[logFC<0, LFC:="NEGLFC"]
  f3De[logFC<0, LFC:="NEGLFC"]
  f4De[logFC<0, LFC:="NEGLFC"]
  p1De[logFC>0, LFC:="POSLFC"]
  p2De[logFC>0, LFC:="POSLFC"]
  p3De[logFC>0, LFC:="POSLFC"]
  p4De[logFC>0, LFC:="POSLFC"]
  d1De[logFC>0, LFC:="POSLFC"]
  d2De[logFC>0, LFC:="POSLFC"]
  d3De[logFC>0, LFC:="POSLFC"]
  d4De[logFC>0, LFC:="POSLFC"]
  f1De[logFC>0, LFC:="POSLFC"]
  f2De[logFC>0, LFC:="POSLFC"]
  f3De[logFC>0, LFC:="POSLFC"]
  f4De[logFC>0, LFC:="POSLFC"]
  myres <- list( "pro1" = p1De, "pro2" = p2De, "pro3" = p3De, "pro4" = p4De, 
                "dex1" = d1De, "dex2" = d2De, "dex3" = d3De, "dex4" = d4De, 
                "fen1" = f1De, "fen2" = f2De, "fen3" = f3De, "fen4" = f4De )

# correlate logFC values
  lCor <- c()
  iter <-  t(combn(names(myres), 2))
  for (i in 1:nrow(iter)){
      n1 <- iter[i,1]
      n2 <- iter[i,2]
      x0 <- lbp$livpmDE[,.(gene,logFC0=logFC)]
      x1 <- myres[[n1]][,.(gene, logFC1=logFC)] 
      x2 <- myres[[n2]][,.(gene, logFC2=logFC)] 
      mx <- merge(merge(x0, x1), x2)
      v1 <- cor(mx$logFC0, mx$logFC1, method="spearman")
      v2 <- cor(mx$logFC0, mx$logFC2, method="spearman")
      v3 <- cor(mx$logFC1, mx$logFC2, method="spearman")
      add <- data.table(ds1=n1, ds2=n2, ds1_vs_livpm=v1, ds2_vs_livpm=v2, ds1_vs_ds2=v3)
      lCor <- rbind(lCor, add)
  }


  # lCorx <- c()
  # iter <-  t(combn(names(myres), 2))
  # for (i in 1:nrow(iter)){
  #     n1 <- iter[i,1]
  #     n2 <- iter[i,2]
  #     x0 <- lbp$livpmDE[,.(gene,logFC0=logFC)]
  #     x1 <- myres[[n1]][,.(gene, logFC1=logFC)] 
  #     x2 <- myres[[n2]][,.(gene, logFC2=logFC)] 
  #     mx <- merge(merge(x0, x1), x2)
  #     v1 <- cor.test(mx$logFC0, mx$logFC1, method="spearman")$estimate
  #     v2 <- cor.test(mx$logFC0, mx$logFC2, method="spearman")$estimate
  #     v3 <- cor.test(mx$logFC1, mx$logFC2, method="spearman")$estimate
  #     p1 <- cor.test(mx$logFC0, mx$logFC1, method="spearman")$pvalue
  #     p2 <- cor.test(mx$logFC0, mx$logFC2, method="spearman")$pvalue
  #     p3 <- cor.test(mx$logFC1, mx$logFC2, method="spearman")$pvalue
  #     add <- data.table(ds1=n1, ds2=n2, ds1_vs_livpm=v1, ds1_vs_livpm.pval=p1, ds2_vs_livpm=v2, ds2_vs_livpm.pval=p2, ds1_vs_ds2=v3, ds1_vs_ds2.pval=p3)
  #     lCorx <- rbind(lCorx, add)
  # }


#only thing I need for the write-up from lCo is the ds1_vs_livpm corr for pro 1-4 dex 1-4, fen 1-4
lCor[!duplicated(lCor$ds1_vs_livpm),]

#      ds1  ds2 ds1_vs_livpm ds2_vs_livpm ds1_vs_ds2
#  1: pro1 pro2    0.9821339    0.9787190  0.9500085
#  2: pro2 pro3    0.9787190    0.9778465  0.9422922
#  3: pro3 pro4    0.9778465    0.9805065  0.9490005
#  4: pro4 dex1    0.9805065    0.9740595  0.9553877
#  5: dex1 dex2    0.9740595    0.9811381  0.9384219
#  6: dex2 dex3    0.9811381    0.9790153  0.9492573
#  7: dex3 dex4    0.9790153    0.9757209  0.9440981
#  8: dex4 fen1    0.9757209    0.9796996  0.9416855
#  9: fen1 fen2    0.9796996    0.9761962  0.9467593
# 10: fen2 fen3    0.9761962    0.9812104  0.9398807
# 11: fen3 fen4    0.9812104    0.9737455  0.9505678
# 12: fen4 -- see above in ds2_vs_livpm -- 0.9737455


cors <- lCor[!duplicated(lCor$ds1_vs_livpm),]
cors2 <- cors[, c(1, 3)]

cors2 <- rbind(cors2, data.table(ds1="fen4", ds1_vs_livpm=as.numeric(cors[11, 4])))

#     ds1 ds1_vs_livpm
#  1: pro1    0.9821339
#  2: pro2    0.9787190
#  3: pro3    0.9778465
#  4: pro4    0.9805065
#  5: dex1    0.9740595
#  6: dex2    0.9811381
#  7: dex3    0.9790153
#  8: dex4    0.9757209
#  9: fen1    0.9796996
# 10: fen2    0.9761962
# 11: fen3    0.9812104
# 12: fen4    0.9737455

summary(cors2$ds1_vs_livpm)
 #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 # 0.9737  0.9761  0.9789  0.9783  0.9807  0.9821

# save 
  #ggplot(lCor, aes(ds1_vs_livpm)) + geom_histogram(fill="grey", color="black", bins=50) + theme_base()
  #dev.new()
  #ggplot(lCor, aes(ds2_vs_livpm)) + geom_histogram(fill="grey", color="black", bins=50) + theme_base()
  #dev.new()
  #ggplot(lCor, aes(ds1_vs_ds2)) + geom_histogram(fill="grey", color="black", bins=50) + theme_base()
  #dev.off()
  #dev.off()
  #dev.off()
  saveRDS(list("de"=myres, "correlations"=lCor),
          file="/sc/arion/projects/psychgen/lbp/results/lel2021_livpm_by_anesthesia_dosing_11APR2022.RDS")

# plot dose distributions
  dt2[,propquartile:=paste0("Q",propquartile)]
  dt2[,dexquartile:=paste0("Q",dexquartile)]
  dt2[,fenquartile:=paste0("Q",fenquartile)]
  dt3 <- rbind( dt2[,.(sid, drug="propofol", dose=propofol, quartile=propquartile)],
               dt2[,.(sid, drug="dexmedetomidine", dose=dexmedetomidine, quartile=dexquartile)],
               dt2[,.(sid, drug="fentanyl", dose=fentanyl, quartile=fenquartile)])
  lCor2 <- unique(rbind(lCor[,.(ds=ds1, ds_vs_livpm=ds1_vs_livpm)], lCor[,.(ds=ds2, ds_vs_livpm=ds2_vs_livpm)]))
  ggplot(dt3, aes(dose)) + geom_histogram(fill="grey", color="black", bins=50) + facet_wrap(quartile ~ drug, scales="free") + theme_base()
  ggplot(lCor2, aes(ds_vs_livpm)) + geom_histogram(fill="grey", color="black", bins=50) + theme_base()
 
#+END_SRC

