module load R/4.0.3
R

rm(list=ls())
options(stringsAsFactors=F)
##################################################
#.libPaths(c("~/.Rlib", .libPaths())) #  need to do this before loading any packages!!!

library(variancePartition)#, lib.loc="~/.Rlib")
library(data.table)

####################################################################################################################

m1 <- fread("/sc/arion/projects/psychgen/lbp/files/m1TableForLiharska2021_updated30JUNE2021.tsv")

m1[, .(SAMPLE_NAME, DataProcessingStatus != "no_notable_issues")]

m1[DataProcessingStatus != "no_notable_issues", .(SAMPLE_NAME, DataProcessingStatus)]

issues <- as.data.table(m1[DataProcessingStatus != "no_notable_issues", .(SAMPLE_NAME, DataProcessingStatus)])


scp liharl02@chimera.hpc.mssm.edu:/sc/arion/projects/psychgen/lbp/files/m1TableForLiharska2021_updated30JUNE2021.tsv ~/Desktop/pca_plots/

#it looks like 4 samples are completely out due to either very low coverage (2) or the fastq being corrupt (2) -- BRAIN737 is one of these
#25 blood samples were run separately by eric for some unclear reason but do exist
#and then one that exists but fastq_in_two_deliveries|EV_RAPiD_run_on_one_fastq_delivery_because_other_massive
#and another that exists but fastq_in_three_deliveries|EV_RAPiD_run_on_two_fastq_delivery_because_unclear_if_third_is_unique


#qc.data <- readRDS("/sc/arion/projects/psychgen2/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_MYMET_QCMetrics_Merged_Compiled_777samples_NoBRAIN107_NoBRAIN734_5MAY2021.RDS")

dt <- readRDS("/sc/arion/projects/psychgen2/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_featureCounts_Compiled_777samples_NoBRAIN107_NoBRAIN734_5MAY2021.RDS")


#Decision: Proceed with what we have, add FCs and QC metrics for LBPSEMA4BRAIN107, the one that eventually finished but was not included in the original compiling

path <- "/sc/arion/projects/psychgen2/lbp/data/RAW/rna/bulk/fromSema4/LBPSEMA4BRAIN107_reRUN/LBPSEMA4BRAIN107/RAPiD/featureCounts/LBPSEMA4BRAIN107.primary.txt"
fc <- fread(path, header=T, skip=1)
dt <- merge(dt, fc, by=c("Geneid","Chr","Start","End","Strand","Length"), all.x=T, all.y=T)

saveRDS(dt, file="/sc/arion/projects/psychgen/lbp/files/files_from_psychgen2/lbp_allBatches_QC/lbp_allBatches_RAPiD_featureCounts_Compiled_777samples_NoBRAIN107_NoBRAIN734_AddedBRAIN107_04JAN2022.RDS")
####################################################################################################################
#get the qc metrics for LBPSEMA4BRAIN107 -- reference lbp_allBatches_RAPiD_qc1_fastqc_unzip_createQCmetrics_4MAY2021_Correct.sh

# summary files across samples; compile the individual files into one file for each qc metric  

i=LBPSEMA4BRAIN107

OUT_DIR=/sc/arion/projects/psychgen2/lbp/scratch/RAPiD_allBatches_qcMetrics_mergedWith_LBPSEMA4BRAIN107
WORK_DIR=/sc/arion/projects/psychgen2/lbp/data/RAW/rna/bulk/fromSema4/LBPSEMA4BRAIN107_reRUN/LBPSEMA4BRAIN107/RAPiD

fqc=${WORK_DIR}/fastqc/${i}.Aligned.out_fastqc/summary.txt
str=${WORK_DIR}/star/${i}.Log.final.out
ftc=${WORK_DIR}/featureCounts/${i}.primary.txt.summary
qc1=${WORK_DIR}/qc_metrics/${i}.AlignmentSummaryMetrics
qc2=${WORK_DIR}/qc_metrics/${i}.DupMetrics
qc3=${WORK_DIR}/qc_metrics/${i}.GcBiasSummaryMetrics
qc4=${WORK_DIR}/qc_metrics/${i}.InsertSizeMetrics
qc5=${WORK_DIR}/qc_metrics/${i}.RNASeqMetrics
      
cat ${fqc} >> ${OUT_DIR}/RAPiD.fastqc_summary.txt 
grep "|" ${str} | awk 'NR>5' | awk '{$1=$1}1' | sed s/" | "/"|"/g | tr ' ' '_' | awk -F"|" -v OFS="\t" -v I=${i} '{print I, $1, $2}' >> ${OUT_DIR}/RAPiD.star_summary.txt
grep -i assigned ${ftc} | awk -v OFS="\t" -v I=${i} '{print I, $1, $2}' >> ${OUT_DIR}/RAPiD.featureCounts_summary.txt
awk -v OFS="\t" -v I=${i} 'NR>7 {print I, $0}' ${qc1} >> ${OUT_DIR}/RAPiD.qc_metrics.AlignmentSummaryMetrics
      
if [[ -f ${qc2} ]]
then awk -v OFS="\t" 'NR==8 {print $0}' ${qc2} >> ${OUT_DIR}/RAPiD.qc_metrics.DupMetrics #some samples do not have dupmetrics, a handful
fi 
      
awk -v OFS="\t" -v I=${i} 'NR==8 {print I, $0}' ${qc3} >> ${OUT_DIR}/RAPiD.qc_metrics.GcBiasSummaryMetrics
awk -v OFS="\t" -v I=${i} 'NR==8 {print I, $0}' ${qc4} >> ${OUT_DIR}/RAPiD.qc_metrics.InsertSizeMetrics
awk -v OFS="\t" -v I=${i} 'NR==8 {print I, $0}' ${qc5} >> ${OUT_DIR}/RAPiD.qc_metrics.RNASeqMetrics
   

#unix logical aside
if [[ 1 == 1 ]]; then echo yes; fi #can't just do a boolean like in R
if [[ ! -f ${qc2} ]]; then echo no; fi


####################################################################################################################
library(data.table)
  library(readxl)
  library(ggplot2)
  library(ggthemes)
  library(PCAmixdata)
  library(openxlsx)

setwd("/sc/arion/projects/psychgen2/lbp/scratch/scratch_RAPiD/allBatches_qcMetrics/")

# master metadata for lbp
  mymet = fread("/sc/arion/projects/psychgen/lbp/files/sema4_bulk_rna_sample_sheet/Bulk_RNA_Isolation_Mastertable_BRAINANDBLOOD_forSEMA4_awcFormatted.tsv")

  #mymet <- read.xlsx("/sc/arion/projects/psychgen/lbp/files/sema4_bulk_rna_sample_sheet/Bulk_RNA_Isolation_Mastertable_BRAINANDBLOOD.xlsx")  
  
  #mymet <- as.data.table(mymet)
  ms.idmap = mymet[,.(SAMPLE.ISMMS=LBPSEMA4_ID, BARCODE.ISMMS=Barcode, IID=iid)]

  keep = c("LBPSEMA4_ID", "PLATE", "PLATE_ROW", "PLATE_COLUMN", "extraction_batch", "bank", "age", "sex", "race", "ethnicity", "living", 
           "phe", "collection_date", "collection_time", "extraction_date", "extraction_kit", "extraction_conc_ngul", 
           "extraction_amount", "extraction_rin","brain", "timepoint", "Sequencing_Batch", "Depletion_Batch")

  mymet = mymet[,keep,with=F]
  colnames(mymet) = c("SAMPLE.ISMMS", "mymet.PLATE", "mymet.PLATEROW", "mymet.PLATECOL", "mymet.extractionbatch", 
                      "mymet.bank", "mymet.age", "mymet.sex", "mymet.race", "mymet.ethnicity", "mymet.living", "mymet.phe", "mymet.collectiondate", 
                      "mymet.collectiontime", "mymet.extractiondate", "mymet.extractionkit", "mymet.rna_conc_ngul", 
                      "mymet.rna_amount_ng", "mymet.rin", "mymet.brain", "mymet.timepoint", 
                      "mymet.seqbatch", "mymet.depletionbatch")

# spreadsheets from sema4 for first 192 samples -- we do not have this for these batches
  s4qc.alignment = "/sc/arion/projects/psychgen/lbp/files/spreadsheets_from_sema4/16JUL2019/lbp_fromSEMA4_alignmentQC.livingbrain_SB1_SB3_QC.xlsx"
  s4qc.sample = "/sc/arion/projects/psychgen/lbp/files/spreadsheets_from_sema4/25JUL2019/lbp_fromSEMA4_sampleQC.Alex_Update2_25Jul2019.xlsx"
  s4qc.alignment = as.data.table(read_excel(s4qc.alignment))
  #s4qc.sample = as.data.table(read_excel(s4qc.sample))
  s4qc.sample = as.data.table(read_excel(s4qc.sample, col_types=c(rep("text", 7), rep("numeric", 9), rep("text", 4), rep("numeric", 4), "text")))
  s4qc.sample.names = c("RSM.SEMA4","sema4sampleqc.nycinfo.SequencingBatch","sema4sampleqc.nycinfo.DepletionBatch", 
                        "sema4sampleqc.nycinfo.ISM","sema4sampleqc.nycinfo.SAMPLE_NAME","sema4sampleqc.nycinfo.INDIVIDUAL_NAME","sema4sampleqc.nycinfo.LIMS_Name",
                        "sema4sampleqc.rnaqc.AlexLabConc_ngul","sema4sampleqc.rnaqc.VolumeSubmitted_ul","sema4sampleqc.rnaqc.BranfordLabConc_ngul",
                        "sema4sampleqc.rnaqc.Sampleusedforprep_ul","sema4sampleqc.rnaqc.Sampleusedforprep_ng","sema4sampleqc.rnaqc.RemainingVolume_ul",
                        "sema4sampleqc.rnaqc.RemainingSample_ng","sema4sampleqc.rnaqc.dv200","sema4sampleqc.rnaqc.RIN",
                        "sema4sampleqc.libprep.i7Adapter","sema4sampleqc.libprep.AdapterSeqi7","sema4sampleqc.libprep.i5Adapter",
                        "sema4sampleqc.libprep.AdapterSeqi5","sema4sampleqc.postenrichmentqc.QubitConc_ngul",
                        "sema4sampleqc.postenrichmentqc.TapeStationAvgSize_bp","sema4sampleqc.postenrichmentqc.TapeStationMolarity_nM",
                        "sema4sampleqc.postenrichmentqc.CalculatedMolarity_nM","sema4sampleqc.comments")
  s4qc.sample = s4qc.sample[2:nrow(s4qc.sample),]
  colnames(s4qc.sample) = s4qc.sample.names

# id map 
  s4.idmap = s4qc.sample[,.(SAMPLE.ISMMS=sema4sampleqc.nycinfo.SAMPLE_NAME,
                         RSM.SEMA4=RSM.SEMA4,ISM.SEMA4=sema4sampleqc.nycinfo.ISM, 
                         IID=sema4sampleqc.nycinfo.INDIVIDUAL_NAME, 
                         LIMS.SEMA4=sema4sampleqc.nycinfo.LIMS_Name)]
  idmap = merge(s4.idmap, ms.idmap, by="SAMPLE.ISMMS", suffixes=c(".SEMA4", ".ISMMS"))
  idmap[, ISM.SEMA4:=tstrsplit(ISM.SEMA4, split="-", fixed=T, keep=1L)]


# combine metadata with sample qc from sema4
  a = merge(idmap, mymet)
  b = merge(idmap, s4qc.sample, by="RSM.SEMA4")
  master.qc = merge(a, b, by=intersect(colnames(a), colnames(b)))

###################################################################################################################

# read in RAPiD QC metrics
#setwd("/sc/arion/projects/psychgen2/lbp/scratch/RAPiD_allBatches_qcMetrics/")
setwd("/sc/arion/projects/psychgen2/lbp/scratch/RAPiD_allBatches_qcMetrics_mergedWith_LBPSEMA4BRAIN107")


# RAPiD.GcBiasSummaryMetrics -- the header
# RAPiD.InsertSizeMetrics -- the header
# RAPiD.fastqc_summary.txt
# RAPiD.star_summary.txt
# RAPiD.featureCounts_summary.txt
# RAPiD.qc_metrics.AlignmentSummaryMetrics
# RAPiD.qc_metrics.DupMetrics -- some samples do not have these
# RAPiD.qc_metrics.GcBiasSummaryMetrics
# RAPiD.qc_metrics.InsertSizeMetrics
# RAPiD.qc_metrics.RNASeqMetrics

  ##
  ##lbpBatch1 - fastqc
  ##
  lbpBatch1.fqc = fread("RAPiD.fastqc_summary.txt", header=F)
  lbpBatch1.fqc[, V2:=paste("FASTQC.", make.names(V2), sep="")]
  lbpBatch1.fqc[, V3:=gsub(".Aligned.out.bam", "", fixed=T, V3)]
  #lbpBatch1.fqc[, V3:=tstrsplit(gsub("Sample_", "", V3), split="-", fixed=T, keep=1L)]
  lbpBatch1.fqc = dcast( lbpBatch1.fqc, V3 ~ V2, value.var="V1", fill=NA)
  #colnames(lbpBatch1.fqc)[1] = "ISM.SEMA4" 
  colnames(lbpBatch1.fqc)[1] = "SAMPLE.ISMMS"

  length(unique(lbpBatch1.fqc$V3))
  #[1] 0?; previous was 778

  dim(lbpBatch1.fqc)
  #[1] 779  12

  ##
  ##lbpBatch1 - featurecounts
  ##
  lbpBatch1.ftc = fread("RAPiD.featureCounts_summary.txt", header=F)
  lbpBatch1.ftc[, V2:=paste("FEATURECOUNTS.", make.names(V2), sep="")]
  #lbpBatch1.ftc[, V1:=tstrsplit(gsub("Sample_", "", V1), split="-", fixed=T, keep=1L)]
  lbpBatch1.ftc = dcast( lbpBatch1.ftc, V1 ~ V2, value.var="V3", fill=NA)
  #colnames(lbpBatch1.ftc)[1] = "ISM.SEMA4"
  colnames(lbpBatch1.ftc)[1] = "SAMPLE.ISMMS"
  ##
  
  ##lbpBatch1 - star
  ##
  lbpBatch1.str = fread("RAPiD.star_summary.txt", header=F) 
  #lbpBatch1.str[, V1:=tstrsplit(gsub("Sample_", "", V1), split="-", fixed=T, keep=1L)]
  lbpBatch1.str[, V2:=paste("STAR.", gsub("._",".",fixed=T,make.names(gsub("%","pct",V2))), sep="")]
  lbpBatch1.str = dcast( lbpBatch1.str, V1 ~ V2, value.var="V3", fill=NA)
  colnames(lbpBatch1.str)[1] = "SAMPLE.ISMMS"
  haspctsign = c("STAR.Deletion_rate_per_base","STAR.Insertion_rate_per_base","STAR.Mismatch_rate_per_base.pct","STAR.Uniquely_mapped_reads_pct",
                 "STAR.pct_of_chimeric_reads","STAR.pct_of_reads_mapped_to_multiple_loci","STAR.pct_of_reads_mapped_to_too_many_loci",
                 "STAR.pct_of_reads_unmapped.other","STAR.pct_of_reads_unmapped.too_many_mismatches","STAR.pct_of_reads_unmapped.too_short")
  for (i in haspctsign) lbpBatch1.str[[i]] = as.numeric(gsub("%","", lbpBatch1.str[[i]]))/100
  ##
  
  ##lbpBatch1 - AlignmentSummaryMetrics
  ##
  lbpBatch1.qc1 = fread("RAPiD.qc_metrics.AlignmentSummaryMetrics", fill=T, na=c("NA",""))
  lbpBatch1.qc1 = lbpBatch1.qc1[,1:(ncol(lbpBatch1.qc1)-3)] #3 empty columns at the end
  lbpBatch1.qc1 = lbpBatch1.qc1[!is.na(CATEGORY)]
  lbpBatch1.qc1[, SAMPLE.ISMMS:=SAMPLE] 
  lbpBatch1.qc1[, SAMPLE:=NULL]
  tmp0 = lbpBatch1.qc1[CATEGORY=="FIRST_OF_PAIR"]
  tmp1 = lbpBatch1.qc1[CATEGORY=="SECOND_OF_PAIR"]
  tmp2 = lbpBatch1.qc1[CATEGORY=="PAIR"]
  tmp0[,CATEGORY:=NULL]
  tmp1[,CATEGORY:=NULL]
  tmp2[,CATEGORY:=NULL]
  colnames(tmp0)[colnames(tmp0) != "SAMPLE.ISMMS"] = paste("AlignmentSummaryMetrics", paste(colnames(tmp0)[colnames(tmp0) != "SAMPLE.ISMMS"], "FIRST_OF_PAIR", sep="."), sep=".")
  colnames(tmp1)[colnames(tmp1) != "SAMPLE.ISMMS"] = paste("AlignmentSummaryMetrics", paste(colnames(tmp1)[colnames(tmp1) != "SAMPLE.ISMMS"], "SECOND_OF_PAIR", sep="."), sep=".")
  colnames(tmp2)[colnames(tmp2) != "SAMPLE.ISMMS"] = paste("AlignmentSummaryMetrics", paste(colnames(tmp2)[colnames(tmp2) != "SAMPLE.ISMMS"], "PAIR", sep="."), sep=".")
  lbpBatch1.qc1 = merge(merge(tmp2, tmp0), tmp1)
  ##
  

  ##lbpBatch1 - DupMetrics
  ##
  lbpBatch1.qc2 = fread("RAPiD.qc_metrics.DupMetrics")
  #lbpBatch1.qc2[, ISM.SEMA4:=tstrsplit(gsub("Sample_", "", SAMPLE), split="-", fixed=T, keep=1L)]
  lbpBatch1.qc2[, SAMPLE.ISMMS:=SAMPLE]
  lbpBatch1.qc2[, SAMPLE:=NULL]
  colnames(lbpBatch1.qc2)[colnames(lbpBatch1.qc2) != "SAMPLE.ISMMS"] = paste("DupMetrics",colnames(lbpBatch1.qc2)[colnames(lbpBatch1.qc2) != "SAMPLE.ISMMS"], sep=".")

  dim(lbpBatch1.qc2)
  #[1] 368  10
  
  ##
  ##lbpBatch1 - GcBiasSummaryMetrics
  

  header <- fread("RAPiD.GcBiasSummaryMetrics")
  header <- colnames(header)

  lbpBatch1.qc3 = fread("RAPiD.qc_metrics.GcBiasSummaryMetrics") #no header wtf
  lbpBatch1.qc3 = lbpBatch1.qc3[,1:(ncol(lbpBatch1.qc3)-3)]

  colnames(lbpBatch1.qc3) <- header

  lbpBatch1.qc3[, SAMPLE.ISMMS:=SAMPLE]
  #lbpBatch1.qc3[, ISM.SEMA4:=tstrsplit(gsub("Sample_", "", SAMPLE), split="-", fixed=T, keep=1L)]
  lbpBatch1.qc3[, SAMPLE:=NULL]
  colnames(lbpBatch1.qc3)[colnames(lbpBatch1.qc3) != "SAMPLE.ISMMS"] = paste("GcBiasMetrics",colnames(lbpBatch1.qc3)[colnames(lbpBatch1.qc3) != "SAMPLE.ISMMS"], sep=".")
  

  ##
  ##lbpBatch1 - InsertSizeMetrics

  ##
  header <- fread("RAPiD.InsertSizeMetrics")
  header <- colnames(header)

  lbpBatch1.qc4 = fread("RAPiD.qc_metrics.InsertSizeMetrics") #no header
  lbpBatch1.qc4 = lbpBatch1.qc4[,1:(ncol(lbpBatch1.qc4)-3)]
  
  colnames(lbpBatch1.qc4) <- header
  
  #lbpBatch1.qc4[, ISM.SEMA4:=tstrsplit(gsub("Sample_", "", SAMPLE), split="-", fixed=T, keep=1L)]
  lbpBatch1.qc4[, SAMPLE.ISMMS:=SAMPLE]
  lbpBatch1.qc4[, SAMPLE:=NULL]
  colnames(lbpBatch1.qc4)[colnames(lbpBatch1.qc4) != "SAMPLE.ISMMS"] = paste("InsertSizeMetrics",colnames(lbpBatch1.qc4)[colnames(lbpBatch1.qc4) != "SAMPLE.ISMMS"], sep=".")
  
  ##
  ##lbpBatch1 - RNASeqMetrics

  ##
  lbpBatch1.qc5 = fread("RAPiD.qc_metrics.RNASeqMetrics")
  lbpBatch1.qc5 = lbpBatch1.qc5[,1:(ncol(lbpBatch1.qc5)-3)]
  lbpBatch1.qc5[, SAMPLE.ISMMS:=SAMPLE]
  lbpBatch1.qc5[, SAMPLE:=NULL]
  colnames(lbpBatch1.qc5)[colnames(lbpBatch1.qc5) != "SAMPLE.ISMMS"] = paste("RNASeqMetrics",colnames(lbpBatch1.qc5)[colnames(lbpBatch1.qc5) != "SAMPLE.ISMMS"], sep=".")
  
  ##
  ##lbpBatch1 - merge
  ##
  
  #lbpBatch1.all = merge(merge(merge(merge(merge(merge(lbpBatch1.fqc,lbpBatch1.str),lbpBatch1.qc1),lbpBatch1.qc2),lbpBatch1.qc3),lbpBatch1.qc4),lbpBatch1.qc5)
  

  lbpBatch1.all = merge(merge(merge(merge(merge(merge(lbpBatch1.fqc,lbpBatch1.str),lbpBatch1.qc1),lbpBatch1.ftc),lbpBatch1.qc3),lbpBatch1.qc4),lbpBatch1.qc5)

  dim(lbpBatch1.all)
  #[1] 779 179

  #lbpBatch1.qc2 -- do DupMetrics separately

  lbpBatch1.all = merge(lbpBatch1.all, lbpBatch1.qc2, all.x=TRUE)

  dim(lbpBatch1.all)
  #[1] 779 188
  dim(lbpBatch1.qc2)
  #[1] 368  10

#REMOVE THE FAKE BAM SAMPLE LBPSEMA4BRAIN734 FROM THE TABLE BEFORE SAVING!!

grep("LBPSEMA4BRAIN734", lbpBatch1.all$SAMPLE.ISMMS)

lbpBatch1.all <- lbpBatch1.all[!SAMPLE.ISMMS == "LBPSEMA4BRAIN734", ]

dim(lbpBatch1.all)
#[1] 778 188


  old <- readRDS(file="/sc/arion/projects/psychgen2/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_QCMetrics_Compiled_777samples_NoBRAIN107_NoBRAIN734_5MAY2021.RDS") #formatting looks like the current so that's good

  saveRDS(lbpBatch1.all, file="/sc/arion/projects/psychgen/lbp/files/files_from_psychgen2/lbp_allBatches_QC/lbp_allBatches_RAPiD_QCMetrics_Compiled_778samples_NoBRAIN107_NoBRAIN734_AddedBRAIN107_04JAN2022.RDS")

  ##########################################################################

  # RAPiD metadata
  #load("formatted_qcdata.Rdata") #qc.data
  qc.data <- readRDS("/sc/arion/projects/psychgen/lbp/files/files_from_psychgen2/lbp_allBatches_QC/lbp_allBatches_RAPiD_QCMetrics_Compiled_778samples_NoBRAIN107_NoBRAIN734_AddedBRAIN107_04JAN2022.RDS")


  #merge with the ID table
  #ids <- readRDS("/sc/arion/projects/psychgen2/lbp/files/lbp_batch1_reRAPiD_QC3/lbp_batch1_reRAPiD_106brainSamples_IDmapping_ISM-SEMA4-to-SAMPLE-ISMMS_3FEB2021.RDS")

  map <- as.data.table(read_excel("/sc/arion/projects/psychgen/lbp/files/spreadsheets_from_sema4/25JUL2019/lbp_fromSEMA4_sampleQC.Alex_Update2_25Jul2019.xlsx", skip=1))
  ids <- map[,1:7]
  #ids[,ISM:=paste0("Sample_", ISM, "-2")]
  ids[,ISM:=gsub("-1", "", ISM)]

  colnames(ids)[5] <- "SAMPLE.ISMMS"

  qc.data[!qc.data$SAMPLE.ISMMS %in% ids$SAMPLE.ISMMS]$SAMPLE.ISMMS

#  [1] "LBPSEMA4BLOOD082_PT-0057" "LBPSEMA4BLOOD116_PT-0152"
#  [3] "LBPSEMA4BLOOD283_PT-0109" "LBPSEMA4BLOOD387_PT-0114"
#  [5] "LBPSEMA4BLOOD552_PT-0086" "LBPSEMA4BRAIN050_PT-0082"
#  [7] "LBPSEMA4BRAIN152_S14518"  "LBPSEMA4BRAIN192_PT-0145"
#  [9] "LBPSEMA4BRAIN388_PT-0178" "LBPSEMA4BRAIN393_S01745" 
# [11] "LBPSEMA4BRAIN410_PT-0159" "LBPSEMA4BRAIN435_PT-0160"
# [13] "LBPSEMA4BRAIN450_S02474"  "LBPSEMA4BRAIN489_S19433" 
# [15] "LBPSEMA4BRAIN511_PT-0033" "LBPSEMA4BRAIN747_PT-0123"
# [17] "LBPSEMA4BRAIN759_S12545" 


 qc.data$SAMPLE.ISMMS <- gsub("_", ".", qc.data$SAMPLE.ISMMS)
 qc.data$SAMPLE.ISMMS<- gsub("\\..*","", qc.data$SAMPLE.ISMMS)

qc2 <- merge(qc.data, ids, by= "SAMPLE.ISMMS")

  dim(qc.data)
  #[1] 778 188
  dim(qc2)
  #[1] 778 194



# metadata 
mymet = fread("/sc/arion/projects/psychgen/lbp/files/sema4_bulk_rna_sample_sheet/Bulk_RNA_Isolation_Mastertable_BRAINANDBLOOD_forSEMA4_awcFormatted.tsv")

keep = c("LBPSEMA4_ID", "PLATE", "PLATE_ROW", "PLATE_COLUMN", "extraction_batch", "bank", "age", "sex", "race", "ethnicity", "living", "phe", "collection_date", "collection_time", "extraction_date", "extraction_kit", "extraction_conc_ngul", "extraction_amount", "extraction_rin","brain", "timepoint", "Sequencing_Batch", "Depletion_Batch", "iid", "tissue", "preservative", "Barcode")

colnames(mymet)[! colnames(mymet) %in% keep]
# [1] "Position"                   "extraction_conc_bioA"      
# [3] "extraction_conc_bioA_units" "notes" 


mymet = mymet[,keep,with=F]

colnames(mymet) = c("SAMPLE.ISMMS", "mymet.PLATE", "mymet.PLATEROW", "mymet.PLATECOL", "mymet.extractionbatch", "mymet.bank", "mymet.age", "mymet.sex", "mymet.race", "mymet.ethnicity", "mymet.living", "mymet.phe", "mymet.collectiondate", "mymet.collectiontime", "mymet.extractiondate", "mymet.extractionkit", "mymet.rna_conc_ngul", "mymet.rna_amount_ng", "mymet.rin", "mymet.brain", "mymet.timepoint", "mymet.seqbatch", "mymet.depletionbatch", "mymet.iid", "mymet.tissue", "mymet.preservative", "mymet.barcode")

dim(mymet)
#[1] 779  27

  #add mymet. to all column names and replace _ with .
  #colnames(mymet) <- paste("mymet.", colnames(mymet), sep = "_")
  #colnames(mymet) <- paste0("mymet.", colnames(mymet))
  #colnames(mymet) <- gsub("_", ".", colnames(mymet))


  ##
  ## define a postmortem indicator
  ##
  mymet <- as.data.table(mymet)
  # mymet[,mymet.postmortem:=0]
  # mymet[mymet.living=="Postmortem", mymet.postmortem:=1]
  #colnames(mymet)[1] <- "SAMPLE.ISMMS"
  ##

  mymet2 <- merge(mymet, qc2, by="SAMPLE.ISMMS")
  dim(mymet2)
  #[1] 777 220

  mymetx <- mymet
  mymet <- mymet2

  ## categorize columns by source
  ##
  # mymet.colnames.ids <- c( "ISM.SEMA4", "SAMPLE.ISMMS", "RSM.SEMA4", "IID.SEMA4", "LIMS.SEMA4", "BARCODE.ISMMS", "IID.ISMMS", 
  #                          "sema4sampleqc.nycinfo.ISM", "sema4sampleqc.nycinfo.SAMPLE_NAME", "sema4sampleqc.nycinfo.INDIVIDUAL_NAME",
  #                          "sema4sampleqc.nycinfo.LIMS_Name")

  which(! mymet$INDIVIDUAL_NAME %in% mymet$mymet.iid)
  mymet[118, c("mymet.iid", "INDIVIDUAL_NAME")]
   # mymet.iid INDIVIDUAL_NAME
   # PT-0022               0

   mymet$INDIVIDUAL_NAME <- NULL


  identical(mymet$mymet.seqbatch, mymet$"Sequencing Batch")
  #[1] TRUE

  identical(mymet$mymet.depletionbatch, mymet$"Depletion Batch")
  #[1] TRUE

  mymet$"Depletion Batch" <- NULL
  mymet$"Sequencing Batch" <- NULL


  colnames(mymet)[which(colnames(mymet) == "ISM")] <- "ISM.SEMA4"
  colnames(mymet)[which(colnames(mymet) == "LIMS Name")] <- "LIMS.SEMA4"
  colnames(mymet)[which(colnames(mymet) == "RSM")] <- "RSM.SEMA4"
  colnames(mymet)[which(colnames(mymet) == "mymet.barcode")] <- "BARCODE.ISMMS"
  colnames(mymet)[which(colnames(mymet) == "mymet.iid")] <- "IID.ISMMS"

  mymet.colnames.ids <- c("ISM.SEMA4", "SAMPLE.ISMMS", "RSM.SEMA4", "LIMS.SEMA4", "BARCODE.ISMMS", "IID.ISMMS") 

  mymet.colnames.fqc <- grep("FASTQC", colnames(mymet), value=T)
  mymet.colnames.ftc <- grep("FEATURECOUNTS", colnames(mymet), value=T)
  mymet.colnames.str <- grep("STAR", colnames(mymet), value=T)
  mymet.colnames.asm <- grep("AlignmentSummaryMetrics", colnames(mymet), value=T)
  mymet.colnames.dup <- grep("DupMetrics", colnames(mymet), value=T)
  mymet.colnames.gcb <- grep("GcBiasMetrics", colnames(mymet), value=T)
  mymet.colnames.ism <- grep("InsertSizeMetrics", colnames(mymet), value=T)
  mymet.colnames.rsm <- grep("RNASeqMetrics", colnames(mymet), value=T)
  
  # mymet.colnames.sm4 <- c( "sema4sampleqc.nycinfo.DepletionBatch", "sema4sampleqc.rnaqc.AlexLabConc_ngul", 
  #                         "sema4sampleqc.rnaqc.VolumeSubmitted_ul", "sema4sampleqc.rnaqc.BranfordLabConc_ngul", 
  #                         "sema4sampleqc.rnaqc.Sampleusedforprep_ul", "sema4sampleqc.rnaqc.Sampleusedforprep_ng", 
  #                         "sema4sampleqc.rnaqc.RemainingVolume_ul", "sema4sampleqc.rnaqc.RemainingSample_ng", "sema4sampleqc.rnaqc.dv200", 
  #                         "sema4sampleqc.rnaqc.RIN", "sema4sampleqc.postenrichmentqc.QubitConc_ngul", 
  #                         "sema4sampleqc.postenrichmentqc.TapeStationAvgSize_bp", "sema4sampleqc.postenrichmentqc.TapeStationMolarity_nM",
  #                         "sema4sampleqc.postenrichmentqc.CalculatedMolarity_nM")

  # mymet.colnames.awc <- c( "mymet.PLATE", "mymet.extractionbatch", "mymet.bank", "mymet.age", 
  #                         "mymet.sex", "mymet.race", "mymet.ethnicity", "mymet.living", "mymet.postmortem", 
  #                         "mymet.phe", "mymet.collectiondate", "mymet.collectiontime", 
  #                         "mymet.extractiondate", "mymet.rna_conc_ngul", "mymet.rna_amount_ng", "mymet.rin", "mymet.brain", 
  #                         "mymet.timepoint", "mymet.seqbatch", "mymet.depletionbatch")

  mymet.colnames.awc <- c("mymet.PLATE", "mymet.extractionbatch", "mymet.bank", "mymet.age", "mymet.sex", "mymet.race", "mymet.ethnicity", "mymet.living", "mymet.phe", "mymet.collectiondate", "mymet.collectiontime", "mymet.extractiondate", "mymet.rna_conc_ngul", "mymet.rna_amount_ng", "mymet.rin", "mymet.timepoint", "mymet.seqbatch", "mymet.depletionbatch", "mymet.brain", "mymet.tissue", "mymet.preservative", "mymet.extractionkit")


  #colnames(mymet)[colnames(mymet)=="mymet.extraction_batch"] <- "mymet.extractionbatch"
  #colnames(mymet)[colnames(mymet)=="mymet.collection_date"] <- "mymet.collectiondate"
  #colnames(mymet)[colnames(mymet)=="mymet.collection_time"] <- "mymet.collectiontime"
  #colnames(mymet)[colnames(mymet)=="mymet.extraction_date"] <- "mymet.extractiondate"
  # colnames(mymet)[colnames(mymet)=="mymet.rna.conc.ngul"] <- "mymet.rna_conc_ngul"
  # colnames(mymet)[colnames(mymet)=="mymet.rna.amount.ng"] <- "mymet.rna_amount_ng"
  #colnames(mymet)[colnames(mymet)=="mymet.extraction_rin"] <- "mymet.rin"
  #colnames(mymet)[colnames(mymet)=="mymet.Depletion_Batch"] <- "mymet.depletionbatch"
  #colnames(mymet)[colnames(mymet)=="mymet.Sequencing_Batch"] <- "mymet.seqbatch"
  #colnames(mymet)[colnames(mymet)=="mymet.notes.delete"] <- "mymet.notes"

  #colnames(mymet)[! colnames(mymet) %in% mymet.colnames.awc]
  mymet.colnames.awc[! mymet.colnames.awc %in% colnames(mymet)] 
  #character(0)

  #mymet.colnames.awc <- mymet.colnames.awc[-9] #this removed mymet.postmortem but it was already removed from colnames.awc above

  ##
  ## list the columns that have actual covariates (as opposed to IDs)
  ##
  # mymet.datacol <- c( mymet.colnames.fqc, mymet.colnames.ftc, mymet.colnames.str, mymet.colnames.asm, 
  #                    mymet.colnames.dup, mymet.colnames.gcb, mymet.colnames.ism, mymet.colnames.rsm, 
  #                    mymet.colnames.sm4, mymet.colnames.awc)

  mymet.datacol <- c( mymet.colnames.fqc, mymet.colnames.ftc, mymet.colnames.str, mymet.colnames.asm, 
                     mymet.colnames.dup, mymet.colnames.ism, mymet.colnames.rsm, 
                     mymet.colnames.awc, mymet.colnames.gcb)

  mymet.datacol[!mymet.datacol %in% colnames(mymet)]
  colnames(mymet)[!colnames(mymet) %in% c(mymet.datacol, mymet.colnames.ids)]

  # [1] "mymet.PLATEROW"      "mymet.PLATECOL"  -- not needed

  ##
  ## make some binary variables (ie postmortem indicator) into characters
  ##
  for (i in c("mymet.brain", "mymet.timepoint", "mymet.sex", "mymet.living")) mymet[[i]] <- as.character(mymet[[i]])
  

   mymet[mymet.brain==1, mymet.brain:="brain"]
   mymet[mymet.brain==0, mymet.brain:="blood"]

  mymet[mymet.timepoint==1, mymet.timepoint:="left"]
  mymet[mymet.timepoint==0, mymet.timepoint:="right"]
  
  mymet[mymet.sex==1, mymet.sex:="F"]
  mymet[mymet.sex==0, mymet.sex:="M"]

  # mymet[mymet.living=="Living", mymet.postmortem:=0]
  # mymet[mymet.living=="Postmortem", mymet.postmortem:=1]


  ##
  ## remove columns with no variation
  ##
  novar <- names(which( apply(mymet, 2, uniqueN) == 1 ))

  dim(mymet)
  #[1] 778 217

  #mymet <- mymet[, !colnames(mymet) %in% novar, with=FALSE]
  

  mymet.datacol <- mymet.datacol[!mymet.datacol %in% novar]
  
  mymet <- mymet[,c(mymet.colnames.ids,mymet.datacol),with=F]
  
  dim(mymet)
  #[1] 778 164

  mymet$mymet.postmortem <- "postmortem"
  mymet[mymet.living==0, mymet.postmortem:="postmortem"]
  mymet[mymet.living==1, mymet.postmortem:="living"]

  table(mymet$mymet.postmortem)
    # living postmortem 
    #    532        246 

  dim(mymet)
  #[1] 778 165 + platerow and platecol which were removed bc not needed = 167 -- checks out

  
  ##
  ## fix column classes for STAR variables (manual review showed numeric variables were characters)
  ##
  for(i in mymet.colnames.str){if (class(mymet[[i]]) == "character") mymet[[i]] <- as.numeric(mymet[[i]])}
  ##
  
  ## reclass variables
  ##
  mymet.colclasses <- sapply(mymet[,mymet.datacol,with=F], class)
  for(i in mymet.datacol){
      if (mymet.colclasses[i] == "character") mymet[[i]] <- factor(mymet[[i]])
      if (mymet.colclasses[i] == "logical") mymet[[i]] <- factor(mymet[[i]])
      if (mymet.colclasses[i] == "integer64") mymet[[i]] <- as.numeric(mymet[[i]])
      if (mymet.colclasses[i] == "integer") mymet[[i]] <- as.numeric(mymet[[i]])
  }
  mymet.colclasses <- sapply(mymet[,mymet.datacol,with=F], class)

# $mymet.extractiondate
# [1] "POSIXct" "POSIXt"

mymet$mymet.extractiondate <- factor(mymet$mymet.extractiondate)

# $mymet.collectiondate
# [1] "POSIXct" "POSIXt"

mymet[, mymet.collectiondate:=NULL] #lots of NAs

dim(mymet)
#[1] 778 164

  # saveRDS(mymet, file="/sc/arion/projects/psychgen2/lbp/scratch/scratch_RAPiD/LBP_AllBatches_11JAN2021_RAPiD/200708_A00734_0062_BH7YCKDSXY/lbp_batch1Remainder_RAPiD_MYMET-QcMetricsCompiled_MERGED_rough_103BrainSamples_11FEB2021.rds")

saveRDS(mymet, file = "/sc/arion/projects/psychgen/lbp/files/files_from_psychgen2/lbp_allBatches_QC/lbp_allBatches_RAPiD_MYMET_QCMetrics_Merged_Compiled_778samples_NoBRAIN107_NoBRAIN734_AddedBRAIN107_04JAN2022.RDS")


old <- readRDS("/sc/arion/projects/psychgen2/lbp/data/RAW/rna/bulk/fromSema4/CompiledData/lbp_allBatches_RAPiD_MYMET_QCMetrics_Merged_Compiled_777samples_NoBRAIN107_NoBRAIN734_5MAY2021.RDS")

###########################################################################################
#Check the dumb NAs

# df[complete.cases(df), ]
# df[rowSums(is.na(df))==0,]
# na.omit(df)

# mymet[rowSums(is.na(mymet))>0,]

col.has.na <- apply(mymet, 2, function(x){any(is.na(x))})
sum(col.has.na)
#15

na.cols <- which(col.has.na == TRUE)
na.cols <- names(na.cols)

na.cols
# [1] "DupMetrics.UNPAIRED_READS_EXAMINED"        "DupMetrics.READ_PAIRS_EXAMINED"            "DupMetrics.SECONDARY_OR_SUPPLEMENTARY_RDS" "DupMetrics.UNMAPPED_READS"                 "DupMetrics.UNPAIRED_READ_DUPLICATES"       "DupMetrics.READ_PAIR_DUPLICATES"           "DupMetrics.READ_PAIR_OPTICAL_DUPLICATES"   "DupMetrics.PERCENT_DUPLICATION"            "DupMetrics.ESTIMATED_LIBRARY_SIZE"        
# [10] "mymet.age"                                 "mymet.race"                                "mymet.ethnicity"                           "mymet.collectiontime"                      "mymet.rna_amount_ng"                       "mymet.preservative"                       
  

#Will address this in the QC script dated 5MAY2021

lbpcov <- mymet

fc <- dt

##################################################################################################

#colnames(lbpcov) <- gsub("[.]", "_", colnames(lbpcov))


test <- apply(lbpcov, 2, function(x) any(is.na(x))) 
test2 <- test[test==TRUE]

names(test2)
# [1] "mymet_age"            "mymet_race"           "mymet_ethnicity"     
# [4] "mymet_preservative"   "mymet_collectiondate" "mymet_collectiontime"
# [7] "mymet_rna_amount_ng"


NA_counts<-sapply(lbpcov, function(x) sum(is.na(x)))

NAcols <- NA_counts[! NA_counts == 0] #columns that have more than 0 NAs and the number of NAs in each column


length(NAcols <- NA_counts[! NA_counts == 0])
#[1] 15

length(NA_counts[ NA_counts == 0]) #columns that have 0 NAs
# 149

NA_counts[! NA_counts == 0]
#        DupMetrics.UNPAIRED_READS_EXAMINED 
#                                       409 
#            DupMetrics.READ_PAIRS_EXAMINED 
#                                       409 
# DupMetrics.SECONDARY_OR_SUPPLEMENTARY_RDS 
#                                       409 
#                 DupMetrics.UNMAPPED_READS 
#                                       409 
#       DupMetrics.UNPAIRED_READ_DUPLICATES 
#                                       409 
#           DupMetrics.READ_PAIR_DUPLICATES 
#                                       409 
#   DupMetrics.READ_PAIR_OPTICAL_DUPLICATES 
#                                       409 
#            DupMetrics.PERCENT_DUPLICATION 
#                                       409 
#         DupMetrics.ESTIMATED_LIBRARY_SIZE 
#                                       409 
#                                 mymet.age 
#                                        16 
#                                mymet.race 
#                                       148 
#                           mymet.ethnicity 
#                                       157 
#                      mymet.collectiontime 
#                                       397 
#                       mymet.rna_amount_ng 
#                                        14 
#                        mymet.preservative 
#                                       491 


# remove mymet_collectiondate, mymet_collectiontime, sema4sampleqc_comments (already removed above) -- this is what was done in the QC of the previous data

lbpcov2 <- lbpcov[, -c("DupMetrics.UNPAIRED_READS_EXAMINED", "DupMetrics.READ_PAIRS_EXAMINED", "DupMetrics.SECONDARY_OR_SUPPLEMENTARY_RDS", "DupMetrics.UNMAPPED_READS", "DupMetrics.UNPAIRED_READ_DUPLICATES", "DupMetrics.READ_PAIR_DUPLICATES", "DupMetrics.READ_PAIR_OPTICAL_DUPLICATES", "DupMetrics.PERCENT_DUPLICATION", "DupMetrics.ESTIMATED_LIBRARY_SIZE", "mymet.collectiontime", "mymet.preservative")]

lbpcov2$mymet.ethnicity[is.na(lbpcov2$mymet.ethnicity)] <- "Unknown"
lbpcov2$mymet.race[is.na(lbpcov2$mymet.race)] <- "Unknown" 


NA_counts<-sapply(lbpcov2, function(x) sum(is.na(x)))
NA_counts[! NA_counts == 0]

# mymet.age mymet.rna_amount_ng 
#      16                  14 

norna2 <- lbpcov2[is.na(lbpcov2$mymet.rna_amount_ng),]$SAMPLE.ISMMS

# b1cov[which(b1cov$ISM_SEMA4 %in% norna)]$mymet_rna_amount_ng

# any(apply(b1cov, 2, function(x) grepl("RNA_Later", x)))

# which(apply(b1cov, 2, function(x) grepl("ice", x)))
# #  [1] 15691 15693 15694 15701 15710 15719 15728 15734 15739 15745 15748 15762
# # [13] 15765 15771 15774 15785 15791

lbpcov2[, mymet.rna_amount_ng:=NULL] #this was done in previous qcs

#######################################################################

# noage <- lbpcov2[is.na(lbpcov2$mymet.age),]$ISM.SEMA4
# norna <- lbpcov2[is.na(lbpcov2$mymet.rna_amount_ng),]$ISM.SEMA4

#noage2 <- lbpcov2[is.na(lbpcov2$mymet.age),]$SAMPLE.ISMMS

noage2 <- lbpcov2[is.na(lbpcov2$mymet.age), c("SAMPLE.ISMMS", "mymet.age")]
#grep("89+", lbpcov2$mymet.age, value=T)


library(openxlsx)
master_table <- read.xlsx("/sc/arion/projects/psychgen2/lbp/files/lbp_batch1_reRAPiD_QC3/Bulk_RNA_Isolation_Mastertable_BRAINANDBLOOD.xlsx")

master_table <- as.data.table(master_table, keep.rownames = TRUE)

master_table[master_table$LBPSEMA4_ID %in% noage2$SAMPLE.ISMMS]$age 

# [1] "89+" "89+" "89+" "na"  "89+" "89+" "89+" "89+" "89+" "89+" "89+" "89+"
# [13] "89+" "na"  "89+" "89+"

noagemaster <- master_table[LBPSEMA4_ID %in% noage2$SAMPLE.ISMMS, c("LBPSEMA4_ID", "age")] 

#          LBPSEMA4_ID age
#  1: LBPSEMA4BRAIN437 89+
#  2: LBPSEMA4BRAIN272 89+
#  3: LBPSEMA4BRAIN665 89+
#  4: LBPSEMA4BRAIN308  na
#  5: LBPSEMA4BRAIN685 89+
#  6: LBPSEMA4BRAIN420 89+
#  7: LBPSEMA4BRAIN276 89+
#  8: LBPSEMA4BRAIN343 89+
#  9: LBPSEMA4BRAIN736 89+
# 10: LBPSEMA4BRAIN133 89+
# 11: LBPSEMA4BRAIN178 89+
# 12: LBPSEMA4BRAIN413 89+
# 13: LBPSEMA4BRAIN432 89+
# 14: LBPSEMA4BRAIN779  na
# 15: LBPSEMA4BRAIN730 89+
# 16: LBPSEMA4BRAIN699 89+


#nineties <- noagemaster[c(1:3, 5:13, 15:16), ]$LBPSEMA4_ID

nineties <- noagemaster[age=="89+"]$LBPSEMA4_ID

lbpcov2[SAMPLE.ISMMS %in% nineties]$mymet.age <- 90

noage3 <- lbpcov2[is.na(lbpcov2$mymet.age), c("SAMPLE.ISMMS", "mymet.age")]
 #SAMPLE.ISMMS mymet.age
# 1: LBPSEMA4BRAIN308        NA
# 2: LBPSEMA4BRAIN779        NA

#LBPSEMA4BRAIN779, LBPSEMA4BRAIN308 -- set to the mean age of Columbia, their brain bank
# setting missing quantitative values to the mean is a simple form of imputation that is defensible in most of these types of situations where its a small number of samples for a metadata variable

unique(master_table$bank)
#[1] "Sinai"    "Harvard"  "Columbia" "Miami" 

columbia <- master_table[bank == "Columbia", c("LBPSEMA4_ID", "age")] #colnames need to be in quotes or else it does weird stuff

columbia[age=="89+"]$age <- 90
columbia2 <- columbia[!age=="na", ] 
columbia2$age <- as.numeric(columbia2$age)

mean(columbia2$age)
#[1] 75.875; round to 76

lbpcov2[SAMPLE.ISMMS %in% noage3$SAMPLE.ISMMS]$mymet.age <- 76

#######################################################################

#lbpcov3 <- lbpcov2[, !c("SAMPLE_ISMMS", "mymet_iid", "mymet_preservative")]

lbpcov3 <- lbpcov2
colnames(lbpcov3) <- gsub("[.]", "_", colnames(lbpcov3))
#######################################################################
#Remove skin samples -- pretty sure it wasn't actually a skin sample tho

lbpcov3[mymet_tissue=="Skin"]$SAMPLE_ISMMS
#[1] "LBPSEMA4BLOOD657"

# Left side and right side -- keep everyone in for now, look at if it's driving the signal -- if it is driving the variance, include it as cov bc we're not interested in left vs right and are accounting for individual anyway 

table(lbpcov3$mymet_tissue)

# L_Blood L_Brain R_Blood R_Brain    Skin 
#     128     284     115     249       1 

#Recategorize the "skin" sample to a blood sample -- R or L? It is R_Blood because the patient PT-0029 already has L_Brain, R_Brain, and L_Blood


lbpcov3[SAMPLE_ISMMS=="LBPSEMA4BLOOD657"]$mymet_tissue <- "R_Blood"

lbpcov3$mymet_tissue <- droplevels(lbpcov3$mymet_tissue)

table(lbpcov3$mymet_tissue)

# L_Blood L_Brain R_Blood R_Brain 
#     128     285     116     249 


test <- apply(lbpcov3, 2, function(x) any(is.na(x))) 
test2 <- test[test==TRUE]
#named logical(0)

lbpcov4 <- lbpcov3
###################################################################################################################################
#The sequencing batch situation
###################################################################################################################################

# Sequencing batch as a covariate now because we have a lot of batches 
# Think about the relationship between the batches and what we called sequencing batch before 

# Probably get rid of depletion batch because we don't know what they actually did; for us it's sequencing batch/3 but we don't know for them


# Make a variable that should be the batch as the merged batch that the data lives in -- did a batch data delivery reflect a batch that they sequenced; the merged batch folders represent sequencing batches in a way -- working assumption; this will give me trouble for the instances where there's only sample so maybe create an artificial category of other -- put them all in other together?

#Create a new variable with merged sequencing batch

grep("seq", colnames(lbpcov4), value=T)

table(lbpcov4$mymet_seqbatch)

# sb1 sb2 sb3 sb4 sb5 sb6 sb7 sb8 sb9 
#  95  96  96  94  95  95  95  96  16

# dirs <- read.table("/sc/arion/projects/psychgen2/lbp/data/RAW/rna/bulk/fromSema4/samplelists/FolderList_allBatches_Compile_4MAY2021.txt", header = F)


oldm1 <- readRDS(file = "/sc/arion/projects/psychgen2/lbp/data/RAW/rna/bulk/fromSema4/samplelists/m1Table_List_allBatches_4MAY2021.RDS")

oldm1$rapidDir3 <- basename(oldm1$rapidDir2)
oldm2 <- oldm1[SAMPLE_NAME  %in%  old$SAMPLE.ISMMS, ] 


m1 <- fread("/sc/arion/projects/psychgen/lbp/files/m1TableForLiharska2021_updated30JUNE2021.tsv")


m2 <- m1[, c("SAMPLE_NAME", "fqDirForLEL2021")]
m2$rapidDir3 <- dirname(m2$fqDirForLEL2021)
m2$rapidDir4 <- basename(m2$rapidDir3)



#m2 <- m2[! SAMPLE_NAME  %in% c("LBPSEMA4BRAIN107", "LBPSEMA4BRAIN734"), ] 

#remove the bad brain sample BRAIN734

m2 <- m2[SAMPLE_NAME  %in%  lbpcov4$SAMPLE_ISMMS, ] 


dim(m2)
#[1] 778   3
dim(m1)
#[1] 779  28

##################################

table(m2$rapidDir4)

#                                      190703_A00732_0026_AHML5TDSXX 
#                                                                  1 
#                                      190705_A00732_0027_BHMK7LDSXX 
#                                                                 16 
# 200529_A00732_0072_BH7YKMDSXY.200529_A00732_0073_AH7YHYDSXY.Merged 
#                                                                190 
#                                      200708_A00732_0077_AH7WMHDSXY 
#                                                                  2 
# 200708_A00732_0077_AH7WMHDSXY.200708_A00732_0078_BH7YHKDSXY.Merged 
#                                                                172 
#                                      200708_A00734_0062_BH7YCKDSXY 
#                                                                 25 
#                        200708_A00734_0062_BH7YCKDSXY.NPP536.Merged 
#                                                                155 
# 201112_A00734_0092_AHLJH3DSXY.201112_A00734_0093_BHLGNFDSXY.Merged 
#                                                                190 
#                                      210129_A00408_0556_BHMCVWDSXY 
#                                                                 14 
#                                      210205_A00734_0107_AHMFTTDSXY 
#                                                                  1 
#                                      210324_A01115_0091_AHTWMWDRXX 
#                                                                  5 
#                                                   LBPSEMA4BRAIN107 
#                                                                  1 
#                        NPP536.210205_A00734_0107_AHMFTTDSXY.Merged 
#                                                                  5 
#                        NPP536.210324_A01115_0091_AHTWMWDRXX.Merged 
#                                                                  1 

sum(table(m2$rapidDir4))
#[1] 778

#Merge all of the dirs that have <10 into their own "combinedbatch" factor level

#which batch was LBPSEMA4BRAIN107 actually in?
m1[SAMPLE_NAME=="LBPSEMA4BRAIN107", ]
#samplelist_201112_A00734_0093_BHLGNFDSXY and samplelist_201112_A00734_0092_AHLJH3DSXY so this actually should have been part of 201112_A00734_0092_AHLJH3DSXY.201112_A00734_0093_BHLGNFDSXY.Merged but was run separately so it should NOT go into otherbatch


m3 <- m2
m3$s4newbatch <- m3$rapidDir4

m3[m3$rapidDir4 == "NPP536.210324_A01115_0091_AHTWMWDRXX.Merged"]$s4newbatch <- "otherbatch"
m3[m3$rapidDir4 == "NPP536.210205_A00734_0107_AHMFTTDSXY.Merged"]$s4newbatch <- "otherbatch"
m3[m3$rapidDir4 == "210324_A01115_0091_AHTWMWDRXX"]$s4newbatch <- "otherbatch"
m3[m3$rapidDir4 == "210205_A00734_0107_AHMFTTDSXY"]$s4newbatch <- "otherbatch"
m3[m3$rapidDir4 == "200708_A00732_0077_AH7WMHDSXY"]$s4newbatch <- "otherbatch"
m3[m3$rapidDir4 == "190703_A00732_0026_AHML5TDSXX"]$s4newbatch <- "otherbatch"


############################
m3$s4newbatch <- factor(m3$s4newbatch)
m3$rapidDir4 <- factor(m3$rapidDir4)

colnames(m3)[1] <- "SAMPLE_ISMMS"

m3[SAMPLE_ISMMS=="LBPSEMA4BRAIN107"]$rapidDir4 <- "201112_A00734_0092_AHLJH3DSXY.201112_A00734_0093_BHLGNFDSXY.Merged"

m3[SAMPLE_ISMMS=="LBPSEMA4BRAIN107"]$s4newbatch <- "201112_A00734_0092_AHLJH3DSXY.201112_A00734_0093_BHLGNFDSXY.Merged"


m3[SAMPLE_ISMMS=="LBPSEMA4BRAIN107", ]


m3$rapidDir3 <- NULL

colnames(m3)[3] <- "rapidBatch"
colnames(m3)[2] <- "fastqDir"

lbpcov5 <- merge(lbpcov4, m3, by = "SAMPLE_ISMMS")

dim(lbpcov5)
#[1] 778 155


# What is the relationship between sequencing batch and the things that they gave us that we merged into these other things -- maybe make a heatmap with sequencing batch as we defined it by the merged batch things to see if the sequencing was generally done before we batched it; need to determine whether to remove our sequencing batch var bc it is being replaced by the batches by which it was delivered

batches <- lbpcov5[, c("SAMPLE_ISMMS", "mymet_seqbatch", "s4newbatch", "rapidBatch", "mymet_depletionbatch")]

library(variancePartition)
form=as.formula(paste("~",paste(colnames(batches)[colnames(batches) != "SAMPLE_ISMMS"],collapse="+")))
C = canCorPairs(form,batches[,colnames(batches)!= "SAMPLE_ISMMS", with=F])


#                      mymet_seqbatch s4newbatch rapidBatch mymet_depletionbatch
# mymet_seqbatch            1.0000000  0.7620066  0.7196162            1.0000000
# s4newbatch                0.7620066  1.0000000  1.0000000            0.7832940
# rapidBatch                0.7196162  1.0000000  1.0000000            0.6117082
# mymet_depletionbatch      1.0000000  0.7832940  0.6117082            1.0000000



#old, brain only
#                        mymet_seqbatch s4newbatch rapidDir3 mymet_depletionbatch
# mymet_seqbatch            1.0000000  0.6250523 0.5623651            1.0000000
# s4newbatch                0.6250523  1.0000000 1.0000000            0.7026412
# rapidDir3                 0.5623651  1.0000000 1.0000000            0.5019659
# mymet_depletionbatch      1.0000000  0.7026412 0.5019659            1.0000000


#######################################################################
#featureCounts Checks

# Look at the featureCounts to do some checks -- a lot of the qc metrics are derived from the fc matrix or are approximations of it 

# Sum of each individual's total reads in the fc matrix and plot the distributions -- expectation is that there will be some crazy outliers; 
# Manual sums and find the equivalent of that variable in the qc data

# For all of the quantitative traits in the matrix, make a histogram pdf with a lot of plots 
fc.copy <- fc

grep("_", colnames(fc), value=T)

#  [1] "LBPSEMA4BLOOD082_PT-0057"       "LBPSEMA4BLOOD116_PT-0152"      
#  [3] "LBPSEMA4BLOOD283_PT-0109"       "LBPSEMA4BLOOD387_PT-0114"      
#  [5] "LBPSEMA4BLOOD552_PT-0086"       "LBPSEMA4BRAIN050_PT-0082"      
#  [7] "LBPSEMA4BRAIN152_S14518"        "LBPSEMA4BRAIN192_PT-0145"      
#  [9] "LBPSEMA4BRAIN388_PT-0178"       "LBPSEMA4BRAIN393_S01745"       
# [11] "LBPSEMA4BRAIN410_PT-0159"       "LBPSEMA4BRAIN435_PT-0160"      
# [13] "LBPSEMA4BRAIN450_S02474"        "LBPSEMA4BRAIN489_S19433"       
# [15] "LBPSEMA4BRAIN511_PT-0033"       "LBPSEMA4BRAIN759_S12545"       
# [17] "LBPSEMA4BRAIN747_PT-0123"       "LBPSEMA4BRAIN107.no_filter.bam"

colnames(fc) <- gsub("_", ".", colnames(fc))

grep("[.]", colnames(fc), value=T)
#  [1] "LBPSEMA4BLOOD082.PT-0057"       "LBPSEMA4BLOOD116.PT-0152"      
#  [3] "LBPSEMA4BLOOD283.PT-0109"       "LBPSEMA4BLOOD387.PT-0114"      
#  [5] "LBPSEMA4BLOOD552.PT-0086"       "LBPSEMA4BRAIN050.PT-0082"      
#  [7] "LBPSEMA4BRAIN152.S14518"        "LBPSEMA4BRAIN192.PT-0145"      
#  [9] "LBPSEMA4BRAIN388.PT-0178"       "LBPSEMA4BRAIN393.S01745"       
# [11] "LBPSEMA4BRAIN410.PT-0159"       "LBPSEMA4BRAIN435.PT-0160"      
# [13] "LBPSEMA4BRAIN450.S02474"        "LBPSEMA4BRAIN489.S19433"       
# [15] "LBPSEMA4BRAIN511.PT-0033"       "LBPSEMA4BRAIN759.S12545"       
# [17] "LBPSEMA4BRAIN747.PT-0123"       "LBPSEMA4BRAIN107.no.filter.bam"

colnames(fc) <- gsub("\\..*","", colnames(fc))

grep("[.]", colnames(fc), value=T)
#character(0)

grep("_", colnames(fc), value=T)
#character(0)


fc <- data.frame(fc)
rownames(fc) <- fc$Geneid

#fc2 <- fc
#rownames(fc2) <- fc2$Geneid
#fc2 <- fc2[, c(7:783)]

#dim(fc2)
#[1] 58929   784

dim(fc)
#[1] 58929   784

fc2 <- fc[, colnames(fc) %in% lbpcov5$SAMPLE_ISMMS]

dim(fc2)
#[1] 58929   778

fc_counts <- apply(fc2, 2, sum)

summary(fc_counts)

  # Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
  # 4034  30387783  39968210  43049575  50151810 551668284 

#old, brain only
     # Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
     # 4034  32122563  41593469  44289522  51338159 551668284 

counts <- as.data.table(fc_counts, keep.rownames=TRUE)

counts[fc_counts==4034]$rn
#[1] "LBPSEMA4BRAIN722" -- REMOVE THIS SAMPLE



lbpcov5[SAMPLE_ISMMS == "LBPSEMA4BRAIN722", FEATURECOUNTS_Assigned]
#[1] 3648 #this sample 

lbpcov5[SAMPLE_ISMMS== "LBPSEMA4BRAIN341", FEATURECOUNTS_Assigned]
# [1] 47773234

m1[SAMPLE_NAME == "LBPSEMA4BRAIN722", ]
#/sc/arion/projects/psychgen2/lbp/data/RAW/rna/bulk/fromSema4/200708_A00734_0062_BH7YCKDSXY/



summary(lbpcov5$FEATURECOUNTS_Assigned)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 3648  27448031  36487874  39445572  46360211 509197771 

#old, brain only
     # Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
     # 3648  29551759  38647790  41064411  47972111 509197771 


lbpcov6 <- lbpcov5
lbpcov6 <- lbpcov6[match(colnames(fc2), lbpcov6$SAMPLE_ISMMS), ]


identical(lbpcov6$SAMPLE_ISMMS, colnames(fc2))
# [1] TRUE

identical(lbpcov6$SAMPLE_ISMMS, counts$rn)
#[1] TRUE

count.cor <- cbind(counts, lbpcov6$FEATURECOUNTS_Assigned) 


cor.test(count.cor$fc_counts, count.cor$V2, method = "spearman")

#         Spearman's rank correlation rho

# data:  count.cor$fc_counts and count.cor$V2
# S = 255192, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.9967485  

summary(fc_counts) - summary(lbpcov5$FEATURECOUNTS_Assigned)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 386  2939752  3480336  3604004  3791599 42470513 


#old, brain only
    # Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    #  386  2570804  2945679  3225110  3366048 42470513 

summary(fc_counts) - (summary(lbpcov5$FEATURECOUNTS_Assigned) + summary(lbpcov5$FEATURECOUNTS_Unassigned_Unmapped))


 # [13] "FEATURECOUNTS_Assigned"
 # [14] "FEATURECOUNTS_Unassigned_MultiMapping"
 # [15] "FEATURECOUNTS_Unassigned_NoFeatures"
 # [16] "FEATURECOUNTS_Unassigned_Unmapped"



lbpcov5[SAMPLE_ISMMS == "LBPSEMA4BLOOD042", ]
###########################
summary(fc_counts)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 4034  30387783  39968210  43049575  50151810 551668284 

#old, brain only
     # Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
     # 4034  32122563  41593469  44289522  51338159 551668284 


counts.melt <- melt(fc_counts)
rownames(counts.melt)[counts.melt$value==50151810]
#nothing, must be the mean of the two on each side

counts[fc_counts==50151810]$rn
#nothing

#this one was identified as the third q in the brain only 
counts[counts$rn=="LBPSEMA4BRAIN341",]$fc_counts
#[1] 51338159

lbpcov5[SAMPLE_ISMMS== "LBPSEMA4BRAIN341", FEATURECOUNTS_Assigned]
# [1] 47773234



/sc/arion/projects/psychgen2/lbp/data/RAW/rna/bulk/fromSema4/Merged_Batches/201112_A00734_0092_AHLJH3DSXY.201112_A00734_0093_BHLGNFDSXY.Merged

fctest <- read.table("/sc/arion/projects/psychgen2/lbp/data/RAW/rna/bulk/fromSema4/Merged_Batches/201112_A00734_0092_AHLJH3DSXY.201112_A00734_0093_BHLGNFDSXY.Merged/LBPSEMA4BRAIN341/RAPiD/featureCounts/LBPSEMA4BRAIN341.primary.txt", header = T)
fctest2 <- fctest$LBPSEMA4BRAIN341.no_filter.bam
sum(fctest2)
#[1] 51338159


#first quartile of the old brain only
m1[SAMPLE_NAME == counts[fc_counts==32122563]$rn, ]
#"LBPSEMA4BRAIN533"


fctest <- read.table("/sc/arion/projects/psychgen2/lbp/data/RAW/rna/bulk/fromSema4/Merged_Batches/200708_A00732_0077_AH7WMHDSXY.200708_A00732_0078_BH7YHKDSXY.Merged/LBPSEMA4BRAIN533/RAPiD/featureCounts/LBPSEMA4BRAIN533.primary.txt", header = T)
fctest3 <- fctest$LBPSEMA4BRAIN533.no_filter.bam
sum(fctest3)
#[1] 32122563

#Confirmed that the featureCounts have been merged correctly and the numbers are correct
#################################################################################
#REMOVE THE NEWLY DISCOVERED BAD SAMPLE FROM EVERYTHING "LBPSEMA4BRAIN722"

lbpcov5 <- lbpcov5[!SAMPLE_ISMMS == "LBPSEMA4BRAIN722", ]
dim(lbpcov5)
#[1] 777 155

table(lbpcov5$mymet_brain)

# blood brain 
#   244   533

names(which( apply(lbpcov5, 2, uniqueN) == 1 ))
#character(0)
#for the brain only it was [1] "FASTQC_Per_tile_sequence_quality"

#lbpcov5$FASTQC_Per_tile_sequence_quality <- NULL

table(lbpcov5$FASTQC_Per_tile_sequence_quality)
# FAIL PASS WARN 
#    1  775    1

lbpcov5[FASTQC_Per_tile_sequence_quality=="FAIL"]$SAMPLE_ISMMS
#[1] "LBPSEMA4BLOOD582"

counts[counts$rn=="LBPSEMA4BLOOD582",]$fc_counts
#[1] 4757 -- checks out, this is one of the ones identified in M1 as having very low coverage; remove in addition to BRAIN722, which was removed above

lbpcov5 <- lbpcov5[!SAMPLE_ISMMS == "LBPSEMA4BLOOD582", ]

dim(lbpcov5)
#[1] 776 155


##############################################################################
#the 25 samples run by eric separately -- are they in the correct batch folder?
##############################################################################

m3[SAMPLE_ISMMS=="LBPSEMA4BLOOD042",]
eric25 <- m3[rapidBatch=="200708_A00734_0062_BH7YCKDSXY"]$SAMPLE_ISMMS

# [1] "LBPSEMA4BLOOD042" "LBPSEMA4BLOOD061" "LBPSEMA4BLOOD582" "LBPSEMA4BLOOD601" "LBPSEMA4BLOOD670"
#  [6] "LBPSEMA4BLOOD744" "LBPSEMA4BLOOD795" "LBPSEMA4BRAIN007" "LBPSEMA4BRAIN099" "LBPSEMA4BRAIN115"
# [11] "LBPSEMA4BRAIN199" "LBPSEMA4BRAIN259" "LBPSEMA4BRAIN314" "LBPSEMA4BRAIN318" "LBPSEMA4BRAIN390"
# [16] "LBPSEMA4BRAIN415" "LBPSEMA4BRAIN502" "LBPSEMA4BRAIN524" "LBPSEMA4BRAIN561" "LBPSEMA4BRAIN565"
# [21] "LBPSEMA4BRAIN615" "LBPSEMA4BRAIN722" "LBPSEMA4BRAIN727" "LBPSEMA4BRAIN729" "LBPSEMA4BRAIN733"

m1[SAMPLE_NAME %in% eric25, .(SAMPLE_NAME, fqDirForLEL2021, DataProcessingStatus)]

samples <- m1[samplelist_200708_A00734_0062_BH7YCKDSXY==1]$SAMPLE_NAME
length(samples)
#[1] 180 -- this is this folder plus any that were in this folder plus another one and thus and merged

length(which(eric25 %in% samples))
#[1] 25

nrow(m3[rapidBatch=="200708_A00734_0062_BH7YCKDSXY", ])
#[1] 25 -- ok so there are only 25 samples that are just in this folder and not in any others

length(which(m3[rapidBatch=="200708_A00734_0062_BH7YCKDSXY"]$SAMPLE_ISMMS %in% eric25))
#[1] 25 -- ok so these are the 25 samples that are only in this folder

#so yes, the 25 eric samples are present and in the correct folder

length(which(eric25 %in% lbpcov5$SAMPLE_ISMMS))
#[1] 23 -- checks out, two eric25 samples removed for low coverage 

eric25[!eric25 %in% lbpcov5$SAMPLE_ISMMS]
#[1] "LBPSEMA4BLOOD582" "LBPSEMA4BRAIN722"

##############################################################################
#Next issue -- LBPSEMA4BRAIN004, per M1: fastq_in_two_deliveries|EV_RAPiD_run_on_one_fastq_delivery_because_other_corrupt
##############################################################################
summary(fc_counts)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 4034  30387783  39968210  43049575  50151810 551668284 

counts[counts$rn=="LBPSEMA4BRAIN004",]$fc_counts
#[1] 24193426

length(which(counts$fc_counts<24193426))
#[1] 122 with counts less than this one -- keep this sample, it is fine

length(which(counts$fc_counts<30387783))
#[1] 195 samples are below the first quartile, and this is one of them 

firstq <- counts[fc_counts<30387783,]

summary(firstq$fc_counts)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 4034 16475394 21927865 20087646 26700446 30380694 


#Decision: DO NOT REMOVE THIS SAMPLE, it has sufficient coverage despite the other delivery being corrupt

##############################################################################

#fc3 <- fc2[, colnames(fc2) %in% lbpcov5$SAMPLE_ISMMS, with=FALSE]
fc3 <- fc2[, colnames(fc2) %in% lbpcov5$SAMPLE_ISMMS,]

dim(fc3)
#[1] 58929   776


saveRDS(lbpcov5, file="/sc/arion/projects/psychgen/lbp/files/files_from_psychgen2/lbp_allBatches_QC/lbp_allBatches_RAPiD_MYMET_QCMetrics_Merged_Compiled_BLOODandBRAIN_776Samples_NoBLOOD582_NoBRAIN734_NoBRAIN722_05JAN2022.RDS")

#saveRDS(lbpcov5, "/sc/arion/projects/psychgen2/lbp/data/RAW/rna/bulk/fromSema4/lbp_allBatches_RAPiD_MYMET_QCMetrics_Merged_Compiled_onlyBRAIN_532Samples_NoBRAIN107_NoBRAIN734_NoBRAIN722_7MAY2021.RDS")

# fc <- readRDS("/sc/arion/projects/psychgen2/lbp/data/RAW/rna/bulk/fromSema4/200708_A00734_0062_BH7YCKDSXY/lbp_batch-200708_A00734_0062_BH7YCKDSXY_WhichIncludes-batch1Remaining_RAPiD_featureCounts_Compiled_Ordered_11FEB2021.RDS")

fc3$Geneid <- rownames(fc3)

fc3 <- fc3[, c(777, 1:776)]

saveRDS(fc3, "/sc/arion/projects/psychgen/lbp/files/files_from_psychgen2/lbp_allBatches_QC/lbp_allBatches_RAPiD_featureCounts_Compiled_BLOODandBRAIN_776samples_NoBLOOD582_NoBRAIN734_NoBRAIN722_05JAN2022.RDS")


#saveRDS(fc3, "/sc/arion/projects/psychgen2/lbp/data/RAW/rna/bulk/fromSema4/lbp_allBatches_RAPiD_featureCounts_Compiled_onlyBRAIN_532samples_NoBRAIN107_NoBRAIN734_NoBRAIN722_7MAY2021.RDS")
