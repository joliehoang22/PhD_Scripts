##check to see if this file already exclude LBP patients who do not consent to data sharing

## List of people who did not consent: RedCap > Excluded/Withdrawn Participants OR > Data Sharing Reports (post-QC)

/sc/arion/projects/psychgen/lbp/data/runEQTL/20OCT2022/DNA/LIVPM_dna.vcf.gz

/sc/arion/projects/mscic1/results/jolie/LBP/eqtl/LIV_only.vcf.gz


library(VariantAnnotation)

vcf <- readVcf(
  "/sc/arion/projects/psychgen/lbp/data/runEQTL/20OCT2022/DNA/LIVPM_dna.vcf.gz",
  genome = "hg38"   # or hg19 if applicable
)

rowRanges(vcf)   # variants
colData(vcf)     # samples
geno(vcf)        # genotypes
info(vcf)        # INFO fields


library(data.table)

vcf_dt <- fread(
  cmd = "zcat /sc/arion/projects/mscic1/results/jolie/LBP/eqtl/LIV_only.vcf.gz",
  skip = "#CHROM"
)

colnames(vcf_dt)
pt_cols <- grep("^PT-", colnames(vcf_dt), value = TRUE)
length(pt_cols) #160

pt_num <- as.integer(sub("PT-", "", pt_cols))
pt_cols_sorted <- pt_cols[order(pt_num)]
pt_cols_sorted

vcf_dt <- fread(
  cmd = "zcat /sc/arion/projects/psychgen/lbp/data/runEQTL/20OCT2022/DNA/LIVPM_dna.vcf.gz",
  skip = "#CHROM"
)

colnames(vcf_dt)
pt_cols <- grep("^PT-", colnames(vcf_dt), value = TRUE)
length(pt_cols) #160

pt_num <- as.integer(sub("PT-", "", pt_cols))
pt_cols_sorted <- pt_cols[order(pt_num)]
pt_cols_sorted
