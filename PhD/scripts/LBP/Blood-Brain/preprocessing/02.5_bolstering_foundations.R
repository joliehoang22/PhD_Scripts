#######this code is a continuation of 02_building_foundations.R
##started this code to prepare for Alex's lab meeting Nov 13th 2025
library(data.table)
library(dplyr)
library(ggplot2)
library(variancePartition)
library(tidyr)
library(stringr)
library(tibble)
library(limma)
library(glmnet)
library(reshape2)
library(org.Hs.eg.db)
library(AnnotationDbi)

load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood-brain_sample_baseline_with_dictionary_20250815.RData")
#head(blood_only_metadata)
head(dictionary)
#v_blood$E and v_brain$E

### test for same gene enrichment with new null/formula 
######### THIS IS WITH FULL BRAIN #########
form5_no_residID_full <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/comparison_of_null_and_real_data_summaries_blood_form5_no_residID_full_brain.txt", data.table=FALSE)
form5_no_residID_full
#            stat observed_value empirical_p_value
# 9      p99.995     0.27210186             0.143
# 10    p99.9995     0.31611339             0.020 #this still translates to ~2,240 gene pairs 
# 11   p99.99995     0.59640267             0.000

##this null was computed based on absolute value-- would raw value makes a difference?

###code from LBP_blood_brain_QC_20250603.rmd
blood_form5_by_brain_cor <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form5_no_indivdualID_brain_full_spearman_cor_matrix.txt", data.table=FALSE)
row.names(blood_form5_by_brain_cor) <- blood_form5_by_brain_cor$V1
blood_form5_by_brain_cor$V1 <- NULL
dim(blood_form5_by_brain_cor) #21046 21356
# blood_form5_by_brain_cor <- abs(blood_form3_by_brain_cor)
range(as.matrix(blood_form5_by_brain_cor)) #-0.7925811  0.8780499
blood_form5_99.9995 <- quantile(as.matrix(blood_form5_by_brain_cor), probs = 0.999995); blood_form5_99.9995 #0.3069401

# Step 1: Compute the 99.9995th percentile
threshold <- quantile(as.matrix(blood_form5_by_brain_cor), probs = 0.999995)
# Step 2: Identify which gene pairs meet or exceed the threshold
high_corr_indices <- which(blood_form5_by_brain_cor >= threshold, arr.ind = TRUE)
# Step 3: Create a dataframe with the corresponding gene names and correlation values
top_gene_pairs <- data.frame(
  blood_gene = rownames(blood_form5_by_brain_cor)[high_corr_indices[, 1]],
  brain_gene = colnames(blood_form5_by_brain_cor)[high_corr_indices[, 2]],
  correlation = blood_form5_by_brain_cor[high_corr_indices]
)
dim(top_gene_pairs) #2248    3
# Step 4 (optional): sort by correlation strength
top_gene_pairs <- top_gene_pairs[order(-top_gene_pairs$correlation), ]
length(unique(top_gene_pairs$brain_gene)) #1103
length(unique(top_gene_pairs$blood_gene)) #1375

head(top_gene_pairs)
             blood_gene         brain_gene correlation
1309 ENSG00000226259.10 ENSG00000226259.10   0.8780499
1320  ENSG00000226752.9  ENSG00000226752.9   0.8543236
517  ENSG00000145736.14 ENSG00000145736.14   0.8382796
1545  ENSG00000241945.8  ENSG00000241945.8   0.8358639
1514  ENSG00000237541.3  ENSG00000237541.3   0.8304804
1989  ENSG00000274602.5  ENSG00000274602.5   0.8253097

head(sort(table(top_gene_pairs$blood_gene), decreasing = TRUE))
#  ENSG00000228903.7 ENSG00000112576.12 ENSG00000160439.16  ENSG00000168209.5 
#                 37                 14                 12                 12 
#  ENSG00000203663.4  ENSG00000228049.7 
#                 12                 11 

head(sort(table(top_gene_pairs$brain_gene), decreasing = TRUE))
#  ENSG00000230445.4  ENSG00000186326.3 ENSG00000100906.10 ENSG00000090104.12 
#                 59                 46                 32                 27 
# ENSG00000047634.15  ENSG00000277785.1 
#                 26                 24

top_gene_pairs$type <- ifelse(top_gene_pairs$blood_gene == top_gene_pairs$brain_gene,
                              "same genes", "different genes")

# Count how many gene pairs have the same blood and brain gene name
n_same <- sum(top_gene_pairs$blood_gene == top_gene_pairs$brain_gene, na.rm = TRUE)

# Total number of pairs
n_total <- nrow(top_gene_pairs)

# Number of different gene pairs
n_different <- n_total - n_same

# Print the result
# cat("Same gene pairs:", n_same, "\nDifferent gene pairs:", n_different, "\n")
# Same gene pairs: 441 
# Different gene pairs: 1807

##calculating odds ratio
# counts
n_same_top <- 441
n_diff_top <- 1807
n_same_not <- 17092
n_diff_not <- 449439036

# build table: rows = Pair type, cols = Set
tab <- matrix(
  c(n_same_top, n_same_not,   # Same row: [Top, NotTop]
    n_diff_top, n_diff_not),  # Different row: [Top, NotTop]
  nrow = 2, byrow = TRUE,
  dimnames = list(Pair = c("Same","Different"),
                  Set  = c("Top","NotTop"))
)
tab
#           Set
# Pair           Top    NotTop
# Same           441      17092
# Different     1807   449439036

# Fisher’s exact test (enrichment of Same in Top)
ft <- fisher.test(tab, alternative = "greater")
ft$p.value        # p-value
ft$estimate       # odds ratio
log(ft$estimate)  # ln(OR)

#######EXPLORE SOME PLOTTING AND GENE ONTOLOGY PLEASE
#######lollipop
library(dplyr)
library(ggplot2)
library(ggrepel)

df <- top_gene_pairs %>%
  mutate(is_same = blood_gene == brain_gene,
         rank = row_number(desc(correlation)))

dim(df) #2248    5

p <- ggplot(df, aes(x = rank, y = correlation, color = is_same)) +
  geom_segment(aes(xend = rank, y = 0, yend = correlation), linewidth = 0.3) +
  geom_point(size = 1.2) +
  scale_color_manual(values = c(`TRUE` = "#7B1FA2", `FALSE` = "grey50"),
                     labels = c("Different gene", "Same gene"),
                     name = NULL) +
  labs(x = "Rank among top correlations", y = "Spearman Correlation",
       title = "Top blood brain gene expression correlations") +
  theme_classic() +
  theme(legend.position = "top")

ggsave("/hpc/users/hoangd02/www/plots/lbp/top_correlated_gene_pairs_lollipop.pdf", p, width = 6, height = 4)

#### GENE ONTOLOGY
library(clusterProfiler)
library(org.Hs.eg.db)

head(top_gene_pairs)
             blood_gene         brain_gene correlation
1309 ENSG00000226259.10 ENSG00000226259.10   0.8780499
1320  ENSG00000226752.9  ENSG00000226752.9   0.8543236
517  ENSG00000145736.14 ENSG00000145736.14   0.8382796
1545  ENSG00000241945.8  ENSG00000241945.8   0.8358639
1514  ENSG00000237541.3  ENSG00000237541.3   0.8304804
1989  ENSG00000274602.5  ENSG00000274602.5   0.8253097

########different ways to run GO
blood_genes <- unique(top_gene_pairs$blood_gene)
brain_genes <- unique(top_gene_pairs$brain_gene)
same_genes  <- unique(top_gene_pairs$blood_gene[top_gene_pairs$blood_gene == top_gene_pairs$brain_gene])

# remove the dot and everything after it
blood_genes <- gsub("\\..*", "", unique(top_gene_pairs$blood_gene))
brain_genes <- gsub("\\..*", "", unique(top_gene_pairs$brain_gene))
same_genes  <- gsub("\\..*", "", 
                    unique(top_gene_pairs$blood_gene[top_gene_pairs$blood_gene == top_gene_pairs$brain_gene]))

#sanity check
length(blood_genes) #1375
length(brain_genes) #1103
length(same_genes) #441

#1. Blood genes only: 1,375 unique blood genes, answering:
#What biological processes in blood associate with strong cross-tissue correlations?
ego_blood <- enrichGO(gene          = blood_genes,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "ENSEMBL",
                      ont           = "ALL",          # BP = Biological Process
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.1) #nothing found at 0.05
ego_blood #1 enriched terms found

#2. Brain genes only: 1,103 unique brain genes, answering:
#What biological processes in brain associate with strong cross-tissue correlations?
ego_brain <- enrichGO(gene = brain_genes, OrgDb = org.Hs.eg.db,
                      keyType = "ENSEMBL", ont = "ALL", pAdjustMethod = "BH",
                      pvalueCutoff = 0.05) #7 enriched terms found at 0.1; 5 at 0.05
ego_brain

#3. Same genes: Which biological processes are shared across genes with highly correlated blood and brain expression?
ego_same <- enrichGO(gene = same_genes, OrgDb = org.Hs.eg.db,
                     keyType = "ENSEMBL", ont = "ALL", pAdjustMethod = "BH",
                     pvalueCutoff = 0.1)
ego_same #1 enriched terms found at 0.05 and 0.1

pdf("/hpc/users/hoangd02/www/plots/lbp/ego_blood_dotplot.pdf", width = 8, height = 6)
dotplot(ego_blood, showCategory = 15, title = "Blood Genes")
dev.off()

dotplot(ego_brain, showCategory = 15, title = "Brain Genes")
dotplot(ego_same,  showCategory = 15, title = "Shared Genes (Blood–Brain)")

#compare the 3:
compareClusterResult <- compareCluster(
  geneCluster = list("Blood" = blood_genes, "Brain" = brain_genes, "Both Tissues" = same_genes),
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "ALL", pAdjustMethod = "BH", pvalueCutoff= 0.1
)

pdf("/hpc/users/hoangd02/www/plots/lbp/GO_results_all_tissues.pdf", width = 8, height = 6)
dotplot(compareClusterResult, showCategory = 10, title = "GO enrichment across tissues")
dev.off()

#######what are the genes? non-coding? coding? 
library(biomaRt)

# connect to Ensembl human dataset
mart <- biomaRt::useEnsembl(biomart = "ensembl",
                            dataset = "hsapiens_gene_ensembl",
                            mirror = "asia")
# remove version suffixes
genes <- gsub("\\..*", "", unique(c(top_gene_pairs$blood_gene, top_gene_pairs$brain_gene)))

# query Ensembl for gene type
annot <- getBM(
  attributes = c("ensembl_gene_id", "gene_biotype", "hgnc_symbol", "description"),
  filters    = "ensembl_gene_id",
  values     = genes,
  mart       = mart
)

head(annot)

#########################if biomart doesnt work ~ this is less comprehensive
library(org.Hs.eg.db)
library(AnnotationDbi)

blood_genes <- gsub("\\..*", "", top_gene_pairs$blood_gene)
brain_genes <- gsub("\\..*", "", top_gene_pairs$brain_gene)

mapped_blood <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = blood_genes,
  columns = c("ENSEMBL", "SYMBOL", "GENETYPE"),
  keytype = "ENSEMBL"
)

head(mapped_blood)

mapped_brain <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = brain_genes,
  columns = c("ENSEMBL", "SYMBOL", "GENETYPE"),
  keytype = "ENSEMBL"
)

head(mapped_brain)

library(ggalluvial)

# Simplify to one row per gene type for each side
df <- full_join(
  mapped_blood %>% 
    count(type = GENETYPE, name = "blood_n"),
  mapped_brain %>% 
    count(type = GENETYPE, name = "brain_n"),
  by = "type"
) %>%
  mutate(blood_n = replace_na(blood_n, 0),
         brain_n = replace_na(brain_n, 0))

# Long format for ggalluvial
df_long <- df %>%
  pivot_longer(cols = c(blood_n, brain_n),
               names_to = "Dataset",
               values_to = "Count") %>%
  mutate(Dataset = recode(Dataset,
                          blood_n = "Blood",
                          brain_n = "Brain"))

# Example alluvial plot
p<- ggplot(df_long,
       aes(axis1 = Dataset, axis2 = type, y = Count)) +
  geom_alluvium(aes(fill = type), width = 1/8) +
  geom_stratum(width = 1/8, fill = "grey80", color = "grey40") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("Blood", "Brain"), expand = c(.05, .05)) +
  labs(title = "Gene biotype overlap between blood and brain",
       y = "Number of genes") +
  theme_minimal()

ggsave("/hpc/users/hoangd02/www/plots/lbp/top_correlated_gene_pairs_type_of_genes.pdf", p, width = 6, height = 4)

library(dplyr)
library(tidyr)
library(ggplot2)

df <- full_join(
  mapped_blood %>% count(type = GENETYPE, name = "Blood"),
  mapped_brain %>% count(type = GENETYPE, name = "Brain"),
  by = "type"
) %>% replace_na(list(Blood = 0, Brain = 0)) %>%
  pivot_longer(c(Blood, Brain), names_to = "Dataset", values_to = "Count")

p<- ggplot(df, aes(x = Dataset, y = Count, fill = type)) +
  geom_col() +
  labs(x = NULL, y = "Number of genes",
       title = "Gene biotype distribution in Blood vs Brain", fill = "Biotype") +
  theme_minimal()

ggsave("/hpc/users/hoangd02/www/plots/lbp/top_correlated_gene_pairs_type_of_genes.pdf", p, width = 6, height = 4)

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggalluvial)

pairs <- merge(
  mapped_blood[, c("ENSEMBL","GENETYPE")], 
  mapped_brain[, c("ENSEMBL","GENETYPE")],
  by = "ENSEMBL", all = TRUE
)
names(pairs) <- c("ENSEMBL","blood_type","brain_type")

pairs$blood_type[is.na(pairs$blood_type)] <- "NA"
pairs$brain_type[is.na(pairs$brain_type)] <- "NA"

pairs <- as.data.frame(table(pairs$blood_type, pairs$brain_type))
names(pairs) <- c("blood_type","brain_type","n")

p<- ggplot(pairs, aes(axis1 = blood_type, axis2 = brain_type, y = n)) +
  geom_alluvium(aes(fill = blood_type), width = 1/8) +
  geom_stratum(width = 1/8, fill = "grey85", color = "grey40") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("Blood", "Brain")) +
  labs(title = "Biotype flow: Blood to Brain", y = "Number of genes", x = NULL, fill = "Blood biotype") +
  theme_minimal()

ggsave("/hpc/users/hoangd02/www/plots/lbp/top_correlated_gene_pairs_type_of_genes_galluvial.pdf", p, width = 6, height = 8)


library(dplyr)
library(tidyr)
library(ggplot2)
library(ggalluvial)

## 1) Build Blood↔Brain biotype pairs ---------------------------------------
pairs <- full_join(
  mapped_blood %>% dplyr::select(ENSEMBL, blood_type = GENETYPE),
  mapped_brain %>% dplyr::select(ENSEMBL, brain_type = GENETYPE),
  by = "ENSEMBL"
) %>%
  mutate(
    blood_type = tidyr::replace_na(blood_type, "NA"),
    brain_type = tidyr::replace_na(brain_type, "NA")
  ) %>%
  count(blood_type, brain_type, name = "n")

## 2) Put NA at the bottom & set palette -------------------------------------
levs <- c("protein-coding","pseudo","snoRNA","snRNA","ncRNA","other","NA")
pairs <- pairs %>%
  mutate(
    blood_type = factor(blood_type, levels = levs),
    brain_type = factor(brain_type, levels = levs)
  )

pal <- c(
  "protein-coding" = "#C62828",
  "pseudo"         = "#6A3D9A",
  "snoRNA"         = "#1F78B4",
  "snRNA"          = "#33A02C",
  "ncRNA"          = "#FF7F00",
  "other"          = "#B2DF8A",
  "NA"             = "#D9D9D9"   # light grey
)

## 3) Plot: blocks colored by *their own* type, counts inside, neutral flows --
p<- ggplot(pairs, aes(axis1 = blood_type, axis2 = brain_type, y = n)) +
  # neutral flows
  geom_alluvium(fill = "grey80", color = "grey70", width = 0.15, alpha = 0.7) +
  # blocks (both axes) colored by their own stratum label
  geom_stratum(aes(fill = after_stat(stratum)), width = 0.15, color = "grey40") +
  # numbers inside blocks (sum of n per stratum)
  geom_text(stat = "stratum", aes(label = after_stat(n)), size = 3) +
  scale_fill_manual(values = pal, name = "Biotype") +
  scale_x_discrete(limits = c("Blood", "Brain"), expand = c(.03, .03)) +
  labs(title = "Biotype flow: Blood to Brain",
       x = NULL, y = "Number of genes") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.x = element_blank(),
        legend.position = "right")

ggsave("/hpc/users/hoangd02/www/plots/lbp/top_correlated_gene_pairs_type_of_genes_galluvial.pdf", p, width = 6, height = 8)

https://hoangd02.u.hpc.mssm.edu/plots/lbp/top_correlated_gene_pairs_type_of_genes_galluvial.pdf 
######################################## EXPLORING SOME PLOTS ########################################
######################################################################################################
head(top_gene_pairs)
             blood_gene         brain_gene correlation
1309 ENSG00000226259.10 ENSG00000226259.10   0.8780499
1320  ENSG00000226752.9  ENSG00000226752.9   0.8543236
517  ENSG00000145736.14 ENSG00000145736.14   0.8382796
1545  ENSG00000241945.8  ENSG00000241945.8   0.8358639
1514  ENSG00000237541.3  ENSG00000237541.3   0.8304804
1989  ENSG00000274602.5  ENSG00000274602.5   0.8253097

library(dplyr)
library(ggplot2)
library(forcats)

# Degree per blood gene
blood_deg <- top_gene_pairs %>%
  group_by(blood_gene) %>%
  summarise(
    n_pairs   = n(),
    mean_cor  = mean(correlation),
    max_cor   = max(correlation),
    .groups = "drop"
  ) %>%
  mutate(tissue = "blood",
         gene   = blood_gene)

# Degree per brain gene
brain_deg <- top_gene_pairs %>%
  group_by(brain_gene) %>%
  summarise(
    n_pairs   = n(),
    mean_cor  = mean(correlation),
    max_cor   = max(correlation),
    .groups = "drop"
  ) %>%
  mutate(tissue = "brain",
         gene   = brain_gene)

deg_df <- bind_rows(blood_deg, brain_deg)

# Pick top 20 hubs per tissue
top_deg <- deg_df %>%
  group_by(tissue) %>%
  slice_max(n_pairs, n = 20) %>%
  ungroup() %>%
  mutate(gene = fct_reorder(gene, n_pairs))

# Plot
pdf("/hpc/users/hoangd02/www/plots/lbp/top_corr_gene_pairs_degree.pdf",
    width = 8, height = 6)

ggplot(top_deg,
       aes(x = n_pairs, y = gene, color = tissue, size = mean_cor)) +
  geom_point() +
  facet_wrap(~ tissue, scales = "free_y") +
  scale_size_continuous(name = "Mean correlation") +
  labs(
    x = "# of connections",
    y = "Gene",
    title = "Top hub genes by number of blood brain connections"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    strip.background = element_rect(fill = "grey90"),
    axis.text.y = element_text(size = 6)
  )

dev.off()


##########################
library(dplyr)
library(ggplot2)
library(ggalluvial)

# 1) Degree per gene (full table)
deg_blood <- top_gene_pairs %>%
  count(blood_gene, name = "deg_blood")

deg_brain <- top_gene_pairs %>%
  count(brain_gene, name = "deg_brain")

# 2) Join degrees + correlation bins + same-gene flag
pairs_gene <- top_gene_pairs %>%
  left_join(deg_blood, by = "blood_gene") %>%
  left_join(deg_brain, by = "brain_gene") %>%
  mutate(
    same_gene = (blood_gene == brain_gene),
    cor_bin = cut(
      correlation,
      breaks = c(0.30, 0.33, 0.38, 0.45, 0.60, 0.88),
      include.lowest = TRUE,
      right = FALSE,
      labels = c(
        "[0.30-0.33)",
        "[0.33-0.38)",
        "[0.38-0.45)",
        "[0.45-0.60)",
        "[0.60-0.88]"
      )
    )
  )

# 3) Filter out singletons for readability
pairs_plot <- pairs_gene %>%
  filter(deg_blood >= 3 | deg_brain >= 3)

# Recompute degrees in the plotted data
deg_blood_plot <- pairs_plot %>%
  count(blood_gene, name = "deg_blood_plot")

deg_brain_plot <- pairs_plot %>%
  count(brain_gene, name = "deg_brain_plot")

# Levels: low -> high so high-degree genes end up at the TOP
blood_levels <- deg_blood_plot %>%
  arrange(deg_blood_plot) %>%
  pull(blood_gene)

brain_levels <- deg_brain_plot %>%
  arrange(deg_brain_plot) %>%
  pull(brain_gene)

# Create *axis-specific* labels so ggalluvial doesn't merge strata across axes
pairs_plot <- pairs_plot %>%
  mutate(
    blood_label = factor(paste0("Blood: ", blood_gene),
                         levels = paste0("Blood: ", blood_levels)),
    brain_label = factor(paste0("Brain: ", brain_gene),
                         levels = paste0("Brain: ", brain_levels))
  )

p <- ggplot(
  pairs_plot,
  aes(axis1 = blood_label, axis2 = brain_label, y = 1)
) +
  # flows: fill = correlation bin, dashed if same gene
  geom_alluvium(
    aes(fill = cor_bin, linetype = same_gene),
    color = NA,
    alpha = 0.9,
    width = 0.18
  ) +
  # strata: neutral blocks
  geom_stratum(
    width = 0.15,
    fill  = "grey90",
    color = "grey50"
  ) +
  # labels: only show counts when n >= 10
  geom_text(
    stat = "stratum",
    aes(label = ifelse(after_stat(n) >= 10, after_stat(n), "")),
    size = 2
  ) +
  scale_fill_brewer(
    palette = "YlOrRd",
    name    = "Correlation (bin)"
  ) +
  scale_linetype_manual(
    name   = "Same gene?",
    values = c(`FALSE` = "solid", `TRUE` = "dashed")
  ) +
  scale_x_discrete(
    limits = c("Blood gene", "Brain gene"),
    expand = c(.03, .03)
  ) +
  labs(
    title = "Blood Brain gene pairs (filtered):\nflows colored by correlation, same genes dashed",
    x = NULL,
    y = "Number of gene to gene pairs (each = 1)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "right"
  )

ggsave(
  "/hpc/users/hoangd02/www/plots/lbp/top_corr_gene_pairs_alluvial_sorted3.pdf",
  p,
  width = 7,
  height = 10
)

https://hoangd02.u.hpc.mssm.edu/plots/lbp/top_corr_gene_pairs_alluvial_sorted2.pdf

p <- ggplot(
  pairs_plot,
  aes(axis1 = blood_label, axis2 = brain_label, y = 1)
) +
  # flows: correlation color + dashed if same gene
  geom_alluvium(
    aes(fill = cor_bin, linetype = same_gene),
    color = "grey40",   # subtle outline so dashes are visible
    size  = 0.3,        # thin border, not bold
    alpha = 0.9,
    width = 0.18
  ) +
  # strata (blocks)
  geom_stratum(
    width = 0.15,
    fill  = "grey90",
    color = "grey50"
  ) +
  # labels: only when n >= 10
  geom_text(
    stat  = "stratum",
    aes(label = ifelse(after_stat(n) >= 10, after_stat(n), "")),
    size  = 2
  ) +
  scale_fill_brewer(
    palette = "YlOrRd",
    name    = "Correlation (bin)"
  ) +
  scale_linetype_manual(
    name   = "Same gene?",
    values = c(`FALSE` = "solid", `TRUE` = "dashed")
  ) +
  scale_x_discrete(
    limits = c("Blood gene", "Brain gene"),
    expand = c(.03, .03)
  ) +
  labs(
    title = "Blood Brain gene pairs (filtered):\nflows colored by correlation, same genes dashed",
    x = NULL,
    y = "Number of gene to gene pairs (each = 1)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    legend.position    = "right"
  )

ggsave(
  "/hpc/users/hoangd02/www/plots/lbp/top_corr_gene_pairs_alluvial_dashed_clean.pdf",
  p,
  width = 7,
  height = 12
)


############## ############## ############## ############## 
############## SOME CLUSTERING ############## #############
############## ############## ############## ############## 
blood_form5_by_brain_cor <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form5_no_indivdualID_brain_full_spearman_cor_matrix.txt", data.table=FALSE)
row.names(blood_form5_by_brain_cor) <- blood_form5_by_brain_cor$V1
blood_form5_by_brain_cor$V1 <- NULL
dim(blood_form5_by_brain_cor) #21046 21356
# blood_form5_by_brain_cor <- abs(blood_form3_by_brain_cor)
range(as.matrix(blood_form5_by_brain_cor)) #-0.7925811  0.8780499
blood_form5_99.9995 <- quantile(as.matrix(blood_form5_by_brain_cor), probs = 0.999995); blood_form5_99.9995 #0.3069401

library(igraph)

## 1) matrix + threshold
mat <- as.matrix(blood_form5_by_brain_cor)

thr <- quantile(mat, 0.999995, na.rm = TRUE); thr #0.3069401

## 2) subset to edges above threshold
idx <- which(mat >= thr, arr.ind = TRUE)

edges <- data.table(
  from = paste0("B_", rownames(mat)[idx[,1]]),
  to   = paste0("R_", colnames(mat)[idx[,2]]),
  w    = mat[idx]
)

## 3) graph + clustering
g <- graph_from_data_frame(edges, directed = FALSE)
cl <- cluster_louvain(g, weights = E(g)$w)

V(g)$cluster <- membership(cl)

nodes <- data.table(
  node = V(g)$name,
  cluster = V(g)$cluster,
  degree = degree(g)
)

png("/hpc/users/hoangd02/www/plots/lbp/blood_brain_graph.png",
    width = 1400, height = 1100, res = 150)

plot(induced_subgraph(g, names(sort(degree(g), decreasing = TRUE)[1:2000])),
     vertex.size = 2,
     vertex.label = NA)

dev.off()


#############non igraph ##############
library(data.table)
library(irlba)

mat <- as.matrix(blood_form5_by_brain_cor)

## 1) threshold + sparsify
thr <- quantile(mat, 0.999995, na.rm = TRUE)
mat_thr <- mat
mat_thr[mat_thr < thr] <- 0

## 2) low-rank structure
sv <- irlba(mat_thr, nv = 30)

blood_embed <- sv$u %*% diag(sv$d)
brain_embed <- sv$v %*% diag(sv$d)

## 3) cluster
set.seed(1)
k_blood <- 20
k_brain <- 20

blood_cl <- kmeans(blood_embed[, 1:20], centers = k_blood)$cluster
brain_cl <- kmeans(brain_embed[, 1:20], centers = k_brain)$cluster

## 4) save outputs (ABSOLUTE PATHS)
blood_dt <- data.table(gene = rownames(mat), cluster = blood_cl)
brain_dt <- data.table(gene = colnames(mat), cluster = brain_cl)

fwrite(blood_dt,
       "/hpc/users/hoangd02/www/plots/lbp/blood_clusters.tsv.gz",
       sep = "\t")
fwrite(brain_dt,
       "/hpc/users/hoangd02/www/plots/lbp/brain_clusters.tsv.gz",
       sep = "\t")

## 5) block map plot
idx <- which(mat >= thr, arr.ind = TRUE)

dt_edges <- data.table(
  b = blood_cl[idx[, 1]],
  r = brain_cl[idx[, 2]],
  w = mat[idx]
)

block <- dt_edges[, .(mean_w = mean(w), n = .N), by = .(b, r)]
block_mat <- dcast(block, b ~ r, value.var = "mean_w", fill = 0)

bm <- as.matrix(block_mat[, -1])
rownames(bm) <- paste0("B", block_mat$b)
colnames(bm) <- paste0("R", colnames(bm))

png("/hpc/users/hoangd02/www/plots/lbp/cluster_blockmap.png",
    width = 1400, height = 1100, res = 150)
heatmap(bm, Rowv = NA, Colv = NA, scale = "none",
        main = sprintf("Mean cross-tissue association (thr=%.4f)", thr))
dev.off()

