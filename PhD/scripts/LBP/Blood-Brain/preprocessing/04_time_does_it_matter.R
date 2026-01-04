library(data.table)
library(dplyr)
library(ggplot2)
library(effsize)
library(ggsignif)
library(tidyr)
library(stringr)
library(tibble)
library(limma)
library(reshape2)

#### FOR NEWER STUFF, SCROLL DOWN!
####################################
load("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/top_gene_pairs_all_genes.RData")

## extract top entries with gene names
get_top_entries <- function(mat, cuts_named) {
  m_abs <- abs(mat)
  out <- lapply(cuts_named, function(thr) {
    # indices where correlation ≥ threshold
    keep_idx <- which(m_abs >= thr, arr.ind = TRUE)
    data.frame(
      blood_gene = rownames(mat)[keep_idx[, 1]],
      brain_gene = colnames(mat)[keep_idx[, 2]],
      correlation = mat[keep_idx]
    )
  })
  names(out) <- names(cuts_named)
  out
}

topentries_t0    <- get_top_entries(cor_t0_only, cuts_t0)
topentries_t1    <- get_top_entries(cor_t1_only, cuts_t1)
topentries_cross <- get_top_entries(cor_cross, cuts_cross)

topentries_t0$99.95%   topentries_t0$99.995% 

##just looking at 99.95th percentile for now 
t0_top <- topentries_t0[["99.95%"]]
dim(t0_top)
#224870      3
head(t0_top)
          blood_gene         brain_gene correlation
1 ENSG00000149418.11 ENSG00000000003.14   0.4186721
2 ENSG00000154447.15 ENSG00000000003.14   0.4078119
3  ENSG00000240065.8 ENSG00000000003.14  -0.3980624
4  ENSG00000282420.1 ENSG00000000003.14  -0.4005615

## 1. How many rows where blood_gene == brain_gene
same_gene_rows <- sum(t0_top$blood_gene == t0_top$brain_gene, na.rm = TRUE); same_gene_rows
#322

## 2. How many unique blood genes
n_unique_blood <- length(unique(t0_top$blood_gene)); n_unique_blood
#20085

## 3. How many unique brain genes
n_unique_brain <- length(unique(t0_top$brain_gene)); n_unique_brain
#19747

t1_top <- topentries_t1[["99.95%"]]
dim(t1_top)
#224891      3
head(t1_top)
          blood_gene         brain_gene correlation
1 ENSG00000135842.17 ENSG00000000003.14   0.4455757
2 ENSG00000143226.13 ENSG00000000003.14   0.4257991
3 ENSG00000168228.15 ENSG00000000003.14   0.4138899
4 ENSG00000196663.16 ENSG00000000003.14   0.4163581

## 1. How many rows where blood_gene == brain_gene
same_gene_rows <- sum(t1_top$blood_gene == t1_top$brain_gene, na.rm = TRUE); same_gene_rows
# 310

## 2. How many unique blood genes
n_unique_blood <- length(unique(t1_top$blood_gene)); n_unique_blood
# 19447

## 3. How many unique brain genes
n_unique_brain <- length(unique(t1_top$brain_gene)); n_unique_brain
# 18910

#Which gene pairs are in T0 but not T1 and vice versa?
# Add a unique pair ID to each
t0_top$pair_id <- paste(t0_top$blood_gene, t0_top$brain_gene, sep = "_")
t1_top$pair_id <- paste(t1_top$blood_gene, t1_top$brain_gene, sep = "_")

# In T0 but not in T1
t0_only <- t0_top[!t0_top$pair_id %in% t1_top$pair_id, ]
dim(t0_only) #224186      4

# In T1 but not in T0
t1_only <- t1_top[!t1_top$pair_id %in% t0_top$pair_id, ]
dim(t1_only) #224207      4

overlap <- t0_top[t0_top$pair_id %in% t1_top$pair_id, ]
dim(overlap) #684   4

###
hypergeometric problem: drawing two independent random subsets of size 𝑘, k from a universe of 𝑁,
N pairs. The overlap 𝑋, X has X∼Hypergeometric(N,K=k,n=k).

N <- 21046L * 21356L
k <- floor(0.0005 * N)   # 224,729
K <- k

# P(X >= 684) where X ~ Hypergeometric(N, K, k)
pval <- phyper(q = 684 - 1, m = K, n = N - K, k = k, lower.tail = FALSE)

# Also report expectation & SD
mu <- k * (K / N)
sd <- sqrt(k * (K / N) * (1 - K / N) * ((N - k) / (N - 1)))

list(pval = pval, mean = mu, sd = sd)
# $pval
# [1] 2.184709e-291
# $mean
# [1] 112.3644
# $sd
# [1] 10.59491
###
set.seed(42)

N <- 21046L * 21356L                 # total pairs
k <- floor(0.0005 * N)               # top 0.05% count
K <- k                               # same size in both sets
x_obs <- 684

B <- 1e6
overlaps <- rhyper(B, m = K, n = N - K, k = k)

mc_p <- mean(overlaps >= x_obs)
mc_mean <- mean(overlaps)
mc_sd <- sd(overlaps)
quant <- quantile(overlaps, c(.5,.9,.95,.99,.999))

list(p = mc_p, mean = mc_mean, sd = mc_sd, quantiles = quant)
# $p
# [1] 0

# $mean
# [1] 112.3722

# $sd
# [1] 10.5971

# $quantiles
#   50%   90%   95%   99% 99.9% 
#   112   126   130   138   146 

#########################
#########################
out_dir <- "/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/swap"

cor_t0_only <- readRDS(file.path(out_dir, "spearman_same_timepoint_t0_only_20250923.rds"))
cor_t1_only <- readRDS(file.path(out_dir, "spearman_same_timepoint_t1_only_20250923.rds"))
cor_cross <- readRDS(file.path(out_dir, "spearman_cross_timepoint_20250923.rds"))

probs <- c(`99.995th` = 0.99995,`99.9995th` = 0.999995)

## Helper: thresholds on full (rectangular) matrix
get_thresholds <- function(mat, probs) {
  vals <- as.vector(abs(mat))
  quantile(vals, probs = probs, names = TRUE)
}

## Helper: create top-X% matrices (values < cutoff -> NA)
make_top_mats <- function(mat, cuts_named) {
  m_abs <- abs(mat)
  out <- lapply(names(cuts_named), function(nm) {
    thr <- cuts_named[[nm]]
    m <- m_abs
    m[m < thr] <- NA_real_
    m
  })
  names(out) <- names(cuts_named)
  out
}

## Apply to your three matrices
cuts_t0    <- get_thresholds(cor_t0_only, probs); 
cuts_t0
#   99.995%  99.9995% 
# 0.4536308 0.5079292 
cuts_t1    <- get_thresholds(cor_t1_only, probs);
cuts_t1
#  99.995%  99.9995% 
# 0.4730038 0.5282303
cuts_cross <- get_thresholds(cor_cross,   probs)
cuts_cross
#   99.995%  99.9995% 
# 0.3352357 0.3837040  

top_t0    <- make_top_mats(cor_t0_only, cuts_t0)     # list of matrices: $`99.95th`, $`99.995th`
top_t0$99.995%   top_t0$99.9995%
dim(top_t0$`99.995%`) # 21046 21356
dim(top_t0$`99.9995%`) #21046 21356

top_t1    <- make_top_mats(cor_t1_only,  cuts_t1)    # same structure
top_cross <- make_top_mats(cor_cross,    cuts_cross) # same structure

# ## If you also want the raw vectors >= cutoff (for violin/boxplots later):
# get_top_vectors <- function(mat, cuts_named) {
#   v <- as.vector(abs(mat))
#   lapply(cuts_named, function(thr) v[v >= thr])
# }

# topvec_t0    <- get_top_vectors(cor_t0_only, cuts_t0)       
# #topvec_t0$99.995%   topvec_t0$99.9995%
# topvec_t1    <- get_top_vectors(cor_t1_only, cuts_t1)
# topvec_cross <- get_top_vectors(cor_cross, cuts_cross)

# # How many correlations survived each cutoff
# length(topvec_t0[["99.995%"]]) # 22473
# length(topvec_t0[["99.9995%"]]) # 2250

# summary(topvec_t0[["99.995%"]])
# #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# #  0.4537  0.4605  0.4700  0.4785  0.4866  0.8534 

# summary(topvec_t1[["99.995%"]])
# #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# #  0.4730  0.4800  0.4898  0.4986  0.5061  0.8947 

# summary(topvec_cross[["99.995%"]])
# #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# #  0.3352  0.3408  0.3488  0.3599  0.3627  0.8529

get_top_entries <- function(mat, cuts_named) {
  m_abs <- abs(mat)
  out <- lapply(cuts_named, function(thr) {
    # indices where correlation ≥ threshold
    keep_idx <- which(m_abs >= thr, arr.ind = TRUE)
    data.frame(
      blood_gene = rownames(mat)[keep_idx[, 1]],
      brain_gene = colnames(mat)[keep_idx[, 2]],
      correlation = mat[keep_idx]
    )
  })
  names(out) <- names(cuts_named)
  out
}

topentries_t0    <- get_top_entries(cor_t0_only, cuts_t0)
topentries_t1    <- get_top_entries(cor_t1_only, cuts_t1)
topentries_cross <- get_top_entries(cor_cross, cuts_cross)

##just looking at 99.95th percentile for now 
t0_top <- topentries_t0[["99.9995%"]]
dim(t0_top)
#2250    3
head(t0_top)
          blood_gene         brain_gene correlation
1 ENSG00000167460.16 ENSG00000002822.15   0.5924349
2  ENSG00000253910.2 ENSG00000003147.19  -0.5228002
3 ENSG00000155008.14  ENSG00000004799.8  -0.5309762
4 ENSG00000136003.15 ENSG00000004864.13   0.5191596
5  ENSG00000273445.1 ENSG00000005194.15   0.5192830
6 ENSG00000167460.16 ENSG00000005206.17   0.5696656

## 1. How many rows where blood_gene == brain_gene
same_gene_rows <- sum(t0_top$blood_gene == t0_top$brain_gene, na.rm = TRUE); same_gene_rows
#322 | 161

## 2. How many unique blood genes
n_unique_blood <- length(unique(t0_top$blood_gene)); n_unique_blood
#20085 | 1494

## 3. How many unique brain genes
n_unique_brain <- length(unique(t0_top$brain_gene)); n_unique_brain
#19747 | 1250

t1_top <- topentries_t1[["99.9995%"]]
dim(t1_top)
#2250    3
head(t1_top)
          blood_gene         brain_gene correlation
1 ENSG00000072952.18 ENSG00000003989.17   0.5488091
2 ENSG00000132612.16 ENSG00000004866.20   0.5351413
3 ENSG00000101255.11 ENSG00000006194.10  -0.5433790
4 ENSG00000079931.15 ENSG00000006625.18  -0.5372701
5  ENSG00000227986.1 ENSG00000006756.16   0.5295878
6 ENSG00000125812.16 ENSG00000007129.18   0.5668271

##agnostic of rows -- 
length(intersect(t0_top$blood_gene, t1_top$blood_gene)) #202
length(intersect(t0_top$brain_gene, t1_top$brain_gene)) #202
##are these the most reliable genes?

length(intersect(t0_top$blood_gene, t0_top$brain_gene)) #214 - lower than all samples combined (441)
length(intersect(t1_top$blood_gene, t1_top$brain_gene)) #217 - lower than all samples combined (441)

##cross-timepoint
length(intersect(t0_top$blood_gene, t1_top$brain_gene)) #185
length(intersect(t0_top$brain_gene, t1_top$blood_gene)) #177

## 1. How many rows where blood_gene == brain_gene
same_gene_rows <- sum(t1_top$blood_gene == t1_top$brain_gene, na.rm = TRUE); same_gene_rows
# 310 | 170


## how many exact blood–brain gene pairs in T0 also appear in T1?
t0_pairs <- paste(t0_top$blood_gene, t0_top$brain_gene, sep = "__")
t1_pairs <- paste(t1_top$blood_gene, t1_top$brain_gene, sep = "__")

n_overlap <- intersect(t0_pairs, t1_pairs)
n_overlap #324

##324 pairs out of 2250 pairs
df_overlap <- t1_top[t1_pairs %in% n_overlap, ]

df_overlap$rank <- rank(-abs(df_overlap$correlation), ties.method = "first")

df_overlap$is_same <- factor(
  df_overlap$blood_gene == df_overlap$brain_gene,
  levels = c(FALSE, TRUE)
)

p <- ggplot(df_overlap, aes(x = rank, y = correlation, color = is_same)) +
  geom_segment(aes(xend = rank, y = 0, yend = correlation),
               linewidth = 0.3, alpha = 0.7) +
  geom_point(size = 1.2) +
  scale_color_manual(
    values = c("FALSE" = "grey50", "TRUE" = "#7B1FA2"),
    labels = c("Different gene", "Same gene"),
    name = NULL
  ) +
  labs(
    x = "Rank among top correlations (T1 reference)",
    y = "Spearman correlation",
    title = "T0 T1 overlapping blood brain gene pairs"
  ) +
  theme_classic() +
  theme(legend.position = "top")

ggsave("/hpc/users/hoangd02/www/plots/lbp/gene_pairs_in_both_t0andt1_lollipop.pdf", p, width = 6, height = 4)

###T0 as reference
##324 pairs out of 2250 pairs
df_overlap <- t0_top[t0_pairs %in% n_overlap, ]

df_overlap$rank <- rank(-abs(df_overlap$correlation), ties.method = "first")

df_overlap$is_same <- factor(
  df_overlap$blood_gene == df_overlap$brain_gene,
  levels = c(FALSE, TRUE)
)

p <- ggplot(df_overlap, aes(x = rank, y = correlation, color = is_same)) +
  geom_segment(aes(xend = rank, y = 0, yend = correlation),
               linewidth = 0.3, alpha = 0.7) +
  geom_point(size = 1.2) +
  scale_color_manual(
    values = c("FALSE" = "grey50", "TRUE" = "#7B1FA2"),
    labels = c("Different gene", "Same gene"),
    name = NULL
  ) +
  labs(
    x = "Rank among top correlations (T0 reference)",
    y = "Spearman correlation",
    title = "T0 T1 overlapping blood brain gene pairs"
  ) +
  theme_classic() +
  theme(legend.position = "top")

ggsave("/hpc/users/hoangd02/www/plots/lbp/gene_pairs_in_both_t0andt1_lollipop_t0_as_reference.pdf", p, width = 6, height = 4)


## 2. How many unique blood genes
n_unique_blood <- length(unique(t1_top$blood_gene)); n_unique_blood
# 19447 | 1261

## 3. How many unique brain genes
n_unique_brain <- length(unique(t1_top$brain_gene)); n_unique_brain
# 18910 | 1226

#Which gene pairs are in T0 but not T1 and vice versa?
# Add a unique pair ID to each
t0_top$pair_id <- paste(t0_top$blood_gene, t0_top$brain_gene, sep = "_")
t1_top$pair_id <- paste(t1_top$blood_gene, t1_top$brain_gene, sep = "_")

# In T0 but not in T1
t0_only <- t0_top[!t0_top$pair_id %in% t1_top$pair_id, ]
dim(t0_only) #224186      4 | 1926    4

# In T1 but not in T0
t1_only <- t1_top[!t1_top$pair_id %in% t0_top$pair_id, ]
dim(t1_only) #224207      4 | 1926    4

overlap <- t0_top[t0_top$pair_id %in% t1_top$pair_id, ]
dim(overlap) #684   4 | 324   4

###
hypergeometric problem: drawing two independent random subsets of size 𝑘, k from a universe of 𝑁,
N pairs. The overlap 𝑋, X has X∼Hypergeometric(N,K=k,n=k).

N <- 21046L * 21356L
k <- floor(0.0005 * N)   # 224,729
K <- k

# P(X >= 684) where X ~ Hypergeometric(N, K, k)
pval <- phyper(q = 684 - 1, m = K, n = N - K, k = k, lower.tail = FALSE)

# Also report expectation & SD
mu <- k * (K / N)
sd <- sqrt(k * (K / N) * (1 - K / N) * ((N - k) / (N - 1)))

list(pval = pval, mean = mu, sd = sd)
# $pval
# [1] 2.184709e-291
# $mean
# [1] 112.3644
# $sd
# [1] 10.59491
###monte carlo
N <- 21046L * 21356L             # 449,319,576
m <- 17533L                      # total same-gene pairs in population
k <- 2248L                       # size of "top 0.0005%" set
x_obs <- 324L                    # observed same-gene in the top set

# one-sided exact p-value: P[X >= x_obs], X ~ Hypergeo(N, m, k)
p_exact <- phyper(q = x_obs - 1L, m = m, n = N - m, k = k, lower.tail = FALSE)

# expectation & SD under null
mu  <- k * (m / N)
sd0 <- sqrt(k * (m / N) * (1 - m / N) * ((N - k) / (N - 1)))

list(p_exact = p_exact, expected_mean = mu, sd = sd0)

set.seed(42)
B <- 1e6
overlaps <- rhyper(B, m = m, n = N - m, k = k)

mc_p    <- mean(overlaps >= x_obs)
mc_mean <- mean(overlaps)
mc_sd   <- sd(overlaps)
quant   <- quantile(overlaps, c(.5, .9, .95, .99, .999))

list(mc_p = mc_p, mean = mc_mean, sd = mc_sd, quantiles = quant)



set.seed(42)

N <- 21046L * 21356L                 # total pairs
k <- floor(0.000005 * N)               # top 0.05% count
K <- k                               # same size in both sets
x_obs <- 324 #684

B <- 1e6
overlaps <- rhyper(B, m = K, n = N - K, k = k)

mc_p <- mean(overlaps >= x_obs)
mc_mean <- mean(overlaps)
mc_sd <- sd(overlaps)
quant <- quantile(overlaps, c(.5,.9,.95,.99,.999))

list(p = mc_p, mean = mc_mean, sd = mc_sd, quantiles = quant)
# $p
# [1] 0

# $mean
# [1] 112.3722

# $sd
# [1] 10.5971

# $quantiles
#   50%   90%   95%   99% 99.9% 
#   112   126   130   138   146 
