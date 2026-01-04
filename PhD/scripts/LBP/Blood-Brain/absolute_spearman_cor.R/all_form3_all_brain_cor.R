## code derived from LBP_blood_brain_QC_20250603.rmd
rm(list=ls())
library(data.table)
library(ggplot2)
library(readxl)
library(biomaRt)
library(dplyr)
library(edgeR)
library(limma)
library(variancePartition)
library(matrixStats)
library(purrr)
library(caret)
library(tidyr)

set.seed(2025)
##################################################
######### NO INDIVIDUAL | Brain baseline ######### 
##################################################
blood_form3_by_brain_cor <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form3_no_indivdualID_brain_baseline_spearman_cor_matrix.txt", data.table=FALSE)
row.names(blood_form3_by_brain_cor) <- blood_form3_by_brain_cor$V1; blood_form3_by_brain_cor$V1 <- NULL
blood_form3_by_brain_cor <- abs(blood_form3_by_brain_cor)
blood_form3_99.95 <- quantile(as.matrix(blood_form3_by_brain_cor), probs = 0.9995); blood_form3_99.95 #0.2483944
blood_form3_99.995 <- quantile(as.matrix(blood_form3_by_brain_cor), probs = 0.99995); blood_form3_99.995 #0.2906922

# Compute the 99.95th percentile
threshold <- quantile(as.matrix(blood_form3_by_brain_cor), probs = 0.9995)
# Identify which gene pairs meet or exceed the threshold
high_corr_indices <- which(blood_form3_by_brain_cor >= threshold, arr.ind = TRUE) #arr.ind = TRUE returns a matrix with which() otherwise it returns a vector of linear indices

# Create a dataframe with the corresponding gene names and correlation values
top_gene_pairs <- data.frame(
  blood_gene = rownames(blood_form3_by_brain_cor)[high_corr_indices[, 1]],
  brain_gene = colnames(blood_form3_by_brain_cor)[high_corr_indices[, 2]],
  correlation = blood_form3_by_brain_cor[high_corr_indices]
)
# Sort by correlation strength
top_gene_pairs <- top_gene_pairs[order(-top_gene_pairs$correlation), ]
length(unique(top_gene_pairs$brain_gene)) #15695
length(unique(top_gene_pairs$blood_gene)) #15366

head(top_gene_pairs)
#               blood_gene         brain_gene correlation
#167582 ENSG00000226259.10 ENSG00000226259.10   0.8344016
#181336  ENSG00000241945.8  ENSG00000241945.8   0.8261620
#210686  ENSG00000274602.5  ENSG00000274602.5   0.8222672
#208259 ENSG00000233327.10  ENSG00000273018.6   0.8079362
#208262  ENSG00000273018.6  ENSG00000273018.6   0.8060124
#175112 ENSG00000233327.10 ENSG00000233327.10   0.8029794

dim(top_gene_pairs) #224734      3

summary(top_gene_pairs$correlation)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2484  0.2537  0.2610  0.2700  0.2736  0.8344 

# Count how many gene pairs have the same blood and brain gene name
n_same <- sum(top_gene_pairs$blood_gene == top_gene_pairs$brain_gene, na.rm = TRUE); n_same #507
# Total number of pairs
n_total <- nrow(top_gene_pairs)
# Number of different gene pairs
n_different <- n_total - n_same; n_different #224,227

highest_brain_gene <- filter(top_gene_pairs, brain_gene == "ENSG00000226259.10")
#           blood_gene         brain_gene correlation
# 1 ENSG00000226259.10 ENSG00000226259.10   0.8344016
# 2 ENSG00000145736.14 ENSG00000226259.10   0.5498357
# 3  ENSG00000230847.4 ENSG00000226259.10   0.4523894
# 4  ENSG00000215630.6 ENSG00000226259.10   0.3453835
# 5  ENSG00000235558.3 ENSG00000226259.10   0.3189686
# 6 ENSG00000050165.17 ENSG00000226259.10   0.2651001
# 7 ENSG00000104722.14 ENSG00000226259.10   0.2563169

####### how many brain genes does each blood gene have the highest correlation with?
top_gene_pairs <- top_gene_pairs %>%
  group_by(blood_gene) %>%
  mutate(n_predictors = n()) %>%
  ungroup()

top_gene_pairs <- as.data.frame(top_gene_pairs) 
dim(top_gene_pairs) #224734      4

##summarize n_predictor for each blood gene
unique_blood_data <- top_gene_pairs[!duplicated(top_gene_pairs$blood_gene), ]
summary(unique_blood_data$n_predictors)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   1.00    2.00    4.00   14.63   10.00 1397.00

table(unique_blood_data$n_predictors)
   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
3778 2347 1511 1086  869  589  515  406  294  282  249  216  180  161  150  143 
  17   18   19   20   21   22   23   24   25   26   27   28   29   30   31   32 
 119  114  104  107   81   89   88   72   70   68   61   45   55   53   47   38 
  33   34   35   36   37   38   39   40   41   42   43   44   45   46   47   48 
  39   43   31   25   30   31   28   29   28   31   30   26   19   16   22   22 
  49   50   51   52   53   54   55   56   57   58   59   60   61   62   63   64 
  19   13   15   20   18   18   16   14   12   24    7   19   13   12   14   16 
  65   66   67   68   69   70   71   72   73   74   75   76   77   78   79   80 
  15   15   13   12   14   12   15    7    9    7    8    8    4    6    9    4 
  81   82   83   84   85   86   87   88   89   90   91   92   93   94   95   96 
   6    8    7    8    6    4    4    7    3    5    3    6    5   11    5    3 
  97   98   99  100  101  102  103  104  105  106  107  108  109  110  111  112 
   3    5    8    8    5    6    3    5    4    4    8    8    3    3    2    3 
 113  114  115  116  117  118  119  120  121  122  123  124  125  126  127  128 
   4    6    6    6    1    4    3    8    5    2    3    4    6    3    3    1 
 129  130  131  132  133  135  136  137  138  139  140  142  143  144  145  146 
   2    2    1    2    4    2    2    2    4    1    2    1    2    3    2    5 
 147  148  149  150  151  152  153  154  155  156  157  158  159  160  161  162 
   5    1    2    5    3    4    4    3    1    3    2    3    1    2    1    1 
 163  165  166  167  168  170  171  172  173  174  175  176  177  178  179  180 
   1    3    2    2    1    1    2    1    1    4    2    2    1    2    1    3 
 181  182  183  184  185  187  188  189  190  193  195  196  197  198  199  200 
   2    2    3    1    2    1    2    3    1    2    1    1    1    1    3    1 
 201  205  207  210  211  212  213  215  217  220  222  223  224  225  226  227 
   2    4    2    3    1    2    1    1    2    2    1    1    1    1    1    2 
 228  229  231  232  234  239  240  243  244  247  248  249  250  251  252  255 
   1    2    1    1    1    1    1    2    1    1    1    2    1    1    1    3 
 256  261  262  268  270  272  275  276  278  281  284  285  290  291  292  296 
   2    1    1    1    1    2    3    1    1    1    1    1    3    1    1    1 
 301  304  305  307  309  310  311  313  315  323  324  326  333  337  341  342 
   1    1    1    1    1    1    1    2    2    1    1    1    2    1    2    1 
 344  347  348  350  352  353  354  357  360  361  364  369  372  375  385  387 
   1    1    1    1    2    1    2    2    1    1    1    1    2    1    2    1 
 394  399  404  418  425  426  433  443  444  484  487  492  493  495  496  503 
   1    2    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
 519  521  526  530  548  559  566  584  594  595  596  620  635  655  683  722 
   1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
 728  731  801  808  810  955  981 1058 1397 
   1    1    1    1    1    1    1    1    1

dim(unique_blood_data) #15366     4

filter(unique_blood_data, n_predictors == "1397") #
#         blood_gene        brain_gene correlation n_predictors
#1 ENSG00000267532.5 ENSG00000236570.1   0.3265655         1397 MIR497HG| lncRNA that acts as a host for the miR‑497‑195 microRNA cluster
https://useast.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000267532;r=17:6875220-7022133

## removing brain_genes with more than n_predictors
filtered <- unique_blood_data[unique_blood_data$n_predictors <= 100, ]

plot <- ggplot(filtered, aes(x = n_predictors)) +
  geom_histogram(bins = 30, fill = "lightblue", color = "white", alpha = 0.7) +
  labs(title = "Distribution of Number of Connections to Brain Genes per Blood Gene with less than 100 connections",
       x = "Number of Connections to Brain Genes",
       y = "Count") +
  theme_minimal()

ggsave("/hpc/users/hoangd02/www/plots/number_of_brain_gene_connections_per_blood_gene_with_lessthan_100_connections.pdf", plot = plot, width = 9, height = 6) 

## next step: computing the blood gene expression co-expression matrix


################## EACH BRAIN GENE ##################
##adding n_predictors (number of blood genes associated per brain gene)
top_gene_pairs <- top_gene_pairs %>%
  group_by(brain_gene) %>%
  mutate(n_predictors = n()) %>%
  ungroup()

top_gene_pairs <- as.data.frame(top_gene_pairs) #224734      4
##summarize n_predictor for each brain gene
unique_brain_data <- top_gene_pairs[!duplicated(top_gene_pairs$brain_gene), ]
summary(unique_brain_data$n_predictors)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   1.00    1.00    3.00   14.32    9.00 2611.00 

## each brain gene have ~ 14 predictors on average

##one brain gene has A LOT of connection with other blood gene; that brain gene is 
filter(unique_brain_data, n_predictors == "2611") #ENSG00000109906.13 | ZBTB16, zinc finger and BTB domain containing 16
#seems to be present in a lot of tissues
https://gtexportal.org/home/gene/ENSG00000109906 

plot<- ggplot(unique_brain_data, aes(y = n_predictors)) +
  geom_boxplot(fill = "lightblue", alpha = 0.7, width = 0.5) +
  labs(title = "Distribution of Number of Predictors per Brain Gene",
       y = "Number of Predictors") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggsave("/hpc/users/hoangd02/www/plots/test.pdf", plot = plot, width = 7, height = 6) 

table(unique_brain_data$n_predictors)
   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
4281 2615 1542 1057  742  606  502  366  321  266  228  213  157  182  137  113 
  17   18   19   20   21   22   23   24   25   26   27   28   29   30   31   32 
 118  100   89   86   76   88   52   71   57   46   47   37   44   43   34   27 
  33   34   35   36   37   38   39   40   41   42   43   44   45   46   47   48 
  32   37   35   30   40   24   25   24   24   21   24   21   16   17   25   21 
  49   50   51   52   53   54   55   56   57   58   59   60   61   62   63   64 
  20   17   22   19   12   15   12   16   16   11   11   10   13   12   19   10 
  65   66   67   68   69   70   71   72   73   74   75   76   77   78   79   80 
   8   15   10   16    8   10   17    4   10    7    4   11   11    4    4    6 
  81   82   83   84   85   86   87   88   89   90   91   92   93   94   95   96 
  11    7   12    8    3    7   13    4    7    9    6    4    2    5    8    5 
  97   98   99  100  101  102  103  104  105  106  107  108  109  110  111  112 
   4    4    6    3    6    5    2    8    4    6    6    5    6    7    5    4 
 113  115  116  117  118  119  120  121  122  123  124  125  126  127  128  129 
   4    2    3    2    1    5    2    3    1    4    2    3    5    2    2    3 
 130  131  133  134  135  136  137  139  140  141  142  143  144  145  146  148 
   2    2    2    1    3    3    5    3    1    4    4    4    2    7    3    4 
 149  150  151  152  153  154  155  156  157  158  159  160  161  162  163  164 
   2    2    1    3    3    4    6    3    6    3    1    4    1    1    2    1 
 165  166  167  168  169  170  172  173  174  175  176  177  178  179  180  181 
   3    2    1    1    1    4    2    1    2    1    4    2    2    1    2    2 
 184  185  187  188  189  190  191  192  193  194  195  196  197  198  201  204 
   3    2    2    4    2    3    2    1    2    3    2    1    3    2    1    2 
 205  206  207  208  209  211  212  213  214  215  216  217  219  221  224  226 
   1    1    1    3    1    4    1    1    1    2    2    2    1    2    2    1 
 227  228  230  231  233  234  235  238  242  247  248  250  252  255  256  257 
   1    1    4    2    2    2    4    1    1    1    2    2    2    1    1    2 
 258  259  260  262  263  264  265  266  267  268  269  270  273  274  280  282 
   2    2    1    1    1    1    1    1    1    1    1    1    1    1    3    1 
 283  286  287  288  289  290  291  293  294  299  300  304  310  312  314  315 
   1    1    2    1    1    1    1    1    1    2    1    1    1    1    1    1 
 322  324  328  330  335  336  337  341  342  354  357  359  360  361  363  368 
   1    1    1    1    1    2    1    1    1    1    1    1    2    2    1    1 
 369  372  374  379  380  391  393  400  409  417  435  441  448  449  450  453 
   1    1    1    1    2    2    1    1    2    1    1    1    2    1    1    1 
 463  465  466  473  483  485  487  494  511  530  562  598  599  640  667  672 
   1    3    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
 687  700  713  747  905  997  998 1145 1694 2466 2611 
   1    1    1    1    1    1    1    1    1    1    1 

## removing brain_genes with more than n_predictors
filtered <- unique_brain_data[unique_brain_data$n_predictors <= 150, ]

plot <- ggplot(filtered, aes(x = n_predictors)) +
  geom_histogram(bins = 30, fill = "lightblue", color = "white", alpha = 0.7) +
  labs(title = "Distribution of Number of Predictors per Brain Gene for Brain Gene with less than 150 predictors",
       x = "Number of Predictors",
       y = "Count") +
  theme_minimal()

ggsave("/hpc/users/hoangd02/www/plots/test2.pdf", plot = plot, width = 9, height = 6) 


##########################################################
######### Blood WITH INDIVIDUAL | Brain baseline ######### 
##########################################################
blood_form3_by_brain_cor <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form3_brain_baseline_spearman_cor_matrix.txt", data.table=FALSE)
row.names(blood_form3_by_brain_cor) <- blood_form3_by_brain_cor$V1; blood_form3_by_brain_cor$V1 <- NULL
blood_form3_by_brain_cor <- abs(blood_form3_by_brain_cor)
blood_form3_99.95 <- quantile(as.matrix(blood_form3_by_brain_cor), probs = 0.9995); blood_form3_99.95 #0.2300379
blood_form3_99.995 <- quantile(as.matrix(blood_form3_by_brain_cor), probs = 0.99995); blood_form3_99.995 #0.2665065

# Compute the 99.95th percentile
threshold <- quantile(as.matrix(blood_form3_by_brain_cor), probs = 0.9995)
# Identify which gene pairs meet or exceed the threshold
high_corr_indices <- which(blood_form3_by_brain_cor >= threshold, arr.ind = TRUE)

# Create a dataframe with the corresponding gene names and correlation values
top_gene_pairs <- data.frame(
  blood_gene = rownames(blood_form3_by_brain_cor)[high_corr_indices[, 1]],
  brain_gene = colnames(blood_form3_by_brain_cor)[high_corr_indices[, 2]],
  correlation = blood_form3_by_brain_cor[high_corr_indices]
)
# Sort by correlation strength
top_gene_pairs <- top_gene_pairs[order(-top_gene_pairs$correlation), ]
length(unique(top_gene_pairs$brain_gene)) #17843  - A LOT MORE THAN NO RESID ID
length(unique(top_gene_pairs$blood_gene)) #16702 - A LOT MORE THAN NO RESID ID

head(top_gene_pairs)
#               blood_gene         brain_gene correlation
#10261  ENSG00000184368.16 ENSG00000067048.17   0.5660830
#10469  ENSG00000184368.16 ENSG00000067646.12   0.5657606
#208801 ENSG00000184368.16  ENSG00000274655.1   0.5608660
#164383 ENSG00000184368.16 ENSG00000229807.11   0.5426970
#144605 ENSG00000184368.16 ENSG00000198692.10   0.5409555
#58327  ENSG00000184368.16 ENSG00000129824.16   0.5384197

summary(top_gene_pairs$correlation)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2300  0.2348  0.2413  0.2465  0.2522  0.5661 

# Count how many gene pairs have the same blood and brain gene name
n_same <- sum(top_gene_pairs$blood_gene == top_gene_pairs$brain_gene, na.rm = TRUE); n_same #361
# Total number of pairs
n_total <- nrow(top_gene_pairs)
# Number of different gene pairs
n_different <- n_total - n_same; n_different #224,374

##################################################
######### NO INDIVIDUAL | Brain full ############# 
##################################################
blood_form3_by_brain_cor <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form3_no_indivdualID_brain_full_spearman_cor_matrix.txt", data.table=FALSE)
row.names(blood_form3_by_brain_cor) <- blood_form3_by_brain_cor$V1; blood_form3_by_brain_cor$V1 <- NULL
blood_form3_by_brain_cor <- abs(blood_form3_by_brain_cor)
blood_form3_99.95 <- quantile(as.matrix(blood_form3_by_brain_cor), probs = 0.9995); blood_form3_99.95 #0.2414433
blood_form3_99.995 <- quantile(as.matrix(blood_form3_by_brain_cor), probs = 0.99995); blood_form3_99.995 #0.2862253 

# Compute the 99.95th percentile
threshold <- quantile(as.matrix(blood_form3_by_brain_cor), probs = 0.9995)
# Identify which gene pairs meet or exceed the threshold
high_corr_indices <- which(blood_form3_by_brain_cor >= threshold, arr.ind = TRUE)

# Create a dataframe with the corresponding gene names and correlation values
top_gene_pairs <- data.frame(
  blood_gene = rownames(blood_form3_by_brain_cor)[high_corr_indices[, 1]],
  brain_gene = colnames(blood_form3_by_brain_cor)[high_corr_indices[, 2]],
  correlation = blood_form3_by_brain_cor[high_corr_indices]
)
# Sort by correlation strength
top_gene_pairs <- top_gene_pairs[order(-top_gene_pairs$correlation), ]
length(unique(top_gene_pairs$brain_gene)) #16789
length(unique(top_gene_pairs$blood_gene)) #19683

head(top_gene_pairs)
#               blood_gene         brain_gene correlation
#168716 ENSG00000226259.10 ENSG00000226259.10   0.8779214
#169105  ENSG00000226752.9  ENSG00000226752.9   0.8499747
#180166  ENSG00000241945.8  ENSG00000241945.8   0.8441066
#85321  ENSG00000145736.14 ENSG00000145736.14   0.8389223
#178047  ENSG00000237541.3  ENSG00000237541.3   0.8304720
#207509  ENSG00000274602.5  ENSG00000274602.5   0.8243647

summary(top_gene_pairs$correlation)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2414  0.2469  0.2546  0.2619  0.2677  0.8779

# Count how many gene pairs have the same blood and brain gene name
n_same <- sum(top_gene_pairs$blood_gene == top_gene_pairs$brain_gene, na.rm = TRUE); n_same #621
# Total number of pairs
n_total <- nrow(top_gene_pairs)
# Number of different gene pairs
n_different <- n_total - n_same; n_different #224,113

##########################################################
######### Blood WITH INDIVIDUAL | Brain full ############# 
##########################################################
blood_form3_by_brain_cor <- fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/blood_form3_brain_full_spearman_cor_matrix.txt", data.table=FALSE)
row.names(blood_form3_by_brain_cor) <- blood_form3_by_brain_cor$V1; blood_form3_by_brain_cor$V1 <- NULL
blood_form3_by_brain_cor <- abs(blood_form3_by_brain_cor)
blood_form3_99.95 <- quantile(as.matrix(blood_form3_by_brain_cor), probs = 0.9995); blood_form3_99.95 #0.23496
blood_form3_99.995 <- quantile(as.matrix(blood_form3_by_brain_cor), probs = 0.99995); blood_form3_99.995 #0.2741446

# Compute the 99.95th percentile
threshold <- quantile(as.matrix(blood_form3_by_brain_cor), probs = 0.9995)
# Identify which gene pairs meet or exceed the threshold
high_corr_indices <- which(blood_form3_by_brain_cor >= threshold, arr.ind = TRUE)

# Create a dataframe with the corresponding gene names and correlation values
top_gene_pairs <- data.frame(
  blood_gene = rownames(blood_form3_by_brain_cor)[high_corr_indices[, 1]],
  brain_gene = colnames(blood_form3_by_brain_cor)[high_corr_indices[, 2]],
  correlation = blood_form3_by_brain_cor[high_corr_indices]
)
# Sort by correlation strength
top_gene_pairs <- top_gene_pairs[order(-top_gene_pairs$correlation), ]
length(unique(top_gene_pairs$brain_gene)) #17956  - A LOT MORE THAN NO RESID ID
length(unique(top_gene_pairs$blood_gene)) #20177 - A LOT MORE THAN NO RESID ID

head(top_gene_pairs)
#               blood_gene         brain_gene correlation
#19718 ENSG00000169902.15 ENSG00000096060.14   0.5491456
#19072 ENSG00000123836.15 ENSG00000096060.14   0.5439760
#16784  ENSG00000185950.9 ENSG00000090104.12   0.5179783
#19861 ENSG00000180509.12 ENSG00000096060.14   0.5112611
#43016 ENSG00000134463.15 ENSG00000118515.11   0.5032986
#18951 ENSG00000115590.14 ENSG00000096060.14   0.5009966

summary(top_gene_pairs$correlation)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 #0.2350  0.2401  0.2471  0.2523  0.2589  0.5491

# Count how many gene pairs have the same blood and brain gene name
n_same <- sum(top_gene_pairs$blood_gene == top_gene_pairs$brain_gene, na.rm = TRUE); n_same #204
# Total number of pairs
n_total <- nrow(top_gene_pairs)
# Number of different gene pairs
n_different <- n_total - n_same; n_different #224,534

#####################
#####################
# Observed counts
n_same_top <- 204
n_diff_top <- 224534

# Background totals
n_same_total <- 17533  # number of overlapping genes between blood & brain
n_diff_total <- 21046 * 21356 - n_same_total  # total pairs minus overlap

# Not in top
n_same_not <- n_same_total - n_same_top
n_diff_not <- n_diff_total - n_diff_top

# Build contingency table
contingency <- matrix(c(n_same_top, n_same_not,
                        n_diff_top, n_diff_not),
                      nrow = 2,
                      byrow = TRUE)

colnames(contingency) <- c("Top_0.05pct", "Not_Top")
rownames(contingency) <- c("Same_gene", "Different_gene")

contingency
# Fisher's exact test with odds ratio
fisher_result <- fisher.test(contingency, alternative = "greater")

fisher_result$estimate   # odds ratio, 23.55151
fisher_result$p.value    #8.329662e-197

lnOR <- log(fisher_result$estimate); lnOR # 3.15919 


##########################################################
##########################################################
form3_full_brain_no_residID<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_results_form3_full_brain_predicted_from_blood.txt", data.table = FALSE)
quantile(form3_full_brain_no_residID$r2_test, 0.9995, na.rm = TRUE) #0.8400378

form3_full_brain_with_residID<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_results_form3_full_brain_with_residID_predicted_from_blood.txt", data.table = FALSE)
quantile(form3_full_brain_with_residID$r2_test, 0.9995, na.rm = TRUE) #0.4831525

form3_base_brain_with_residID<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_results_form3_base_brain_with_residID_predicted_from_blood.txt", data.table = FALSE)
quantile(form3_base_brain_with_residID$r2_test, 0.9995, na.rm = TRUE) #0.808647

form3_base_brain_no_residID<-fread("/sc/arion/projects/mscic1/results/jolie/LBP/blood-brain/qc2/elastic_net_results_form3_base_brain_predicted_from_blood.txt", data.table = FALSE)
quantile(form3_base_brain_no_residID$r2_test, 0.9995, na.rm = TRUE) #0.8350979

