##archived

gse_id <- "GSE7621" #see GSE7621_sup doc
gse <- getGEO(gse_id, GSEMatrix = TRUE) # Download the GEO dataset
expression_set <- gse[[1]] # Extract the expression set
metadata <- pData(expression_set) # Extract metadata
gene_expression <- exprs(expression_set) # Extract gene expression matrix; dim:  54339  x 25
sum(is.na(gene_expression)) #there are 525 NAs in the gene expression matrix
gene_expression <- na.omit(gene_expression) #drop rows with NAs; dim: 54318 x 25 (delta = 21 rows)

metadata[,c("characteristics_ch1.1","source_name_ch1")]
#write.table(metadata, file = "metadata_GSE7621.txt", sep = "\t", row.names = TRUE, col.names = NA)
#write.table(gene_expression, file = "gene_expression_GSE7621.txt", sep = "\t", row.names = TRUE, col.names = NA)

##disease status
metadata$characteristics_ch1 <- ifelse(metadata$characteristics_ch1 == "Parkinson's Disease", "PD", metadata$characteristics_ch1)
metadata$characteristics_ch1 <- ifelse(metadata$characteristics_ch1 == "Old Control", "control", metadata$characteristics_ch1)
metadata$disease_status <- as.factor(metadata$characteristics_ch1)
table(metadata$disease_status)

#covariates
metadata$gender<-metadata$characteristics_ch1.1
metadata$gender <- ifelse(metadata$gender == "female", "F", metadata$gender)
metadata$gender <- ifelse(metadata$gender == "male", "M", metadata$gender)
metadata$gender<-as.factor(metadata$gender); table(metadata$gender)

##DEA: log2 transformation, quantile normalization, covariate correction and limma
# Log2 transformation
gene_expression <- log2(gene_expression + 0.5) #add 0.5 to avoid taking log of 0

# Quantile normalization
gene_expression2<- normalize.quantiles(as.matrix(gene_expression))

##run until here
###base DEA, no SVA
# Design matrix for limma
design <- model.matrix(~ 0 + disease_status +  gender, data = metadata); head(design)
colnames(design) <- c("control", "PD", "gender")

# Fit the linear model
fit <- lmFit(gene_expression, design)

# Create a contrast matrix to compare PD vs control
contrast.matrix <- makeContrasts(PD_vs_control = PD - control, levels = design)

# Fit the contrasts
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Summarize the results
results <- topTable(fit2, coef = "PD_vs_control", adjust.method = "fdr", number = Inf)
results$probeID <- rownames(results)
results<- results[, c(7,1,2,3,4,5,6)]
rownames(results) <- seq_len(nrow(results))
head(results)

write.table(results, file = "/sc/arion/projects/mscic1/results/jolie/GEO/eli/GSE7621_DEA_results_without_sva_final.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

GSE7621_no_sva<-fread("/sc/arion/projects/mscic1/results/jolie/GEO/eli/GSE7621_DEA_results_without_sva_final.txt",data.table = FALSE)
#results$X<-rownames(results)





# Quantile normalization
gene_ids <- rownames(gene_expression)
gene_expression <- normalize.quantiles(as.matrix(gene_expression))
rownames(gene_expression) <- gene_ids


results$probeID <- rownames(results)
results<- results[, c(7,1,2,3,4,5,6)]
rownames(results) <- seq_len(nrow(results))
head(results)


write.table(results, file = "/sc/arion/projects/mscic1/results/jolie/GEO/eli/GSE20141_DEA_results_without_sva_final.txt", sep = "\t", row.names = FALSE, col.names = TRUE)






















#######
res_0_sv_sim <- res %>% 
  filter(param.n_sv == 0)

res0_sv_sim_sub <- res_0_sv_sim[,c("method","nsv")] #there are some bugs here, @ryan, @noam

table(res0_sv_sim_sub$method, res0_sv_sim_sub$nsv)

be_0sv_sim <- res0_sv_sim_sub %>% 
  filter(method == "be")
table(be_0sv_sim$nsv)
#                 0 0.111111111111111             0.125 0.142857142857143 
#              1193                 6                15                 4 
# 0.166666666666667               0.2 0.333333333333333 
#                 7                24                 1 

leek_0sv_sim <- res0_sv_sim_sub %>% 
  filter(method == "leek")
table(leek_0sv_sim$nsv)
#    0 
# 1250

length(unique(res$nsv)) #1749 unique values in nsv (while we simulated 11 values of nsv)

# Get unique whole numbers
unique_whole <- unique(res$nsv[res$nsv %% 1 == 0])
whole_count <- length(unique_whole) #132 of which are whole numbers 

# Get unique decimal numbers
unique_decimal <- unique(res$nsv[res$nsv %% 1 != 0])
decimal_count <- length(unique_decimal) # 1617 of which are decimal numbers 

nsv_0.2 <- res %>% 
  filter(nsv == 0.2)  

table(nsv_0.2$param.n_sv)
#insane
#  0   10   20   30   40   50   60   70   80   90  100 
# 24  497  406  835  742 1128  441  589  472  629  298 

table(nsv_0.2$method)
  be leek 
3099 2962

#pick 1 sim for be, row 12 when do head(nsv_0.2)
nsv_0.2[12, ]

sim1 <- grid_table %>%
  filter(
    abs(total_sv_sd - 0) < 1e-6,
    abs(fraction_degs - 0.10) < 1e-6,
    abs(status_sd - 0.1) < 1e-6,
    abs(resid_sd - 0.3) < 1e-6,
    abs(fraction_sv - 0.3) < 1e-6,
    seed == 12345
  )

all_nsv_0.2_rows <- nsv_0.2 %>%
  filter(
    #abs(total_sv_sd - 0) < 1e-6,
    abs(param.fraction_degs - 0.10) < 1e-6,
    abs(param.status_sd - 0.1) < 1e-6,
    abs(param.resid_sd - 0.3) < 1e-6,
    abs(param.fraction_sv - 0.3) < 1e-6,
    #seed == 12345
  ); all_nsv_0.2_rows #3 BE's and 1 leek

#is it an average across seed problem?? YES!


#for the rest of the simulation, see sim_bug.R
#be = 0.02 with the folloing param 
row=1
n_sv = 50 #0 to 100 (all = 0.2 SVs for be in res)
total_sv_sd = 0
grid_sv = as.data.frame(expand.grid(n_sv=n_sv,total_sv_sd=total_sv_sd))
grid_sv$sv_sd = sqrt(grid_sv$total_sv_sd^2 / grid_sv$n_sv)
grid_sv$sv_sd[!is.finite(grid_sv$sv_sd)] = 0
grid_list <- list(grid_sv,
                  data.frame(fraction_degs = 0.1),
                  data.frame(status_sd = 0.1),
                  data.frame(resid_sd = 0.3),
                  data.frame(fraction_sv = 0.3),
                  data.frame(seed = 12345))
gridsubset=as.data.frame(t(as.data.frame(unlist(grid_list))))

n_sv = 0 
             TN FP  FN   TP    method nsv param.n_sv param.fraction_degs
be        17932 68 609 1391        be   0          0                 0.1
leek      17932 68 609 1391      leek   0          0                 0.1
none      17932 68 609 1391      none   0          0                 0.1
known_nsv 17932 68 607 1393 known_nsv   1          0                 0.1
oracle    17932 68 609 1391    oracle   1          0                 0.1
          param.fraction_sv param.status_sd param.sv_sd param.resid_sd
be                      0.3             0.1           0            0.3
leek                    0.3             0.1           0            0.3
none                    0.3             0.1           0            0.3
known_nsv               0.3             0.1           0            0.3
oracle                  0.3             0.1           0            0.3
          param.seed
be             12345
leek           12345
none           12345
known_nsv      12345
oracle         12345

n_sv = 50 
#              TN FP  FN   TP    method nsv param.n_sv param.fraction_degs
# be        17926 74 567 1433        be   0         50                 0.1
# leek      17926 74 567 1433      leek   0         50                 0.1
# none      17926 74 567 1433      none   0         50                 0.1
# known_nsv 17913 87 575 1425 known_nsv  50         50                 0.1
# oracle    17926 74 567 1433    oracle  50         50                 0.1
#           param.fraction_sv param.status_sd param.sv_sd param.resid_sd
# be                      0.3             0.1           0            0.3
# leek                    0.3             0.1           0            0.3
# none                    0.3             0.1           0            0.3
# known_nsv               0.3             0.1           0            0.3
# oracle                  0.3             0.1           0            0.3
#           param.seed
# be             12345
# leek           12345
# none           12345
# known_nsv      12345
# oracle         12345

n_sv = 100 
#              TN FP  FN   TP    method nsv param.n_sv param.fraction_degs
# be        17926 74 596 1404        be   0        100                 0.1
# leek      17926 74 596 1404      leek   0        100                 0.1
# none      17926 74 596 1404      none   0        100                 0.1
# known_nsv 17940 60 605 1395 known_nsv 100        100                 0.1
# oracle    17926 74 596 1404    oracle 100        100                 0.1
#           param.fraction_sv param.status_sd param.sv_sd param.resid_sd
# be                      0.3             0.1           0            0.3
# leek                    0.3             0.1           0            0.3
# none                    0.3             0.1           0            0.3
# known_nsv               0.3             0.1           0            0.3
# oracle                  0.3             0.1           0            0.3
#           param.seed
# be             12345
# leek           12345
# none           12345
# known_nsv      12345
# oracle         12345

nsv_0.3 <- res %>% 
  filter(abs(nsv - 0.333333333333333) < 1e-9)

row=1
n_sv = 100 #0 to 100 (all = 0.2 SVs for be in res)
total_sv_sd = 0.5 #to get sv_sd = 0.05
grid_sv = as.data.frame(expand.grid(n_sv=n_sv,total_sv_sd=total_sv_sd))
grid_sv$sv_sd = sqrt(grid_sv$total_sv_sd^2 / grid_sv$n_sv)
grid_sv$sv_sd[!is.finite(grid_sv$sv_sd)] = 0
grid_list <- list(grid_sv,
                  data.frame(fraction_degs = 0.1),
                  data.frame(status_sd = 0.4),
                  data.frame(resid_sd = 0.4),
                  data.frame(fraction_sv = 0.3),
                  data.frame(seed = 12350))
gridsubset=as.data.frame(t(as.data.frame(unlist(grid_list))))
#              TN FP  FN   TP    method nsv param.n_sv param.fraction_degs
# be        17909 91 227 1773        be   1        100                 0.1
# leek      17909 91 227 1773      leek   0        100                 0.1
# none      17909 91 227 1773      none   0        100                 0.1
# known_nsv 17911 89 231 1769 known_nsv 100        100                 0.1
# oracle    17921 79 235 1765    oracle 100        100                 0.1
#           param.fraction_sv param.status_sd param.sv_sd param.resid_sd
# be                      0.3             0.4        0.05            0.4
# leek                    0.3             0.4        0.05            0.4
# none                    0.3             0.4        0.05            0.4
# known_nsv               0.3             0.4        0.05            0.4
# oracle                  0.3             0.4        0.05            0.4
#           param.seed
# be             12350
# leek           12350
# none           12350
# known_nsv      12350
# oracle         12350

##find nsv 
##see if nsv rounds up or down

res_sub <- res[res$nsv < 20 & > 10]


#pick 1 sim for leek
