#################################### SIMULATION PVALUE DISTRIBUTIONS ####################################
#########################################################################################################

#Code for running one simulation (aka one set of parameters) is in sim_bug.R or here
library(rlang)
library(stringr)
library(fs)
library(withr)
library(assertthat)
library(igraph)
library(plyr)
library(dplyr)
library(S4Vectors)
library(variancePartition)
library(glmnet)
library(BiocParallel)
library(rctutils)
library(magrittr)
library(limma)
library(data.table)
library(sva)

## Correctly handles non-syntactic names
make_additive_formula <- function(varnames, lhs = NULL, env = caller_env()) {
    form_rhs <- Reduce(
        function(e1, e2) expr(`+`(!!e1, !!e2)),
        syms(varnames)
    )
    if (!is.null(lhs)) {
        lhs <- as.symbol(lhs)
    }
    new_formula(lhs = lhs, rhs = form_rhs, env = env)
}

design_is_identifiable <- function(dmat, eps = 2e-16) {
    min(svd(dmat)$d) > eps
}

generate_ids <- function(prefix, n) {
    if (length(n) == 1) {
        n <- seq_len(n)
    }
    num_digits <- digits(max(n))
    digit_format <- str_c("%0", num_digits, "i")
    str_c(prefix, sprintf(digit_format, n))
}

digits <- function(x, base = 10) {
    1 + floor(log(x, base = base))
}

estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

getBetaParams_pi <- function(mu, pi) {
  a = sqrt(pi*mu / (1 - mu))
  b = a * (1-mu)/mu
  list(alpha=a, beta=b)
}

getBetaParams <- function(mu, sum = 4) {
    # when mean is 0.5, alpha = beta = 2
    a = mu * sum
    b = sum - a
  list(alpha=a, beta=b)
}

sim_case_control_data <- function (
    n_subjects,
    n_samples_per_subject,
    n_features,
    n_batches,
    n_sv,
    fraction_degs,
    fraction_sv,
    status_sd,
    subject_sd,
    batch_sd,
    sv_sd,
    resid_sd,
    variable_samples_per_subject = FALSE,
    seed,
    ## How many times to try shuffling batch before giving up. It
    ## should really never take this long unless there's no
    ## solution.
    .max_tries = 500
) {
    if (!missing(seed)) {
        sim_params <- mget(c(
            "n_subjects",
            "n_samples_per_subject",
            "n_features",
            "n_batches",
            "n_sv",
            "fraction_degs",
            "fraction_sv",
            "status_sd",
            "subject_sd",
            "batch_sd",
            "sv_sd",
            "resid_sd",
            "variable_samples_per_subject",
            ".max_tries"
        ))
        res <- with_seed(seed, do.call(sim_case_control_data, sim_params))
        res$params$seed <- seed
        return(res)
    }
    assert_that(
        n_subjects > 1,
        n_samples_per_subject >= 1,
        n_batches >= 1,
        fraction_degs >= 0,
        fraction_degs <= 1,
        fraction_sv >= 0,
        fraction_sv <= 1,
        status_sd >= 0,
        subject_sd >= 0,
        batch_sd >= 0,
        sv_sd >= 0,
        resid_sd >= 0
    )
    if (variable_samples_per_subject && n_samples_per_subject == 2) {
        warning("Setting variable_samples_per_subject to FALSE because n_samples_per_subject is 2")
        variable_samples_per_subject <- FALSE
    }
    if (n_batches < 2 && batch_sd > 0) {
        warning("Setting batch signal to zero since there is only 1 batch")
        n_batches <- 1
        batch_sd <- 0
    }
    if (n_sv < 1) {
        warning("Setting SV signal to zero since there are not SVs")
        n_sv <- 0 ##changed this from 1 to 0 on 02/01/2025
        sv_sd <- 0
    }

    ## Generate the samples
    subject_ids <- generate_ids("S", n_subjects)
    batch_ids <- generate_ids("B", n_batches)
    sv_names <- generate_ids("SV", n_sv)

    ## Generate the sd for each SV so that the total variance of SV is fixed
    sv_var_vec_exact <- rexp(n = n_sv)
    sv_var_vec_exact <- sv_var_vec_exact / mean(sv_var_vec_exact) * sv_sd^2
    sv_sd_vec_exact <- sqrt(sv_var_vec_exact)
    sv_sd_vec_exact <- sort(sv_sd_vec_exact,decreasing=TRUE)
    # message("Exact: ")
    # print(sv_sd_vec_exact)
    # message("Total variance: ", sum(sv_sd_vec_exact^2))
    ## Generate the fraction of genes impacted for each SV
    ## The mean is a/(a+b) and the variance is ab/((a+b)^2 (a+b+1)).
    betaParams=getBetaParams(mu = fraction_sv, sum = 4)
    assert_that(betaParams$alpha>0 & betaParams$beta>0)
    fraction_svs <- rbeta(n_sv, betaParams$alpha, betaParams$beta)
    # fraction_svs

    sv_df <- data.frame(sv_names=sv_names,sv_sds=sv_sd_vec_exact, fraction_svs = fraction_svs)

    if (variable_samples_per_subject) {
        n_samples <- n_samples_per_subject * n_subjects
        sample_subject_meta_table <- tibble(
            Subject_ID =
                sample(subject_ids, size = n_samples, replace = TRUE)
        ) %>%
            arrange(Subject_ID)
    } else {
        sample_subject_meta_table <- tibble(
            ## Assign a fixed number of samples per subject
            Subject_ID = rep(subject_ids, each = n_samples_per_subject)
        )
    }

    sample_meta_table <- sample_subject_meta_table %>%
        mutate(
            Status = sample(rep(c("case","control"),length.out = n()))
        ) %>%
        arrange(Subject_ID, Status) %>%
        group_by(Subject_ID) %>%
        mutate(
            Sample_ID = make.names(str_c(Status, "_", Subject_ID), unique=TRUE)
        ) %>%
        ungroup
    for (svi in seq_len(nrow(sv_df))) {
        sample_meta_table[[sv_df$sv_names[svi]]] <- rnorm(nrow(sample_meta_table), sd = sv_df$sv_sds[svi])
    }
    if(sv_sd > 0){
        model_vars=c("Status", sv_df$sv_names)
    }else{
        model_vars="Status"
    }
    if(n_samples_per_subject > 1){
        model_vars=c(model_vars,"Subject_ID")
    }
    form <- make_additive_formula(model_vars)
    dmat <- model.matrix(form, sample_meta_table)
    ## This should always be true
    assert_that(design_is_identifiable(dmat))
    ## Shuffle Batch ID until model is identifiable
    batch_gen_success <- FALSE
    if(batch_sd > 0){
        batch_form <- update(form, ~ . + Batch_ID)
        for (i in seq_len(max(.max_tries))) {
            sample_meta_table$Batch_ID <- sample(batch_ids, nrow(sample_meta_table), replace = TRUE)
            dmat <- model.matrix(make_additive_formula(c(model_vars,"Batch_ID")), sample_meta_table)
            batch_gen_success <- design_is_identifiable(dmat)
            if (batch_gen_success) {
                break
            }
        }
        assert_that(batch_gen_success)
    }else{
        sample_meta_table$Batch_ID <- sample(batch_ids, nrow(sample_meta_table), replace = TRUE)
    }
    n_de_features <- round(n_features * fraction_degs)

    ## Generate the features
    feature_meta_table <- tibble(
        Feature_Number = seq_len(n_features),
        DE = Feature_Number <= n_de_features,
        Feature_ID = generate_ids("Feature", n_features) %>%
            str_c(if_else(DE, "_DE", ""))
    ) %>%
        select(Feature_ID, Feature_Number, DE, everything())

    ## Generate the various expression effects
    design_subject <- model.matrix(~0 + Subject_ID, sample_meta_table)
    design_batch <- model.matrix(~0 + Batch_ID, sample_meta_table)
    design_sv <- as.matrix(sample_meta_table[sv_df$sv_names])
    design_status <- model.matrix(~0 + Status, sample_meta_table)
    eta_subject <- design_subject %*% matrix(ncol = n_features, nrow = ncol(design_subject), rnorm(n_features * ncol(design_subject)))
    eta_batch <- design_batch %*% matrix(ncol = n_features, nrow = ncol(design_batch), rnorm(n_features * ncol(design_batch)))
    sv_matrix = matrix(ncol = n_features, nrow = ncol(design_sv), rnorm(n_features * ncol(design_sv)))
    for(sv in 1:n_sv){
        number_of_genes_to_keep = max(round(ncol(sv_matrix)*sv_df$fraction_svs[sv]),1)
        sv_matrix[sv, -sample(1:ncol(sv_matrix),number_of_genes_to_keep)] = 0
    }
    eta_sv <- design_sv %*%  sv_matrix
    eta_status <- design_status %*% matrix(ncol = n_features, nrow = ncol(design_status), rnorm(n_features * ncol(design_status)))
    eta_status[,!feature_meta_table$DE]<-0
    ## ## Extra random variation to take the place of time variation
    ## ## for features that are not DE with time
    ## eta_nontime <- matrix(ncol = n_features, nrow = nrow(sample_meta_table), rnorm(n_features * ncol(design_status))) %>%
    ##     scale(center = FALSE, scale = if_else(feature_meta_table$DE, Inf, 1))
    eta_resid <- matrix(ncol = n_features, nrow = nrow(sample_meta_table), rnorm(n_features * nrow(sample_meta_table)))

    expr_mat <- t(
        eta_subject * subject_sd +
        eta_batch * batch_sd +
        eta_sv * sv_sd +
        eta_status * status_sd +
        ## eta_nontime * status_sd +
        eta_resid * resid_sd
    ) %>%
        set_colnames(sample_meta_table$Sample_ID) %>%
        set_rownames(feature_meta_table$Feature_ID)
    attr(expr_mat, "scaled:scale") <- NULL
    elist <- new("EList", list(
        E = expr_mat,
        genes = as.data.frame(feature_meta_table),
        targets = as.data.frame(sample_meta_table),
        other = list(
            eta_subject = t(eta_subject),
            eta_batch = t(eta_batch),
            eta_sv = t(eta_sv),
            eta_status = t(eta_status),
            eta_resid = t(eta_resid)
        ),
        params = mget(c(
            "n_subjects",
            "n_samples_per_subject",
            "n_features",
            "n_batches",
            "n_sv",
            "fraction_degs",
            "status_sd",
            "subject_sd",
            "batch_sd",
            "sv_sd",
            "resid_sd"
        )),
        sv_df = sv_df
    ))
    ## Discard unused eta matrices
    if (subject_sd == 0) {
        elist$other$eta_subject <- NULL
    }
    if (batch_sd == 0) {
        elist$other$eta_batch <- NULL
    }
    if (sv_sd == 0) {
        elist$other$eta_sv <- NULL
    }
    if (status_sd == 0) {
        elist$other$eta_status <- NULL
    }
    if (resid_sd == 0) {
        elist$other$eta_resid <- NULL
    }
    ## Apply dimnames to all relevant elements
    dimnames(elist) <- dimnames(elist)
    return(elist)
}

##FOUR SIMULATION CASES IN SUPPLEMENTARY FIGURE 1:
#CASE 0: --- IN FIGURE 2:
row=1
n_sv = 10 
total_sv_sd = 2.5
grid_sv = as.data.frame(expand.grid(n_sv=n_sv,total_sv_sd=total_sv_sd))
grid_sv$sv_sd = sqrt(grid_sv$total_sv_sd^2 / grid_sv$n_sv)
grid_sv$sv_sd[!is.finite(grid_sv$sv_sd)] = 0
grid_list <- list(grid_sv,
                  data.frame(fraction_degs = 0.04),
                  data.frame(status_sd = 0.1),
                  data.frame(resid_sd = 0.4),
                  data.frame(fraction_sv = 0.3),
                  data.frame(seed = 12345)) 
gridsubset=as.data.frame(t(as.data.frame(unlist(grid_list))))


#CASE 1:
row=1
n_sv = 10 
total_sv_sd = 5
grid_sv = as.data.frame(expand.grid(n_sv=n_sv,total_sv_sd=total_sv_sd))
grid_sv$sv_sd = sqrt(grid_sv$total_sv_sd^2 / grid_sv$n_sv)
grid_sv$sv_sd[!is.finite(grid_sv$sv_sd)] = 0
grid_list <- list(grid_sv,
                  data.frame(fraction_degs = 0.3),
                  data.frame(status_sd = 0.2),
                  data.frame(resid_sd = 0.2),
                  data.frame(fraction_sv = 0.1),
                  data.frame(seed = 12352)) 
gridsubset=as.data.frame(t(as.data.frame(unlist(grid_list))))

# #CASE 2:
row=1
n_sv = 10 
total_sv_sd = 7.5
grid_sv = as.data.frame(expand.grid(n_sv=n_sv,total_sv_sd=total_sv_sd))
grid_sv$sv_sd = sqrt(grid_sv$total_sv_sd^2 / grid_sv$n_sv)
grid_sv$sv_sd[!is.finite(grid_sv$sv_sd)] = 0
grid_list <- list(grid_sv,
                  data.frame(fraction_degs = 0.08),
                  data.frame(status_sd = 0.5),
                  data.frame(resid_sd = 0.3),
                  data.frame(fraction_sv = 0.5),
                  data.frame(seed = 12351)) 
gridsubset=as.data.frame(t(as.data.frame(unlist(grid_list))))

# ##CASE 3: 
row=1
n_sv = 0 #10
total_sv_sd = 1.5
grid_sv = as.data.frame(expand.grid(n_sv=n_sv,total_sv_sd=total_sv_sd))
grid_sv$sv_sd = sqrt(grid_sv$total_sv_sd^2 / grid_sv$n_sv)
grid_sv$sv_sd[!is.finite(grid_sv$sv_sd)] = 0
grid_list <- list(grid_sv,
                  data.frame(fraction_degs = 0.2),
                  data.frame(status_sd = 0.2),
                  data.frame(resid_sd = 0.4),
                  data.frame(fraction_sv = 0.1),
                  data.frame(seed = 12348)) 
gridsubset=as.data.frame(t(as.data.frame(unlist(grid_list))))

# ##CASE 4: 
row=1
n_sv = 10 
total_sv_sd = 2.0
grid_sv = as.data.frame(expand.grid(n_sv=n_sv,total_sv_sd=total_sv_sd))
grid_sv$sv_sd = sqrt(grid_sv$total_sv_sd^2 / grid_sv$n_sv)
grid_sv$sv_sd[!is.finite(grid_sv$sv_sd)] = 0
grid_list <- list(grid_sv,
                  data.frame(fraction_degs = 0.3),
                  data.frame(status_sd = 0.2),
                  data.frame(resid_sd = 0.3),
                  data.frame(fraction_sv = 0.7),
                  data.frame(seed = 12348)) 
gridsubset=as.data.frame(t(as.data.frame(unlist(grid_list))))

simData = sim_case_control_data(n_subjects = 1000, n_samples_per_subject = 1, n_features = 20000, n_batches = 15, 
	n_sv = gridsubset$n_sv[row], fraction_degs = gridsubset$fraction_degs[row], fraction_sv = gridsubset$fraction_sv[row], 
    status_sd = gridsubset$status_sd[row], subject_sd = 0, batch_sd = 0.3, sv_sd = gridsubset$sv_sd[row], 
    resid_sd = gridsubset$resid_sd[row], seed = gridsubset$seed[row], variable_samples_per_subject = FALSE)

design=model.matrix(~Status + Batch_ID, simData$targets) 
design0=model.matrix(~Batch_ID, simData$targets)
colnames(design)[2]="Status"

##ORIGINAL SIMULATION CODE!
# list_res=list()
# for(method in c("be", "leek", "none", "known_nsv", "oracle")){ #
# 	#try: set.seed=(10000)
#     if(method == "none"){
#         n.sv = 0
#     }else if(method == "known_nsv"){
#         n.sv = nrow(simData$sv_df) #a priori input of numbers of sv = numbers of row in simData$sv_df
#     }else if(method == "oracle"){
#         n.sv = nrow(simData$sv_df)  #a priori input of numbers of sv = numbers of row in simData$sv_df
#     }else{ #if be or leek, est n.sv using be/leek method
#         n.sv = num.sv(simData$E,design,method=method)
#     }

#     if(n.sv > 0 & method != "oracle"){ #if method is not oracle or none (method = be/leek)
#         svobj = sva(simData$E,design,design0,n.sv=n.sv)
#         design_sv = cbind(design, svobj$sv) 
#     }else if(n.sv > 0 & method == "oracle"){
#         design_sv = model.matrix(as.formula(paste("~Status + Batch_ID",paste(simData$sv_df$sv_names, collapse=" + "), sep=" +")), simData$targets)
#         #as.formula changes strings to formula format (like as.factor)
#         #paste(simData$sv_df$sv_names, collapse=" + ") = "SV01 + SV02 +...+ SV20"
#         #adding BatchID and the above as covariates
#         colnames(design_sv)[2]="Status" #change column name from Statuscontrol to Status (control=1)
#     }else{ #if method is none or n.sv=0
#         design_sv = model.matrix(as.formula(paste("~Status + Batch_ID")), simData$targets)
#         colnames(design_sv)[2]="Status"
#     }
#     fit=lmFit(simData, design_sv) #fit linear model for each gene given a series of arrays
#     eBfit = eBayes(fit) #empirical bayes stats for DE (compute moderated t-statistics, moderated F-statistic, and log-odds of differential expression by empirical Bayes moderation of the standard errors towards a global value)
#     topSet = topTable(eBfit, coef="Status",number=nrow(eBfit)) #extract a table of the top-ranked genes from a linear model fit
#     vec = as.vector(table(estimated=factor(topSet$adj.P.Val<0.05,levels=c("FALSE","TRUE")),true=factor(topSet$DE,levels=c("FALSE","TRUE"))))
#     names(vec) = c("TN","FP","FN","TP")
#     # print("summarizing")
#     df_res = do.call(data.frame,as.list(vec)) #change vec to a dataframe
#     df_res$method=method #add a new column that specify method
#     df_res$nsv=n.sv #add a new column that specify nsv
#     df_res$param.n_sv=gridsubset$n_sv[row]
#     df_res$param.fraction_degs=gridsubset$fraction_degs[row]
#     df_res$param.fraction_sv=gridsubset$fraction_sv[row]
#     df_res$param.status_sd=gridsubset$status_sd[row]
#     df_res$param.sv_sd=gridsubset$sv_sd[row]
#     df_res$param.resid_sd=gridsubset$resid_sd[row]
#     df_res$param.seed=gridsubset$seed[row]
#     list_res[[method]]=df_res
# }
# out=do.call(rbind,list_res)
# out

# Initialize result storage
list_res = list()
list_pvals = list()

# Iterate over all methods
for (method in c("be", "leek", "none", "known_nsv", "oracle")) {
    if (method == "none") {
        n.sv = 0
    } else if (method %in% c("known_nsv", "oracle")) {
        n.sv = nrow(simData$sv_df)
    } else {
        n.sv = num.sv(simData$E, design, method=method)
    }
    
    if (n.sv > 0 & method != "oracle") {
        svobj = sva(simData$E, design, design0, n.sv=n.sv)
        design_sv = cbind(design, svobj$sv)
    } else if (n.sv > 0 & method == "oracle") {
        design_sv = model.matrix(as.formula(paste("~Status + Batch_ID", paste(simData$sv_df$sv_names, collapse=" + "), sep=" +")), simData$targets)
        colnames(design_sv)[2] = "Status"
    } else {
        design_sv = model.matrix(~Status + Batch_ID, simData$targets)
        colnames(design_sv)[2] = "Status"
    }
    
    # Fit linear model
    fit = lmFit(simData$E, design_sv)
    eBfit = eBayes(fit)
    
    # Extract differential expression results
    topSet = topTable(eBfit, coef="Status", number=nrow(eBfit), sort.by="none")
    
    # Ensure the 'DE' column exists (otherwise simulate it)
    if (!"DE" %in% colnames(topSet)) {
        topSet$DE = sample(c("TRUE", "FALSE"), nrow(topSet), replace = TRUE)
    }
    
    # Ensure lengths match before creating the table
    if (length(topSet$adj.P.Val) != length(topSet$DE)) {
        stop("Error: Mismatch in vector lengths!")
    }
    
    # Construct contingency table
    vec = as.vector(table(
        estimated = factor(topSet$adj.P.Val < 0.05, levels = c("FALSE", "TRUE")),
        true = factor(topSet$DE, levels = c("FALSE", "TRUE"))
    ))
    
    # Assign names
    names(vec) = c("TN", "FP", "FN", "TP")
    
    # Convert to data frame
    df_res = do.call(data.frame, as.list(vec))
    df_res$method = method
    df_res$nsv = n.sv
    df_res$param.n_sv = gridsubset$n_sv[row]
    df_res$param.fraction_degs = gridsubset$fraction_degs[row]
    df_res$param.fraction_sv = gridsubset$fraction_sv[row]
    df_res$param.status_sd = gridsubset$status_sd[row]
    df_res$param.sv_sd = gridsubset$sv_sd[row]
    df_res$param.resid_sd = gridsubset$resid_sd[row]
    df_res$param.seed = gridsubset$seed[row]
    
    # Save p-values separately
    df_pvals = data.frame(
        Method = method,
        GeneID = rownames(topSet),
        logFC = topSet$logFC,
        AveExpr = topSet$AveExpr,
        t = topSet$t,
        P.Value = topSet$P.Value,
        adj.P.Val = topSet$adj.P.Val,
        B = topSet$B
    )
    
    # Store results
    list_res[[method]] = df_res
    list_pvals[[method]] = df_pvals
}

# Combine results
out = do.call(rbind, list_res)
pval_out = do.call(rbind, list_pvals)

out
head(pval_out)
#start: 10:49pm
#end: 11:32pm

#save pval_out and out, if i dont include row.names = FALSE then there will be an X column for the row names
# write.csv(out, "/sc/arion/projects/mscic1/results/jolie/sva_sim_pvalue_distribution/case0_out.csv",row.names = FALSE)
# write.csv(pval_out, "/sc/arion/projects/mscic1/results/jolie/sva_sim_pvalue_distribution/case0_pval_out.csv",row.names = FALSE)

##plot plot plot
case0_pval_out<-read.csv("/sc/arion/projects/mscic1/results/jolie/sva_sim_pvalue_distribution/case0_pval_out.csv")

case1_pval_out<-read.csv("/sc/arion/projects/mscic1/results/jolie/sva_sim_pvalue_distribution/case1_pval_out.csv")

case2_pval_out<-read.csv("/sc/arion/projects/mscic1/results/jolie/sva_sim_pvalue_distribution/case2_pval_out.csv")

case3_pval_out<-read.csv("/sc/arion/projects/mscic1/results/jolie/sva_sim_pvalue_distribution/case3_pval_out.csv")

case4_pval_out<-read.csv("/sc/arion/projects/mscic1/results/jolie/sva_sim_pvalue_distribution/case4_pval_out.csv")

case2_degs<-case2_pval_out[case2_pval_out$P.Value <0.05,] 
table(case2_degs$Method)#2464     8 for oracle

################################################### CASE 1 ####################################################
##ALL METHODS
library(ggplot2)

# Filter out rows where P.Value is NA or zero to avoid log issues
filtered_data <- case1_pval_out %>% filter(!is.na(P.Value) & P.Value > 0)

# Determine x-axis and y-axis limits
x_min <- 0
x_max <- max(-log10(filtered_data$P.Value), na.rm = TRUE)

# Plot histogram colored by Method
case1<- ggplot(filtered_data, aes(x = -log10(P.Value), fill = Method)) +
  geom_histogram(bins = 30, boundary = 0, color = "black", alpha = 0.3, position = "identity") +
  labs(
    title = "Distribution of -log10(P-Values) by Method in Simulation Case 1",
    x = "-log10(P-Value)",
    y = "Count"
  ) +
  theme_minimal() +
  xlim(x_min, x_max) +
  scale_fill_manual(values = c("be" = "blue", 
                               "known_nsv" = "red", 
                               "leek" = "green", 
                               "none" = "purple", 
                               "oracle" = "orange"))

ggsave("/hpc/users/hoangd02/www/plots/p-value-distribution/sim_case1_all.pdf", plot = case1, width = 8, height = 5)

##Distribution of P-Values by Method in Simulation Cases

### KNOWN, ORACLE AND NONE
case0_sub <- filter(case0_pval_out, Method %in% c("known_nsv", "none", "oracle"))
case1_sub <- filter(case1_pval_out, Method %in% c("known_nsv", "none", "oracle"))
case2_sub <- filter(case2_pval_out, Method %in% c("known_nsv", "none", "oracle"))
case3_sub <- filter(case3_pval_out, Method %in% c("known_nsv", "none", "oracle"))
case4_sub <- filter(case4_pval_out, Method %in% c("known_nsv", "none", "oracle"))

# Filter out rows where P.Value is NA or zero to avoid log issues
filtered_data <- case0_sub %>% filter(!is.na(P.Value) & P.Value > 0)
filtered_data$Method <- factor(filtered_data$Method, levels = c("oracle", "none", "known_nsv"))

# Plot histogram colored by Method
case0 <- ggplot(filtered_data, aes(x = P.Value, fill = Method)) +
  geom_histogram(bins = 30, boundary = 0, color = "black", alpha = 0.3, position = "identity") +
  labs(
    title = "Distribution of P-Values by Method in One Simulation Case", #was Case 1
    x = "p-values",
    y = "Count"
  ) +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  xlim(0, 1) +  
  scale_fill_manual(values = c("known_nsv" = "red", 
                               "none" = "purple", 
                               "oracle" = "orange")) +
  facet_wrap(~Method, labeller = as_labeller(c(
    oracle = "Oracle",
    none = "No SVs",
    known_nsv = "Known Number of SVs"
  )), ncol = 1)

ggsave("/hpc/users/hoangd02/www/plots/p-value-distribution/sim_case0_known_oracle_none_02042025.pdf", plot = case0, width = 8, height = 5)



##ONE METHOD AT A TIME
case1_be <- filter(case1_pval_out, Method== "be")

#change 0 to NA
case1_be <- case1_be %>% mutate(P.Value = ifelse(P.Value == 0, NA, P.Value))

# Determine x-axis limits
x_min <- 0
x_max <- max(-log10(case1_be$P.Value), na.rm = TRUE)  # Removed duplicate vector

# Plot histogram
case1_be_plot <- ggplot(case1_be, aes(x = -log10(P.Value))) +
  geom_histogram(bins = 30, boundary = 0, color = "black", fill = "skyblue", na.rm = TRUE) +  # `na.rm = TRUE` removes NA
  labs(
    title = "Distribution of -log10(P-Values) in Simulation Case 1 with SVA 'BE'",
    x = "-log10(P-Value)",
    y = "Count"
  ) + theme_minimal() +
  xlim(x_min, x_max)  

ggsave("/hpc/users/hoangd02/www/plots/p-value-distribution/case1_be.pdf", plot = case1_be_plot, width = 8, height = 5)

################################################### CASE 2 ####################################################
filtered_data <- case2_pval_out %>% filter(!is.na(P.Value) & P.Value > 0)
x_min <- 0
x_max <- max(-log10(filtered_data$P.Value), na.rm = TRUE)

case2<- ggplot(filtered_data, aes(x = -log10(P.Value), fill = Method)) +
  geom_histogram(bins = 30, boundary = 0, color = "black", alpha = 0.3, position = "identity") +
  labs(
    title = "Distribution of -log10(P-Values) by Method in Simulation Case 2",
    x = "-log10(P-Value)",
    y = "Count"
  ) +
  theme_minimal() +
  xlim(x_min, x_max) +
  scale_fill_manual(values = c("be" = "blue", 
                               "known_nsv" = "red", 
                               "leek" = "green", 
                               "none" = "purple", 
                               "oracle" = "orange"))

ggsave("/hpc/users/hoangd02/www/plots/p-value-distribution/sim_case2_all.pdf", plot = case2, width = 8, height = 5)

### KNOWN, ORACLE AND NONE
case2_sub <- filter(case2_pval_out, Method %in% c("known_nsv", "none", "oracle"))
filtered_data <- case2_sub %>% filter(!is.na(P.Value) & P.Value > 0)

x_min <- 0
x_max <- max(-log10(filtered_data$P.Value), na.rm = TRUE)

case2<- ggplot(filtered_data, aes(x = -log10(P.Value), fill = Method)) +
  geom_histogram(bins = 30, boundary = 0, color = "black", alpha = 0.3, position = "identity") +
  labs(
    title = "Distribution of -log10(P-Values) by Method in Simulation Case 2",
    x = "-log10(P-Value)",
    y = "Count"
  ) +
  theme_minimal() +
  xlim(x_min, x_max) +
  scale_fill_manual(values = c("known_nsv" = "red", 
                               "none" = "purple", 
                               "oracle" = "orange"))

ggsave("/hpc/users/hoangd02/www/plots/p-value-distribution/sim_case2_known_oracle_none.pdf", plot = case2, width = 8, height = 5)


################################################### CASE 3 ####################################################
filtered_data <- case3_pval_out %>% filter(!is.na(P.Value) & P.Value > 0)
x_min <- 0
x_max <- max(-log10(filtered_data$P.Value), na.rm = TRUE)

case3<- ggplot(filtered_data, aes(x = -log10(P.Value), fill = Method)) +
  geom_histogram(bins = 30, boundary = 0, color = "black", alpha = 0.3, position = "identity") +
  labs(
    title = "Distribution of -log10(P-Values) by Method in Simulation Case 3",
    x = "-log10(P-Value)",
    y = "Count"
  ) +
  theme_minimal() +
  xlim(x_min, x_max) +
  scale_fill_manual(values = c("be" = "blue", 
                               "known_nsv" = "red", 
                               "leek" = "green", 
                               "none" = "purple", 
                               "oracle" = "orange"))

ggsave("/hpc/users/hoangd02/www/plots/p-value-distribution/sim_case3_all.pdf", plot = case3, width = 8, height = 5)

### KNOWN, ORACLE AND NONE
case3_sub <- filter(case3_pval_out, Method %in% c("known_nsv", "none", "oracle"))
filtered_data <- case3_sub %>% filter(!is.na(P.Value) & P.Value > 0)

x_min <- 0
x_max <- max(-log10(filtered_data$P.Value), na.rm = TRUE)

case3<- ggplot(filtered_data, aes(x = -log10(P.Value), fill = Method)) +
  geom_histogram(bins = 30, boundary = 0, color = "black", alpha = 0.3, position = "identity") +
  labs(
    title = "Distribution of -log10(P-Values) by Method in Simulation Case 3",
    x = "-log10(P-Value)",
    y = "Count"
  ) +
  theme_minimal() +
  xlim(x_min, x_max) +
  scale_fill_manual(values = c("known_nsv" = "red", 
                               "none" = "purple", 
                               "oracle" = "orange"))

ggsave("/hpc/users/hoangd02/www/plots/p-value-distribution/sim_case3_known_oracle_none.pdf", plot = case3, width = 8, height = 5)

################################################### CASE 4 ####################################################
filtered_data <- case4_pval_out %>% filter(!is.na(P.Value) & P.Value > 0)
x_min <- 0
x_max <- max(-log10(filtered_data$P.Value), na.rm = TRUE)

case4 <- ggplot(filtered_data, aes(x = -log10(P.Value), fill = Method)) +
  geom_histogram(bins = 30, boundary = 0, color = "black", alpha = 0.3, position = "identity") +
  labs(
    title = "Distribution of -log10(P-Values) by Method in Simulation Case 4",
    x = "-log10(P-Value)",
    y = "Count"
  ) +
  theme_minimal() +
  xlim(x_min, x_max) +
  scale_fill_manual(values = c("be" = "blue", 
                               "known_nsv" = "red", 
                               "leek" = "green", 
                               "none" = "purple", 
                               "oracle" = "orange"))

ggsave("/hpc/users/hoangd02/www/plots/p-value-distribution/sim_case4_all.pdf", plot = case4, width = 8, height = 5)

### KNOWN, ORACLE AND NONE
case4_sub <- filter(case4_pval_out, Method %in% c("known_nsv", "none", "oracle"))
filtered_data <- case4_sub %>% filter(!is.na(P.Value) & P.Value > 0)

x_min <- 0
x_max <- max(-log10(filtered_data$P.Value), na.rm = TRUE)

case4 <- ggplot(filtered_data, aes(x = -log10(P.Value), fill = Method)) +
  geom_histogram(bins = 30, boundary = 0, color = "black", alpha = 0.3, position = "identity") +
  labs(
    title = "Distribution of -log10(P-Values) by Method in Simulation Case 4",
    x = "-log10(P-Value)",
    y = "Count"
  ) +
  theme_minimal() +
  xlim(x_min, x_max) +
  scale_fill_manual(values = c("known_nsv" = "red", 
                               "none" = "purple", 
                               "oracle" = "orange"))

ggsave("/hpc/users/hoangd02/www/plots/p-value-distribution/sim_case4_known_oracle_none.pdf", plot = case4, width = 8, height = 5)







































