##from /sc/arion/projects/mscic1/results/jolie/sva_sim/scripts/sva_simulation_analysis_parallel.R
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
        n_sv <- 1
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

##tried a problematic set: set10[[1]] three times with different total_sv_sd
row=1
n_sv = 0
total_sv_sd = 5 ##try 0, 5, and 10
grid_sv = as.data.frame(expand.grid(n_sv=n_sv,total_sv_sd=total_sv_sd))
grid_sv$sv_sd = sqrt(grid_sv$total_sv_sd^2 / grid_sv$n_sv)
grid_sv$sv_sd[!is.finite(grid_sv$sv_sd)] = 0
grid_list <- list(grid_sv, 
                   data.frame(fraction_degs = 0), 
                   data.frame(status_sd = 0.2),
                   data.frame(resid_sd = 0.4), 
                   data.frame(fraction_sv = 0.9), 
                   #data.frame(seed = 20240))
                   data.frame(seed = 12350))


##tried a problematic set: set10[[4]]
row=1
n_sv = 0
total_sv_sd = 5 #try 0 and 5
grid_sv = as.data.frame(expand.grid(n_sv=n_sv,total_sv_sd=total_sv_sd))
grid_sv$sv_sd = sqrt(grid_sv$total_sv_sd^2 / grid_sv$n_sv)
grid_sv$sv_sd[!is.finite(grid_sv$sv_sd)] = 0
grid_list <- list(grid_sv, 
                   data.frame(fraction_degs = 0.04), 
                   data.frame(status_sd = 0.1),
                   data.frame(resid_sd = 0.3), 
                   data.frame(fraction_sv = 0.9), 
                   data.frame(seed = 20240))
                   #data.frame(seed = 12350))

##tried a problematic set: set6[[3000]]; set6_3000_grid
row=1
n_sv = 90
total_sv_sd = 0.5 
grid_sv = as.data.frame(expand.grid(n_sv=n_sv,total_sv_sd=total_sv_sd))
grid_sv$sv_sd = sqrt(grid_sv$total_sv_sd^2 / grid_sv$n_sv)
grid_sv$sv_sd[!is.finite(grid_sv$sv_sd)] = 0
grid_list <- list(grid_sv, 
                   data.frame(fraction_degs = 0.4), 
                   data.frame(status_sd = 0.5),
                   data.frame(resid_sd = 0.2), 
                   data.frame(fraction_sv = 0.3),
                   data.frame(seed = 12345))
                   #data.frame(seed = 20240)) 
                   #data.frame(seed = 99876)) 

##tried a problematic set: set7[[720]]; set7_720_grid
row=1
n_sv = 100
total_sv_sd = 9
grid_sv = as.data.frame(expand.grid(n_sv=n_sv,total_sv_sd=total_sv_sd))
grid_sv$sv_sd = sqrt(grid_sv$total_sv_sd^2 / grid_sv$n_sv)
grid_sv$sv_sd[!is.finite(grid_sv$sv_sd)] = 0
grid_list <- list(grid_sv,
                   data.frame(fraction_degs = 0.1),
                   data.frame(status_sd = 0.5),
                   data.frame(resid_sd = 0.4),
                   data.frame(fraction_sv = 0.9),
                   data.frame(seed = 12346))
gridsubset=as.data.frame(t(as.data.frame(unlist(grid_list))))

#be nsv = 0 when seed = 12345, 12347, 12351, 12354
#be nsv = 1 when seed = 12352

#test[1,]
row=1
n_sv = 0 #0 to 100
total_sv_sd = 0
grid_sv = as.data.frame(expand.grid(n_sv=n_sv,total_sv_sd=total_sv_sd))
grid_sv$sv_sd = sqrt(grid_sv$total_sv_sd^2 / grid_sv$n_sv)
grid_sv$sv_sd[!is.finite(grid_sv$sv_sd)] = 0
grid_list <- list(grid_sv,
                  data.frame(fraction_degs = 0),
                  data.frame(status_sd = 0.1),
                  data.frame(resid_sd = 0.1),
                  data.frame(fraction_sv = 0.3),
                  data.frame(seed = 12352)) #12345,12347,12351,12354,12352
gridsubset=as.data.frame(t(as.data.frame(unlist(grid_list))))

simData = sim_case_control_data(n_subjects = 1000, n_samples_per_subject = 1, n_features = 20000, n_batches = 15, 
	n_sv = gridsubset$n_sv[row], fraction_degs = gridsubset$fraction_degs[row], fraction_sv = gridsubset$fraction_sv[row], 
    status_sd = gridsubset$status_sd[row], subject_sd = 0, batch_sd = 0.3, sv_sd = gridsubset$sv_sd[row], 
    resid_sd = gridsubset$resid_sd[row], seed = gridsubset$seed[row], variable_samples_per_subject = FALSE)

design=model.matrix(~Status + Batch_ID, simData$targets) 
design0=model.matrix(~Batch_ID, simData$targets)
colnames(design)[2]="Status"

list_res=list()
for(method in c("be", "leek", "none", "known_nsv", "oracle")){ #
	#try: set.seed=(10000)
    if(method == "none"){
        n.sv = 0
    }else if(method == "known_nsv"){
        n.sv = nrow(simData$sv_df) #a priori input of numbers of sv = numbers of row in simData$sv_df
    }else if(method == "oracle"){
        n.sv = nrow(simData$sv_df)  #a priori input of numbers of sv = numbers of row in simData$sv_df
    }else{ #if be or leek, est n.sv using be/leek method
        n.sv = num.sv(simData$E,design,method=method)
    }

    if(n.sv > 0 & method != "oracle"){ #if method is not oracle or none (method = be/leek)
        svobj = sva(simData$E,design,design0,n.sv=n.sv)
        design_sv = cbind(design, svobj$sv) 
    }else if(n.sv > 0 & method == "oracle"){
        design_sv = model.matrix(as.formula(paste("~Status + Batch_ID",paste(simData$sv_df$sv_names, collapse=" + "), sep=" +")), simData$targets)
        #as.formula changes strings to formula format (like as.factor)
        #paste(simData$sv_df$sv_names, collapse=" + ") = "SV01 + SV02 +...+ SV20"
        #adding BatchID and the above as covariates
        colnames(design_sv)[2]="Status" #change column name from Statuscontrol to Status (control=1)
    }else{ #if method is none or n.sv=0
        design_sv = model.matrix(as.formula(paste("~Status + Batch_ID")), simData$targets)
        colnames(design_sv)[2]="Status"
    }
    fit=lmFit(simData, design_sv) #fit linear model for each gene given a series of arrays
    eBfit = eBayes(fit) #empirical bayes stats for DE (compute moderated t-statistics, moderated F-statistic, and log-odds of differential expression by empirical Bayes moderation of the standard errors towards a global value)
    topSet = topTable(eBfit, coef="Status",number=nrow(eBfit)) #extract a table of the top-ranked genes from a linear model fit
    vec = as.vector(table(estimated=factor(topSet$adj.P.Val<0.05,levels=c("FALSE","TRUE")),true=factor(topSet$DE,levels=c("FALSE","TRUE"))))
    names(vec) = c("TN","FP","FN","TP")
    # print("summarizing")
    df_res = do.call(data.frame,as.list(vec)) #change vec to a dataframe
    df_res$method=method #add a new column that specify method
    df_res$nsv=n.sv #add a new column that specify nsv
    df_res$param.n_sv=gridsubset$n_sv[row]
    df_res$param.fraction_degs=gridsubset$fraction_degs[row]
    df_res$param.fraction_sv=gridsubset$fraction_sv[row]
    df_res$param.status_sd=gridsubset$status_sd[row]
    df_res$param.sv_sd=gridsubset$sv_sd[row]
    df_res$param.resid_sd=gridsubset$resid_sd[row]
    df_res$param.seed=gridsubset$seed[row]
    list_res[[method]]=df_res
}
out=do.call(rbind,list_res)
out

###IT'S A TOTAL_SV_SD PROBLEM
#running it once at a time seems fine: out of set10[[1]] & total_sv_sd = 0 
#Number of significant surrogate variables is:  1 
#Iteration (out of 5 ):1  2  3  4  5  Coefficients not estimable: SV1 
#Warning messages:
#1: Partial NA coefficients for 20000 probe(s) 
#2: In .ebayes(fit = fit, proportion = proportion, stdev.coef.lim = stdev.coef.lim,  :
#Estimation of var.prior failed - set to default value
             TN FP FN TP    method nsv param.n_sv param.fraction_degs
be        20000  0  0  0        be   0          0                   0
leek      20000  0  0  0      leek   0          0                   0
none      20000  0  0  0      none   0          0                   0
known_nsv 20000  0  0  0 known_nsv   1          0                   0
oracle    20000  0  0  0    oracle   1          0                   0
          param.fraction_sv param.status_sd param.sv_sd param.resid_sd
be                      0.9             0.2           0            0.4
leek                    0.9             0.2           0            0.4
none                    0.9             0.2           0            0.4
known_nsv               0.9             0.2           0            0.4
oracle                  0.9             0.2           0            0.4
          param.seed
be             12350
leek           12350
none           12350
known_nsv      12350
oracle         12350

##same result as a different seed: 20240
be        20000  0  0  0        be   0          0                   0
leek      20000  0  0  0      leek   0          0                   0
none      20000  0  0  0      none   0          0                   0
known_nsv 20000  0  0  0 known_nsv   1          0                   0
oracle    20000  0  0  0    oracle   1          0                   0
          param.fraction_sv param.status_sd param.sv_sd param.resid_sd
be                      0.9             0.2           0            0.4
leek                    0.9             0.2           0            0.4
none                    0.9             0.2           0            0.4
known_nsv               0.9             0.2           0            0.4
oracle                  0.9             0.2           0            0.4
          param.seed
be             20240
leek           20240
none           20240
known_nsv      20240
oracle         20240

#out of set10[[1]] & total_sv_sd = 5 
#Warning message (occurs after initializing simData)
#In (function (n_subjects, n_samples_per_subject, n_features, n_batches,  :
# Setting SV signal to zero since there are not SVs
             TN FP FN TP    method nsv param.n_sv param.fraction_degs
be        20000  0  0  0        be  10          0                   0
leek      20000  0  0  0      leek   0          0                   0
none      20000  0  0  0      none   0          0                   0
known_nsv 20000  0  0  0 known_nsv   1          0                   0
oracle    20000  0  0  0    oracle   1          0                   0
          param.fraction_sv param.status_sd param.sv_sd param.resid_sd
be                      0.9             0.2           0            0.4
leek                    0.9             0.2           0            0.4
none                    0.9             0.2           0            0.4
known_nsv               0.9             0.2           0            0.4
oracle                  0.9             0.2           0            0.4
          param.seed
be             12350
leek           12350
none           12350
known_nsv      12350
oracle         12350

#out of set10[[1]] & total_sv_sd = 10 
#Number of significant surrogate variables is:  9 
#Iteration (out of 5 ):1  2  3  4  5  Number of significant surrogate variables is:  1 
#Iteration (out of 5 ):1  2  3  4  5  Coefficients not estimable: SV1 
#Warning messages:
#1: Partial NA coefficients for 20000 probe(s) 
#2: In .ebayes(fit = fit, proportion = proportion, stdev.coef.lim = stdev.coef.lim,  :
#Estimation of var.prior failed - set to default value
             TN FP FN TP    method nsv param.n_sv param.fraction_degs
be        20000  0  0  0        be   9          0                   0
leek      20000  0  0  0      leek   0          0                   0
none      20000  0  0  0      none   0          0                   0
known_nsv 20000  0  0  0 known_nsv   1          0                   0
oracle    20000  0  0  0    oracle   1          0                   0
          param.fraction_sv param.status_sd param.sv_sd param.resid_sd
be                      0.9             0.2           0            0.4
leek                    0.9             0.2           0            0.4
none                    0.9             0.2           0            0.4
known_nsv               0.9             0.2           0            0.4
oracle                  0.9             0.2           0            0.4
          param.seed
be             12350
leek           12350
none           12350
known_nsv      12350
oracle         12350

#running it once at a time seems fine: out of set10[[4]] & total_sv_sd=0
             TN FP  FN  TP    method nsv param.n_sv param.fraction_degs
be        19169 31 257 543        be   9          0                0.04
leek      19168 32 256 544      leek   0          0                0.04
none      19168 32 256 544      none   0          0                0.04
known_nsv 19168 32 257 543 known_nsv   1          0                0.04
oracle    19168 32 256 544    oracle   1          0                0.04
          param.fraction_sv param.status_sd param.sv_sd param.resid_sd
be                      0.9             0.1           0            0.3
leek                    0.9             0.1           0            0.3
none                    0.9             0.1           0            0.3
known_nsv               0.9             0.1           0            0.3
oracle                  0.9             0.1           0            0.3
          param.seed
be             12350
leek           12350
none           12350
known_nsv      12350
oracle         12350


#running it once at a time seems fine: out of set10[[4]] & total_sv_sd=5
#Warning message:
#In (function (n_subjects, n_samples_per_subject, n_features, n_batches,  :
#  Setting SV signal to zero since there are not SVs

#Number of significant surrogate variables is:  1 
#Iteration (out of 5 ):1  2  3  4  5  Coefficients not estimable: SV1 
#Warning messages:
#1: Partial NA coefficients for 20000 probe(s) 
#2: In .ebayes(fit = fit, proportion = proportion, stdev.coef.lim = stdev.coef.lim,  :
#Estimation of var.prior failed - set to default value
             TN FP  FN  TP    method nsv param.n_sv param.fraction_degs
be        19168 32 256 544        be   0          0                0.04
leek      19168 32 256 544      leek   0          0                0.04
none      19168 32 256 544      none   0          0                0.04
known_nsv 19168 32 257 543 known_nsv   1          0                0.04
oracle    19168 32 256 544    oracle   1          0                0.04
          param.fraction_sv param.status_sd param.sv_sd param.resid_sd
be                      0.9             0.1           0            0.3
leek                    0.9             0.1           0            0.3
none                    0.9             0.1           0            0.3
known_nsv               0.9             0.1           0            0.3
oracle                  0.9             0.1           0            0.3
          param.seed
be             12350
leek           12350
none           12350
known_nsv      12350
oracle         12350

##change seed and it gives a different TN, FP, etc.
             TN FP  FN  TP    method nsv param.n_sv param.fraction_degs
be        19171 29 286 514        be   0          0                0.04
leek      19171 29 286 514      leek   0          0                0.04
none      19171 29 286 514      none   0          0                0.04
known_nsv 19170 30 287 513 known_nsv   1          0                0.04
oracle    19171 29 286 514    oracle   1          0                0.04
          param.fraction_sv param.status_sd param.sv_sd param.resid_sd
be                      0.9             0.1           0            0.3
leek                    0.9             0.1           0            0.3
none                    0.9             0.1           0            0.3
known_nsv               0.9             0.1           0            0.3
oracle                  0.9             0.1           0            0.3
          param.seed
be             20240
leek           20240
none           20240
known_nsv      20240
oracle         20240


############PROBLEM 4: SAME SET OF PARAMETERS (SAME SEED) BUT (at least) 3 DIFFERENT RESULTS
#out of set6[[3000]] & total_sv_sd = 0.5 (have only one total_sv_sd), seed = 12345
#Number of significant surrogate variables is:  111 
#Iteration (out of 5 ):1  2  3  4  5  Number of significant surrogate variables is:  90 
#Iteration (out of 5 ):1  2  3  4  5  
             TN  FP  FN   TP    method nsv param.n_sv param.fraction_degs
be        11735 265 278 7722        be 111         90                 0.4
leek      11756 244 262 7738      leek   0         90                 0.4
none      11756 244 262 7738      none   0         90                 0.4
known_nsv 11741 259 278 7722 known_nsv  90         90                 0.4
oracle    11765 235 273 7727    oracle  90         90                 0.4
          param.fraction_sv param.status_sd param.sv_sd param.resid_sd
be                      0.3             0.5  0.05270463            0.2
leek                    0.3             0.5  0.05270463            0.2
none                    0.3             0.5  0.05270463            0.2
known_nsv               0.3             0.5  0.05270463            0.2
oracle                  0.3             0.5  0.05270463            0.2
          param.seed
be             12345
leek           12345
none           12345
known_nsv      12345
oracle         12345

#out of set6[[3000]] & seed = 12345, try 2 -- STILL DIFFERENT ANSWER
             TN  FP  FN   TP    method nsv param.n_sv param.fraction_degs
be        11749 251 277 7723        be 101         90                 0.4
leek      11756 244 262 7738      leek   0         90                 0.4
none      11756 244 262 7738      none   0         90                 0.4
known_nsv 11741 259 278 7722 known_nsv  90         90                 0.4
oracle    11765 235 273 7727    oracle  90         90                 0.4
          param.fraction_sv param.status_sd param.sv_sd param.resid_sd
be                      0.3             0.5  0.05270463            0.2
leek                    0.3             0.5  0.05270463            0.2
none                    0.3             0.5  0.05270463            0.2
known_nsv               0.3             0.5  0.05270463            0.2
oracle                  0.3             0.5  0.05270463            0.2
          param.seed
be             12345
leek           12345
none           12345
known_nsv      12345
oracle         12345

##this is very problematic because this yields a different of nsv=111 compared to nsv=112 in set6[[3000]]
set6[[3000]]
     TN  FP  FN   TP    method nsv param.n_sv param.fraction_degs
1 11747 253 274 7726        be  79         90                 0.4
2 11756 244 262 7738      leek   0         90                 0.4
3 11756 244 262 7738      none   0         90                 0.4
4 11741 259 278 7722 known_nsv  90         90                 0.4
5 11765 235 273 7727    oracle  90         90                 0.4
6 11740 260 273 7727        be 112         90                 0.4
  param.fraction_sv param.status_sd param.sv_sd param.resid_sd param.seed
1               0.3             0.5  0.05270463            0.2      12345
2               0.3             0.5  0.05270463            0.2      12345
3               0.3             0.5  0.05270463            0.2      12345
4               0.3             0.5  0.05270463            0.2      12345
5               0.3             0.5  0.05270463            0.2      12345
6               0.3             0.5  0.05270463            0.2      12345

#out of set6[[3000]] & seed = 20240 to see if the outcome is deterministic, it should be, but it's not - SO CONFUSED, and it gives a different nsv
             TN  FP  FN   TP    method nsv param.n_sv param.fraction_degs
be        11750 250 270 7730        be 102         90                 0.4
leek      11770 230 266 7734      leek   0         90                 0.4
none      11770 230 266 7734      none   0         90                 0.4
known_nsv 11761 239 271 7729 known_nsv  90         90                 0.4
oracle    11757 243 289 7711    oracle  90         90                 0.4
          param.fraction_sv param.status_sd param.sv_sd param.resid_sd
be                      0.3             0.5  0.05270463            0.2
leek                    0.3             0.5  0.05270463            0.2
none                    0.3             0.5  0.05270463            0.2
known_nsv               0.3             0.5  0.05270463            0.2
oracle                  0.3             0.5  0.05270463            0.2
          param.seed
be             20240
leek           20240
none           20240
known_nsv      20240
oracle         20240

#out of set6[[3000]] & seed = 99876
             TN  FP  FN   TP    method nsv param.n_sv param.fraction_degs
be        11731 269 263 7737        be  61         90                 0.4
leek      11733 267 264 7736      leek   0         90                 0.4
none      11733 267 264 7736      none   0         90                 0.4
known_nsv 11726 274 274 7726 known_nsv  90         90                 0.4
oracle    11754 246 281 7719    oracle  90         90                 0.4
          param.fraction_sv param.status_sd param.sv_sd param.resid_sd
be                      0.3             0.5  0.05270463            0.2
leek                    0.3             0.5  0.05270463            0.2
none                    0.3             0.5  0.05270463            0.2
known_nsv               0.3             0.5  0.05270463            0.2
oracle                  0.3             0.5  0.05270463            0.2
          param.seed
be             99876
leek           99876
none           99876
known_nsv      99876
oracle         99876
#random-number-dependent

#out of set7[[720]] manually ran 10/13 vs  set7[[720]]
             TN    FP   FN   TP    method nsv param.n_sv param.fraction_degs
be         4014 13986  264 1736        be  83        100                 0.1
leek       1901 16099  108 1892      leek  99        100                 0.1
none      17997     3 1942   58      none   0        100                 0.1
known_nsv  1933 16067  111 1889 known_nsv 100        100                 0.1
oracle    17920    80  164 1836    oracle 100        100                 0.1
          param.fraction_sv param.status_sd param.sv_sd param.resid_sd
be                      0.9             0.5         0.9            0.4
leek                    0.9             0.5         0.9            0.4
none                    0.9             0.5         0.9            0.4
known_nsv               0.9             0.5         0.9            0.4
oracle                  0.9             0.5         0.9            0.4
          param.seed
be             12346
leek           12346
none           12346
known_nsv      12346
oracle         12346

##set7[[720]] from res
     TN    FP   FN   TP    method nsv param.n_sv param.fraction_degs
1  4014 13986  264 1736        be  83        100                 0.1
2  1932 16068  111 1889      leek 100        100                 0.1
3 17997     3 1942   58      none   0        100                 0.1
4  1932 16068  111 1889 known_nsv 100        100                 0.1
5 17920    80  164 1836    oracle 100        100                 0.1
6  1901 16099  108 1892      leek  99        100                 0.1
7  1933 16067  111 1889 known_nsv 100        100                 0.1
  param.fraction_sv param.status_sd param.sv_sd param.resid_sd param.seed
1               0.9             0.5         0.9            0.4      12346
2               0.9             0.5         0.9            0.4      12346
3               0.9             0.5         0.9            0.4      12346
4               0.9             0.5         0.9            0.4      12346
5               0.9             0.5         0.9            0.4      12346
6               0.9             0.5         0.9            0.4      12346
7               0.9             0.5         0.9            0.4      12346











