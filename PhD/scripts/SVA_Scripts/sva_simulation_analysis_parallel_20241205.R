#for job submission
args <- commandArgs(trailingOnly=TRUE) #take arguments 
start = as.numeric(args[[1]]) #start line of grid
end = as.numeric(args[[2]])
Sys.setenv(OMP_THREAD_LIMIT = 1, OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, 
    MKL_NUM_THREADS = 1, GOTO_NUM_THREADS = 1, RCPP_PARALLEL_NUM_THREADS= 1)

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

# n_sv <- 10
# sv_sd <- 0.5
# target_var <- n_sv * sv_sd^2
# message("Target total variance: ", target_var)

# set.seed(1986)
# sv_sd_vec_approx <- sqrt(rexp(n = n_sv, rate = 1/sv_sd^2))
# message("Approx: ")
# print(sv_sd_vec_approx)
# message("Total variance: ", sum(sv_sd_vec_approx^2))

# set.seed(1986)
# sv_var_vec_exact <- rexp(n = n_sv)
# sv_var_vec_exact <- sv_var_vec_exact / mean(sv_var_vec_exact) * sv_sd^2
# sv_sd_vec_exact <- sqrt(sv_var_vec_exact)
# message("Exact: ")
# print(sv_sd_vec_exact)
# message("Total variance: ", sum(sv_sd_vec_exact^2))

   # ## Example values for testing
   #  local({
   #      n_subjects = 200
   #      n_samples_per_subject = 1
   #      n_features = 10000
   #      n_batches = 15
   #      n_sv = 20
   #      fraction_degs = 0.3
   #      status_sd = 0.1
   #      subject_sd = 0.1
   #      batch_sd = 0.1
   #      sv_sd = 0.1
   #      resid_sd = 0.1
   #      seed = 1986
   #      variable_samples_per_subject = TRUE
   #      .max_tries = 500
   #  })

n_subjects = 1000; n_samples_per_subject = 1; n_features = 20000; n_batches = 15; n_sv = 20; fraction_degs = 0.3; fraction_sv = 0.3; status_sd = 0.1; subject_sd = 0; batch_sd = 0.3; sv_sd = 1.58113883; resid_sd = 0.1; seed = 1986; variable_samples_per_subject = FALSE

sim_case_control_data <- function (
    n_subjects,
    n_samples_per_subject,
    n_features,
    n_batches,
    n_sv,
    fraction_degs,
    fraction_sv, #how to distribute this to each EIVs 
    #add fraction_de_sv_overlap (check if this is less or equal to fraction_degs 
    #AND fraction_sv; greater than fraction_degs +fraction_sv -1; min = 0 if they're disjoint, add error message)
    #default = NULL --> it happens at random (current code)
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
    n_sv >= 0, ##changed on 2024-12-05
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
    warning("Setting SV signal to zero since there are no SVs")
    n_sv <- 0  ##changed on 2024-12-05
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
      Status = sample(rep(c("case","control"),length.out = n())) #500 cases and 500 controls
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


# defining grid on which we want to run the analysis
n_sv = seq(from = 0,to = 100, by=10)
total_sv_sd = seq(from = 0,to = 10, by=0.5)
grid_sv = as.data.frame(expand.grid(n_sv=n_sv,total_sv_sd=total_sv_sd))
grid_sv$sv_sd = sqrt(grid_sv$total_sv_sd^2 / grid_sv$n_sv)
grid_sv$sv_sd[!is.finite(grid_sv$sv_sd)] = 0
grid_list <- list( grid_sv, 
    data.frame(fraction_degs = c(seq(from = 0,to = 0.099, by=0.02),seq(from = 0.1,to = 0.5, by=0.1))), 
    data.frame(status_sd = seq(from = 0.1,to = 0.5, by=0.1)), 
    data.frame(resid_sd = seq(from = 0.1,to = 0.5, by=0.1)), 
    data.frame(fraction_sv = seq(from = 0.1,to = 0.9, by=0.2)), 
    data.frame(seed = seq(from = 12345, length.out = 10))
    )

grid_table = Reduce(f = dplyr::cross_join, x = grid_list)
dim(grid_table)#2887500 x 8
# grid_table is the grid

library(foreach)
library(doParallel)
registerDoParallel(cores=10) 
gridsubset = grid_table[start:end,] #start and end are external inputs
#run 
res_final=foreach(row=1:nrow(gridsubset)) %dopar%{ 
    # running the simulation
    # print("pre sim")
    print(paste("##########################\n##########################\n##########################\n##########################\n##########################\n##########################\n##########################\n##########################\n##########################\n##########################\n##########################\n##########################\n",row))
    simData = sim_case_control_data(n_subjects = 1000, n_samples_per_subject = 1, n_features = 20000, n_batches = 15, n_sv = gridsubset$n_sv[row], 
                                    fraction_degs = gridsubset$fraction_degs[row], fraction_sv = gridsubset$fraction_sv[row], status_sd = gridsubset$status_sd[row], 
                                    subject_sd = 0, batch_sd = 0.3, sv_sd = gridsubset$sv_sd[row], 
                                    resid_sd = gridsubset$resid_sd[row], seed = gridsubset$seed[row], variable_samples_per_subject = FALSE)
    # print("post sim")
    # design=model.matrix(as.formula(paste("~Status + Batch_ID",paste(simData$sv_df$sv_names, collapse=" + "), sep=" +")), simData$targets)
    # from limma package, running DE on the simulated data example
    design=model.matrix(~Status + Batch_ID, simData$targets) #create a design matrix, or a matrix of values that represents the variables in a regression model in a format suitable for linear modeling analysis
    #only have predictor variables, no outcome variables
    design0=model.matrix(~Batch_ID, simData$targets)
    colnames(design)[2]="Status"
    # fit=lmFit(simData, design)
    # eBfit = eBayes(fit)
    # topSet=topTable(eBfit, coef="Status",number=nrow(eBfit))
    # table(estimated=topSet$adj.P.Val<0.05,true=topSet$DE)
    #   
    # actual analysis - 1 iteration: lines 377-410
    # print("main loop")
    list_res=list()
    for(method in c("be", "leek", "none", "known_nsv", "oracle")){ #
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
}

final = do.call(rbind,res_final)
fwrite(final,file=paste0("results_of_sva_simulation_start_",start,"_end_",end,".txt"),row.names = FALSE,col.names = TRUE,sep="\t",quote = FALSE)
# design=model.matrix(as.formula(paste("~Status + Batch_ID",paste(simData$sv_df$sv_names, collapse=" + "), sep=" +")), simData$targets)
# colnames(design)[2]="Status"
# fit=lmFit(simData, design)
# eBfit = eBayes(fit)
# topSet=topTable(eBfit, coef="Status",number=nrow(eBfit))
# table(estimated=topSet$adj.P.Val<0.05,true=topSet$DE)

