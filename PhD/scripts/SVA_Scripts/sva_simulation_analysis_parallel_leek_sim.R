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

#n_subjects = 1000; n_samples_per_subject = 1; n_features = 20000; n_batches = 15; n_sv = 20; fraction_degs = 0.3; fraction_sv = 0.3; status_sd = 0.1; subject_sd = 0; batch_sd = 0.3; sv_sd = 1.58113883; resid_sd = 0.1; seed = 1986; variable_samples_per_subject = FALSE

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

# Define grid -- weird error so lets just try it once at a time
n_sv = seq(from = 1, to = 10, by = 1)  # Generates 11 values: 0 to 10
total_sv_sd = 2
grid_sv = expand.grid(n_sv = n_sv, total_sv_sd = total_sv_sd)
grid_sv$sv_sd = sqrt(grid_sv$total_sv_sd^2 / grid_sv$n_sv)
grid_sv$sv_sd[!is.finite(grid_sv$sv_sd)] = 0  # Handle division by zero
grid_list <- list(grid_sv,
                  data.frame(fraction_degs = 0.3),
                  data.frame(status_sd = 0.3),
                  data.frame(resid_sd = 0.3),
                  data.frame(fraction_sv = 0.5),
                  data.frame(seed = 12348))
gridsubset = Reduce(f = dplyr::cross_join, x = grid_list)

## TESTING IT ONE SIM AT A TIME 
row=1
n_sv = 1 #change this from 1 to 10, run either 0 or 10 and compare to make sure the output is correct
total_sv_sd = 2
grid_sv = as.data.frame(expand.grid(n_sv=n_sv,total_sv_sd=total_sv_sd))
grid_sv$sv_sd = sqrt(grid_sv$total_sv_sd^2 / grid_sv$n_sv)
grid_sv$sv_sd[!is.finite(grid_sv$sv_sd)] = 0
grid_list <- list( grid_sv, 
                   data.frame(fraction_degs = 0.3), 
                   data.frame(status_sd = 0.3), 
                   data.frame(resid_sd = 0.3), 
                   data.frame(fraction_sv = 0.5), 
                   data.frame(seed = 12348))
gridsubset=as.data.frame(t(as.data.frame(unlist(grid_list))))
simData = sim_case_control_data(n_subjects = 1000, n_samples_per_subject = 1, n_features = 20000, n_batches = 15, n_sv = gridsubset$n_sv[row], 
                                fraction_degs = gridsubset$fraction_degs[row], fraction_sv = gridsubset$fraction_sv[row], status_sd = gridsubset$status_sd[row], 
                                subject_sd = 0, batch_sd = 0.3, sv_sd = gridsubset$sv_sd[row], 
                                resid_sd = gridsubset$resid_sd[row], seed = gridsubset$seed[row], variable_samples_per_subject = FALSE)

###NEW 03/21/2025
## TESTING IT ONE SIM AT A TIME 
row=1
n_sv = 1 #change this from 1 to 10, run either 0 or 10 and compare to make sure the output is correct
total_sv_sd = 2
grid_sv = as.data.frame(expand.grid(n_sv=n_sv,total_sv_sd=total_sv_sd))
grid_sv$sv_sd = sqrt(grid_sv$total_sv_sd^2 / grid_sv$n_sv)
grid_sv$sv_sd[!is.finite(grid_sv$sv_sd)] = 0
grid_list <- list( grid_sv, 
                   data.frame(fraction_degs = 0.4), 
                   data.frame(status_sd = 0.3), 
                   data.frame(resid_sd = 0.3), 
                   data.frame(fraction_sv = 0.4), 
                   data.frame(seed = 12350))
gridsubset=as.data.frame(t(as.data.frame(unlist(grid_list))))
simData = sim_case_control_data(n_subjects = 1000, n_samples_per_subject = 1, n_features = 20000, n_batches = 1, n_sv = gridsubset$n_sv[row], 
                                fraction_degs = gridsubset$fraction_degs[row], fraction_sv = gridsubset$fraction_sv[row], status_sd = gridsubset$status_sd[row], 
                                subject_sd = 0, batch_sd = 0, sv_sd = gridsubset$sv_sd[row], 
                                resid_sd = gridsubset$resid_sd[row], seed = gridsubset$seed[row], variable_samples_per_subject = FALSE)

design = model.matrix(~ Status + Batch_ID, simData$targets)
design0 = model.matrix(~ Batch_ID, simData$targets)
colnames(design)[2] = "Status"

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
##remember to change nsv = 0 for oracle and known_nsv
write.csv(out, file="/sc/arion/projects/mscic1/results/jolie/sva_testing_leek_sim/nobatch_totalsvsd2_fractiondegs0.4_fractionsv0.4_seed12350.csv", row.names = FALSE)


write.csv(out, file="/sc/arion/projects/mscic1/results/jolie/sva_testing_leek_sim/nobatch_totalsvsd2_fractiondegs0.4_fractionsv0.4_seed12349.csv", row.names = FALSE)
write.csv(out, file="/sc/arion/projects/mscic1/results/jolie/sva_testing_leek_sim/nobatch_totalsvsd2_fractiondegs0.4_fractionsv0.4_seed12348.csv", row.names = FALSE)

##
seeds <- c(12349, 12350)

for (current_seed in seeds) {
  
  row = 1
  n_sv = 1
  total_sv_sd = 2
  
  # Construct SV grid and sv_sd
  grid_sv = as.data.frame(expand.grid(n_sv = n_sv, total_sv_sd = total_sv_sd))
  grid_sv$sv_sd = sqrt(grid_sv$total_sv_sd^2 / grid_sv$n_sv)
  grid_sv$sv_sd[!is.finite(grid_sv$sv_sd)] = 0
  
  # Create grid_list with current seed
  grid_list <- list(
    grid_sv,
    data.frame(fraction_degs = 0.4),
    data.frame(status_sd = 0.3),
    data.frame(resid_sd = 0.3),
    data.frame(fraction_sv = 0.4),
    data.frame(seed = current_seed)
  )
  
  gridsubset = as.data.frame(t(as.data.frame(unlist(grid_list))))
  
  # Simulate data
  simData = sim_case_control_data(
    n_subjects = 1000, 
    n_samples_per_subject = 1, 
    n_features = 20000, 
    n_batches = 0, 
    n_sv = gridsubset$n_sv[row],
    fraction_degs = gridsubset$fraction_degs[row], 
    fraction_sv = gridsubset$fraction_sv[row], 
    status_sd = gridsubset$status_sd[row],
    subject_sd = 0, 
    batch_sd = 0, 
    sv_sd = gridsubset$sv_sd[row],
    resid_sd = gridsubset$resid_sd[row], 
    seed = gridsubset$seed[row], 
    variable_samples_per_subject = FALSE
  )
  
  # Build design matrices
  design = model.matrix(~ Status + Batch_ID, simData$targets)
  design0 = model.matrix(~ Batch_ID, simData$targets)
  colnames(design)[2] = "Status"
  
  list_res = list()
  
  for (method in c("be", "leek", "none", "known_nsv", "oracle")) {
    
    if (method == "none") {
      n.sv = 0
    } else if (method %in% c("known_nsv", "oracle")) {
      n.sv = nrow(simData$sv_df)
    } else {
      n.sv = num.sv(simData$E, design, method = method)
    }
    
    if (n.sv > 0 & method != "oracle") {
      svobj = sva(simData$E, design, design0, n.sv = n.sv)
      design_sv = cbind(design, svobj$sv)
    } else if (n.sv > 0 & method == "oracle") {
      design_sv = model.matrix(as.formula(paste("~Status + Batch_ID", paste(simData$sv_df$sv_names, collapse = " + "), sep = " +")), simData$targets)
      colnames(design_sv)[2] = "Status"
    } else {
      design_sv = model.matrix(~ Status + Batch_ID, simData$targets)
      colnames(design_sv)[2] = "Status"
    }
    
    # Fit model
    fit = lmFit(simData, design_sv)
    eBfit = eBayes(fit)
    topSet = topTable(eBfit, coef = "Status", number = nrow(eBfit))
    
    # Confusion matrix
    vec = as.vector(table(
      estimated = factor(topSet$adj.P.Val < 0.05, levels = c("FALSE", "TRUE")),
      true = factor(topSet$DE, levels = c("FALSE", "TRUE"))
    ))
    names(vec) = c("TN", "FP", "FN", "TP")
    
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
    
    list_res[[method]] = df_res
  }
  
  out = do.call(rbind, list_res)
  
  # Save to CSV with seed in filename
  outfile = paste0("/sc/arion/projects/mscic1/results/jolie/sva_testing_leek_sim/nobatch_totalsvsd2_fractiondegs0.4_fractionsv0.4_seed", current_seed, ".csv")
  write.csv(out, file = outfile, row.names = FALSE)
}


#sanity check for 10sv - compare the new 10sv with the old sim 10sv and the results matched!
##they're the SAME!!! but in different order so i thought i was tripping LOL okay good to run for the rest!
write.csv(out, file="/sc/arion/projects/mscic1/results/jolie/sva_testing_leek_sim/sim_convergence_with_leek_9sv.csv", row.names = FALSE)
write.csv(out, file="/sc/arion/projects/mscic1/results/jolie/sva_testing_leek_sim/sim_convergence_with_leek_8sv.csv", row.names = FALSE)
write.csv(out, file="/sc/arion/projects/mscic1/results/jolie/sva_testing_leek_sim/sim_convergence_with_leek_7sv.csv", row.names = FALSE)
write.csv(out, file="/sc/arion/projects/mscic1/results/jolie/sva_testing_leek_sim/sim_convergence_with_leek_6sv.csv", row.names = FALSE)
write.csv(out, file="/sc/arion/projects/mscic1/results/jolie/sva_testing_leek_sim/sim_convergence_with_leek_5sv.csv", row.names = FALSE)
write.csv(out, file="/sc/arion/projects/mscic1/results/jolie/sva_testing_leek_sim/sim_convergence_with_leek_4sv.csv", row.names = FALSE)
#write.csv(out, file="/sc/arion/projects/mscic1/results/jolie/sva_testing_leek_sim/sim_convergence_with_leek_3sv.csv", row.names = FALSE)
#write.csv(out, file="/sc/arion/projects/mscic1/results/jolie/sva_testing_leek_sim/sim_convergence_with_leek_2sv.csv", row.names = FALSE)
#write.csv(out, file="/sc/arion/projects/mscic1/results/jolie/sva_testing_leek_sim/sim_convergence_with_leek_1sv.csv", row.names = FALSE)
#write.csv(out, file="/sc/arion/projects/mscic1/results/jolie/sva_testing_leek_sim/sim_convergence_with_leek_10sv.csv", row.names = FALSE)

out_1sv<-fread("/sc/arion/projects/mscic1/results/jolie/sva_testing_leek_sim/sim_convergence_with_leek_1sv.csv",data.table = FALSE)
out_2sv<-fread("/sc/arion/projects/mscic1/results/jolie/sva_testing_leek_sim/sim_convergence_with_leek_2sv.csv",data.table = FALSE)
out_3sv<-fread("/sc/arion/projects/mscic1/results/jolie/sva_testing_leek_sim/sim_convergence_with_leek_3sv.csv",data.table = FALSE)
out_4sv<-fread("/sc/arion/projects/mscic1/results/jolie/sva_testing_leek_sim/sim_convergence_with_leek_4sv.csv",data.table = FALSE)
out_5sv<-fread("/sc/arion/projects/mscic1/results/jolie/sva_testing_leek_sim/sim_convergence_with_leek_5sv.csv",data.table = FALSE)
out_6sv<-fread("/sc/arion/projects/mscic1/results/jolie/sva_testing_leek_sim/sim_convergence_with_leek_6sv.csv",data.table = FALSE)
out_7sv<-fread("/sc/arion/projects/mscic1/results/jolie/sva_testing_leek_sim/sim_convergence_with_leek_7sv.csv",data.table = FALSE)
out_8sv<-fread("/sc/arion/projects/mscic1/results/jolie/sva_testing_leek_sim/sim_convergence_with_leek_8sv.csv",data.table = FALSE)
out_9sv<-fread("/sc/arion/projects/mscic1/results/jolie/sva_testing_leek_sim/sim_convergence_with_leek_9sv.csv",data.table = FALSE)
out_10sv<-fread("/sc/arion/projects/mscic1/results/jolie/sva_testing_leek_sim/sim_convergence_with_leek_10sv.csv",data.table = FALSE)

##check if the datasets have the same colnames 
# Put all dataframes into a list
dfs <- list(out_1sv, out_2sv, out_3sv, out_4sv, out_5sv, out_6sv, out_7sv, out_8sv, out_9sv, out_10sv)
# Use the first dataframe's colnames as reference
ref_cols <- colnames(dfs[[1]])
# Check if all other dataframes have the same colnames
sapply(dfs, function(df) identical(colnames(df), ref_cols))

res1 <- readRDS("/sc/arion/projects/mscic1/results/jolie/sim_bug/all_feature_res_for_set5.rds") #91% of simulations 
#none<-subset(res1, method == "none")
out_1sv

test <- res1 %>% ##the saved output used to be filtered_res1
  filter(param.n_sv == 10, 
         param.fraction_degs == 0.3, 
         param.status_sd == 0.3, 
         param.resid_sd == 0.3, 
         param.fraction_sv == 0.4, 
         param.seed == 12348)

# Print the filtered result
head(filtered_res1)
#subset for what i need:
closest_value <- filtered_res1$param.sv_sd[which.min(abs(filtered_res1$param.sv_sd - 0.6324555))]
test <- filtered_res1 %>%
  filter(param.sv_sd == closest_value)
head(test)

out_0sv$Status_mean <- NULL
out_0sv$Status_sd <- NULL
out_0sv$Status_median <- NULL
out_0sv$Feature <- NULL
out_0sv$Ratio_DE_non_DE <- NULL

#reorder columns of out_0sv so that it has the same order as the rest of the out_1sv
out_0sv <- out_0sv[, colnames(out_1sv)]

all<-rbind(out_0sv, out_1sv, out_2sv, out_3sv, out_4sv, out_5sv, out_6sv, out_7sv, out_8sv,out_9sv,out_10sv)

##change nsv of known_nsv and oracle from 0 to 1 when simulated sv = 0 
all$nsv[all$param.n_sv == 0 & all$method %in% c("oracle", "known_nsv")] <- 0

out_1sv <- out_1sv %>%
  mutate(
    TPR = TP / (TP + FN),  # True Positive Rate: TP / (TP + FN) truth: all positive
    FPR = FP / (FP + TN),  # False Positive Rate: FP / (FP + TN) truth: all negative
    TNR = TN / (TN + FP),  # True Negative Rate: TN / (TN + FP) truth: all negative
    FNR = FN / (FN + TP)   # False Negative Rate: FN / (FN + TP) truth: all positive
  )

#all$param.n_sv2 <- factor(all$param.n_sv, levels = 0:8)
all$method <- factor(all$method, levels = c("oracle", "none", "known_nsv", "be", "leek"))

library(ggplot2)
plot<-ggplot(all, aes(x = param.n_sv, y = FPR)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~ method) +
  labs(
    title = "False Positive Rate by Number of Simulated SVs",
    x = "Number of Simulated SVs",
    y = "False Positive Rate"
  ) +
  scale_x_continuous(breaks = 0:10, limits = c(0, 10)) +
  theme_minimal()

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/FPR_by_simulatedSV_0_to_10.pdf", 
       plot = plot, width = 8, height = 5) 


plot<-ggplot(all, aes(x = param.n_sv, y = nsv)) +
  geom_point(alpha = 0.6) +
  geom_smooth(se = FALSE, method = "loess", color = "blue") +
  facet_wrap(~ method) +
  labs(
    title = "Number of Estimated SVs by Number of Simulated SVs",
    x = "Number of Simulated SVs",
    y = "Number of Estimated SVs"
  ) +
  scale_x_continuous(breaks = 0:10, limits = c(0, 10)) +
  theme_minimal()

ggsave("/hpc/users/hoangd02/www/plots/sva/sim_results_oct2024/estimatedSV_simulatedSV_0_to_10.pdf", 
       plot = plot, width = 8, height = 5) 


#save pval_out and out, if i dont include row.names = FALSE then there will be an X column for the row names
#write.csv(out, "/sc/arion/projects/mscic1/results/jolie/sva_testing_leek_sim/our_sim_convergence_with_leek2.csv",row.names = FALSE)
#write.csv(pval_out, "/sc/arion/projects/mscic1/results/jolie/sva_testing_leek_sim/our_sim_convergence_with_leek2.csv",row.names = FALSE)

#start 8:57pm march 11
#/sc/arion/projects/mscic1/results/jolie/
#there's one screen section outside of cbipm
#write.csv(out, "/sc/arion/projects/mscic1/results/jolie/sva_testing_leek_sim/our_sim_convergence_with_leek_out.csv",row.names = FALSE)
#write.csv(pval_out, "/sc/arion/projects/mscic1/results/jolie/sva_testing_leek_sim/our_sim_convergence_with_leek.csv",row.names = FALSE)

#there's one in cbipm01-1 
write.csv(out, "/sc/arion/projects/mscic1/results/jolie/sva_testing_leek_sim/our_sim_convergence_with_leek_out.csv",row.names = FALSE)
#write.csv(pval_out, "/sc/arion/projects/mscic1/results/jolie/sva_testing_leek_sim/our_sim_convergence_with_leek2.csv",row.names = FALSE)

out<-read.csv("/sc/arion/projects/mscic1/results/jolie/sva_testing_leek_sim/our_sim_convergence_with_leek_out.csv")
p_value<-read.csv("/sc/arion/projects/mscic1/results/jolie/sva_testing_leek_sim/our_sim_convergence_with_leek.csv")
dim(out) #50x13
#10 unique simulation parameter sets (n_sv = 1 to 11) x 5 methods
#plot FPR with lower SVs
res<-out
##get the n_sv = 0
res1 <- readRDS("/sc/arion/projects/mscic1/results/jolie/sim_bug/all_feature_res_for_set5.rds") #91% of simulations 
filtered_res1 <- res1 %>%
  filter(param.n_sv == 0, 
         param.fraction_degs == 0.3, 
         param.status_sd == 0.3, 
         param.resid_sd == 0.3, 
         param.fraction_sv == 0.5, 
         param.seed == 12348)

# Print the filtered result
head(filtered_res1)
filtered_res1$Ratio_DE_non_DE <- NULL
filtered_res1$Feature <- NULL
filtered_res1$Status_median <- NULL
filtered_res1$Status_sd <- NULL
filtered_res1$Status_mean <- NULL
dim(filtered_res1) #5 x 13

### START
filtered_res10 <- res1 %>%
  filter(param.n_sv == 10, 
         param.fraction_degs == 0.3, 
         param.status_sd == 0.3, 
         param.resid_sd == 0.3, 
         param.fraction_sv == 0.5, 
         param.seed == 12348)

filtered_res10 <- filtered_res10 %>%
  mutate(
    TPR = TP / (TP + FN),  # True Positive Rate: TP / (TP + FN) truth: all positive
    FPR = FP / (FP + TN),  # False Positive Rate: FP / (FP + TN) truth: all negative
    TNR = TN / (TN + FP),  # True Negative Rate: TN / (TN + FP) truth: all negative
    FNR = FN / (FN + TP)   # False Negative Rate: FN / (FN + TP) truth: all positive
  )

oracle <- subset(filtered_res10, method == "oracle")

oracle_in_new_sim <- res %>%
  filter(param.n_sv == 10, 
         param.fraction_degs == 0.3, 
         param.status_sd == 0.3, 
         param.resid_sd == 0.3, 
         param.fraction_sv == 0.5, 
         param.seed == 12348)
#param.sv_sd = 0.6324555

old_sim <- res1 %>%
  filter(param.n_sv == 10, 
         param.fraction_degs == 0.3, 
         param.status_sd == 0.3, 
         param.resid_sd == 0.3, 
         param.fraction_sv == 0.5, 
         param.seed == 12348)

oracle_in_old_sim <- subset(old_sim, method == "oracle")

#raw numbers
#nsv   TP    TN  FP  FN method param.n_sv param.fraction_degs
#10  5460 13794 206 540 oracle         10                 0.3

### END
res <- res %>%
  mutate(
    TPR = TP / (TP + FN),  # True Positive Rate: TP / (TP + FN) truth: all positive
    FPR = FP / (FP + TN),  # False Positive Rate: FP / (FP + TN) truth: all negative
    TNR = TN / (TN + FP),  # True Negative Rate: TN / (TN + FP) truth: all negative
    FNR = FN / (FN + TP)   # False Negative Rate: FN / (FN + TP) truth: all positive
  )

filtered_res1 <- filtered_res1[, colnames(res)]

dim(res) #50x13
res <- rbind(res, filtered_res1)
dim(res) #55x13

res <- res %>%
  mutate(
    TPR = TP / (TP + FN),  # True Positive Rate: TP / (TP + FN) truth: all positive
    FPR = FP / (FP + TN),  # False Positive Rate: FP / (FP + TN) truth: all negative
    TNR = TN / (TN + FP),  # True Negative Rate: TN / (TN + FP) truth: all negative
    FNR = FN / (FN + TP)   # False Negative Rate: FN / (FN + TP) truth: all positive
  )

#res1$difference_of_nsv_simulated_sv <- res1$nsv - res1$param.n_sv 
res$true_degs<-res$param.fraction_degs * 20000
res$all_positive=res$TP + res$FP
res$all_negative=res$TN + res$FN
res$accuracy=(res$TP + res$TN)/(res$TP + res$TN + res$FP + res$FN)
res$precision=res$TP/(res$TP + res$FP) 
res$specificity=res$TN/(res$TN + res$FP) #ability to correctly identify true negatives; It indicates how well a test can avoid false positives.
#res$falsePositiveRate=res$FP/(res$FP + res$TN) #type 1 error = false positive #SVA commits this!
res$recall=res$TP/(res$TP + res$FN) #true positive rate; measure of a test's ability to correctly identify true positives
res$negativePredictiveValue=res$TN/(res$TN + res$FN)
#res$falseNegativeRate=res$FN/(res$FN + res$TP) #type 2 error = false negative
res$totalDEGsignal=res$param.fraction_degs*res$param.status_sd
dim(res) #55 26

# 
# plot <- ggplot(res, aes(x = param.n_sv, y = FPR, color = as.factor(nsv))) +
#   geom_point() +
#   labs(
#     title = "FPR vs. Simulated and Detected Surrogate Variables",
#     x = "Simulated Number of SVs (param.n_sv)",
#     y = "False Positive Rate (FPR)",
#     color = "Detected SVs (nsv)"
#   ) +
#   theme_minimal()
# 
# ggsave("/hpc/users/hoangd02/www/plots/sva_leek_convergence.pdf", 
#        plot = plot, width = 8, height = 5) 
# 
# ggplot(res, aes(x = param.n_sv, y = FPR, color = as.factor(nsv))) +
#   geom_point() +
#   geom_line() +
#   labs(
#     title = "FPR vs. Simulated and Detected Surrogate Variables by Method",
#     x = "Simulated Number of SVs (param.n_sv)",
#     y = "False Positive Rate (FPR)",
#     color = "Detected SVs (nsv)"
#   ) +
#   theme_minimal() +
#   scale_x_continuous(limits = c(0, 1)) + 
#   facet_wrap(~ method)
# 
# ggsave("/hpc/users/hoangd02/www/plots/sva_leek_convergence2.pdf", 
#        plot = plot, width = 8, height = 5) 

res$method <- factor(res$method, levels = c("oracle", "none", "known_nsv", "be", "leek"))
#oracle <- subset(res, method == "oracle")
res$nsv[res$param.n_sv == 0 & res$method %in% c("oracle", "known_nsv")] <- 0

plot<-ggplot(res, aes(x = param.n_sv, y = FPR, color = as.factor(nsv))) +
  geom_point(size = 3) +  # Larger points for better visibility
  geom_text(aes(label = nsv), vjust = -1, size = 4) +  # Adds text labels slightly above points
  scale_x_continuous(breaks = 0:10, limits = c(0, 10)) +  # Sets x-axis range and breaks
  facet_wrap(~ method) + # Separate histograms by method
  labs(
    title = "FPR vs. Simulated and Detected Surrogate Variables",
    x = "Simulated Number of SVs (param.n_sv)",
    y = "False Positive Rate (FPR)",
    color = "Detected SVs (nsv)"
  ) +
  theme_minimal()

#ggsave("/hpc/users/hoangd02/www/plots/sva_leek_convergence3.pdf", 
#       plot = plot, width = 8, height = 10) 
ggsave("/hpc/users/hoangd02/www/plots/sva_leek_convergence_with_0.pdf", 
       plot = plot, width = 8, height = 10) 



