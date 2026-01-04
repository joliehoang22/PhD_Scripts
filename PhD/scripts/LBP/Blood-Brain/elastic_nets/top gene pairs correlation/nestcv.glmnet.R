nestcv.glmnet(
       y,
       x,
       family = c("gaussian", "binomial", "poisson", "multinomial", "cox", "mgaussian"),
       filterFUN = NULL,
       filter_options = NULL,
       balance = NULL,
       balance_options = NULL,
       modifyX = NULL,
       modifyX_useY = FALSE,
       modifyX_options = NULL,
       outer_method = c("cv", "LOOCV"),
       n_outer_folds = 10,
       n_inner_folds = 10,
       outer_folds = NULL,
       pass_outer_folds = FALSE,
       alphaSet = seq(0.1, 1, 0.1),
       min_1se = 0,
       keep = TRUE,
       outer_train_predict = FALSE,
       weights = NULL,
       penalty.factor = rep(1, ncol(x)),
       parallel_mode = NULL,
       cv.cores = 1,
       finalCV = TRUE,
       na.option = "omit",
       verbose = FALSE,
       ...
     )

Arguments:

       y: Response vector or matrix. Matrix is only used for ‘family =
          'mgaussian'’ or ‘'cox'’.

       x: Matrix of predictors. Dataframes will be coerced to a matrix
          as is necessary for glmnet.

  family: Either a character string representing one of the built-in
          families, or else a ‘glm()’ family object. Passed to
          glmnet::cv.glmnet and glmnet::glmnet

filterFUN: Filter function, e.g. ttest_filter or relieff_filter. Any
          function can be provided and is passed ‘y’ and ‘x’. Must
          return a numeric vector with indices of filtered predictors.

filter_options: List of additional arguments passed to the filter
          function specified by ‘filterFUN’.

 balance: Specifies method for dealing with imbalanced class data.
          Current options are ‘"randomsample"’ or ‘"smote"’. See
          ‘randomsample()’ and ‘smote()’

balance_options: List of additional arguments passed to the balancing
          function

modifyX: Character string specifying the name of a function to modify
          ‘x’. This can be an imputation function for replacing missing
          values, or a more complex function which alters or even adds
          columns to ‘x’. The required return value of this function
          depends on the ‘modifyX_useY’ setting.

modifyX_useY: Logical value whether the ‘x’ modifying function makes
          use of response training data from ‘y’. If ‘FALSE’ then the
          ‘modifyX’ function simply needs to return a modified ‘x’
          object, which will be coerced to a matrix as required by
          ‘glmnet’. If ‘TRUE’ then the ‘modifyX’ function must return a
          model type object on which ‘predict()’ can be called, so that
          train and test partitions of ‘x’ can be modified
          independently.

modifyX_options: List of additional arguments passed to the ‘x’
          modifying function

outer_method: String of either ‘"cv"’ or ‘"LOOCV"’ specifying whether
          to do k-fold CV or leave one out CV (LOOCV) for the outer
          folds

n_outer_folds: Number of outer CV folds

n_inner_folds: Number of inner CV folds

outer_folds: Optional list containing indices of test folds for outer
          CV. If supplied, ‘n_outer_folds’ is ignored.

pass_outer_folds: Logical indicating whether the same outer folds are
          used for fitting of the final model when final CV is applied.
          Note this can only be applied when ‘n_outer_folds’ and
          ‘n_inner_folds’ are the same and no balancing is applied.

alphaSet: Vector of alphas to be tuned

min_1se: Value from 0 to 1 specifying choice of optimal lambda from
          0=lambda.min to 1=lambda.1se

    keep: Logical indicating whether inner CV predictions are retained
          for calculating left-out inner CV fold accuracy etc. See
          argument ‘keep’ in glmnet::cv.glmnet.

outer_train_predict: Logical whether to save predictions on outer
          training folds to calculate performance on outer training
          folds.

 weights: Weights applied to each sample. Note ‘weights’ and ‘balance’
          cannot be used at the same time. Weights are only applied in
          glmnet and not in filters.

penalty.factor: Separate penalty factors can be applied to each
          coefficient. Can be 0 for some variables, which implies no
          shrinkage, and that variable is always included in the model.
          Default is 1 for all variables. See glmnet::glmnet. Note this
          works separately from filtering. For some ‘nestedcv’ filter
          functions you might need to set ‘force_vars’ to avoid
          filtering out features.

parallel_mode: Either "mclapply", "parLapply" or "future". This
          determines which parallel backend to use. The default is
          ‘parallel::mclapply’ on unix/mac and ‘parallel::parLapply’ on
          windows.

cv.cores: Number of cores for parallel processing of the outer loops.
          Ignored if ‘parallel_mode = "future"’.

finalCV: Logical whether to perform one last round of CV on the whole
          dataset to determine the final model parameters. If set to
          ‘FALSE’, the median of hyperparameters from outer CV folds
          are used for the final model. Performance metrics are
          independent of this last step. If set to ‘NA’, final model
          fitting is skipped altogether, which gives a useful speed
          boost if performance metrics are all that is needed.

na.option: Character value specifying how ‘NA’s are dealt with.
          ‘"omit"’ (the default) is equivalent to ‘na.action =
          na.omit’. ‘"omitcol"’ removes cases if there are ‘NA’ in 'y',
          but columns (predictors) containing ‘NA’ are removed from 'x'
          to preserve cases. Any other value means that ‘NA’ are
          ignored (a message is given).

 verbose: Logical whether to print messages and show progress

     ...: Optional arguments passed to glmnet::cv.glmnet