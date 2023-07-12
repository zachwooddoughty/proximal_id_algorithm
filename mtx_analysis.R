# packages that require installation
suppressMessages(library(R.utils, warn.conflicts=FALSE))  # for arg parsing
suppressMessages(library(jsonlite))
library(ipw)  # for reweighing GMM when M is continuous
library(cubature)  # for computing integrals for ground-truth effect
library(nnet)  # for multiclass logistic regression
library(numDeriv)

library(data.table)
library(MASS)
library(stats)
library(stringr)

source("utils.R")


GMM_func  <- function(mrf, para, data, weights=NA, vec1vars=NA, vec2vars=NA){
    g0 <- mrf(para=para ,data=data, weights=weights, vec1vars=vec1vars, vec2vars=vec2vars)
    g <- apply(g0, 2, mean)
    gmmf  <- sum(g^2)
    return(gmmf)
}

naive_frontdoor <- function(data, K=100, a_ints=FALSE, a_vals=NULL){

  N <- dim(data)[1]

  formula <- as.formula(paste("y~", paste(c(dfvars$m, dfvars$z, dfvars$x), collapse="+")))
  eym <- lm(formula, data=data)

  pmaz_funcs <- list()
  for (m_var in dfvars$m) {
    m_predictors <- c(dfvars$x, dfvars$z)
    if (a_ints){
        m_predictors <- c(m_predictors, dfvars$a_ints)
    }
    formula <- as.formula(paste(m_var, "~a+", paste(m_predictors, collapse="+")))
    pmaz <- glm(formula, family=binomial(), data=data)
    pmaz_funcs[[m_var]] <- pmaz
  }

  if (is.null(a_vals)){
    a_vals <- sort(unique(data$a))
  }

  counterfactuals <- c()
  for (a in a_vals){
      counterfactual <- 0
      for (i in 1:N){
        z_row <- data[i, dfvars$z]
        x_row <- data[i, dfvars$x]

        m0_preds <- c()
        m1_preds <- c()
        for (m_var in dfvars$m){
          pmaz_inputs <- cbind(1, a, z_row, x_row, row.names=NULL)
          if (a_ints) {
              a_int_cols <- cbind(a * x_row, a*a, row.names=NULL)
              names(a_int_cols) <- dfvars$a_ints
              pmaz_inputs <- cbind(pmaz_inputs, a_int_cols)
          }
          coefs <- as.matrix(pmaz_funcs[[m_var]]$coefficients)
          # Hack in case we have linearly-dependent data
          coefs[is.na(coefs)] <- 0
          m1_prob <- plogis(as.matrix(pmaz_inputs) %*% coefs)
          m1_pred <- rbinom(n=K, size=1, p=m1_prob)
          m1_preds <- cbind(m1_preds, m1_pred)

        }
        m1 <- m1_preds
        colnames(m1) <- dfvars$m

        m_effect <- sum(as.matrix(cbind(1, m1, z_row, x_row, row.names=NULL))
                        %*% as.matrix(eym$coef))
        counterfactual <- counterfactual + m_effect / K
      }
      counterfactuals <- c(counterfactuals, counterfactual / N)
  }

  return(counterfactuals)
}

simple_mrf <- function(para, data1, weights=NA, vec1vars, vec2vars){
  # Moment restriction function for simple proximal

  intercept <- rep(1, dim(data1)[1])
  vec1 <- as.matrix(cbind(intercept, data1[, c(vec1vars)]))
  vec2 <- as.matrix(cbind(intercept, data1[, c(vec2vars)]))

  hlink <- vec1 %*% para
  g <- (as.vector(data1[, "y"] - hlink)) * vec2
  return(g)
}

simple_proximal <- function(data1, a_ints=FALSE, a_vals=NULL) {

  vec1vars <- c("a", dfvars$w, dfvars$x)
  vec2vars <- c("a", dfvars$z, dfvars$x)
  if (a_ints){
      vec1vars <- c(vec1vars, dfvars$a_ints)
      vec2vars <- c(vec2vars, dfvars$a_ints)
  }
  inioptim <- runif(1 + length(vec1vars), -2, 2)

  bpar_am  <- optim(
      par = inioptim, fn = GMM_func, mrf = simple_mrf, data = data1,
      weights = all_weights, vec1vars = vec1vars, vec2vars = vec2vars,
      method = "BFGS", hessian = FALSE
      )$par

  # gmm_fit <- GMM_func(simple_mrf, bpar_am, data=data1, weights=all_weights,
  #                 vec1vars=vec1vars, vec2vars=vec2vars)
  # print(paste("gmm fit :", gmm_fit, sep=" "))

  w <- data1[, dfvars$w]
  x <- data1[, dfvars$x]
  
  if (is.null(a_vals)){
    a_vals <- sort(unique(data1$a))
  }

  counterfactuals <- c()
  for (a in a_vals){
      b_mat <- cbind(1, a, w, x, row.names=NULL)

      if (a_ints){
          a_int_cols <- cbind(a * x, a * a, row.names=NULL)
          names(a_int_cols) <- dfvars$a_ints
          b_mat <- cbind(b_mat, a_int_cols)
      }
      counterfactual <- mean(as.matrix(b_mat) %*% as.matrix(bpar_am))
      counterfactuals <- c(counterfactuals, counterfactual)
  }
  return(counterfactuals)
}

mtx_mrf <- function(para, data1, weights, vec1vars, vec2vars){
  # Moment restriction function for proximal front-door

  intercept <- rep(1, dim(data1)[1])
  vec1 <- as.matrix(cbind(intercept, data1[, c(vec1vars)]))
  vec2 <- as.matrix(cbind(intercept, data1[, c(vec2vars)]))

  hlink <- (vec1 / weights) %*% para
  g <- (as.vector(data1[, "y"] / weights - hlink)) * vec2

  return(g)
}

# Proximal frontdoor set up for MTX dataset
proximal_frontdoor <- function(data1, a_ints=FALSE, a_vals=NULL) {
    N <- dim(data1)[1]
    all_weights <- rep(1, N)
    pmaz_funcs <- list()
    for (m_var in dfvars$m) {
        formula <- as.formula(paste(m_var, "~ a+", paste(c(dfvars$x, dfvars$z), collapse="+")))
        pmaz <- glm(formula, family=binomial(), data=data1)
        pmaz_funcs[[m_var]] <- pmaz
        weights <- plogis(predict(pmaz, newdata=data1))
        weights <- data1[, m_var] * weights + (1 - data1[, m_var]) * (1 - weights)
        all_weights <- all_weights * weights
    }

    all_weights <- clip_weights(all_weights)
    data1$all_weights <- all_weights

    vec1vars <- c("a", dfvars$m, dfvars$w, dfvars$x)
    vec2vars <- c("a", dfvars$m, dfvars$z, dfvars$x)
    # If we are considering treatment interactions, add those in
    if (a_ints){
        vec1vars <- c(vec1vars, dfvars$a_ints)
        vec2vars <- c(vec2vars, dfvars$a_ints)
    }

    # Fit the bridge function
    inioptim <- runif(1 + length(vec1vars), -2, 2)
    bpar_am  <- optim(
        par = inioptim, fn = GMM_func, mrf = mtx_mrf, data = data1,
        weights = all_weights, vec1vars = vec1vars, vec2vars = vec2vars,
        method = "BFGS", hessian = FALSE
        )$par

    # gmm_fit <- GMM_func(mtx_mrf, bpar_am, data=data1, weights=all_weights, vec1vars=vec1vars, vec2vars=vec2vars)
    # print(paste(c("gmm fit :", gmm_fit), sep=" "))

    # Fit p(A | X, Z)
    formula <- as.formula(paste("a~", paste(c(dfvars$x, dfvars$z), collapse="+")))
    paz <- nnet::multinom(formula, data=data1, trace=FALSE)

    # Fit p(W | A, M, Z, X) for each W variable
    pwamz_funcs <- list()
    for (w_var in dfvars$w) {
        formula <- as.formula(paste(w_var, "~a+", paste(c(dfvars$m, dfvars$z, dfvars$x), collapse="+"), sep=''))
        pwamz <- nnet::multinom(formula, data=data1, weights=1 / all_weights, trace=FALSE)

        pwamz_funcs[[w_var]] <- pwamz
    }

    if (is.null(a_vals)){
      a_vals <- sort(unique(data1$a))
    }

    # E(Y(a)) = 
    # \sum_{m,z} {\sum_w b(w, m, a) \sum_{a'} p(w | m, a', z)p(a' | z)}
    #            p(m | a, z)p(z)
    K <- 100  # How many trajectories to compute
    counterfactuals <- c()
    for (a in a_vals){  # For each possible treatment value
        counterfactual <- 0
        for (i in 1:N){
            z_row <- data1[i, dfvars$z]
            x_row <- data1[i, dfvars$x]

            # Sample M ~ A, Z, X for each M variable
            m_preds <- c()
            for (m_var in dfvars$m){
                m_prob <- plogis(as.matrix(cbind(1, a, z_row, x_row, row.names=NULL))
                                 %*% as.matrix(pmaz_funcs[[m_var]]$coefficients))
                m_pred <- rbinom(n=K, size=1, p=m_prob)
                m_preds <- cbind(m_preds, m_pred)
            }
            m <- m_preds
            colnames(m) <- dfvars$m

            # Compute a'
            a_prime_probs <- predict(paz, newdata=cbind(1, z_row, x_row, row.names=NULL), type="probs")
            a_sample_array <- rmultinom(K, 1, a_prime_probs)
            a_prime <- apply(a_sample_array, 2, function(row){row %*% c(1:dim(a_sample_array)[1])})

            # Sample W ~ A', M, Z, X for each W variable
            w_preds <- c()
            for (w_var in dfvars$w){
                w_dim <- length(unique(data1[, w_var]))
                w_pred_data <- cbind(1, a_prime, m, z_row, x_row, row.names=NULL)
                names(w_pred_data)[names(w_pred_data) == 'a_prime'] <- 'a'

                # This outputs a (K, w_dim) matrix
                w_probs <- predict(pwamz_funcs[[w_var]], newdata=w_pred_data, type="probs")

                if (w_dim == 2){
                    w_pred <- rbinom(K, 1, w_probs)
                } else {
                    # For each row of w_probs, sample a one-hot array
                    w_sample_array <- apply(w_probs, 1, function(row){rmultinom(1, 1, row)})

                    # Convert those one-hot arrays to actual values
                    # Note: rmultinom returns a (w_shape x K) matrix, so we apply to columns
                    w_pred <- apply(w_sample_array, 2,
                                    function(row){row %*% c(1:dim(w_sample_array)[1])})
                }

                w_preds <- cbind(w_preds, w_pred)
            }
            w <- w_preds
            colnames(w) <- dfvars$w

            b_mat <- cbind(1, a, m, w, x_row, row.names=NULL)

            # Optionally consider treatment interactions
            if (a_ints){
                a_int_cols <- cbind(a * x_row, a*a, row.names=NULL)
                names(a_int_cols) <- dfvars$a_ints
                b_mat <- cbind(b_mat, a_int_cols)
            }

            # sum the causal effects for K trajectories sampled from the i'th
            # row and add it to the running sum for this treatment assignment
            b <- as.matrix(b_mat) %*% as.matrix(bpar_am)
            counterfactual <- counterfactual + sum(b)
          }
      counterfactuals <- c(counterfactuals, counterfactual / (N * K))
    }
    return(counterfactuals)
}

# For the given arguments, find the correct dataset
#   These files need to be created first by preprocess_mtx_data.py
get_dataset_fn <- function(args){
    if (args$binary_treat){
        base_str <- "mtx_data/mtx_bin"
    } else {
        base_str <- "mtx_data/mtx_cont"
    }
    timestep_str <- paste(c("m", args$mtime, "_w", args$wtime, "_y", args$ytime),
                          sep="", collapse="")

    infn_cols <- c(base_str, timestep_str)
    if (args$no_mtx_after > -1){
        infn_cols <- c(infn_cols, paste(c("mtx_to", args$no_mtx_after),
                                        sep="", collapse=""))
    }
    if (args$no_mtx_before > -1){
        infn_cols <- c(infn_cols, paste(c("mtx_from", args$no_mtx_before),
                                        sep="", collapse=""))
    }
    if (args$delta_outcome){
        infn_cols <- c(infn_cols, "ydelta")
    }
    if (args$a_ints) {
        infn_cols <- c(infn_cols, "a_int")
    } 
    infn_cols <- c(infn_cols, c("csv", "gz"))
    infn <- paste(infn_cols, sep="", collapse=".")
    return(infn)
}

# Which mtx data columns we use for which variables
dfvars = list(
  x=c('age_0', 'sex', 'year_0', 'onprd2_0', 'duration_0',
      'edu_0', 'smoke_0', 'dmrd_0', 'rapos_0'),
  m=c("dmrd_m", "onprd2_m"),
  z=c('haqc_z', 'esrc_z', 'jc_z', 'gsc_z'),
  w=c('haqc_w', 'esrc_w', 'jc_w', 'gsc_w'),
  a_ints=c('a.age_0', 'a.sex', 'a.year_0', 'a.onprd2_0', 'a.duration_0',
      'a.edu_0', 'a.smoke_0', 'a.dmrd_0', 'a.rapos_0', 'a.a')
)

# Argument defaults; see preprocess_mtx_data.py for corresponding args 
defaults <- list(
    models="naive_frontdoor,simple_proximal,proximal_frontdoor",
    n_bootstrap=0,      # How many bootstrap datasets to resample?
    binary_treat=TRUE,  # Are we using a binary or categorical treatment?
    a_ints=FALSE,       # Are we consider treatment interactions?
    delta_outcome=TRUE, # Is our outcome the difference in jc from baseline or just jc?
    mtime=6,            # What month timestep do we use for mediators?
    wtime=7,            # What month timestep do we use for outcome-inducing proxy?
    ytime=12,           # What month timestep do we use for outcome?
    no_mtx_before=-1,   # Ignore patients who started mtx before this month
    no_mtx_after=-1     # Ignore patients who started mtx after this month
)

main <- function(){
    args <- commandArgs(trailingOnly=TRUE, asValues=TRUE, defaults=defaults)

    # Check args are choosing allowed models
    for (model in unlist(strsplit(args$models, ","))){
        bad <- TRUE
        for (allowed in c("naive_frontdoor", "simple_proximal", "proximal_frontdoor")){
            if (model == allowed){ bad <- FALSE }
        }
        if (bad){ stop(paste("Unknown model: ", model, sep="", collapse="")) }
    }

    df <- read.csv(get_dataset_fn(args))
    rownames(df) <- NULL

    # Printout unique A values, distribution of observed treatments
    a_vals <- sort(unique(df$a))
    a_val_counts <- c()
    for (a_val in a_vals){
        count_a_val <- sum(df$a == a_val)
        a_val_counts <- c(a_val_counts, paste(a_val, "=", count_a_val))
    }
    print(paste("Treatment distribution", a_val_counts, collapse=", "))

    ci_probs <- c(0.025, 0.975)
    causal_effect <- function(row) { row - row[1] }

    # Run naive frontdoor
    if (grepl("naive_frontdoor", args$models, fixed=TRUE)){
        set.seed(42)
        est <- naive_frontdoor(df, a_ints=args$a_ints)
        print(paste("naive front:", paste(round(est, 3), collapse=", "), sep=" "))
        print(paste("causal eff: ", paste(round(causal_effect(est), 3), collapse=", "), sep=" "))
        if (args$n_bootstrap > 0){
            boot <- bootstrap_ci(df, args$n_bootstrap,
                                 naive_frontdoor, a_ints=args$a_ints, a_vals=a_vals)
            ci <- apply(boot, 2, quantile, probs=ci_probs)
            colnames(ci) <- a_vals
            rownames(ci) <- paste(ci_probs * 100, "%", sep="")
            print(round(ci, 3))
            effects <- t(apply(boot, 1, causal_effect))
            ci <- apply(effects, 2, quantile, probs=ci_probs)
            colnames(ci) <- a_vals
            rownames(ci) <- paste(ci_probs * 100, "%", sep="")
            print(round(ci, 3))
        }
    }

    # Run simple proximal
    if (grepl("simple_proximal", args$models, fixed=TRUE)){
        set.seed(42)
        est <- simple_proximal(df, a_ints=args$a_ints)
        print(paste("simple prox:", paste(round(est, 3), collapse=", "), sep=" "))
        print(paste("causal eff: ", paste(round(causal_effect(est), 3), collapse=", "), sep=" "))
        if (args$n_bootstrap > 0){
            boot <- bootstrap_ci(df, args$n_bootstrap,
                                 simple_proximal, a_ints=args$a_ints, a_vals=a_vals)
            ci <- apply(boot, 2, quantile, probs=ci_probs)
            colnames(ci) <- a_vals
            rownames(ci) <- paste(ci_probs * 100, "%", sep="")
            print(round(ci, 3))

            effects <- t(apply(boot, 1, causal_effect))
            ci <- apply(effects, 2, quantile, probs=ci_probs)
            colnames(ci) <- a_vals
            rownames(ci) <- paste(ci_probs * 100, "%", sep="")
            print(round(ci, 3))
        }
    }

    # Run proximal frontdoor
    if (grepl("proximal_frontdoor", args$models, fixed=TRUE)){
        set.seed(42)
        est <- proximal_frontdoor(df, a_ints=args$a_ints)
        print(paste("prox front:", paste(round(est, 3), collapse=", "), sep=" "))
        print(paste("causal eff: ", paste(round(causal_effect(est), 3), collapse=", "), sep=" "))
        if (args$n_bootstrap > 0){
            boot <- bootstrap_ci(df, args$n_bootstrap,
                                 proximal_frontdoor, a_ints=args$a_ints, a_vals=a_vals)
            ci <- apply(boot, 2, quantile, probs=ci_probs)
            colnames(ci) <- a_vals
            rownames(ci) <- paste(ci_probs * 100, "%", sep="")
            print(round(ci, 3))

            effects <- t(apply(boot, 1, causal_effect))
            ci <- apply(effects, 2, quantile, probs=ci_probs)
            colnames(ci) <- a_vals
            rownames(ci) <- paste(ci_probs * 100, "%", sep="")
            print(round(ci, 3))
        }
    }
}

main()
warnings()
