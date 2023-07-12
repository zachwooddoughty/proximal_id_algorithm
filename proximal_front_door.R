# packages that require installation
suppressMessages(library(R.utils, warn.conflicts=FALSE))  # for arg parsing
suppressMessages(library(jsonlite))
library(ipw)  # for reweighing GMM when M is continuous
library(cubature)  # for computing integrals for ground-truth effect

library(MASS)
library(numDeriv)
library(stats)
library(stringr)
library(jsonlite)

source("utils.R")


model_sd <- function(model, N) {
    # Compute standard deviation of model for sampling
    denom <- N - (1 + (length(model$coef) - 1))
    return(sqrt(sum(model$residuals**2) / denom))
}

GMM_func <- function(mrf, para, data, weights=NA){
  # General method of moments function
  g0 <- mrf(para=para ,data=data, weights=weights)
  g <- apply(g0, 2, mean)
  gmmf  <- sum(g^2)
  return(gmmf)
}

simple_mrf <- function(para, data1, weights=NA){
  # Moment restriction function for simple proximal
  W <- data1[, "W"]
  Y <- data1[, "Y"]
  Z <- data1[, "Z"]
  A <- data1[, "A"]
  CC <- data1[, "CC"]

  hlink  <- cbind(1, A, W, CC) %*% para     # b(W, a)
  # E[. | A, Z] for both E[b(W, a) | A, Z] and E[Y | A, Z]
  g0 <- cbind(1, A, Z, CC)
  g <- (as.vector(Y - hlink )) * g0  # E[Y - b(W, a) | A, Z]
  return(g)
}

proximal_frontdoor_mrf <- function(para, data1, weights){
  # Moment restriction function for proximal front-door

  M <- data1[, "M"]
  W <- data1[, "W"]
  Y <- data1[, "Y"]
  Z <- data1[, "Z"]
  A <- data1[, "A"]
  CC <- data1[, "CC"]

  if (any(is.na(weights))){
    weights <- 1
  }
  # We only need the expectation for our simulation studies
  #   if we wanted p(Y | ...), we'd need Y in this function
  hlink  <- (cbind(1, A, M, W, CC) / weights) %*% para
  g0 <- cbind(1, A, M, Z, CC)
  
  g <- (as.vector(Y / weights - hlink)) * g0

  return(g)
}

simulate <- function(sim_parameters, N, z_bin=TRUE, m_bin=TRUE,
                     a_val=NA){
  # Generate synthetic data 

  mbeta <- sim_parameters$mbeta
  wbeta <- sim_parameters$wbeta
  abeta <- sim_parameters$abeta
  ybeta <- sim_parameters$ybeta
  zbeta <- sim_parameters$zbeta
  ubeta <- sim_parameters$ubeta
  cbeta <- sim_parameters$cbeta
  sd_ <- sim_parameters$sd_

  # printout("ybeta", ybeta, c("int", "A", "U", "M", "W"))

  # Generate unobserved confounder
  U0 <- rnorm(N, mean=ubeta, sd=sd_)
  # Generate observed confounder
  C0 <- cbind(1, U0) %*% cbeta + rnorm(N, mean=0, sd=sd_)

  # Generate control (treatment-inducing) proxy Z0
  if (z_bin) {
    Z_prob <- cbind(1, U0, C0) %*% zbeta
    Z0 <- rbinom(n=N, size=1, prob=plogis(Z_prob))
  } else {
    Z0 <- cbind(1, U0, C0) %*% zbeta + rnorm(N, mean=0, sd=sd_)
  }

  # Generate the exposure/treatment
  if (is.na(a_val)){
    A_prob <- cbind(1, U0, Z0, C0) %*% abeta
    A0 <- rbinom(n=N, size=1, prob=plogis(A_prob))
  } else {
    A0 <- rep(a_val, N)
  }

  # Generate the mediator
  if (m_bin) {
    M_prob <- cbind(1, A0, Z0, C0) %*% mbeta
    M0 <- rbinom(n=N, size=1, prob=plogis(M_prob))
  } else {
    M0 <- cbind(1, A0, Z0, C0) %*% mbeta + rnorm(N, mean=0, sd=sd_)
  }

  # Generate outcome-inducing proxy W0
  W0 <- cbind(1, U0, M0, C0) %*% wbeta + rnorm(N,mean=0,sd=sd_)

  # Generate outcome
  Y0 <- cbind(1, A0, U0, M0, W0, C0) %*% ybeta + rnorm(N,mean=0,sd=sd_)

  data1 <- list(A=A0,M=M0,Y=Y0,Z=Z0,W=W0,U=U0,CC=C0)
  return(data1)
}

# Simple Proximal method, assuming no Z->M->W path
simple_prox <- function(data1, inioptim) {

  hpar  <- optim(
      par = inioptim, fn = GMM_func, mrf = simple_mrf, data = data1,
      method = "BFGS", hessian = FALSE)$par

  W <- data1[, "W"]
  Y <- data1[, "Y"]
  Z <- data1[, "Z"]
  A <- data1[, "A"]
  U <- data1[, "U"]
  CC <- data1[, "CC"]

  # Useful debugging statements
  # gmm_fit <- GMM_func(simple_mrf, hpar, data=data1)
  # print(paste("gmm_fit", gmm_fit))
  # ewaz <- lm(W ~ A + Z, data=data1)
  # printout("ewaz", ewaz$coefficients, cols=c("int", "a", "z"))
  # euaz <- lm(U ~ A + Z, data=data1)
  # printout("euaz", euaz$coefficients, cols=c("int", "a", "z"))
  # B <- cbind(1, A, W, CC) %*% hpar
  # printout("compare", sum(apply(Y - B, 2, mean)^2))

  B1 <- cbind(1, 1, W, CC) %*% hpar
  B0 <- cbind(1, 0, W, CC) %*% hpar
  return(mean(B1 - B0))
}

proximal_frontdoor <- function(data1, inioptim, m_bin=TRUE, z_bin=TRUE) {
    # Code for the proximal frontdoor estimator

    N <- dim(data1)[1]
    Y <- data1[, "Y"]
    A <- data1[, "A"]
    M <- data1[, "M"]
    W <- data1[, "W"]
    Z <- data1[, "Z"]
    U <- data1[, "U"]
    CC <- data1[, "CC"]

    # Fit p(m | a, z, c)
    if (m_bin) {
        pmaz <- glm(M ~ A + Z + CC, family=binomial(), data=data1)
        weights <- plogis(predict(pmaz, newdata=data1))
        weights <- M * weights + (1 - M) * (1 - weights)
    } else {
        # This code benefitted greatly from the following guide
        # https://www.andrewheiss.com/blog/2020/12/01/ipw-binary-continuous/
        pmaz <- lm(M ~ A + Z + CC, data=data1)
        weights <- ipwpoint(
            exposure = M,
            family = "gaussian",
            numerator = ~ 1,
            denominator = ~ A + Z + CC,
            data = data1
        )
        weights <- weights$ipw.weights
        pmaz_sd <- model_sd(pmaz, N)
    }

    # Fit the bridge function
    weights <- clip_weights(weights)
    bpar_am  <- optim(
        par = inioptim, fn = GMM_func, mrf = proximal_frontdoor_mrf, data = data1,
        weights = weights, method = "BFGS", hessian = FALSE
        )$par

    # Fit p(a | z, c) and p(w | m, a, z, c)
    paz <- glm(A ~ Z + CC, family=binomial(), data=data1)
    pwmaz <- lm(W ~ M + A + Z + CC, data=data1, weights=1/weights)
    pwmaz_sd <- model_sd(pwmaz, N)

    # Useful debugging statements
    # gmm_fit <- GMM_func(proximal_frontdoor_mrf, bpar_am, data=data1, weights=weights)
    # print(paste("gmm fit :", gmm_fit, sep=" "))
    # printout("gamma", bpar_am, cols=c("int", "A", "M", "W"))
    # euamz <- lm(U ~ A + M + Z, data=data1, weights=weights)
    # printout("euamz", euamz$coef, cols=c("int", "a", "m", "z"))
    # ewamz <- lm(W ~ A + M + Z, data=data1, weights=weights)
    # printout("ewamz", ewamz$coef, cols=c("int", "a", "m", "z"))
    # printout("paz", paz$coef, c("int", "z"))
    # printout("pwmaz", pwmaz$coef, c("int", "m", "a", "z"))


    K <- 100  # How many trajectories of m, a', and w to sample
    effect <- 0  # Start an empty sum for effect estimate

    # Equation (26) in the paper, except C is omitted from the following,
    # and we're using the expectation E[Y(a)] instead of p(Y(a))
    # E[Y(a)] = \sum_{m,z} p(m | a, z) p(z) {\sum_w b(w, m, a) \sum_{a'} p(w | m, a', z) p(a' | z)} 
    for (i in 1:N){
        for (a in 0:1) {
            # Multiply effect by +1 or -1 depending on A
            a_effect_dir <- 2 * a - 1
            z <- Z[i]
            cc <- CC[i]

            # Sample M for outer summation
            if (m_bin) {
              m_prob <- plogis(cbind(1, a, z, cc) %*% pmaz$coef)
              m <- rbinom(n=K, size=1, p=m_prob)
            } else {
              m_mean <- cbind(1, a, z, cc) %*% pmaz$coef
              m <- rnorm(n=K, mean=m_mean, sd=pmaz_sd)
            }

            # Sample a' and w for inner summation
            a_prime_prob <- plogis(c(1, z, cc) %*% paz$coef)
            a_prime <- rbinom(n=K, size=1, p=a_prime_prob)
            w_mean <- cbind(1, m, a_prime, z, cc) %*% pwmaz$coef
            w <- rnorm(n=K, mean=w_mean, sd=pwmaz_sd)

            # plug a, m, w, and c into bridge function
            b <- cbind(1, a, m, w, cc) %*% bpar_am

            # sum the causal effects for K trajectories sampled from the i'th
            # row and add it or subtract it from running sum
            effect <- effect + a_effect_dir * sum(b)
        }
    }
    effect <- effect / (N * K)

    return(effect)
}


frontdoor <- function(data, m_bin=TRUE) {
  # Naive frontdoor estimator
  pym <- lm(Y ~ M, data=data)
  if (m_bin){
    pma <- glm(M ~ A + Z, data=data, family=binomial())
  } else {
    pma <- lm(M ~ A + Z, data=data)
  }

  am <- pma$coef[2]
  my <- pym$coef[2]
  val <- as.matrix(am * my)
  return(val)
}


# \sum_{m,z} p(m | a,z) \sum_{a’} p(Y | a’,m,z) p(a’|z) p(z)
#  subtract the same thing but with p(m | a_2,z)

sample_frontdoor <- function(data, K=100, m_bin=TRUE){
  N <- dim(data)[1]
  eym <- lm(Y ~ M + Z + CC, data=data)
  if (m_bin){
    pma <- glm(M ~ A + Z + CC, data=data, family=binomial())
  } else {
    pma <- lm(M ~ A + Z + CC, data=data)
  }
  # print(paste("sfd:", "am", pma$coef[2], "my", eym$coef[2]))

  effect <- 0
  for (i in 1:N){
    zi <- data$Z[i]
    cci <- data$CC[i]
    if (m_bin){ 
      m1_prob <- plogis(cbind(1, 1, zi, cci) %*% pma$coef)
      m1 <- rbinom(n=K, size=1, prob=m1_prob)
      m0_prob <- plogis(cbind(1, 0, zi, cci) %*% pma$coef)
      m0 <- rbinom(n=K, size=1, prob=m0_prob)
    } else {
      sd_ <- model_sd(pma, N)
      m1 <- rnorm(n=K, mean=cbind(1, 1, zi, cci) %*% pma$coef, sd=sd_)
      m0 <- rnorm(n=K, mean=cbind(1, 0, zi, cci) %*% pma$coef, sd=sd_)
    }
    # TODO why do we use z = 0 instead of z=zi?
    m1_effect <- sum(cbind(1, m1, zi, cci) %*% eym$coef)
    m0_effect <- sum(cbind(1, m0, zi, cci) %*% eym$coef)
    effect <- effect + (m1_effect - m0_effect) / K
  }

  return(effect / N)
}

backdoor <- function(data, m_bin=FALSE) {
  if (!m_bin) {
    pya <- lm(Y ~ A + U + CC + Z, data=data)
    return(as.matrix(pya$coef[2]))
  } else {
    pyamzu <- lm(Y ~ A + M + Z + U + CC, data=data)
    pmazu <- glm(M ~ A + Z + U + CC, data=data, family=binomial())

    Z <- data$Z
    U <- data$U
    CC <- data$CC
    m1 <- plogis(cbind(1, 1, Z, U, CC) %*% pmazu$coef)
    m0 <- plogis(cbind(1, 0, Z, U, CC) %*% pmazu$coef)
    y1 <- cbind(1, 1, m1, Z, U, CC) %*% pyamzu$coef
    y0 <- cbind(1, 0, m0, Z, U, CC) %*% pyamzu$coef

    return(mean(y1 - y0))
  }
}

get_simulation_params <- function(
    m_bin=FALSE, z_bin=FALSE, ay_effect=-1,
    uz_effect=NA, uw_effect=NA, zmw_effect=NA){

  min_ <- -2
  max_ <- 2
  sd_ <- 1

  ubeta <- runif(1, min_, max_)
  cbeta <- runif(2, min_, max_)
  zbeta <- runif(3, min_, max_)
  wbeta <- runif(4, min_, max_)
  abeta <- runif(4, min_, max_)
  mbeta <- runif(4, min_, max_)
  ybeta <- runif(6, min_, max_)

  if (!is.na(uz_effect)){
    zbeta[2] <- uz_effect
  }
  if (!is.na(uw_effect)){
    wbeta[2] <- uw_effect
  }
  if (!is.na(ay_effect)){
    ybeta[2] <- ay_effect
  }
  if (!is.na(zmw_effect)) {
    mbeta[3] <- zmw_effect
    wbeta[3] <- zmw_effect
  }

  # Y0's parameter order is (int, A0, U0, M0, W0, C0)
  y0 <- ybeta[1]
  ay <- ybeta[2] 
  uy <- ybeta[3]
  my <- ybeta[4]
  wy <- ybeta[5]

  # M parameter order is Int, A, Z, CC
  m0 <- mbeta[1]
  am <- mbeta[2]
  zm <- mbeta[3]

  # W parameter order is Int, U, M, CC
  mw <- wbeta[3]

  u0 <- ubeta[1]

  c0 <- cbeta[1]
  uc <- cbeta[2]

  # Z0 <- cbind(1, U0, C0) %*% zbeta
  z0 <- zbeta[1]
  uz <- zbeta[2]
  cz <- zbeta[3]

  # Compute the true effect so we can evaluate all methods
  if (!m_bin) {
      # If M is gaussian, it's easy
      true_effect <- ay + am * my + am * mw * wy
  } else if (z_bin) {
      # If M is binary and Z is binary, use hcubature to integrate
      #   over Z
      ec <- c0 + u0 * uc
      pz <- z0 + u0 * uz + ec * cz
      z_func <- function(x) {
          c <- x[1]
          u <- x[2]
          (dnorm(c, mean=ec, sd=sd_) * dnorm(u, mean=u0, sd=sd_)
           * plogis(z0 + uz * u + cz * c))
      }
      int <- hcubature(z_func, lower=c(ec - 10 * sd_, u0 - 10 * sd_),
                               upper=c(ec + 10 * sd_, u0 + 10 * sd_))
      ez <- int$integral
      
      ma1 <- plogis(m0 + am + zm) * ez + plogis(m0 + am) * (1 - ez)
      ma0 <- plogis(m0 + zm) * ez + plogis(m0) * (1 - ez)
      am_effect <- ma1 - ma0
      true_effect <- ay + (mw * wy + my) * am_effect
      # print(paste("eu", u0, "ec", ec, "pz", pz, "ez", ez))
      # print(paste("m_bin, z_bin true effect", true_effect))
  } else {
      # If M is binary and Z is gaussian, instead integrate over M
      m_func <- function(z) {
          dnorm(z, mean=ez, sd=sd_) * (plogis(am + m0 + z * zm)
                                       - plogis(m0 + z * zm))}
      am_int <- integrate(m_func, lower=ez - 10 * sd_, upper=ez + 10 * sd_)
      am_effect <- am_int$value
      true_effect <- ay + (mw * wy + my) * am_effect
      # print(paste("m_bin, z_gau true effect", true_effect))
  }

  # Debugging printouts
  # printout("ybeta", c(y0, ay, my, wy, uy),
  #          c("int", "ay", "my", "wy", "uy"))
  # printout("mbeta", mbeta, c("int", "a", "z"))
  # printout("wbeta", wbeta, c("int", "u", "m"))
  # print(paste("true", "am", am, "my+mwy", my + mw * wy))
  # print(paste("!m_bin true_effect", am * my + am * mw * wy))

  sim_parameters  <- list(
      ubeta=ubeta, mbeta=mbeta, wbeta=wbeta,
      abeta=abeta, ybeta=ybeta, zbeta=zbeta,
      cbeta=cbeta, sd_=sd_, true_effect=true_effect)
  return(sim_parameters)
}


run <- function(N, n_bootstrap, sim_parameters,
                m_bin=TRUE, z_bin=TRUE, outfn=NULL){

    # Init parameters for bridge function optimizations
    min_ = -2
    max_ = 2
    proximal_frontdoor_init <- runif(5, min_, max_)
    simple_prox_init <- runif(4, min_, max_)

    # Actually generate the dataset
    dataset <- simulate(sim_parameters, N, m_bin=m_bin, z_bin=z_bin)
    truth <- sim_parameters$true_effect
    df <- lst_to_df(dataset)

    # Run all methods without bootstrap
    results <- list(backdoor=list(), simple=list(), pfd=list(), naive=list())
    results$backdoor$estimate <- backdoor(df, m_bin=m_bin) - truth
    results$naive$estimate <- sample_frontdoor(df, m_bin=m_bin) - truth
    results$pfd$estimate <- proximal_frontdoor(df, proximal_frontdoor_init, m_bin=m_bin) - truth
    results$simple$estimate <- simple_prox(df, simple_prox_init) - truth

    if (n_bootstrap <= 1 & !is.null(outfn)) {
        # If not bootstrapping, we're done
        headers <- paste(c("method", "estimate"), collapse=',')
        write(headers, file=outfn, append=TRUE)
        write(paste(c("truth", truth), collapse=","), file=outfn, append=TRUE)
        for (key in c("backdoor", "naive", "simple", "pfd")){
            row <- paste(c(key, results[[key]][["estimate"]]), collapse=',')
            write(row, file=outfn, append=TRUE)
        }
    } else {
        # Run all methods with bootstrap 
        quantiles <- c(0.025, 0.975)

        # Oracle backdoor
        backdoor_bootstrap <- bootstrap_ci(df, n_bootstrap, backdoor)
        results$backdoor$mean <- mean(backdoor_bootstrap) - truth
        i_stats <- interval_stats(backdoor_bootstrap, truth)
        results$backdoor$coverage <- i_stats[1]
        results$backdoor$width <- i_stats[2]

        # Proximal frontdoor
        pfd_func <- function(dataset){
            return(proximal_frontdoor(dataset, proximal_frontdoor_init, m_bin=m_bin))
        }
        nc_bootstrap <- bootstrap_ci(df, n_bootstrap, pfd_func)
        results$pfd$mean <- mean(nc_bootstrap) - truth
        i_stats <- interval_stats(nc_bootstrap, truth)
        results$pfd$coverage <- i_stats[1]
        results$pfd$width <- i_stats[2]

        # Simple proximal
        simple_prox_func <- function(dataset){
            return(simple_prox(dataset, simple_prox_init))
        }
        simple_bootstrap <- bootstrap_ci(df, n_bootstrap, simple_prox_func)
        results$simple$mean <- mean(simple_bootstrap) - truth
        i_stats <- interval_stats(simple_bootstrap, truth)
        results$simple$coverage <- i_stats[1]
        results$simple$width <- i_stats[2]

        # Naive frontdoor
        naive_func <- function(dataset){
            return(sample_frontdoor(dataset, m_bin=m_bin))
        }
        naive_bootstrap <- bootstrap_ci(df, n_bootstrap, naive_func)
        results$naive$mean <- mean(naive_bootstrap) - truth
        i_stats <- interval_stats(naive_bootstrap, truth)
        results$naive$coverage <- i_stats[1]
        results$naive$width <- i_stats[2]

        # Write results to file
        if (!is.null(outfn)){
            headers <- c("method", "estimate", "mean", "width", "coverage")
            write(paste(headers, collapse=','), file=outfn, append=TRUE)
            write(paste(c("truth", truth, 0, 0, 1), collapse=","),
                  file=outfn, append=TRUE)
            for (key in c("backdoor", "naive", "simple", "pfd")){
                row <- c(key)
                for (metric in headers[2:length(headers)]){ 
                    row <- append(row, results[[key]][[metric]])
                }
                row <- paste(row, collapse=',')
                write(row, file=outfn, append=TRUE)
            }
        }
    }

    return(results)
}

format_results <- function(results, name, n_samples, metrics){
    options(scipen = 999)
    row <- c()
    for (metric in metrics){
        row <- append(row, results[[metric]])
    }
    row <- data.frame(rbind(row))
    rownames(row) <- paste(sprintf("%20s @ %7s", name, n_samples))
    colnames(row) <- metrics
    
    return(row)
}

# Default values for the argument. These can be overwritten
#     with command line arguments (e.g., "--n_bootstrap=2").
# NA values indicate that the corresponding effect is just
#     set randomly by `get_simulation_params`
defaults <- list(
    n_bootstrap=1,            # How many bootstrap samples?
    n_dgps=1,                 # How many different DGPs to sample datasets from?
    dgp_start=1,              # What DGP seed to start with (enables parallelization)
    n_datasets=1,             # How many datasets to sample from each DGP?
    dataset_start=1,          # What dataset seed to start with (enables parallelization)
    m_bin=FALSE,              # Is M binary (or Gaussian)?
    z_bin=FALSE,              # Is Z binary (or Gaussian)?
    ay_effect="NA",           # What should be the direct effect of A on Y?
    uz_effect="NA",           # What should be the direct effect of U on Z?
    uw_effect="NA",           # What should be the direct effect of U on W?
    zmw_effect="NA",          # What should be the path-specific-effect of Z->M->W?
    sample_sizes="1000,4000", # How many samples? (single int or comma-separated list)
    outdir="sim_results/",    # Where to save results?
    save=FALSE                # Should we save results?
)

main <- function() {
  
  # Command line argument processing
  args <- commandArgs(trailingOnly=TRUE, asValues=TRUE,
                      defaults=defaults)

  for (name in names(args)){
      if (name != "" & length(defaults[[name]]) == 0) {
         stop(paste("Unknown arg `", name, "`; should have default value in `defaults`."))
      }
  }
  n_dgps <- as.numeric(args$n_dgps)
  n_datasets <- as.numeric(args$n_datasets)
  n_bootstrap <- as.numeric(args$n_bootstrap)
  dgp_start <- as.numeric(args$dgp_start)
  dataset_start <- as.numeric(args$dataset_start)
  m_bin <- as.logical(args$m_bin)
  z_bin <- as.logical(args$z_bin)
  ay_effect <- suppressWarnings(as.numeric(args$ay_effect))
  uz_effect <- suppressWarnings(as.numeric(args$uz_effect))
  uw_effect <- suppressWarnings(as.numeric(args$uw_effect))
  zmw_effect <- suppressWarnings(as.numeric(args$zmw_effect))

  # sample_sizes should either be a single integer (e.g., "4000")
  #   or a string of sample sizes (e.g., "1000,4000")
  if(typeof(args$sample_sizes) == "character") {
      tmp <- sapply(strsplit(args$sample_sizes, ",")[[1]], as.numeric)
      sample_sizes <- as.vector(tmp)
  } else {
      sample_sizes <- args$sample_sizes
  }

  # Printout a summary of all experiment arguments
  dgp_str <- paste(c(dgp_start, dgp_start + n_dgps - 1), collapse=":")
  data_str <- paste(c(dataset_start, dataset_start + n_datasets - 1), collapse=":")
  print(paste("dgps", dgp_str, "data", data_str, "n_boot", n_bootstrap,
              "samples", paste(sample_sizes, collapse=",")))
  print(paste("z_bin", z_bin, "m_bin", m_bin,
              "ay", ay_effect, "uz", uz_effect, "uw", uw_effect,
              "zmw", zmw_effect))

  # If we're going to save these experiment results to be used in compile_experiments.py
  #   we need to pick the right filename
  if (args$save) {
    fn_params <- c()
    for (arg in c("z_bin", "m_bin", "ay_effect", "uz_effect", "uw_effect", "zmw_effect", "n_bootstrap")){ 
      short_key <- substr(arg, 0, 4)
      short_val <- args[[arg]]
      if (typeof(short_val) == "logical") {
        short_val <- as.integer(args[[arg]])
      }
      param <- paste(c(short_key, short_val), collapse="")
      fn_params <- append(fn_params, param)
    }
    outfn <- paste(c(args$outdir, paste(fn_params, collapse="-")), collapse="/")
  } else {
    outfn <- NULL
  }

  # Run the experiment(s)
  metrics <- c("estimate", "mean", "coverage", "width")
  first_row <- TRUE
  for (n_samples in sample_sizes) {

      results <- list()
      # Iterate through the data generating process seeds
      dgp_end <- dgp_start + n_dgps - 1
      for (dgp_seed in dgp_start:dgp_end) {
          dgp_results <- NULL
          # For the specified dgp seed, generate the parameters we will use
          #   to generate the synthetic data; but don't generate data yet.
          set.seed(dgp_seed)
          sim_parameters <- get_simulation_params(
              m_bin=m_bin, z_bin=z_bin, ay_effect=ay_effect,
              uz_effect=uz_effect, uw_effect=uw_effect, zmw_effect=zmw_effect)

          # Iterate through the dataset seeds
          dataset_end <- dataset_start + n_datasets - 1
          for (dataset_seed in dataset_start:dataset_end) {

              # If we're saving results, write headers
              if (args$save){
                  exp_params <- c(paste(c("dgp", dgp_seed), collapse=""),
                                  paste(c("data", dataset_seed), collapse=""),
                                  paste(c("samples", n_samples), collapse=""))
                  exp_outfn <- paste(c(outfn, exp_params), collapse="-")
                  exp_outfn <- paste(c(exp_outfn, "csv"), collapse=".")

                  if (file.exists(exp_outfn)){
                      print(paste("file exists:", exp_outfn))
                      next
                  }

                  # Write experiment args to the top of the results csv file
                  #   for post-hoc sanity checking
                  write(paste(c("//", toJSON(args)), collapse=" "),
                        file=exp_outfn, append=TRUE)
                  write(paste(c("//", toJSON(sim_parameters)), collapse=" "),
                        file=exp_outfn, append=TRUE)
              } else {
                  exp_outfn <- NULL
              }

              # With the specified dataset seed, actually generate the dataset
              #   of the desired sample size and run the experiment.
              set.seed(dataset_seed)
              dataset_results <- run(n_samples, n_bootstrap, sim_parameters,
                                     m_bin=m_bin, z_bin=z_bin, outfn=exp_outfn)

              # Compile results from one experiment
              if (is.null(dgp_results)){
                  dgp_results <- list()
                  for (name in names(dataset_results)) {
                      dgp_results[[name]] <- list()
                  }
              }
              for (name in names(dataset_results)) {
                  for (metric in metrics) {
                      dgp_results[[name]][[metric]] <- append(
                          dgp_results[[name]][[metric]],
                          dataset_results[[name]][[metric]])
            }}} # end for (dataset_seed ...)

          if (is.null(results)){
            results <- list()
          }

          # Compile results from all experiments with this dgp
          for (name in names(dgp_results)){
              if (is.null(results[[name]])){
                  results[[name]] <- list()
              }
              for (metric in metrics) {
                  mean_ <- suppressWarnings(mean(dgp_results[[name]][[metric]]))
                  if (metric == "estimate" | metric == "mean") {
                      pab_key <- paste(c("pab", metric), collapse=".")
                      pab_val <- (abs(mean_)) / abs(sim_parameters$true_effect)

                      results[[name]][[pab_key]] <- append(
                          results[[name]][[pab_key]], pab_val)
                  }
                  results[[name]][[metric]] <- append(
                      results[[name]][[metric]], mean_)
              }
          }
      }  # end for (dgp_seed in 1:n_dgps)

      if (is.null(results) | sum(dim(as.data.frame(results))) == 0){
          print(paste("No results to print for n_samples=", n_samples))
          next
      }

      # Aggregate results across dgp seeds and dataset seeds,
      #   all for this one sample size
      all_metrics <- append(metrics, c("pab.estimate", "pab.mean"))
      for (name in names(results)){
          for (metric in all_metrics) {
              results[[name]][[metric]] <- suppressWarnings(mean(results[[name]][[metric]]))
          }
      }

      fancy_names <- list(backdoor="backdoor", naive="naive frontdoor",
                          simple="simple prox", pfd="frontdoor prox")
      row <- data.frame(c())
      for (name in names(fancy_names)){
          cell <- format_results(results[[name]], fancy_names[[name]],
                                 n_samples, all_metrics)
          row <- rbind(row, cell)
      }

      # Print out aggregated results for this sample size
      if (first_row) {
          first_row <- FALSE
      } else {
          # Suppress column names for subsequent sample sizes
          colnames(row) <- NULL
      }
      print(format(row, digits=3, width=7, nsmall=4))
    }  # end for (n_samples in sample_sizes)
}  # end main()

main()
warnings()
