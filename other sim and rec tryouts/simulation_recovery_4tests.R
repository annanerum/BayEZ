
# Model 1: alpha=cond+slope, nu=cov, tau=intercept
# Model 2: alpha=cond+slope, nu=intercept, tau=intercept
# Model 3: alpha=cond, nu=cov, tau=intercept (no slope)
# Model 4: alpha=cond, nu=cond, tau=intercept (no cov, no slope)

library(rstan)
library(bayesplot)
library(ggplot2)
library(dplyr)
library(tidyr)
library(knitr)
library(parallel)

rstan_options(auto_write = TRUE)
options(mc.cores = 1)

if (!dir.exists("output_models")) dir.create("output_models")

## simulation grid
I_vals <- c(80, 160)
J_vals <- c(80, 160)
N_sim  <- 10
K      <- 2

## true parameters 

true_m1 <- list(
  # Model 1: condition on alpha + random slope, covariate on nu, tau intercept only
  mu_nu    = 1.0,  # drift intercept
  b1_nu    = 0.3,  # covariate effect on drift
  mu_alpha = 1.5,  # boundary intercept
  b_alpha  = 0.2,  # condition effect on boundary
  mu_tau   = 0.3,  # NDT intercept
  sigma_nu    = 0.3,
  sigma_alpha = 0.25,
  sigma_tau   = 0.15,
  sigma_v     = 0.2  # random slope SD
)

true_m2 <- list(
  # Model 2: condition on alpha + random slope, nu intercept only, tau intercept only
  mu_nu    = 1.0,
  mu_alpha = 1.5,
  b_alpha  = 0.2,
  mu_tau   = 0.3,
  sigma_nu    = 0.3,
  sigma_alpha = 0.25,
  sigma_tau   = 0.15,
  sigma_v     = 0.2
)

true_m3 <- list(
  # Model 3: condition on alpha (no slope), covariate on nu, tau intercept only
  mu_nu    = 1.0,
  b1_nu    = 0.3,
  mu_alpha = 1.5,
  b_alpha  = 0.2,
  mu_tau   = 0.3,
  sigma_nu    = 0.3,
  sigma_alpha = 0.25,
  sigma_tau   = 0.15
)

true_m4 <- list(
  # Model 4: condition on alpha and nu, tau intercept only, no cov, no slope
  mu_nu    = 1.0,
  b2_nu    = 0.5,  # condition effect on drift
  mu_alpha = 1.5,
  b_alpha  = 0.2,
  mu_tau   = 0.3,
  sigma_nu    = 0.3,
  sigma_alpha = 0.25,
  sigma_tau   = 0.15
)


## EZ point estimates
ez_point_estimates <- function(C_arr, J_arr, MRT_mat, VRT_mat, I, K) {
  nu_hat <- matrix(NA, I, K); alpha_hat <- matrix(NA, I, K); tau_hat <- matrix(NA, I, K)
  for (i in seq_len(I)) for (k in seq_len(K)) {
    Pc  <- C_arr[i,k] / J_arr[i,k]
    Pc  <- max(min(Pc, 1 - 1/(2*J_arr[i,k])), 1/(2*J_arr[i,k]))
    mrt <- MRT_mat[i,k]; vrt <- VRT_mat[i,k]
    if (is.na(mrt) || is.na(vrt) || vrt <= 0) next
    L   <- log(Pc/(1-Pc))
    nu4 <- L*(Pc^2*L - Pc*L + Pc - 0.5)/vrt
    if (is.na(nu4) || nu4 <= 0) next
    nu_hat[i,k]    <- sign(Pc-0.5)*nu4^(1/4)
    alpha_hat[i,k] <- L/nu_hat[i,k]
    tau_hat[i,k]   <- mrt - (alpha_hat[i,k]/(2*nu_hat[i,k]))*(2*Pc-1)
  }
  list(nu=nu_hat, alpha=alpha_hat, tau=tau_hat)
}


## derive priors 
derive_priors <- function(ez, P_nu, P_alpha, P_tau, has_sigma_v) {
  nu_vals    <- as.vector(ez$nu[!is.na(ez$nu)])
  alpha_vals <- as.vector(ez$alpha[!is.na(ez$alpha) & ez$alpha > 0])
  tau_vals   <- as.vector(ez$tau[!is.na(ez$tau) & ez$tau > 0.05])

  nu_sd    <- min(max(sd(nu_vals),    0.3), 1.5)
  alpha_sd <- min(max(sd(alpha_vals), 0.2), 1.0)
  tau_sd   <- min(max(sd(tau_vals),   0.1), 0.5)

  # beta_nu: 2 x P_nu matrix
  # row 1 = prior means: intercept gets pilot mean, others get 0
  # row 2 = prior SDs: all get nu_sd
  beta_nu_means <- c(mean(nu_vals), rep(0, P_nu - 1))
  beta_nu_sds   <- rep(nu_sd, P_nu)
  beta_nu       <- matrix(c(beta_nu_means, beta_nu_sds), nrow=2, byrow=TRUE)

  beta_alpha_means <- c(mean(alpha_vals), rep(0, P_alpha - 1))
  beta_alpha_sds   <- c(alpha_sd*2, rep(alpha_sd*2, P_alpha - 1))
  beta_alpha       <- matrix(c(beta_alpha_means, beta_alpha_sds), nrow=2, byrow=TRUE)

  beta_tau_means <- c(mean(tau_vals), rep(0, P_tau - 1))
  beta_tau_sds   <- c(tau_sd*2, rep(tau_sd*2, P_tau - 1))
  beta_tau       <- matrix(c(beta_tau_means, beta_tau_sds), nrow=2, byrow=TRUE)

  out <- list(
    beta_nu     = beta_nu,
    beta_alpha  = beta_alpha,
    beta_tau    = beta_tau,
    sigma_nu    = c(nu_sd,    nu_sd/2),
    sigma_alpha = c(alpha_sd, alpha_sd/2),
    sigma_tau   = c(tau_sd,   tau_sd/2)
  )
  if (has_sigma_v)
    out$sigma_v <- c(nu_sd*0.3, nu_sd)
  out
}


## pilot prior derivation 
derive_priors_from_pilot <- function(J_per, true, model_id, I_pilot=20) {
  K <- 2; cond_dummy <- c(0,1)
  J_p   <- matrix(J_per, I_pilot, K)
  C_p   <- matrix(NA_integer_, I_pilot, K)
  MRT_p <- matrix(NA_real_,    I_pilot, K)
  VRT_p <- matrix(NA_real_,    I_pilot, K)

  x_p      <- if (model_id %in% c(1,3)) rnorm(I_pilot, 0, 1) else NULL
  true_v_p <- if (model_id %in% c(1,2)) rnorm(I_pilot, 0, true$sigma_v) else NULL

  for (i in seq_len(I_pilot)) for (k in seq_len(K)) {
    mean_nu_ik <- switch(as.character(model_id),
      "1" = true$mu_nu + true$b1_nu * x_p[i],            # covariate only
      "2" = true$mu_nu,                                    # intercept only
      "3" = true$mu_nu + true$b1_nu * x_p[i],            # covariate only
      "4" = true$mu_nu + true$b2_nu * cond_dummy[k]       # condition only
    )
    mean_alpha_ik <- switch(as.character(model_id),
      "1" = true$mu_alpha + (true$b_alpha + true_v_p[i]) * cond_dummy[k],
      "2" = true$mu_alpha + (true$b_alpha + true_v_p[i]) * cond_dummy[k],
      "3" = true$mu_alpha + true$b_alpha * cond_dummy[k],
      "4" = true$mu_alpha + true$b_alpha * cond_dummy[k]
    )
    mean_tau_ik <- true$mu_tau  # always intercept only

    repeat { a <- rnorm(1, mean_alpha_ik, true$sigma_alpha)
             if (a > 0.1 && a < 3.0) { ap <- a; break } }
    np <- rnorm(1, mean_nu_ik, true$sigma_nu)
    repeat { t <- rnorm(1, mean_tau_ik, true$sigma_tau)
             if (t > 0.05) { tp <- t; break } }

    e      <- exp(-ap*np); pi_c <- 1/(1+e)
    mu_rt  <- tp + (ap/(2*np))*((1-e)/(1+e))
    sig2rt <- (ap/(2*np^3))*((1-2*ap*np*e-e^2)/(1+e)^2)
    C_ik   <- max(1L, min(rbinom(1,J_per,pi_c), J_per-1L))
    MRT_p[i,k] <- rnorm(1, mu_rt, sqrt(sig2rt/C_ik))
    VRT_p[i,k] <- max(rnorm(1, sig2rt, sqrt(2*sig2rt^2/max(C_ik-1,1))), 1e-6)
    C_p[i,k]   <- C_ik
  }

  ez <- ez_point_estimates(C_p, J_p, MRT_p, VRT_p, I_pilot, K)
  P_nu    <- switch(as.character(model_id), "1"=2, "2"=1, "3"=2, "4"=2)
  P_alpha <- 2
  P_tau   <- 1
  has_v   <- model_id %in% c(1, 2)
  derive_priors(ez, P_nu, P_alpha, P_tau, has_v)
}


## run one simulation and fit
run_one <- function(I, J_per, true, model, priors, model_id,
                    rhat_threshold=1.05, hopeless_threshold=10,
                    iter1=1000, warmup1=300, iter2=3000, warmup2=1000) {

  K <- 2; cond_dummy <- c(0,1)
  prs <- priors

  # draw covariate and random slopes depending on model
  x      <- if (model_id %in% c(1,3)) rnorm(I, 0, 1) else NULL
  true_v <- if (model_id %in% c(1,2)) rnorm(I, 0, true$sigma_v) else NULL

  # design matrices 
  P_nu    <- switch(as.character(model_id), "1"=2, "2"=1, "3"=2, "4"=2)
  P_alpha <- 2
  P_tau   <- 1

  X_nu    <- array(NA_real_, dim=c(I, K, P_nu))
  X_alpha <- array(NA_real_, dim=c(I, K, P_alpha))
  X_tau   <- array(NA_real_, dim=c(I, K, P_tau))

  for (i in seq_len(I)) for (k in seq_len(K)) {
    X_nu[i,k,] <- switch(as.character(model_id),
      "1" = c(1, x[i]),           # intercept + covariate
      "2" = c(1),                  # intercept only
      "3" = c(1, x[i]),           # intercept + covariate
      "4" = c(1, cond_dummy[k])   # intercept + condition
    )
    X_alpha[i,k,] <- c(1, cond_dummy[k])  # always intercept + condition
    X_tau[i,k,]   <- c(1)                 # always intercept only
  }

  # true individual parameter generation 
  alpha_true <- matrix(NA_real_, I, K)
  nu_true    <- matrix(NA_real_, I, K)
  tau_true   <- matrix(NA_real_, I, K)

  for (i in seq_len(I)) for (k in seq_len(K)) {
    mean_nu_ik <- switch(as.character(model_id),
      "1" = true$mu_nu + true$b1_nu * x[i],
      "2" = true$mu_nu,
      "3" = true$mu_nu + true$b1_nu * x[i],
      "4" = true$mu_nu + true$b2_nu * cond_dummy[k]
    )
    mean_alpha_ik <- switch(as.character(model_id),
      "1" = true$mu_alpha + (true$b_alpha + true_v[i]) * cond_dummy[k],
      "2" = true$mu_alpha + (true$b_alpha + true_v[i]) * cond_dummy[k],
      "3" = true$mu_alpha + true$b_alpha * cond_dummy[k],
      "4" = true$mu_alpha + true$b_alpha * cond_dummy[k]
    )
    mean_tau_ik <- true$mu_tau

    repeat { a <- rnorm(1, mean_alpha_ik, true$sigma_alpha)
             if (a > 0.1 && a < 3.0) { alpha_true[i,k] <- a; break } }
    nu_true[i,k] <- rnorm(1, mean_nu_ik, true$sigma_nu)
    repeat { t <- rnorm(1, mean_tau_ik, true$sigma_tau)
             if (t > 0.05) { tau_true[i,k] <- t; break } }
  }

  # observed data
  J_mat <- matrix(J_per, I, K)
  C_mat <- matrix(NA_integer_, I, K)
  MRT_mat <- matrix(NA_real_, I, K)
  VRT_mat <- matrix(NA_real_, I, K)

  for (i in seq_len(I)) for (k in seq_len(K)) {
    a <- alpha_true[i,k]; v <- nu_true[i,k]; tau_ik <- tau_true[i,k]
    e      <- exp(-a*v); pi_c <- 1/(1+e)
    mu_rt  <- tau_ik + (a/(2*v))*((1-e)/(1+e))
    sig2rt <- (a/(2*v^3))*((1-2*a*v*e-e^2)/(1+e)^2)
    C_ik   <- max(1L, min(rbinom(1,J_per,pi_c), J_per-1L))
    MRT_mat[i,k] <- rnorm(1, mu_rt, sqrt(sig2rt/C_ik))
    VRT_mat[i,k] <- max(rnorm(1, sig2rt, sqrt(2*sig2rt^2/max(C_ik-1,1))), 1e-6)
    C_mat[i,k]   <- C_ik
  }

  # stan data 
  stan_data <- list(
    I=I, K=K, P_nu=P_nu, P_alpha=P_alpha, P_tau=P_tau,
    J=J_mat, C=C_mat, MRT=MRT_mat, VRT=VRT_mat,
    X_nu=X_nu, X_alpha=X_alpha, X_tau=X_tau,
    prior_beta_nu     = prs$beta_nu,
    prior_beta_alpha  = prs$beta_alpha,
    prior_beta_tau    = prs$beta_tau,
    prior_sigma_nu    = prs$sigma_nu,
    prior_sigma_alpha = prs$sigma_alpha,
    prior_sigma_tau   = prs$sigma_tau
  )
  if (model_id %in% c(1,2))
    stan_data$prior_sigma_v <- prs$sigma_v

  # group parameters to monitor for Rhat
  group_pars <- c(
    paste0("beta_nu[",    seq_len(P_nu),    "]"),
    paste0("beta_alpha[", seq_len(P_alpha), "]"),
    paste0("beta_tau[",   seq_len(P_tau),   "]"),
    "sigma_nu", "sigma_alpha", "sigma_tau"
  )
  if (model_id %in% c(1,2)) group_pars <- c(group_pars, "sigma_v")

  hopeless <- function(fit) {
    rh <- tryCatch(summary(fit)$summary[group_pars, "Rhat"], error=function(e) NA_real_)
    any(is.na(rh)) || any(is.nan(rh)) || any(is.infinite(rh)) ||
      any(rh > hopeless_threshold, na.rm=TRUE)
  }

  fit <- tryCatch(
    sampling(model, data=stan_data,
             iter=iter1, warmup=warmup1, chains=4, cores=1, refresh=0,
             control=list(adapt_delta=0.85, max_treedepth=10)),
    error=function(e) NULL)

  if (is.null(fit)) return(list(status="error", I=I, J=J_per, model_id=model_id))
  if (hopeless(fit)) return(list(status="excluded", reason="degenerate_fit1",
                                  I=I, J=J_per, model_id=model_id))

  rhat_max  <- max(summary(fit)$summary[group_pars, "Rhat"], na.rm=TRUE)
  converged <- "first_fit"

  if (rhat_max > rhat_threshold) {
    message(sprintf("Refit triggered (model=%d, I=%d, J=%d) -- Rhat=%.2f",
                    model_id, I, J_per, rhat_max))
    fit <- tryCatch(
      sampling(model, data=stan_data,
               iter=iter2, warmup=warmup2, chains=4, cores=1, refresh=0,
               control=list(adapt_delta=0.95, max_treedepth=12)),
      error=function(e) NULL)

    if (is.null(fit)) return(list(status="error", I=I, J=J_per, model_id=model_id))
    if (hopeless(fit)) return(list(status="excluded", reason="degenerate_refit",
                                    I=I, J=J_per, model_id=model_id))
    rhat_max  <- max(summary(fit)$summary[group_pars, "Rhat"], na.rm=TRUE)
    converged <- "refit"
    if (rhat_max > rhat_threshold)
      return(list(status="excluded", reason="non_converged",
                  I=I, J=J_per, model_id=model_id))
  }

  s <- summary(fit)$summary

  # extract mean_ik from generated quantities --> design matrix verification
  # what the model thinks the predicted mean was for each ppn x condition
  mean_nu_est    <- matrix(NA_real_, I, K)
  mean_alpha_est <- matrix(NA_real_, I, K)
  mean_tau_est   <- matrix(NA_real_, I, K)
  for (i in seq_len(I)) for (k in seq_len(K)) {
    mean_nu_est[i,k]    <- s[paste0("mean_nu_ik[",    i, ",", k, "]"), "mean"]
    mean_alpha_est[i,k] <- s[paste0("mean_alpha_ik[", i, ",", k, "]"), "mean"]
    mean_tau_est[i,k]   <- s[paste0("mean_tau_ik[",   i, ",", k, "]"), "mean"]
  }

  # group-level results
  param_names <- switch(as.character(model_id),
    "1" = c("mu_nu", "b1_nu", "mu_alpha", "b_alpha", "mu_tau",
            "sigma_nu", "sigma_alpha", "sigma_tau", "sigma_v"),
    "2" = c("mu_nu", "mu_alpha", "b_alpha", "mu_tau",
            "sigma_nu", "sigma_alpha", "sigma_tau", "sigma_v"),
    "3" = c("mu_nu", "b1_nu", "mu_alpha", "b_alpha", "mu_tau",
            "sigma_nu", "sigma_alpha", "sigma_tau"),
    "4" = c("mu_nu", "b2_nu", "mu_alpha", "b_alpha", "mu_tau",
            "sigma_nu", "sigma_alpha", "sigma_tau")
  )

  true_vals <- switch(as.character(model_id),
    "1" = c(true$mu_nu, true$b1_nu, true$mu_alpha, true$b_alpha, true$mu_tau,
            true$sigma_nu, true$sigma_alpha, true$sigma_tau, true$sigma_v),
    "2" = c(true$mu_nu, true$mu_alpha, true$b_alpha, true$mu_tau,
            true$sigma_nu, true$sigma_alpha, true$sigma_tau, true$sigma_v),
    "3" = c(true$mu_nu, true$b1_nu, true$mu_alpha, true$b_alpha, true$mu_tau,
            true$sigma_nu, true$sigma_alpha, true$sigma_tau),
    "4" = c(true$mu_nu, true$b2_nu, true$mu_alpha, true$b_alpha, true$mu_tau,
            true$sigma_nu, true$sigma_alpha, true$sigma_tau)
  )

  stan_par_names <- switch(as.character(model_id),
    "1" = c("beta_nu[1]", "beta_nu[2]", "beta_alpha[1]", "beta_alpha[2]",
            "beta_tau[1]", "sigma_nu", "sigma_alpha", "sigma_tau", "sigma_v"),
    "2" = c("beta_nu[1]", "beta_alpha[1]", "beta_alpha[2]",
            "beta_tau[1]", "sigma_nu", "sigma_alpha", "sigma_tau", "sigma_v"),
    "3" = c("beta_nu[1]", "beta_nu[2]", "beta_alpha[1]", "beta_alpha[2]",
            "beta_tau[1]", "sigma_nu", "sigma_alpha", "sigma_tau"),
    "4" = c("beta_nu[1]", "beta_nu[2]", "beta_alpha[1]", "beta_alpha[2]",
            "beta_tau[1]", "sigma_nu", "sigma_alpha", "sigma_tau")
  )

  group_results <- data.frame(
    model_id  = model_id,
    I = I, J = J_per,
    param     = param_names,
    true_val  = true_vals,
    estimate  = s[stan_par_names, "mean"],
    converged = converged,
    rhat_max  = rhat_max,
    stringsAsFactors = FALSE
  )

  # individual results
  alpha_est <- matrix(NA_real_, I, K)
  nu_est    <- matrix(NA_real_, I, K)
  tau_est   <- matrix(NA_real_, I, K)
  for (i in seq_len(I)) for (k in seq_len(K)) {
    alpha_est[i,k] <- s[paste0("alpha[",i,",",k,"]"), "mean"]
    nu_est[i,k]    <- s[paste0("nu[",   i,",",k,"]"), "mean"]
    tau_est[i,k]   <- s[paste0("tau[",  i,",",k,"]"), "mean"]
  }

  indiv_results <- data.frame(
    model_id  = model_id,
    I = I, J = J_per,
    condition  = rep(seq_len(K), each=I),
    true_alpha = as.vector(alpha_true), est_alpha = as.vector(alpha_est),
    true_nu    = as.vector(nu_true),    est_nu    = as.vector(nu_est),
    true_tau   = as.vector(tau_true),   est_tau   = as.vector(tau_est)
  )

  # mean_ik results for design matrix verification
  mean_ik_results <- data.frame(
    model_id  = model_id,
    I = I, J = J_per,
    participant = rep(seq_len(I), times=K),
    condition   = rep(seq_len(K), each=I),
    # estimated predicted means from generated quantities
    est_mean_nu    = as.vector(mean_nu_est),
    est_mean_alpha = as.vector(mean_alpha_est),
    est_mean_tau   = as.vector(mean_tau_est),
    # true predicted means for comparison
    true_mean_nu = as.vector(matrix(
      sapply(seq_len(K), function(k)
        sapply(seq_len(I), function(i) {
          switch(as.character(model_id),
            "1" = true$mu_nu + true$b1_nu * x[i],
            "2" = true$mu_nu,
            "3" = true$mu_nu + true$b1_nu * x[i],
            "4" = true$mu_nu + true$b2_nu * cond_dummy[k]
          )
        })), nrow=I)),
    true_mean_alpha = as.vector(matrix(
      sapply(seq_len(K), function(k)
        sapply(seq_len(I), function(i) {
          switch(as.character(model_id),
            "1" = true$mu_alpha + (true$b_alpha + true_v[i]) * cond_dummy[k],
            "2" = true$mu_alpha + (true$b_alpha + true_v[i]) * cond_dummy[k],
            "3" = true$mu_alpha + true$b_alpha * cond_dummy[k],
            "4" = true$mu_alpha + true$b_alpha * cond_dummy[k]
          )
        })), nrow=I)),
    true_mean_tau = true$mu_tau
  )

  # random slope results if applicable
  v_results <- if (model_id %in% c(1,2)) {
    data.frame(
      model_id    = model_id,
      I = I, J = J_per,
      participant = seq_len(I),
      true_v      = true_v,
      est_v       = s[paste0("v[", seq_len(I), "]"), "mean"]
    )
  } else NULL

  list(status="ok", group=group_results, indiv=indiv_results,
       mean_ik=mean_ik_results, v=v_results,
       rhat_max=rhat_max, converged=converged,
       I=I, J=J_per, model_id=model_id)
}


## compile all four Stan models
stan_files <- c(
  "model1_cond_alpha_cov_nu.stan",
  "model2_cond_alpha_no_cov.stan",
  "model3_cond_alpha_cov_nu_no_slope.stan",
  "model4_cond_alpha_nu_no_cov_no_slope.stan"
)
for (f in stan_files) {
  cache <- paste0(tools::file_path_sans_ext(f), ".rds")
  if (file.exists(cache)) file.remove(cache)
}

true_list <- list(true_m1, true_m2, true_m3, true_m4)

## derive priors per model per J
set.seed(42)
priors_by_model_J <- lapply(seq_along(stan_files), function(mid) {
  setNames(
    lapply(J_vals, function(j)
      derive_priors_from_pilot(j, true_list[[mid]], mid, I_pilot=20)),
    as.character(J_vals)
  )
})

## build job list across all 4 models
jobs <- do.call(rbind, lapply(seq_along(stan_files), function(mid) {
  g <- expand.grid(I=I_vals, J=J_vals, sim=seq_len(N_sim))
  g$model_id <- mid
  g
}))
jobs_list <- split(jobs, seq_len(nrow(jobs)))

## parallel run
cl <- makeCluster(detectCores() - 1)
clusterSetRNGStream(cl, 42)

clusterExport(cl, varlist=c(
  "run_one", "derive_priors_from_pilot", "ez_point_estimates", "derive_priors",
  "true_m1", "true_m2", "true_m3", "true_m4", "true_list",
  "priors_by_model_J", "stan_files", "K"
))

clusterEvalQ(cl, {
  library(rstan)
  library(dplyr)
  rstan_options(auto_write=FALSE)
  options(mc.cores=1)
  stan_models <<- lapply(stan_files, stan_model)
})

results_list <- parLapply(cl, jobs_list, function(row) {
  mid <- row$model_id
  run_one(
    I        = row$I,
    J_per    = row$J,
    true     = true_list[[mid]],
    model    = stan_models[[mid]],
    priors   = priors_by_model_J[[mid]][[as.character(row$J)]],
    model_id = mid
  )
})

stopCluster(cl)
saveRDS(results_list, "output_models/simulation_results.rds")


## post-processing
ok_idx     <- sapply(results_list, function(r) !is.null(r) && r$status=="ok")
ok_results <- results_list[ok_idx]

results_group  <- bind_rows(lapply(ok_results, `[[`, "group"))
results_indiv  <- bind_rows(lapply(ok_results, `[[`, "indiv"))
results_mean_ik <- bind_rows(lapply(ok_results, `[[`, "mean_ik"))
results_v      <- bind_rows(Filter(Negate(is.null),
                                   lapply(ok_results, `[[`, "v")))

# exclusion report per model and cell
n_total    <- length(results_list)
n_ok       <- sum(ok_idx)
n_excluded <- sum(sapply(results_list, function(r) !is.null(r) && r$status=="excluded"))
n_error    <- sum(sapply(results_list, function(r) is.null(r) || r$status=="error"))
cat("\n OVERALL \n")
print(data.frame(Total=n_total, Converged=n_ok, Excluded=n_excluded, Error=n_error,
                 `Excluded (%)`=round(100*n_excluded/n_total,1), check.names=FALSE))

status_df <- bind_rows(lapply(seq_along(results_list), function(i) {
  r <- results_list[[i]]; row <- jobs[i,]
  data.frame(model_id=row$model_id, I=row$I, J=row$J, sim=row$sim,
             status = if (is.null(r)) "error" else r$status,
             reason = if (is.null(r) || r$status=="ok") NA else r$reason)
}))

cat("\n EXCLUSIONS BY MODEL AND CELL \n")
status_df %>%
  group_by(model_id, I, J) %>%
  summarise(n=n(), excluded=sum(status=="excluded"),
            pct=round(100*excluded/n,1),
            reasons=paste(table(reason[status=="excluded"]), collapse=", "),
            .groups="drop") %>%
  print(n=Inf)

# mean_ik verification -- check that estimated means match true means
cat("Correlation between true and estimated predicted means:\n")
for (mid in 1:4) {
  df <- results_mean_ik %>% filter(model_id == mid)
  if (nrow(df) == 0) next
  cat(sprintf("  Model %d | nu: r=%.3f | alpha: r=%.3f | tau: r=%.3f\n",
              mid,
              cor(df$true_mean_nu,    df$est_mean_nu,    use="complete.obs"),
              cor(df$true_mean_alpha, df$est_mean_alpha, use="complete.obs"),
              cor(df$true_mean_tau,   df$est_mean_tau,   use="complete.obs")))
}

# factor labels for plots
model_labels <- c(
  "1" = "M1: alpha=cond+slope, nu=cov",
  "2" = "M2: alpha=cond+slope, nu=intercept",
  "3" = "M3: alpha=cond, nu=cov",
  "4" = "M4: alpha=cond, nu=cond"
)

results_group <- results_group %>%
  mutate(
    model_label = model_labels[as.character(model_id)],
    I_label = factor(paste0("I = ", I), levels=paste0("I = ", I_vals)),
    J_label = factor(paste0("J = ", J), levels=paste0("J = ", J_vals))
  )

results_indiv <- results_indiv %>%
  mutate(
    model_label = model_labels[as.character(model_id)],
    I_label = factor(paste0("I = ", I), levels=paste0("I = ", I_vals)),
    J_label = factor(paste0("J = ", J), levels=paste0("J = ", J_vals)),
    condition = factor(condition, labels=c("Condition 1", "Condition 2"))
  )

dark_theme <- theme_minimal(base_size=12) +
  theme(
    plot.background   = element_rect(fill="#1e1e1e", colour=NA),
    panel.background  = element_rect(fill="#1e1e1e", colour=NA),
    panel.grid.major  = element_line(colour="#333333"),
    panel.grid.minor  = element_blank(),
    text              = element_text(colour="white"),
    axis.text         = element_text(colour="white"),
    strip.text        = element_text(colour="white"),
    legend.background = element_rect(fill="#1e1e1e"),
    legend.text       = element_text(colour="white")
  )

save_plot <- function(p, filename, width=10, height=8) {
  ggsave(file.path("output_models", filename), p,
         width=width, height=height, dpi=300, bg="#1e1e1e")
}

# individual recovery plots per model
for (mid in 1:4) {
  df <- results_indiv %>% filter(model_id == mid)
  if (nrow(df) == 0) next
  label <- model_labels[as.character(mid)]

  p_alpha <- ggplot(df, aes(x=true_alpha, y=est_alpha, colour=condition)) +
    geom_abline(intercept=0, slope=1, linetype="dashed", colour="white", linewidth=0.5) +
    geom_point(alpha=0.35, size=0.9) +
    scale_colour_manual(values=c("Condition 1"="goldenrod","Condition 2"="steelblue")) +
    facet_grid(I_label ~ J_label) +
    labs(title=paste0(label, " | alpha"), colour="Condition") +
    dark_theme + theme(legend.position="bottom")

  p_nu <- ggplot(df, aes(x=true_nu, y=est_nu, colour=condition)) +
    geom_abline(intercept=0, slope=1, linetype="dashed", colour="white", linewidth=0.5) +
    geom_point(alpha=0.35, size=0.9) +
    scale_colour_manual(values=c("Condition 1"="goldenrod","Condition 2"="steelblue")) +
    facet_grid(I_label ~ J_label) +
    coord_cartesian(xlim=c(-1, NA), ylim=c(-1, NA)) +
    labs(title=paste0(label, " | nu"), colour="Condition") +
    dark_theme + theme(legend.position="bottom")

  p_tau <- ggplot(df, aes(x=true_tau, y=est_tau, colour=condition)) +
    geom_abline(intercept=0, slope=1, linetype="dashed", colour="white", linewidth=0.5) +
    geom_point(alpha=0.35, size=0.9) +
    scale_colour_manual(values=c("Condition 1"="goldenrod","Condition 2"="steelblue")) +
    facet_grid(I_label ~ J_label) +
    labs(title=paste0(label, " | tau"), colour="Condition") +
    dark_theme + theme(legend.position="bottom")

  save_plot(p_alpha, sprintf("m%d_recovery_alpha.png", mid))
  save_plot(p_nu,    sprintf("m%d_recovery_nu.png",    mid))
  save_plot(p_tau,   sprintf("m%d_recovery_tau.png",   mid))
}

# random slope recovery for models 1 and 2
if (nrow(results_v) > 0) {
  results_v <- results_v %>%
    mutate(model_label = model_labels[as.character(model_id)],
           I_label = factor(paste0("I = ", I), levels=paste0("I = ", I_vals)),
           J_label = factor(paste0("J = ", J), levels=paste0("J = ", J_vals)))

  p_v <- ggplot(results_v, aes(x=true_v, y=est_v)) +
    geom_abline(intercept=0, slope=1, linetype="dashed", colour="white", linewidth=0.5) +
    geom_point(colour="goldenrod", alpha=0.4, size=1.2) +
    facet_grid(I_label + model_label ~ J_label) +
    labs(title="Random slope recovery v[i]", x="True v", y="Estimated v") +
    dark_theme
  save_plot(p_v, "recovery_v_models1_2.png", width=10, height=12)
}

cat("\nDone. Results saved to output_models/\n")
