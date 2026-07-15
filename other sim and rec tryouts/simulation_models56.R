# Six model comparison simulation
# Model 5: alpha=cond, nu=cond+covariate, tau=intercept, no random slope
# Model 6: alpha=cond, nu=cond+random slope on nu, tau=intercept, no covariate

library(rstan)
library(bayesplot)
library(ggplot2)
library(dplyr)
library(tidyr)
library(knitr)
library(parallel)

rstan_options(auto_write = TRUE)
options(mc.cores = 1)

if (!dir.exists("output_models56")) dir.create("output_models56")

## simulation grid
I_vals <- c(80, 160)
J_vals <- c(80, 160)
N_sim  <- 10
K      <- 2

## true parameters -- one set per model

true_m5 <- list(
  # Model 5: condition + covariate on nu, condition on alpha, tau intercept only
  # alpha higher in cond 2 (b_alpha positive), nu lower in cond 2 (b2_nu negative)
  mu_nu    = 1.0,   # drift intercept (condition 1)
  b1_nu    = 0.3,   # covariate effect on nu
  b2_nu    = -0.3,  # condition effect on nu -- NEGATIVE: drift lower in cond 2
  mu_alpha = 1.5,   # boundary intercept (condition 1)
  b_alpha  = 0.2,   # condition effect on alpha -- POSITIVE: boundary wider in cond 2
  mu_tau   = 0.3,   # NDT intercept
  sigma_nu    = 0.3,
  sigma_alpha = 0.25,
  sigma_tau   = 0.15
  # no sigma_v -- no random slope
)

true_m6 <- list(
  # Model 6: condition on nu + random slope on nu condition effect,
  #          condition on alpha (fixed), tau intercept only, no covariate
  # same direction: alpha higher in cond 2, nu lower in cond 2
  mu_nu    = 1.0,   # drift intercept (condition 1)
  b2_nu    = -0.3,  # condition effect on nu -- NEGATIVE: drift lower in cond 2
  mu_alpha = 1.5,   # boundary intercept (condition 1)
  b_alpha  = 0.2,   # condition effect on alpha -- POSITIVE: boundary wider in cond 2
  mu_tau   = 0.3,   # NDT intercept
  sigma_nu    = 0.3,
  sigma_alpha = 0.25,
  sigma_tau   = 0.15,
  sigma_v     = 0.2  # SD of random slopes on nu condition effect
  # v[i] shifts how much each person's drift changes in condition 2
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


## derive priors -- flexible for different P sizes
derive_priors <- function(ez, P_nu, P_alpha, P_tau, has_sigma_v) {
  nu_vals    <- as.vector(ez$nu[!is.na(ez$nu)])
  alpha_vals <- as.vector(ez$alpha[!is.na(ez$alpha) & ez$alpha > 0])
  tau_vals   <- as.vector(ez$tau[!is.na(ez$tau) & ez$tau > 0.05])

  nu_sd    <- min(max(sd(nu_vals),    0.3), 1.5)
  alpha_sd <- min(max(sd(alpha_vals), 0.2), 1.0)
  tau_sd   <- min(max(sd(tau_vals),   0.1), 0.5)

  beta_nu_means <- c(mean(nu_vals), rep(0, P_nu - 1))
  beta_nu_sds   <- rep(nu_sd, P_nu)
  beta_nu       <- matrix(c(beta_nu_means, beta_nu_sds), nrow=2, byrow=TRUE)

  beta_alpha_means <- c(mean(alpha_vals), rep(0, P_alpha - 1))
  beta_alpha_sds   <- rep(alpha_sd*2, P_alpha)
  beta_alpha       <- matrix(c(beta_alpha_means, beta_alpha_sds), nrow=2, byrow=TRUE)

  beta_tau_means <- c(mean(tau_vals), rep(0, P_tau - 1))
  beta_tau_sds   <- rep(tau_sd*2, P_tau)
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


## pilot prior derivation -- model-specific data generation
derive_priors_from_pilot <- function(J_per, true, model_id, I_pilot=20) {
  K <- 2; cond_dummy <- c(0,1)
  J_p   <- matrix(J_per, I_pilot, K)
  C_p   <- matrix(NA_integer_, I_pilot, K)
  MRT_p <- matrix(NA_real_,    I_pilot, K)
  VRT_p <- matrix(NA_real_,    I_pilot, K)

  # model 5 has covariate, model 6 has random slope on nu
  x_p      <- if (model_id == 5) rnorm(I_pilot, 0, 1) else NULL
  true_v_p <- if (model_id == 6) rnorm(I_pilot, 0, true$sigma_v) else NULL

  for (i in seq_len(I_pilot)) for (k in seq_len(K)) {
    mean_nu_ik <- switch(as.character(model_id),
      "5" = true$mu_nu + true$b1_nu * x_p[i] + true$b2_nu * cond_dummy[k],
      # intercept + covariate + condition (negative condition effect)
      "6" = true$mu_nu + (true$b2_nu + true_v_p[i]) * cond_dummy[k]
      # intercept + condition + random slope (in cond 1 random slope is zero)
    )
    # alpha: fixed condition effect only for both models
    mean_alpha_ik <- true$mu_alpha + true$b_alpha * cond_dummy[k]
    mean_tau_ik   <- true$mu_tau

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
  P_nu    <- switch(as.character(model_id), "5"=3, "6"=2)
  P_alpha <- 2
  P_tau   <- 1
  has_v   <- model_id == 6
  derive_priors(ez, P_nu, P_alpha, P_tau, has_v)
}


## run one simulation and fit Stan
run_one <- function(I, J_per, true, model, priors, model_id,
                    rhat_threshold=1.05, hopeless_threshold=10,
                    iter1=1000, warmup1=300, iter2=3000, warmup2=1000) {

  K <- 2; cond_dummy <- c(0,1)
  prs <- priors

  x      <- if (model_id == 5) rnorm(I, 0, 1) else NULL
  true_v <- if (model_id == 6) rnorm(I, 0, true$sigma_v) else NULL

  P_nu    <- switch(as.character(model_id), "5"=3, "6"=2)
  P_alpha <- 2
  P_tau   <- 1

  X_nu    <- array(NA_real_, dim=c(I, K, P_nu))
  X_alpha <- array(NA_real_, dim=c(I, K, P_alpha))
  X_tau   <- array(NA_real_, dim=c(I, K, P_tau))

  for (i in seq_len(I)) for (k in seq_len(K)) {
    X_nu[i,k,] <- switch(as.character(model_id),
      "5" = c(1, x[i], cond_dummy[k]),  # intercept + covariate + condition
      "6" = c(1, cond_dummy[k])          # intercept + condition (no covariate)
    )
    X_alpha[i,k,] <- c(1, cond_dummy[k])  # intercept + condition for both
    X_tau[i,k,]   <- c(1)                 # intercept only for both
  }

  # true individual parameters
  alpha_true <- matrix(NA_real_, I, K)
  nu_true    <- matrix(NA_real_, I, K)
  tau_true   <- matrix(NA_real_, I, K)

  for (i in seq_len(I)) for (k in seq_len(K)) {
    mean_nu_ik <- switch(as.character(model_id),
      "5" = true$mu_nu + true$b1_nu * x[i] + true$b2_nu * cond_dummy[k],
      "6" = true$mu_nu + (true$b2_nu + true_v[i]) * cond_dummy[k]
    )
    mean_alpha_ik <- true$mu_alpha + true$b_alpha * cond_dummy[k]
    mean_tau_ik   <- true$mu_tau

    repeat { a <- rnorm(1, mean_alpha_ik, true$sigma_alpha)
             if (a > 0.1 && a < 3.0) { alpha_true[i,k] <- a; break } }
    nu_true[i,k] <- rnorm(1, mean_nu_ik, true$sigma_nu)
    repeat { t <- rnorm(1, mean_tau_ik, true$sigma_tau)
             if (t > 0.05) { tau_true[i,k] <- t; break } }
  }

  # observed data
  J_mat   <- matrix(J_per, I, K)
  C_mat   <- matrix(NA_integer_, I, K)
  MRT_mat <- matrix(NA_real_,    I, K)
  VRT_mat <- matrix(NA_real_,    I, K)

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
  if (model_id == 6)
    stan_data$prior_sigma_v <- prs$sigma_v

  group_pars <- c(
    paste0("beta_nu[",    seq_len(P_nu),    "]"),
    paste0("beta_alpha[", seq_len(P_alpha), "]"),
    paste0("beta_tau[",   seq_len(P_tau),   "]"),
    "sigma_nu", "sigma_alpha", "sigma_tau"
  )
  if (model_id == 6) group_pars <- c(group_pars, "sigma_v")

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

  # extract mean_ik for design matrix verification
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
    "5" = c("mu_nu", "b1_nu", "b2_nu", "mu_alpha", "b_alpha", "mu_tau",
            "sigma_nu", "sigma_alpha", "sigma_tau"),
    "6" = c("mu_nu", "b2_nu", "mu_alpha", "b_alpha", "mu_tau",
            "sigma_nu", "sigma_alpha", "sigma_tau", "sigma_v")
  )

  true_vals <- switch(as.character(model_id),
    "5" = c(true$mu_nu, true$b1_nu, true$b2_nu, true$mu_alpha, true$b_alpha,
            true$mu_tau, true$sigma_nu, true$sigma_alpha, true$sigma_tau),
    "6" = c(true$mu_nu, true$b2_nu, true$mu_alpha, true$b_alpha, true$mu_tau,
            true$sigma_nu, true$sigma_alpha, true$sigma_tau, true$sigma_v)
  )

  stan_par_names <- switch(as.character(model_id),
    "5" = c("beta_nu[1]", "beta_nu[2]", "beta_nu[3]",
            "beta_alpha[1]", "beta_alpha[2]",
            "beta_tau[1]", "sigma_nu", "sigma_alpha", "sigma_tau"),
    "6" = c("beta_nu[1]", "beta_nu[2]",
            "beta_alpha[1]", "beta_alpha[2]",
            "beta_tau[1]", "sigma_nu", "sigma_alpha", "sigma_tau", "sigma_v")
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

  # mean_ik for design matrix verification
  # true predicted means computed from true parameters
  true_mean_nu_vec <- as.vector(matrix(
    sapply(seq_len(K), function(k)
      sapply(seq_len(I), function(i) {
        switch(as.character(model_id),
          "5" = true$mu_nu + true$b1_nu * x[i] + true$b2_nu * cond_dummy[k],
          "6" = true$mu_nu + true$b2_nu * cond_dummy[k]
          # note: true_v[i] not included here because it's the fixed-effects mean
          # the random slope adds individual variation around this mean
        )
      })), nrow=I))

  true_mean_alpha_vec <- as.vector(matrix(
    sapply(seq_len(K), function(k)
      sapply(seq_len(I), function(i)
        true$mu_alpha + true$b_alpha * cond_dummy[k]
      )), nrow=I))

  mean_ik_results <- data.frame(
    model_id       = model_id,
    I = I, J = J_per,
    participant    = rep(seq_len(I), times=K),
    condition      = rep(seq_len(K), each=I),
    est_mean_nu    = as.vector(mean_nu_est),
    est_mean_alpha = as.vector(mean_alpha_est),
    est_mean_tau   = as.vector(mean_tau_est),
    true_mean_nu   = true_mean_nu_vec,
    true_mean_alpha = true_mean_alpha_vec,
    true_mean_tau  = true$mu_tau
  )

  # random slope results for model 6
  v_results <- if (model_id == 6) {
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


## compile Stan models
stan_files <- c(
  "model5_cond_alpha_cond_cov_nu_no_slope.stan",
  "model6_cond_alpha_cond_nu_slope_nu.stan"
)
for (f in stan_files) {
  cache <- paste0(tools::file_path_sans_ext(f), ".rds")
  if (file.exists(cache)) file.remove(cache)
}

true_list <- list(true_m5, true_m6)
# index 1 = model 5, index 2 = model 6
# model_id in run_one is 5 or 6, list index is model_id - 4

## derive priors per model per J
set.seed(42)
priors_by_model_J <- lapply(c(5,6), function(mid) {
  setNames(
    lapply(J_vals, function(j)
      derive_priors_from_pilot(j, true_list[[mid-4]], mid, I_pilot=20)),
    as.character(J_vals)
  )
})
# priors_by_model_J[[1]] = model 5, priors_by_model_J[[2]] = model 6

## build job list
jobs <- do.call(rbind, lapply(c(5,6), function(mid) {
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
  "true_m5", "true_m6", "true_list",
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
  mid       <- row$model_id
  list_idx  <- mid - 4  # model 5 -> index 1, model 6 -> index 2
  run_one(
    I        = row$I,
    J_per    = row$J,
    true     = true_list[[list_idx]],
    model    = stan_models[[list_idx]],
    priors   = priors_by_model_J[[list_idx]][[as.character(row$J)]],
    model_id = mid
  )
})

stopCluster(cl)
saveRDS(results_list, "output_models56/simulation_results.rds")


## post-processing
ok_idx     <- sapply(results_list, function(r) !is.null(r) && r$status=="ok")
ok_results <- results_list[ok_idx]

results_group   <- bind_rows(lapply(ok_results, `[[`, "group"))
results_indiv   <- bind_rows(lapply(ok_results, `[[`, "indiv"))
results_mean_ik <- bind_rows(lapply(ok_results, `[[`, "mean_ik"))
results_v       <- bind_rows(Filter(Negate(is.null),
                                    lapply(ok_results, `[[`, "v")))

# exclusion report
n_total    <- length(results_list)
n_ok       <- sum(ok_idx)
n_excluded <- sum(sapply(results_list, function(r) !is.null(r) && r$status=="excluded"))
n_error    <- sum(sapply(results_list, function(r) is.null(r) || r$status=="error"))
cat("\n=== OVERALL ===\n")
print(data.frame(Total=n_total, Converged=n_ok, Excluded=n_excluded, Error=n_error,
                 `Excluded (%)`=round(100*n_excluded/n_total,1), check.names=FALSE))

status_df <- bind_rows(lapply(seq_along(results_list), function(i) {
  r <- results_list[[i]]; row <- jobs[i,]
  data.frame(model_id=row$model_id, I=row$I, J=row$J, sim=row$sim,
             status = if (is.null(r)) "error" else r$status,
             reason = if (is.null(r) || r$status=="ok") NA else r$reason)
}))

cat("\n=== EXCLUSIONS BY MODEL AND CELL ===\n")
status_df %>%
  group_by(model_id, I, J) %>%
  summarise(n=n(), excluded=sum(status=="excluded"),
            pct=round(100*excluded/n,1),
            reasons=paste(table(reason[status=="excluded"]), collapse=", "),
            .groups="drop") %>%
  print(n=Inf)

# design matrix verification
cat("\n=== DESIGN MATRIX VERIFICATION (mean_ik) ===\n")
cat("Correlation between true and estimated predicted means:\n")
for (mid in c(5,6)) {
  df <- results_mean_ik %>% filter(model_id == mid)
  if (nrow(df) == 0) next
  cat(sprintf("  Model %d | nu: r=%.3f | alpha: r=%.3f | tau: r=%.3f\n",
              mid,
              cor(df$true_mean_nu,    df$est_mean_nu,    use="complete.obs"),
              cor(df$true_mean_alpha, df$est_mean_alpha, use="complete.obs"),
              cor(df$true_mean_tau,   df$est_mean_tau,   use="complete.obs")))
}

# factor labels
model_labels <- c(
  "5" = "M5: alpha=cond, nu=cond+cov",
  "6" = "M6: alpha=cond, nu=cond+slope"
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
  ggsave(file.path("output_models56", filename), p,
         width=width, height=height, dpi=300, bg="#1e1e1e")
}

# individual recovery plots
for (mid in c(5,6)) {
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

# random slope recovery for model 6
if (nrow(results_v) > 0) {
  results_v <- results_v %>%
    mutate(
      model_label = model_labels[as.character(model_id)],
      I_label = factor(paste0("I = ", I), levels=paste0("I = ", I_vals)),
      J_label = factor(paste0("J = ", J), levels=paste0("J = ", J_vals))
    )

  p_v <- ggplot(results_v, aes(x=true_v, y=est_v)) +
    geom_abline(intercept=0, slope=1, linetype="dashed", colour="white", linewidth=0.5) +
    geom_point(colour="goldenrod", alpha=0.4, size=1.2) +
    facet_grid(I_label ~ J_label) +
    labs(title="M6: Random slope recovery v[i] on nu condition effect",
         x="True v", y="Estimated v") +
    dark_theme
  save_plot(p_v, "m6_recovery_v.png", width=10, height=8)
}

cat("\nDone. Results saved to output_models56/\n")
