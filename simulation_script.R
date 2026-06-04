# check for correct stan file!!
library(rstan)
library(bayesplot)
library(ggplot2)
library(dplyr)
library(tidyr)
library(knitr)
library(parallel)

rstan_options(auto_write = TRUE)
# save compiled Stan model automatically
# compilation is slow --> need to do it only once and is reused on future runs

options(mc.cores = 1)
# use 1 CPU core per chain, overrode later when we do parallel mode
# --> each worker runs on 1 core --> no problem with shared cores

if (!dir.exists("output_new")) dir.create("output_new")
# make output folder or store in that folder if already there

## simulation grid
I_vals <- c(40, 160) # number of participants to test
J_vals <- c(40, 160) # trials per participant per condition to test
N_sim  <- 5          # independent simulation replications per I x J combination
K      <- 2          # how many conditions

## true parameters -- simulation generates data from these
true <- list(
  mu_nu    = 0.4,  # mean/intercept for drift (raw scale)
  b1_nu    = 0.1,  # covariate shift on drift in cond 1
  b2_nu    = 0.3,  # condition 2 vs condition 1 difference in drift
  mu_alpha = 1.5,  # average boundary separation in condition 1
  b_alpha  = 0.2,  # condition 2 boundary is on average 0.2 wider
  mu_tau   = 0.3,  # average non-decision time in condition 1 (seconds)
  b_tau    = 0.25, # condition 2 NDT is 0.25 seconds longer
  sigma_nu    = 0.3,  # between-participant SD of drift (raw scale)
  sigma_alpha = 0.25, # between-participant SD of boundary separation
  sigma_tau   = 0.15, # between-participant SD of NDT
  sigma_v     = 0.2   # SD of random slopes on condition effect for boundary
)

# nu safety threshold: matches eps in Stan model
# nu values with |nu| < NU_EPS are moved to +/-NU_EPS before EZ formulas
# prevents mu_rt diverging to Inf when nu is near zero
# used in both data generation and pilot
NU_EPS <- 0.05


## function takes observed data and computes EZ point estimates of nu, alpha, tau
# used only in pilot to derive priors
ez_point_estimates <- function(C_arr, J_arr, MRT_mat, VRT_mat, I, K) {
# inputs:
# C_arr:   matrix of correct responses, I rows x K columns
# J_arr:   matrix of trials (same dimensions)
# MRT_mat: matrix of mean RTs of correct responses
# VRT_mat: matrix of RT variances
# I: number of participants
# K: number of conditions

  nu_hat    <- matrix(NA, I, K)
  alpha_hat <- matrix(NA, I, K)
  tau_hat   <- matrix(NA, I, K)
  # three empty I x K matrices

  for (i in seq_len(I)) {
    for (k in seq_len(K)) {

      Pc <- C_arr[i, k] / J_arr[i, k]
      # proportion correct for participant i in condition k

      Pc <- max(min(Pc, 1 - 1 / (2 * J_arr[i, k])), 1 / (2 * J_arr[i, k]))
      # keep Pc away from 0 and 1 because log(0) is undefined in EZ

      mrt <- MRT_mat[i, k]
      vrt <- VRT_mat[i, k]

      if (is.na(mrt) || is.na(vrt) || vrt <= 0) next
      # skip cell if mean RT or variance is missing or non-positive

      # EZ inversion formula was derived for positive drift
      # below-chance estimates (Pc < 0.5) are not valid EZ solutions
      # skip them so derive_priors() gets clean pilot estimates
      if (Pc < 0.5) next

      L   <- log(Pc / (1 - Pc))
      # logit of proportion correct = log-odds of being correct

      nu4 <- L * (Pc^2 * L - Pc * L + Pc - 0.5) / vrt
      # nu^4 from the EZ inversion formula

      if (is.na(nu4) || nu4 <= 0) next

      nu_hat[i, k]    <- nu4^(1/4)
      # drift rate estimate -- always positive because we skipped Pc < 0.5 above

      alpha_hat[i, k] <- L / nu_hat[i, k]
      # boundary separation estimate

      tau_hat[i, k]   <- mrt - (alpha_hat[i, k] / (2 * nu_hat[i, k])) * (2 * Pc - 1)
      # NDT estimate
    }
  }
  list(nu = nu_hat, alpha = alpha_hat, tau = tau_hat)
}


## function takes EZ estimates from pilot and turns them into prior hyperparameters
derive_priors <- function(ez) {
  # ez is the list returned by ez_point_estimates

  nu_vals <- as.vector(ez$nu[!is.na(ez$nu) & ez$nu > NU_EPS])
  # only keep valid positive estimates (ez_point_estimates already skipped below-chance)
  # additionally filter < NU_EPS for safety

  alpha_vals <- as.vector(ez$alpha[!is.na(ez$alpha) & ez$alpha > 0])

  tau_vals   <- as.vector(ez$tau[!is.na(ez$tau) & ez$tau > 0.05])

  # nu is on the raw scale (no log transform)
  nu_sd    <- min(max(sd(nu_vals), 0.3), 1.5)
  alpha_sd <- min(max(sd(alpha_vals), 0.2), 1.0)
  tau_sd   <- min(max(sd(tau_vals),   0.1), 0.5)

  list(
    beta_nu    = matrix(c(mean(nu_vals), 0,       0,
                          nu_sd,         nu_sd,   nu_sd),
                        nrow = 2, byrow = TRUE),
    # 2 x 3 matrix
    # row 1 = prior means for [intercept, covariate slope, condition slope]
    # row 2 = prior SDs
    # intercept centred on pilot mean drift; slopes centred at 0 (unknown direction)

    beta_alpha = matrix(c(mean(alpha_vals), 0,
                          alpha_sd * 2,     alpha_sd * 2),
                        nrow = 2, byrow = TRUE),

    beta_tau   = matrix(c(mean(tau_vals), 0,
                          tau_sd * 2,     tau_sd * 2),
                        nrow = 2, byrow = TRUE),

    sigma_nu    = c(nu_sd,       nu_sd / 2),
    sigma_alpha = c(alpha_sd,    alpha_sd / 2),
    sigma_tau   = c(tau_sd,      tau_sd / 2),
    sigma_v     = c(nu_sd * 0.3, nu_sd)
  )
}


## simulates pilot participants, runs EZ on them, returns prior hyperparameters
derive_priors_from_pilot <- function(J_per, true, I_pilot = 20) {
  K          <- 2
  cond_dummy <- c(0, 1)

  x_p      <- rnorm(I_pilot, 0, 1)
  true_v_p <- rnorm(I_pilot, 0, true$sigma_v)

  J_p   <- matrix(J_per, I_pilot, K)
  C_p   <- matrix(NA_integer_, I_pilot, K)
  MRT_p <- matrix(NA_real_,    I_pilot, K)
  VRT_p <- matrix(NA_real_,    I_pilot, K)

  for (i in seq_len(I_pilot)) {
    for (k in seq_len(K)) {

      mean_nu_ik    <- true$mu_nu + true$b1_nu * x_p[i] + true$b2_nu * cond_dummy[k]
      mean_alpha_ik <- true$mu_alpha + (true$b_alpha + true_v_p[i]) * cond_dummy[k]
      mean_tau_ik   <- true$mu_tau + true$b_tau * cond_dummy[k]

      repeat { a <- rnorm(1, mean_alpha_ik, true$sigma_alpha)
               if (a > 0.1 && a < 3.0) { ap <- a; break } }

      # nu drawn from normal on raw scale
      # repeat until |nu| >= NU_EPS so the EZ formulas stay finite
      # matches the nu_safe guard in Stan: Stan will never see |nu| < NU_EPS in the data
      repeat { v <- rnorm(1, mean_nu_ik, true$sigma_nu)
               if (abs(v) >= NU_EPS) { np <- v; break } }

      repeat { t <- rnorm(1, mean_tau_ik, true$sigma_tau)
               if (t > 0.05) { tp <- t; break } }

      e      <- exp(-ap * np)
      pi_c   <- 1 / (1 + e)
      mu_rt  <- tp + (ap / (2 * np)) * ((1 - e) / (1 + e))
      sig2rt <- (ap / (2 * np^3)) * ((1 - 2*ap*np*e - e^2) / (1 + e)^2)
      # sig2rt is always positive when |np| >= NU_EPS

      C_ik <- rbinom(1, J_per, pi_c)
      C_ik <- max(1L, min(C_ik, J_per - 1L))

      MRT_p[i, k] <- rnorm(1, mu_rt, sqrt(sig2rt / C_ik))
      VRT_p[i, k] <- max(rnorm(1, sig2rt, sqrt(2 * sig2rt^2 / max(C_ik-1, 1))), 1e-6)
      C_p[i, k]   <- C_ik
    }
  }

  ez <- ez_point_estimates(C_p, J_p, MRT_p, VRT_p, I_pilot, K)
  derive_priors(ez)
}


## simulate one dataset from true parameters, fit Stan, return results
run_one <- function(I, J_per, true, model, priors,
                    rhat_threshold     = 1.05,
                    hopeless_threshold = 10,
                    iter1 = 1000, warmup1 = 300,
                    iter2 = 3000, warmup2 = 1000) {
  # I:     number of participants
  # J_per: trials per participant per condition
  # true:  true parameter list
  # model: compiled Stan model
  # priors: prior hyperparameters for this J value

  K          <- 2
  cond_dummy <- c(0, 1)
  prs <- priors

  x      <- rnorm(I, 0, 1)
  true_v <- rnorm(I, 0, true$sigma_v)

  # design matrices: I x K x predictors
  X_nu    <- array(NA_real_, dim = c(I, K, 3))
  X_alpha <- array(NA_real_, dim = c(I, K, 2))
  X_tau   <- array(NA_real_, dim = c(I, K, 2))

  for (i in seq_len(I)) {
    for (k in seq_len(K)) {
      X_nu[i, k, ]    <- c(1, x[i], cond_dummy[k])
      X_alpha[i, k, ] <- c(1, cond_dummy[k])
      X_tau[i, k, ]   <- c(1, cond_dummy[k])
    }
  }

  # true individual parameters
  alpha_true <- matrix(NA_real_, I, K)
  nu_true    <- matrix(NA_real_, I, K)
  tau_true   <- matrix(NA_real_, I, K)

  for (i in seq_len(I)) {
    for (k in seq_len(K)) {

      mean_nu_ik    <- true$mu_nu + true$b1_nu * x[i] + true$b2_nu * cond_dummy[k]
      mean_alpha_ik <- true$mu_alpha + (true$b_alpha + true_v[i]) * cond_dummy[k]
      mean_tau_ik   <- true$mu_tau + true$b_tau * cond_dummy[k]

      repeat { a <- rnorm(1, mean_alpha_ik, true$sigma_alpha)
               if (a > 0.1 && a < 3.0) { alpha_true[i, k] <- a; break } }

      # nu: repeat until |nu| >= NU_EPS
      # ensures data generation and Stan model treat the same region as off-limits
      # without this, nu_true near 0 produces implausible MRT values (mu_rt --> Inf)
      repeat { v <- rnorm(1, mean_nu_ik, true$sigma_nu)
               if (abs(v) >= NU_EPS) { nu_true[i, k] <- v; break } }

      repeat { t <- rnorm(1, mean_tau_ik, true$sigma_tau)
               if (t > 0.05) { tau_true[i, k] <- t; break } }
    }
  }

  # generate observed data from true parameters via EZ forward equations
  J_mat   <- matrix(J_per, I, K)
  C_mat   <- matrix(NA_integer_, I, K)
  MRT_mat <- matrix(NA_real_,    I, K)
  VRT_mat <- matrix(NA_real_,    I, K)

  for (i in seq_len(I)) {
    for (k in seq_len(K)) {
      a <- alpha_true[i, k]
      v <- nu_true[i, k]    # |v| >= NU_EPS guaranteed by repeat-loop above
      tau_ik <- tau_true[i, k]

      e      <- exp(-a * v)
      pi_c   <- 1 / (1 + e)
      mu_rt  <- tau_ik + (a / (2 * v)) * ((1 - e) / (1 + e))
      sig2rt <- (a / (2 * v^3)) * ((1 - 2*a*v*e - e^2) / (1 + e)^2)
      # both mu_rt and sig2rt are finite and sig2rt > 0 because |v| >= NU_EPS

      C_ik <- rbinom(1, J_per, pi_c)
      C_ik <- max(1L, min(C_ik, J_per - 1L))

      MRT_mat[i, k] <- rnorm(1, mu_rt, sqrt(sig2rt / C_ik))
      VRT_mat[i, k] <- max(rnorm(1, sig2rt, sqrt(2 * sig2rt^2 / max(C_ik-1, 1))), 1e-6)
      C_mat[i, k]   <- C_ik
    }
  }

  stan_data <- list(
    I = I, K = K,
    P_nu = 3L, P_alpha = 2L, P_tau = 2L,
    J = J_mat, C = C_mat, MRT = MRT_mat, VRT = VRT_mat,
    X_nu = X_nu, X_alpha = X_alpha, X_tau = X_tau,
    prior_beta_nu     = prs$beta_nu,
    prior_beta_alpha  = prs$beta_alpha,
    prior_beta_tau    = prs$beta_tau,
    prior_sigma_nu    = prs$sigma_nu,
    prior_sigma_alpha = prs$sigma_alpha,
    prior_sigma_tau   = prs$sigma_tau,
    prior_sigma_v     = prs$sigma_v
  )

  group_pars <- c(
    paste0("beta_nu[",    seq_len(3), "]"),
    paste0("beta_alpha[", seq_len(2), "]"),
    paste0("beta_tau[",   seq_len(2), "]"),
    "sigma_nu", "sigma_alpha", "sigma_tau", "sigma_v"
  )

  hopeless <- function(fit) {
    rh <- tryCatch(summary(fit)$summary[group_pars, "Rhat"], error = function(e) NA_real_)
    any(is.na(rh)) || any(is.nan(rh)) || any(is.infinite(rh)) ||
      any(rh > hopeless_threshold, na.rm = TRUE)
  }

  fit <- tryCatch(
    sampling(model, data = stan_data,
             iter = iter1, warmup = warmup1, chains = 4, cores = 1, refresh = 0,
             control = list(adapt_delta = 0.85, max_treedepth = 10)),
    error = function(e) NULL
  )

  if (is.null(fit))  return(list(status = "error"))
  if (hopeless(fit)) return(list(status = "excluded", reason = "degenerate_fit1"))

  rhat_max  <- max(summary(fit)$summary[group_pars, "Rhat"], na.rm = TRUE)
  converged <- "first_fit"

  if (rhat_max > rhat_threshold) {
    message(sprintf("Refit triggered (I=%d, J=%d) -- Rhat = %.2f", I, J_per, rhat_max))

    fit <- tryCatch(
      sampling(model, data = stan_data,
               iter = iter2, warmup = warmup2, chains = 4, cores = 1, refresh = 0,
               control = list(adapt_delta = 0.95, max_treedepth = 12)),
      error = function(e) NULL
    )

    if (is.null(fit))  return(list(status = "error"))
    if (hopeless(fit)) return(list(status = "excluded", reason = "degenerate_refit"))
    rhat_max  <- max(summary(fit)$summary[group_pars, "Rhat"], na.rm = TRUE)
    converged <- "refit"
    if (rhat_max > rhat_threshold)
      return(list(status = "excluded", reason = "non_converged"))
  }

  s <- summary(fit)$summary

  group_results <- data.frame(
    I = I, J = J_per,
    param = c("mu_nu", "b1_nu", "b2_nu",
              "mu_alpha", "b_alpha",
              "mu_tau", "b_tau",
              "sigma_nu", "sigma_alpha", "sigma_tau", "sigma_v"),
    true_val = c(true$mu_nu,    true$b1_nu,    true$b2_nu,
                 true$mu_alpha, true$b_alpha,
                 true$mu_tau,   true$b_tau,
                 true$sigma_nu, true$sigma_alpha, true$sigma_tau, true$sigma_v),
    estimate = c(
      s["beta_nu[1]",    "mean"], s["beta_nu[2]",    "mean"], s["beta_nu[3]",  "mean"],
      s["beta_alpha[1]", "mean"], s["beta_alpha[2]", "mean"],
      s["beta_tau[1]",   "mean"], s["beta_tau[2]",   "mean"],
      s["sigma_nu",      "mean"], s["sigma_alpha",   "mean"],
      s["sigma_tau",     "mean"], s["sigma_v",       "mean"]
    ),
    prior_mean = c(
      prs$beta_nu[1,1],    prs$beta_nu[1,2],    prs$beta_nu[1,3],
      prs$beta_alpha[1,1], prs$beta_alpha[1,2],
      prs$beta_tau[1,1],   prs$beta_tau[1,2],
      prs$sigma_nu[1], prs$sigma_alpha[1], prs$sigma_tau[1], prs$sigma_v[1]
    ),
    prior_sd = c(
      prs$beta_nu[2,1],    prs$beta_nu[2,2],    prs$beta_nu[2,3],
      prs$beta_alpha[2,1], prs$beta_alpha[2,2],
      prs$beta_tau[2,1],   prs$beta_tau[2,2],
      prs$sigma_nu[2], prs$sigma_alpha[2], prs$sigma_tau[2], prs$sigma_v[2]
    ),
    converged = converged,
    rhat_max  = rhat_max,
    stringsAsFactors = FALSE
  )

  v_results <- data.frame(
    I = I, J = J_per,
    participant = seq_len(I),
    true_v = true_v,
    est_v  = s[paste0("v[", seq_len(I), "]"), "mean"]
  )

  alpha_est <- matrix(NA_real_, I, K)
  nu_est    <- matrix(NA_real_, I, K)
  tau_est   <- matrix(NA_real_, I, K)

  for (i in seq_len(I)) {
    for (k in seq_len(K)) {
      alpha_est[i, k] <- s[paste0("alpha[", i, ",", k, "]"), "mean"]
      nu_est[i, k]    <- s[paste0("nu[",    i, ",", k, "]"), "mean"]
      tau_est[i, k]   <- s[paste0("tau[",   i, ",", k, "]"), "mean"]
    }
  }

  indiv_results <- data.frame(
    I = I, J = J_per,
    condition  = rep(seq_len(K), each = I),
    true_alpha = as.vector(alpha_true), est_alpha = as.vector(alpha_est),
    true_nu    = as.vector(nu_true),    est_nu    = as.vector(nu_est),
    true_tau   = as.vector(tau_true),   est_tau   = as.vector(tau_est)
  )

  list(status = "ok", group = group_results, v = v_results, indiv = indiv_results,
       priors_used = prs, rhat_max = rhat_max, converged = converged)
}


## compile Stan model
cache_file <- paste0(tools::file_path_sans_ext(
  "model.stan"), ".rds")
if (file.exists(cache_file)) file.remove(cache_file)
# force fresh recompile so we never run on a stale binary

## derive priors once per J value
set.seed(42)
priors_by_J <- setNames(
  lapply(J_vals, function(j) derive_priors_from_pilot(j, true, I_pilot = 20)),
  as.character(J_vals)
)

## parallel simulation
jobs      <- expand.grid(I = I_vals, J = J_vals, sim = seq_len(N_sim))
jobs_list <- split(jobs, seq_len(nrow(jobs)))

cl <- makeCluster(detectCores() - 1)
clusterSetRNGStream(cl, 42)

clusterExport(cl, varlist = c("run_one", "derive_priors_from_pilot",
                              "ez_point_estimates", "derive_priors",
                              "true", "priors_by_J", "NU_EPS"))
# NU_EPS exported so workers use the same threshold as the main session

clusterEvalQ(cl, {
  library(rstan)
  library(dplyr)
  rstan_options(auto_write = FALSE)
  options(mc.cores = 1)
  ez_cond <<- stan_model("model.stan")
})

results_list <- parLapply(cl, jobs_list, function(row) {
  run_one(I      = row$I,
          J_per  = row$J,
          true   = true,
          model  = ez_cond,
          priors = priors_by_J[[as.character(row$J)]])
})

stopCluster(cl)

saveRDS(results_list, file = "output_new/simulation_results.rds")


## plotting

ok_idx        <- sapply(results_list, function(r) !is.null(r) && r$status == "ok")
ok_results    <- results_list[ok_idx]
results_group <- bind_rows(lapply(ok_results, `[[`, "group"))
results_v     <- bind_rows(lapply(ok_results, `[[`, "v"))
results_indiv <- bind_rows(lapply(ok_results, `[[`, "indiv"))

# exclusion report
n_total    <- length(results_list)
n_ok       <- sum(ok_idx)
n_excluded <- sum(sapply(results_list, function(r) !is.null(r) && r$status == "excluded"))
n_error    <- sum(sapply(results_list, function(r) is.null(r)  || r$status == "error"))
print(data.frame(
  Total = n_total, Converged = n_ok, Excluded = n_excluded, Error = n_error,
  `Excluded (%)` = round(100 * n_excluded / n_total, 1), check.names = FALSE
))

# exclusion breakdown by cell
status_df <- bind_rows(lapply(seq_along(results_list), function(i) {
  r <- results_list[[i]]; row <- jobs[i, ]
  data.frame(I = row$I, J = row$J, sim = row$sim,
             status = if (is.null(r)) "error" else r$status,
             reason = if (is.null(r) || r$status == "ok") NA else r$reason)
}))

status_df %>%
  group_by(I, J) %>%
  summarise(n        = n(),
            excluded = sum(status == "excluded"),
            pct      = round(100 * excluded / n, 1),
            reasons  = paste(table(reason[status == "excluded"]), collapse = ", "),
            .groups  = "drop") %>%
  print(n = Inf)

# factor labels for plots
add_labels <- function(df) {
  df %>% mutate(
    I_label = factor(paste0("I = ", I), levels = paste0("I = ", I_vals)),
    J_label = factor(paste0("J = ", J), levels = paste0("J = ", J_vals))
  )
}
results_group <- add_labels(results_group)
results_v     <- add_labels(results_v)
results_indiv <- results_indiv %>%
  add_labels() %>%
  mutate(condition = factor(condition, labels = c("Condition 1 (k=1)", "Condition 2 (k=2)")))

# plot theme
dark_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.background   = element_rect(fill = "#1e1e1e", colour = NA),
    panel.background  = element_rect(fill = "#1e1e1e", colour = NA),
    panel.grid.major  = element_line(colour = "#333333"),
    panel.grid.minor  = element_blank(),
    text              = element_text(colour = "white"),
    axis.text         = element_text(colour = "white"),
    strip.text        = element_text(colour = "white"),
    legend.background = element_rect(fill = "#1e1e1e"),
    legend.text       = element_text(colour = "white")
  )

save_plot <- function(plot, filename, width = 8, height = 6) {
  ggsave(file.path("output_new", filename), plot, width = width, height = height,
         dpi = 300, bg = "#1e1e1e")
}

# group-level bias plots
plot_group_bias <- function(param_name, label) {
  df <- results_group %>%
    filter(param == param_name) %>%
    mutate(bias = estimate - true_val)
  ggplot(df, aes(x = J_label, y = bias)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "white", linewidth = 0.5) +
    geom_violin(fill = "goldenrod", colour = NA, alpha = 0.35, trim = TRUE, scale = "width") +
    geom_jitter(colour = "goldenrod", alpha = 0.7, size = 1.4, width = 0.08) +
    stat_summary(fun = mean, geom = "crossbar", colour = "white",
                 linewidth = 0.6, width = 0.4, fatten = 1) +
    facet_grid(I_label ~ .) +
    labs(title = label, x = "Trials per condition", y = "Bias (estimate - true)") +
    dark_theme
}

p_mu_nu    <- plot_group_bias("mu_nu",    expression(mu^(nu)~"(intercept, raw scale)"))
p_mu_alpha <- plot_group_bias("mu_alpha", expression(mu^(alpha)))
p_mu_tau   <- plot_group_bias("mu_tau",   expression(mu^(tau)))
p_b1_nu    <- plot_group_bias("b1_nu",    expression(beta[1]^(nu)~"(covariate effect on drift)"))
p_b2_nu    <- plot_group_bias("b2_nu",    expression(beta[2]^(nu)~"(condition effect on drift)"))
p_b_alpha  <- plot_group_bias("b_alpha",  expression(beta^(alpha)~"(avg condition effect)"))
p_b_tau    <- plot_group_bias("b_tau",    expression(beta^(tau)~"(condition)"))
p_sigma_v  <- plot_group_bias("sigma_v",  expression(sigma[v]~"(random slope SD)"))

for (p in list(p_mu_nu, p_mu_alpha, p_mu_tau,
               p_b1_nu, p_b2_nu, p_b_alpha, p_b_tau, p_sigma_v)) print(p)

save_plot(p_mu_nu,    "bias_mu_nu.png");    save_plot(p_mu_alpha, "bias_mu_alpha.png")
save_plot(p_mu_tau,   "bias_mu_tau.png");   save_plot(p_b1_nu,    "bias_b1_nu.png")
save_plot(p_b2_nu,    "bias_b2_nu.png");    save_plot(p_b_alpha,  "bias_b_alpha.png")
save_plot(p_b_tau,    "bias_b_tau.png");    save_plot(p_sigma_v,  "bias_sigma_v.png")

# individual-level recovery plots
plot_indiv_recovery <- function(true_col, est_col, label) {
  ggplot(results_indiv, aes_string(x = true_col, y = est_col, colour = "condition")) +
    geom_abline(intercept = 0, slope = 1,
                linetype = "dashed", colour = "white", linewidth = 0.5) +
    geom_point(alpha = 0.35, size = 0.9) +
    scale_colour_manual(values = c("Condition 1 (k=1)" = "goldenrod",
                                   "Condition 2 (k=2)" = "steelblue")) +
    facet_grid(I_label ~ J_label) +
    labs(title = label, x = paste("True", label), y = paste("Estimated", label),
         colour = "Condition") +
    dark_theme + theme(legend.position = "bottom")
}

p_indiv_alpha <- plot_indiv_recovery("true_alpha", "est_alpha", expression(alpha[ik]))
p_indiv_nu    <- plot_indiv_recovery("true_nu",    "est_nu",    expression(nu[ik]))
p_indiv_tau   <- plot_indiv_recovery("true_tau",   "est_tau",   expression(tau[ik]))

for (p in list(p_indiv_alpha, p_indiv_nu, p_indiv_tau)) print(p)
save_plot(p_indiv_alpha, "recovery_alpha.png", width = 10, height = 8)
save_plot(p_indiv_nu,    "recovery_nu.png",    width = 10, height = 8)
save_plot(p_indiv_tau,   "recovery_tau.png",   width = 10, height = 8)

# random slope recovery plot
p_v_recovery <- ggplot(results_v, aes(x = true_v, y = est_v)) +
  geom_abline(intercept = 0, slope = 1,
              linetype = "dashed", colour = "white", linewidth = 0.5) +
  geom_point(colour = "goldenrod", alpha = 0.4, size = 1.2) +
  facet_grid(I_label ~ J_label) +
  labs(title = expression("Recovery of by-participant random slopes " * v[i]),
       x = "True v", y = "Estimated v") +
  dark_theme
print(p_v_recovery)
save_plot(p_v_recovery, "recovery_v.png", width = 10, height = 8)

# prior summary table
prior_summary <- results_group %>%
  group_by(param, I, J) %>%
  summarise(`Prior mean (avg)` = round(mean(prior_mean, na.rm = TRUE), 3),
            `Prior SD (avg)`   = round(mean(prior_sd,   na.rm = TRUE), 3),
            .groups = "drop") %>%
  rename(Parameter = param)
write.csv(prior_summary, file.path("output_new", "prior_summary.csv"), row.names = FALSE)
kable(prior_summary,
      caption = "Average derived prior hyperparameters across simulations,
                 by parameter, number of participants, and trials per condition.")
