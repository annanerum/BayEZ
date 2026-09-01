options(error = NULL)

library(cmdstanr)
library(posterior)
library(dplyr)

TEST_RUN <- FALSE
n_reps   <- 10
DATA_DIR <- "data_full7"
FITS_DIR <- "fits/M2 sv hier"

m2sv_file <- "stan/M2_sv_hier.stan"

prior_w       <- c(0.5, 0.15)
prior_beta_sv  <- c(log(0.6), 0.5)   
prior_sigma_sv <- c(2, 4)            

m2sv_extra <- list(prior_w = prior_w, prior_beta_sv = prior_beta_sv, prior_sigma_sv = prior_sigma_sv)

mod_m2sv <- cmdstan_model(m2sv_file, cpp_options = list(stan_threads = TRUE), force_recompile = TRUE)

CHAINS          <- 2
PARALLEL_CHAINS <- 2
N_THREADS         <- max(1, parallel::detectCores() - 1)
THREADS_PER_CHAIN <- max(1, floor(N_THREADS / PARALLEL_CHAINS))

message(sprintf("detectCores()=%d --> N_THREADS=%d, THREADS_PER_CHAIN=%d",
                parallel::detectCores(), N_THREADS, THREADS_PER_CHAIN))

make_init_fun <- function(I) {
  function() {
    list(
      beta_nu     = c(0, 0, 0, 0) + rnorm(4, 0, 0.05),
      beta_alpha1 = log(1.8) + rnorm(1, 0, 0.05),
      beta_alpha2 = 0 + rnorm(1, 0, 0.02),
      beta_alpha3 = 0 + rnorm(1, 0, 0.02),
      beta_tau    = 0.3 + rnorm(1, 0, 0.05),
      sigma_nu_c  = 0.5 + runif(1, 0, 0.1),
      sigma_nu_r  = 0.3 + runif(1, 0, 0.05),
      sigma_nu_cr = 0.3 + runif(1, 0, 0.05),
      sigma_alpha_intercept = 0.3 + runif(1, 0, 0.05),
      sigma_tau_intercept   = 0.3 + runif(1, 0, 0.05),
      z_nu_c  = rnorm(I, 0, 0.1),
      z_nu_r  = rnorm(I, 0, 0.1),
      z_nu_cr = rnorm(I, 0, 0.1),
      z_alpha_intercept = rnorm(I, 0, 0.1),
      z_tau_intercept   = rnorm(I, 0, 0.1),
      w        = 0.5 + runif(1, -0.02, 0.02),
      beta_sv  = log(0.6) + rnorm(1, 0, 0.05),
      sigma_sv = 0.3 + runif(1, 0, 0.1),
      z_sv     = rnorm(I, 0, 0.1)
    )
  }
}

build_stan_data <- function(rep_id) {
  rep_dir <- file.path(DATA_DIR, sprintf("rep%02d", rep_id))
  sim_data   <- read.csv(file.path(rep_dir, "sim_data_full_ddm.csv"))
  min_rt_obs <- readRDS(file.path(rep_dir, "min_rt_observed.rds"))
  
  subjects_used <- sort(unique(sim_data$subj))
  min_rt_used   <- as.numeric(min_rt_obs[as.character(subjects_used)])
  sim_data$pid  <- as.integer(factor(sim_data$subj, levels = subjects_used))
  
  stan_data_base <- list(
    N = nrow(sim_data), I = length(subjects_used), pid = sim_data$pid,
    rt = sim_data$rt, acc = sim_data$acc,
    condition = sim_data$condition, resp_type = sim_data$resp_type,
    min_rt = min_rt_used, grainsize = 1,
    prior_beta_nu = matrix(c(0, 0, 0, 0, 3, 3, 3, 3), nrow = 2, byrow = TRUE),
    prior_beta_alpha1 = c(log(2), 1), prior_beta_alpha2 = c(0, 0.5),
    prior_beta_alpha3 = c(0, 0.5), prior_beta_tau = c(0, 1),
    prior_sigma_nu_c = c(2, 1), prior_sigma_nu_r = c(2, 1), prior_sigma_nu_cr = c(2, 1),
    prior_sigma_alpha_intercept = c(2, 1), prior_sigma_tau_intercept = c(2, 1)
  )
  c(stan_data_base, m2sv_extra)
}

fit_m2sv_rep <- function(rep_id) {
  cat(sprintf("\n Fitting M2 (hierarchical sv), rep %02d / %02d \n", rep_id, n_reps))
  
  stan_data <- build_stan_data(rep_id)
  init_fun  <- make_init_fun(stan_data$I)
  
  out_dir <- file.path(FITS_DIR, sprintf("rep%02d_fullscale", rep_id))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  basename_r <- sprintf("M2sv_rep%02d", rep_id)
  
  fit <- mod_m2sv$sample(
    data = stan_data, chains = CHAINS, parallel_chains = PARALLEL_CHAINS,
    threads_per_chain = THREADS_PER_CHAIN, iter_warmup = 1000, iter_sampling = 1000,
    adapt_delta = 0.9, max_treedepth = 8, refresh = 50, init = init_fun,
    output_dir = out_dir, output_basename = paste0(basename_r, "_stufe1"),
    save_cmdstan_config = TRUE
  )
  
  n_div <- sum(fit$diagnostic_summary(quiet = TRUE)$num_divergent)
  if (n_div > 0) {
    cat(sprintf("  %d Divergences --> Stufe 2 \n", n_div))
    fit <- mod_m2sv$sample(
      data = stan_data, chains = CHAINS, parallel_chains = PARALLEL_CHAINS,
      threads_per_chain = THREADS_PER_CHAIN, iter_warmup = 1000, iter_sampling = 1000,
      adapt_delta = 0.95, max_treedepth = 10, refresh = 50, seed = 2026, init = init_fun,
      output_dir = out_dir, output_basename = paste0(basename_r, "_stufe2"),
      save_cmdstan_config = TRUE
    )
  }
  
  list(rep = rep_id, fit = fit, status = "full_scale_ok", out_dir = out_dir, stan_data = stan_data)
}

dir.create("fits", showWarnings = FALSE)
dir.create(FITS_DIR, showWarnings = FALSE, recursive = TRUE)

results <- list()
for (r in seq_len(n_reps)) {
  res <- fit_m2sv_rep(r)
  results[[r]] <- res
  if (!is.null(res$fit)) {
    out_dir <- res$out_dir
    res$fit$save_object(file.path(out_dir, "fit.rds"))
    write.csv(res$fit$summary(), file.path(out_dir, "summary.csv"), row.names = FALSE)
    cat("  gespeichert:", out_dir, "(fit.rds + summary.csv)\n")
  }
  saveRDS(results, file.path(FITS_DIR, "all_results.rds"))
}

cat("\n\nZusammenfassung:\n")
for (r in seq_len(n_reps)) cat(sprintf("  rep%02d: %s\n", r, results[[r]]$status))

