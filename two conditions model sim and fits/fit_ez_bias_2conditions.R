library(cmdstanr)
library(dplyr)
library(posterior)

PROJECT_DIR <- "/Users/annegiacobello/Desktop/EZ Paper "   
DATA_DIR    <- file.path(PROJECT_DIR, "data_sim_fullddm_2condtions")
FIT_DIR     <- file.path(PROJECT_DIR, "fits_ez_bias_2conditions")          
STAN_FILE   <- file.path(PROJECT_DIR, "ez_bias_ddm_2conditions.stan")      

# safety check: stop if the data are not from the current simulation (with person-specific drift intercept)
stopifnot("nu_intercept" %in% names(readRDS(file.path(DATA_DIR, "rep01", "true_subject_parameters.rds"))))
n_reps    <- 10

# TEST_MODE = TRUE --> only rep01, 2 chains, few iterations
TEST_MODE <- FALSE

dir.create(FIT_DIR, showWarnings = FALSE, recursive = TRUE)

mod <- cmdstan_model(STAN_FILE)


#### summary statistics per person and cell
# cell order must match the Stan model: k = 1 --> cond +1, k = 2 --> cond -1
COND_LEVELS <- c(1, -1)

to_matrices <- function(cells) {
  subj_ids <- sort(unique(cells$subj))
  I <- length(subj_ids)
  K <- length(COND_LEVELS)
  to_mat <- function(x) matrix(x, nrow = I, ncol = K, byrow = TRUE)  # rows = persons, cols = cells
  list(I = I, K = K, subj_ids = subj_ids, to_mat = to_mat)
}

# make EZ statistics (correct RTs only) 
make_ez_classic_stats <- function(dat) {
  cells <- dat %>%
    mutate(k = match(condition, COND_LEVELS)) %>%
    group_by(subj, k) %>%
    summarise(J   = n(),
              C   = sum(acc == 1),
              MRT = mean(rt[acc == 1]),
              VRT = var(rt[acc == 1]),
              .groups = "drop") %>%
    arrange(subj, k)
  if (any(cells$C < 2)) stop("some cells have fewer than 2 correct responses")
  m <- to_matrices(cells)
  list(I = m$I, K = m$K, J = m$to_mat(cells$J), C = m$to_mat(cells$C),
       MRT = m$to_mat(cells$MRT), VRT = m$to_mat(cells$VRT))
}

# make EZ+bias
# U = number of upper responses (upper = correct, acc = +1)
# MRT_U = mean RT of upper responses
# MRT_L = mean RT of lower responses (errors)
# VRT   = variance of ALL RTs in the cell (both boundaries)
# cells without responses at one boundary get mean RT 0 there
# Stan skips that term (U = 0 or L = 0)
make_stan_data <- function(dat) {
  cells <- dat %>%
    mutate(k = match(condition, COND_LEVELS)) %>%
    group_by(subj, k) %>%
    summarise(J     = n(),
              U     = sum(acc == 1),
              MRT_U = if (any(acc == 1))  mean(rt[acc == 1])  else 0,
              MRT_L = if (any(acc == -1)) mean(rt[acc == -1]) else 0,
              VRT   = var(rt),
              .groups = "drop") %>%
    arrange(subj, k)

  n_no_err <- sum(cells$U == cells$J)
  if (n_no_err > 0)
    cat(n_no_err, "cells without any error --> error-RT term skipped there\n")

  m <- to_matrices(cells)
  list(
    I = m$I, K = m$K, cond = COND_LEVELS,
    J     = m$to_mat(cells$J),
    U     = m$to_mat(cells$U),
    MRT_U = m$to_mat(cells$MRT_U),
    MRT_L = m$to_mat(cells$MRT_L),
    VRT   = m$to_mat(cells$VRT),
    subj_ids = m$subj_ids          # for matching persons later
  )
}


#### priors (as in the EZ fit) + bias priors
# EZ point estimates (closed-form equations) for a pilot of 20 persons
# prior means from the average estimates, prior SDs from their spread
# pilot persons are taken from a different rep than the one being fitted

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

derive_priors <- function(ez) {
  nu_vals    <- as.vector(ez$nu[!is.na(ez$nu)])
  alpha_vals <- as.vector(ez$alpha[!is.na(ez$alpha) & ez$alpha > 0])
  tau_vals   <- as.vector(ez$tau[!is.na(ez$tau) & ez$tau > 0.05])

  log_alpha_vals <- log(alpha_vals)

  nu_sd        <- min(max(sd(nu_vals),        0.3), 1.5)
  log_alpha_sd <- min(max(sd(log_alpha_vals), 0.1), 0.7)
  tau_sd       <- min(max(sd(tau_vals),       0.1), 0.5)

  log_alpha_scale <- min(max(sd(log_alpha_vals), 0.05), 0.7)
  tau_scale       <- min(max(sd(tau_vals),       0.03), 0.5)

  list(
    prior_beta_nu_intercept     = c(mean(nu_vals),        nu_sd),
    prior_beta_nu_c             = c(0,                    nu_sd),
    prior_beta_alpha            = c(mean(log_alpha_vals), 2 * log_alpha_sd),
    prior_beta_tau              = c(mean(tau_vals),       2 * tau_sd),
    # bias prior: EZ cannot estimate w, so it can't come from the pilot
    # like Nina but translated it to the logit scale
    # mu_beta    ~ N(0.5, 0.125) --> beta_w  ~ N(0, 0.5)  (logit --> 95% of w between about 0.27 and 0.73)
    # sigma_beta ~ N(0.1, 0.05) T(0.02,0.15) --> about 0.4 on the logit scale
    # used as scale in the same N(0.3 * s, s) form as the other SD priors
    prior_beta_w                = c(0,                    0.5),
    prior_sigma_nu_intercept    = c(0.3 * nu_sd,          nu_sd),
    prior_sigma_nu_c            = c(0.3 * nu_sd,          nu_sd),
    prior_sigma_alpha_intercept = c(0.3 * log_alpha_scale, log_alpha_scale),
    prior_sigma_tau_intercept   = c(0.3 * tau_scale,       tau_scale),
    prior_sigma_w               = c(0.3 * 0.4,            0.4)          # NEW
  )
}

derive_priors_from_pilot <- function(r, I_pilot = 20) {
  pilot_rep <- if (r == n_reps) 1 else r + 1
  dat_p     <- read.csv(file.path(DATA_DIR, sprintf("rep%02d", pilot_rep), "sim_data_full_ddm.csv"))

  set.seed(42 + r)
  pilot_subj <- sample(unique(dat_p$subj), I_pilot)
  sd_p <- make_ez_classic_stats(dat_p[dat_p$subj %in% pilot_subj, ])   

  ez <- ez_point_estimates(sd_p$C, sd_p$J, sd_p$MRT, sd_p$VRT, sd_p$I, sd_p$K)
  pr <- derive_priors(ez)
  attr(pr, "pilot_rep") <- pilot_rep
  pr
}


#### fit all reps

reps_to_fit <- if (TEST_MODE) 1 else seq_len(n_reps)

# refit only selected reps (e.g. a non-converged one), leave NULL to fit all
REFIT_REPS <- NULL          
if (!is.null(REFIT_REPS)) reps_to_fit <- REFIT_REPS

# starting values near the prior means
make_inits <- function(priors, I, n_chains) {
  lapply(seq_len(n_chains), function(ch) list(
    beta_nu_intercept     = priors$prior_beta_nu_intercept[1] + rnorm(1, 0, 0.05),
    beta_nu_c             = 0.3 + rnorm(1, 0, 0.05),
    beta_alpha            = priors$prior_beta_alpha[1] + rnorm(1, 0, 0.05),
    beta_tau              = min(max(priors$prior_beta_tau[1], 0.10), 0.40) + rnorm(1, 0, 0.01),
    beta_w                = rnorm(1, 0, 0.05),                       # w near 0.5
    sigma_nu_intercept    = 0.2,
    sigma_nu_c            = 0.2,
    sigma_alpha_intercept = 0.1,
    sigma_tau_intercept   = 0.03,
    sigma_w               = 0.2,                                    
    z_nu_intercept        = rnorm(I, 0, 0.1),
    z_nu_c                = rnorm(I, 0, 0.1),
    z_alpha_intercept     = rnorm(I, 0, 0.1),
    z_tau_intercept       = rnorm(I, 0, 0.1),
    z_w                   = rnorm(I, 0, 0.1)                         
  ))
}

sampler_args <- if (TEST_MODE) {
  list(chains = 2, parallel_chains = 2, iter_warmup = 300, iter_sampling = 300)
} else {
  list(chains = 4, parallel_chains = 4, iter_warmup = 1000, iter_sampling = 1000)
}

pop_pars    <- c("beta_nu_intercept", "beta_nu_c", "beta_alpha", "beta_tau", "beta_w",           
                 "sigma_nu_intercept", "sigma_nu_c", "sigma_alpha_intercept", "sigma_tau_intercept",
                 "sigma_w", "w_pop_median")                                                      
person_pars <- c("nu_intercept", "nu_c", "alpha_intercept", "tau_intercept", "w_person")        

for (r in reps_to_fit) {
  rep_name <- sprintf("rep%02d", r)
  cat("\n ", rep_name, " \n")

  dat        <- read.csv(file.path(DATA_DIR, rep_name, "sim_data_full_ddm.csv"))
  stan_input <- make_stan_data(dat)
  priors     <- derive_priors_from_pilot(r)
  cat("priors derived from 20 persons of rep", attr(priors, "pilot_rep"), "\n")
  print(round(do.call(rbind, priors), 3))                           # rows = priors, cols = [mean, SD]
  stan_data  <- c(stan_input[setdiff(names(stan_input), "subj_ids")], priors)

  set.seed(2026 + r)
  fit <- do.call(mod$sample, c(list(
    data          = stan_data,
    init          = make_inits(priors, stan_input$I, sampler_args$chains),
    seed          = 2026 + r,
    adapt_delta   = 0.95,
    max_treedepth = 12,
    refresh       = 500
  ), sampler_args))

  out_dir <- file.path(FIT_DIR, rep_name)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  diag <- fit$diagnostic_summary(quiet = TRUE)
  summ <- fit$summary(variables = c(pop_pars, person_pars))
  diag_row <- data.frame(
    rep           = r,
    divergences   = sum(diag$num_divergent),
    treedepth_hit = sum(diag$num_max_treedepth),
    max_rhat      = max(summ$rhat, na.rm = TRUE),
    min_ess_bulk  = min(summ$ess_bulk, na.rm = TRUE),
    min_ess_tail  = min(summ$ess_tail, na.rm = TRUE)
  )
  print(diag_row)

  fit$save_object(file.path(out_dir, "fit_ez_bias.rds"))
  write.csv(summ,     file.path(out_dir, "summary_ez_bias.csv"),     row.names = FALSE)
  write.csv(diag_row, file.path(out_dir, "diagnostics_ez_bias.csv"), row.names = FALSE)
  saveRDS(stan_input$subj_ids, file.path(out_dir, "subj_ids.rds"))
  saveRDS(priors,              file.path(out_dir, "priors_ez_bias.rds"))
}


#### parameter recovery

all_reps <- if (TEST_MODE) 1 else seq_len(n_reps)
diag_all <- bind_rows(lapply(all_reps, function(r)
  read.csv(file.path(FIT_DIR, sprintf("rep%02d", r), "diagnostics_ez_bias.csv"))))
diag_all$converged <- diag_all$max_rhat < 1.05 & diag_all$divergences == 0 &
                      diag_all$treedepth_hit == 0
print(diag_all)
ok_reps <- diag_all$rep[diag_all$converged]
if (any(!diag_all$converged))
  cat("not converged, excluded from recovery:", diag_all$rep[!diag_all$converged],
      "--> refit with REFIT_REPS\n")

recovery_pop    <- list()
recovery_person <- list()

for (r in ok_reps) {
  rep_name <- sprintf("rep%02d", r)
  true_pop <- readRDS(file.path(DATA_DIR, rep_name, "true_population_parameters.rds"))
  true_sub <- readRDS(file.path(DATA_DIR, rep_name, "true_subject_parameters.rds"))
  summ     <- read.csv(file.path(FIT_DIR, rep_name, "summary_ez_bias.csv"))
  subj_ids <- readRDS(file.path(FIT_DIR, rep_name, "subj_ids.rds"))

  # true between-person SDs are not saved by the simulation script --> taken from it here
  true_sd <- c(sigma_nu_intercept = 0.20, sigma_nu_c = 0.30,
               sigma_alpha_intercept = 0.08, sigma_tau_intercept = 0.05,
               sigma_w = 0.25)                                                      

  mean_ndt_factor <- 1 + true_pop$st0_raw

  # "true" for tau = mean NDT (what EZ estimates), true_t0 = t0 for reference
  truth <- data.frame(
    variable = pop_pars,
    true     = c(true_pop$beta_nu1, true_pop$nu_c, true_pop$beta_alpha1,
                 true_pop$beta_tau * mean_ndt_factor,
                 qlogis(true_pop$w),                                                
                 true_sd[["sigma_nu_intercept"]], true_sd[["sigma_nu_c"]],
                 true_sd[["sigma_alpha_intercept"]],
                 true_sd[["sigma_tau_intercept"]] * mean_ndt_factor,
                 true_sd[["sigma_w"]],                                              
                 true_pop$w),                                                       
    true_t0  = c(NA, NA, NA, true_pop$beta_tau, NA,
                 NA, NA, NA, true_sd[["sigma_tau_intercept"]], NA, NA)
  )

  recovery_pop[[r]] <- summ %>%
    filter(variable %in% pop_pars) %>%
    select(variable, mean, q5, q95, rhat) %>%
    left_join(truth, by = "variable") %>%
    mutate(rep = r, covered_90 = true >= q5 & true <= q95)

  person_est <- summ %>%
    filter(grepl("^(nu_intercept|nu_c|alpha_intercept|tau_intercept|w_person)\\[", variable)) %>%   
    mutate(par  = sub("\\[.*", "", variable),
           i    = as.integer(sub(".*\\[(\\d+)\\]", "\\1", variable)),
           subj = subj_ids[i]) %>%
    select(par, subj, est = mean, q5, q95)

  true_long <- bind_rows(
    data.frame(par = "nu_intercept",    subj = true_sub$subj, true = true_sub$nu_intercept),
    data.frame(par = "nu_c",            subj = true_sub$subj, true = true_sub$nu_c),
    data.frame(par = "alpha_intercept", subj = true_sub$subj, true = true_sub$alpha_intercept),  
    data.frame(par = "tau_intercept",   subj = true_sub$subj, true = true_sub$t0 * mean_ndt_factor),
    data.frame(par = "w_person",        subj = true_sub$subj, true = true_sub$w)             
  )

  recovery_person[[r]] <- person_est %>%
    left_join(true_long, by = c("par", "subj")) %>%
    mutate(rep = r)
}

recovery_pop    <- bind_rows(recovery_pop)
recovery_person <- bind_rows(recovery_person)

write.csv(recovery_pop,    file.path(FIT_DIR, "recovery_population.csv"), row.names = FALSE)
write.csv(recovery_person, file.path(FIT_DIR, "recovery_person.csv"),     row.names = FALSE)

recovery_pop_summary <- recovery_pop %>%
  group_by(variable) %>%
  summarise(true          = first(true),
            true_t0       = first(true_t0),
            est_mean      = mean(mean),
            est_min       = min(mean),
            est_max       = max(mean),
            coverage_90   = mean(covered_90),
            .groups = "drop") %>%
  print()

recovery_person_summary <- recovery_person %>%
  group_by(par, rep) %>%
  summarise(r    = cor(est, true),
            bias = mean(est - true),
            rmse = sqrt(mean((est - true)^2)),
            .groups = "drop") %>%
  group_by(par) %>%
  summarise(across(c(r, bias, rmse), mean), .groups = "drop") %>%
  print()

write.csv(recovery_pop_summary,    file.path(FIT_DIR, "recovery_population_summary.csv"), row.names = FALSE)
write.csv(recovery_person_summary, file.path(FIT_DIR, "recovery_person_summary.csv"),     row.names = FALSE)
write.csv(diag_all,                file.path(FIT_DIR, "diagnostics_all_reps.csv"),        row.names = FALSE)

# scatter plots true vs estimated 
draw_scatter <- function() {
  par(mfrow = c(2, 3))
  for (p in person_pars) {
    d <- recovery_person[recovery_person$par == p, ]
    plot(d$true, d$est, pch = 16, cex = 0.5,
         main = sprintf("%s (r = %.2f)", p, cor(d$true, d$est)),
         xlab = switch(p, tau_intercept = "true mean NDT (1.15 * t0)",
                          w_person      = "true w (0-1)", "true"),
         ylab = "EZ+bias estimate (posterior mean)")
    abline(0, 1, lty = 2)
  }
  par(mfrow = c(1, 1))
}
draw_scatter()
png(file.path(FIT_DIR, "recovery_person_scatter.png"), width = 2400, height = 1600, res = 200)
draw_scatter()
dev.off()

cat("\nSaved in", FIT_DIR, ":\n",
    " recovery_population_summary.csv, recovery_person_summary.csv,\n",
    " diagnostics_all_reps.csv, recovery_person_scatter.png\n")
