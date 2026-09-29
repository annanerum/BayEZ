library(cmdstanr)
library(dplyr)
library(posterior)

PROJECT_DIR <- "/Users/annegiacobello/Desktop/EZ Paper "   # same as in the simulation script (note the trailing space)
DATA_DIR    <- file.path(PROJECT_DIR, "data_sim_fullddm_2condtions")
FIT_DIR     <- file.path(PROJECT_DIR, "fits_ez_2conditions")
STAN_FILE   <- file.path(PROJECT_DIR, "ez_ddm_2conditions.stan")

# safety check, stop if the data are not from the current simulation (with person-specific drift intercept)
stopifnot("nu_intercept" %in% names(readRDS(file.path(DATA_DIR, "rep01", "true_subject_parameters.rds"))))
n_reps    <- 10

# TEST_MODE = TRUE --> only rep01, 2 chains, few iterations 
TEST_MODE <- TRUE

dir.create(FIT_DIR, showWarnings = FALSE, recursive = TRUE)

mod <- cmdstan_model(STAN_FILE)


#### summary statistics per person and cell 
# cell order must match the Stan model: k = 1 --> cond +1, k = 2 --> cond -1
COND_LEVELS <- c(1, -1)

make_stan_data <- function(dat) {
  cells <- dat %>%
    mutate(k = match(condition, COND_LEVELS)) %>%
    group_by(subj, k) %>%
    summarise(J   = n(),
              C   = sum(acc == 1),
              MRT = mean(rt[acc == 1]),   
              VRT = var(rt[acc == 1]),
              .groups = "drop") %>%
    arrange(subj, k)
  
  # EZ needs at least 2 correct RTs per cell for a variance
  if (any(cells$C < 2)) stop("some cells have fewer than 2 correct responses")
  
  subj_ids <- sort(unique(cells$subj))
  I <- length(subj_ids)
  K <- length(COND_LEVELS)
  
  to_mat <- function(x) matrix(x, nrow = I, ncol = K, byrow = TRUE)  # rows = persons, cols = cells
  
  list(
    I = I, K = K, cond = COND_LEVELS,
    J   = to_mat(cells$J),
    C   = to_mat(cells$C),
    MRT = to_mat(cells$MRT),
    VRT = to_mat(cells$VRT),
    subj_ids = subj_ids          # for matching persons later
  )
}


#### priors (as in the first simulation, model2_v2_ncp) 
# same as derive_priors_from_pilot() / derive_priors():
# EZ point estimates (closed-form equations) for a pilot of 20 persons
# prior means from the average estimates, prior SDs from their spread
# pilot persons are taken from a different rep than the one being fitted
# so the data are not used twice

# EZ point estimates 
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

# mapping new parameters (first simulation --> this model):
# beta_nu[1]  N(mean nu, nu_sd)            --> beta_nu_intercept (same)
# beta_nu[2]  N(0, nu_sd)                  --> beta_nu_c (same)
# beta_alpha  N(mean alpha, 2*alpha_sd)    --> beta_alpha (same, on log(alpha))
# beta_tau    N(mean tau, 2*tau_sd)        --> beta_tau (same)
# sigma_v     N(0.3*nu_sd, nu_sd)          --> sigma_nu_c (same)
#                                          --> sigma_nu_intercept (new, same form as sigma_v)
# sigma_alpha N(alpha_sd, alpha_sd/2)      --> sigma_alpha_intercept N(0.3*s, s), same form as sigma_v (log scale)
# sigma_tau   N(tau_sd, tau_sd/2)          --> sigma_tau_intercept   N(0.3*s, s), same form as sigma_v
# all SD priors are normal (half-normal via lower=0 in Stan)

derive_priors <- function(ez) {
  nu_vals    <- as.vector(ez$nu[!is.na(ez$nu)])
  alpha_vals <- as.vector(ez$alpha[!is.na(ez$alpha) & ez$alpha > 0])
  tau_vals   <- as.vector(ez$tau[!is.na(ez$tau) & ez$tau > 0.05])
  
  log_alpha_vals <- log(alpha_vals)   # alpha is on the log scale here
  
  nu_sd        <- min(max(sd(nu_vals),        0.3), 1.5)   # same clamps as before
  log_alpha_sd <- min(max(sd(log_alpha_vals), 0.1), 0.7)   # new clamps: before 0.2-1.0 on the natural scale,
  # relative to alpha around 1.5 --> about 0.1 - 0.7
  # on the log scale
  tau_sd       <- min(max(sd(tau_vals),       0.1), 0.5)   # same clamps as before
  
  # scales for the between-person SD priors of alpha and tau
  # pilot SD of EZ point estimates = true between-person SD + noise
  # --> upper-end guess, used as prior scale and not as prior mean
  # small floors, so the prior is not pushed above small true SDs (0.08, 0.05)
  log_alpha_scale <- min(max(sd(log_alpha_vals), 0.05), 0.7)
  tau_scale       <- min(max(sd(tau_vals),       0.03), 0.5)
  
  list(
    prior_beta_nu_intercept     = c(mean(nu_vals),        nu_sd),
    prior_beta_nu_c             = c(0,                    nu_sd),
    prior_beta_alpha            = c(mean(log_alpha_vals), 2 * log_alpha_sd),
    prior_beta_tau              = c(mean(tau_vals),       2 * tau_sd),
    prior_sigma_nu_intercept    = c(0.3 * nu_sd,          nu_sd),
    prior_sigma_nu_c            = c(0.3 * nu_sd,          nu_sd),
    # same form as sigma_v or sigma_nu_c: N(0.3 * scale, scale)
    prior_sigma_alpha_intercept = c(0.3 * log_alpha_scale, log_alpha_scale),
    prior_sigma_tau_intercept   = c(0.3 * tau_scale,       tau_scale)
  )
}

# pilot: I_pilot persons from another rep (the next one, rep10 --> rep01)
derive_priors_from_pilot <- function(r, I_pilot = 20) {
  pilot_rep <- if (r == n_reps) 1 else r + 1
  dat_p     <- read.csv(file.path(DATA_DIR, sprintf("rep%02d", pilot_rep), "sim_data_full_ddm.csv"))
  
  set.seed(42 + r)
  pilot_subj <- sample(unique(dat_p$subj), I_pilot)
  sd_p <- make_stan_data(dat_p[dat_p$subj %in% pilot_subj, ])
  
  ez <- ez_point_estimates(sd_p$C, sd_p$J, sd_p$MRT, sd_p$VRT, sd_p$I, sd_p$K)
  pr <- derive_priors(ez)
  attr(pr, "pilot_rep") <- pilot_rep
  pr
}


#### fit all reps 

reps_to_fit <- if (TEST_MODE) 1 else seq_len(n_reps)

# refit only selected reps (e.g. a non-converged one)
# NULL to fit all
# the other reps' saved fits are kept and still used in the recovery part
REFIT_REPS <- NULL          
if (!is.null(REFIT_REPS)) reps_to_fit <- REFIT_REPS

make_inits <- function(priors, I, n_chains) {
  lapply(seq_len(n_chains), function(ch) list(
    beta_nu_intercept     = priors$prior_beta_nu_intercept[1] + rnorm(1, 0, 0.05),
    beta_nu_c             = 0.3 + rnorm(1, 0, 0.05),
    beta_alpha            = priors$prior_beta_alpha[1] + rnorm(1, 0, 0.05),
    beta_tau              = min(max(priors$prior_beta_tau[1], 0.10), 0.40) + rnorm(1, 0, 0.01),
    sigma_nu_intercept    = 0.2,
    sigma_nu_c            = 0.2,
    sigma_alpha_intercept = 0.1,
    sigma_tau_intercept   = 0.03,
    z_nu_intercept        = rnorm(I, 0, 0.1),
    z_nu_c                = rnorm(I, 0, 0.1),
    z_alpha_intercept     = rnorm(I, 0, 0.1),
    z_tau_intercept       = rnorm(I, 0, 0.1)
  ))
}

sampler_args <- if (TEST_MODE) {
  list(chains = 2, parallel_chains = 2, iter_warmup = 300, iter_sampling = 300)
} else {
  list(chains = 4, parallel_chains = 4, iter_warmup = 1000, iter_sampling = 1000)
}

pop_pars    <- c("beta_nu_intercept", "beta_nu_c", "beta_alpha", "beta_tau",
                 "sigma_nu_intercept", "sigma_nu_c", "sigma_alpha_intercept", "sigma_tau_intercept")
person_pars <- c("nu_intercept", "nu_c", "alpha_intercept", "tau_intercept")

for (r in reps_to_fit) {
  rep_name <- sprintf("rep%02d", r)
  cat("\n ", rep_name, " \n")
  
  dat        <- read.csv(file.path(DATA_DIR, rep_name, "sim_data_full_ddm.csv"))
  stan_input <- make_stan_data(dat)
  priors     <- derive_priors_from_pilot(r)                         #  priors per rep, from pilot
  cat("priors derived from 20 persons of rep", attr(priors, "pilot_rep"), "\n")
  print(round(do.call(rbind, priors), 3))                           # rows = priors, cols = [mean, SD]
  stan_data  <- c(stan_input[setdiff(names(stan_input), "subj_ids")], priors)
  
  set.seed(2026 + r)
  fit <- do.call(mod$sample, c(list(
    data          = stan_data,
    init          = make_inits(priors, stan_input$I, sampler_args$chains),   # NEW
    seed          = 2026 + r,
    adapt_delta   = 0.95,
    max_treedepth = 12,
    refresh       = 500
  ), sampler_args))
  
  out_dir <- file.path(FIT_DIR, rep_name)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # convergence diagnostics
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
  
  # save
  fit$save_object(file.path(out_dir, "fit_ez.rds"))       # full(all draws)
  write.csv(summ,     file.path(out_dir, "summary_ez.csv"),     row.names = FALSE)
  write.csv(diag_row, file.path(out_dir, "diagnostics_ez.csv"), row.names = FALSE)
  saveRDS(stan_input$subj_ids, file.path(out_dir, "subj_ids.rds"))
  saveRDS(priors,              file.path(out_dir, "priors_ez.rds"))   # priors used for this rep
}


#### parameter recovery

# use all reps that have a saved fit, and only the converged ones
all_reps <- if (TEST_MODE) 1 else seq_len(n_reps)
diag_all <- bind_rows(lapply(all_reps, function(r)
  read.csv(file.path(FIT_DIR, sprintf("rep%02d", r), "diagnostics_ez.csv"))))
diag_all$converged <- diag_all$max_rhat < 1.05 & diag_all$divergences == 0 &
  diag_all$treedepth_hit == 0
print(diag_all)
ok_reps <- diag_all$rep[diag_all$converged]
if (any(!diag_all$converged))
  cat("NOT converged, excluded from recovery:", diag_all$rep[!diag_all$converged],
      "--> refit with REFIT_REPS\n")

recovery_pop    <- list()
recovery_person <- list()

for (r in ok_reps) {
  rep_name <- sprintf("rep%02d", r)
  true_pop <- readRDS(file.path(DATA_DIR, rep_name, "true_population_parameters.rds"))
  true_sub <- readRDS(file.path(DATA_DIR, rep_name, "true_subject_parameters.rds"))
  summ     <- read.csv(file.path(FIT_DIR, rep_name, "summary_ez.csv"))
  subj_ids <- readRDS(file.path(FIT_DIR, rep_name, "subj_ids.rds"))
  
  # true between-person SDs are not saved by the simulation script --> taken from it here
  true_sd <- c(sigma_nu_intercept = 0.20, sigma_nu_c = 0.30,
               sigma_alpha_intercept = 0.08, sigma_tau_intercept = 0.05)
  
  mean_ndt_factor <- 1 + true_pop$st0_raw   # mean NDT = t0 + st0/2 = (1 + st0_raw) * t0
  
  # "true" for tau = mean NDT (what EZ estimates), true_t0 = t0 for reference
  truth <- data.frame(
    variable = pop_pars,
    true     = c(true_pop$beta_nu1, true_pop$nu_c, true_pop$beta_alpha1,
                 true_pop$beta_tau * mean_ndt_factor,
                 true_sd[["sigma_nu_intercept"]], true_sd[["sigma_nu_c"]],
                 true_sd[["sigma_alpha_intercept"]],
                 true_sd[["sigma_tau_intercept"]] * mean_ndt_factor),
    true_t0  = c(NA, NA, NA, true_pop$beta_tau,
                 NA, NA, NA, true_sd[["sigma_tau_intercept"]])
  )
  
  recovery_pop[[r]] <- summ %>%
    filter(variable %in% pop_pars) %>%
    select(variable, mean, q5, q95, rhat) %>%
    left_join(truth, by = "variable") %>%
    mutate(rep = r, covered_90 = true >= q5 & true <= q95)
  
  # person level: Stan index i = position in subj_ids
  person_est <- summ %>%
    filter(grepl("^(nu_intercept|nu_c|alpha_intercept|tau_intercept)\\[", variable)) %>%
    mutate(par  = sub("\\[.*", "", variable),
           i    = as.integer(sub(".*\\[(\\d+)\\]", "\\1", variable)),
           subj = subj_ids[i]) %>%
    select(par, subj, est = mean, q5, q95)
  
  true_long <- bind_rows(
    data.frame(par = "nu_intercept",    subj = true_sub$subj, true = true_sub$nu_intercept),
    data.frame(par = "nu_c",            subj = true_sub$subj, true = true_sub$nu_c),
    data.frame(par = "alpha_intercept", subj = true_sub$subj, true = true_sub$alpha_intercept),  # log scale
    # compared to each person's mean NDT (= (1 + st0_raw) * t0), what EZ estimates
    data.frame(par = "tau_intercept",   subj = true_sub$subj, true = true_sub$t0 * mean_ndt_factor)
  )
  
  recovery_person[[r]] <- person_est %>%
    left_join(true_long, by = c("par", "subj")) %>%
    mutate(rep = r)
}

recovery_pop    <- bind_rows(recovery_pop)
recovery_person <- bind_rows(recovery_person)

write.csv(recovery_pop,    file.path(FIT_DIR, "recovery_population.csv"), row.names = FALSE)
write.csv(recovery_person, file.path(FIT_DIR, "recovery_person.csv"),     row.names = FALSE)

# population level: estimate vs truth, averaged over reps
recovery_pop_summary <- recovery_pop %>%     
  group_by(variable) %>%
  summarise(true          = first(true),      # for tau: mean NDT
            true_t0       = first(true_t0),   # for tau: t0, reference only
            est_mean      = mean(mean),
            est_min       = min(mean),
            est_max       = max(mean),
            coverage_90   = mean(covered_90),
            .groups = "drop") %>%
  print()

# person level: correlation, bias and RMSE per parameter (per rep, then averaged)
recovery_person_summary <- recovery_person %>%    
  group_by(par, rep) %>%
  summarise(r    = cor(est, true),
            bias = mean(est - true),
            rmse = sqrt(mean((est - true)^2)),
            .groups = "drop") %>%
  group_by(par) %>%
  summarise(across(c(r, bias, rmse), mean), .groups = "drop") %>%
  print()

# save the summary tables and the convergence table
write.csv(recovery_pop_summary,    file.path(FIT_DIR, "recovery_population_summary.csv"), row.names = FALSE)
write.csv(recovery_person_summary, file.path(FIT_DIR, "recovery_person_summary.csv"),     row.names = FALSE)
write.csv(diag_all,                file.path(FIT_DIR, "diagnostics_all_reps.csv"),        row.names = FALSE)

# scatter plots true vs estimated (person level, all fitted reps pooled)
# drawn twice 
draw_scatter <- function() {
  par(mfrow = c(2, 2))
  for (p in person_pars) {
    d <- recovery_person[recovery_person$par == p, ]
    plot(d$true, d$est, pch = 16, cex = 0.5,
         main = sprintf("%s (r = %.2f)", p, cor(d$true, d$est)),
         xlab = if (p == "tau_intercept") "true mean NDT (1.15 * t0)" else "true",
         ylab = "EZ estimate (posterior mean)")
    abline(0, 1, lty = 2)
  }
  par(mfrow = c(1, 1))
}
draw_scatter()
png(file.path(FIT_DIR, "recovery_person_scatter.png"), width = 2000, height = 1600, res = 200)
draw_scatter()
dev.off()

cat("\nSaved in", FIT_DIR, ":\n",
    " recovery_population_summary.csv, recovery_person_summary.csv,\n",
    " diagnostics_all_reps.csv, recovery_person_scatter.png\n")