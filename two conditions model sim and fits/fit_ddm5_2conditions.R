library(cmdstanr)
library(dplyr)
library(posterior)

PROJECT_DIR <- "/Users/annegiacobello/Desktop/EZ Paper "  
DATA_DIR    <- file.path(PROJECT_DIR, "data_sim_fullddm_2condtions")
FIT_DIR     <- file.path(PROJECT_DIR, "fits_ddm5_2conditions")
STAN_FILE   <- file.path(PROJECT_DIR, "ddm5_2conditions.stan")

stopifnot("nu_intercept" %in% names(readRDS(file.path(DATA_DIR, "rep01", "true_subject_parameters.rds"))))
n_reps <- 10


TEST_MODE     <- FALSE
TEST_N_PERSON <- 20

dir.create(FIT_DIR, showWarnings = FALSE, recursive = TRUE)

# within-chain threading 
mod <- cmdstan_model(STAN_FILE, cpp_options = list(stan_threads = TRUE))


#### data 

COND_LEVELS <- c(1, -1)

# trial-level data for Stan 
make_stan_data <- function(dat) {
  subj_ids <- sort(unique(dat$subj))
  dat$pid  <- match(dat$subj, subj_ids)
  min_rt   <- tapply(dat$rt, dat$pid, min)            # fastest RT per person (t0 must stay below this)
  list(
    N = nrow(dat), I = length(subj_ids), pid = dat$pid,
    rt = dat$rt, acc = as.integer(dat$acc), condition = dat$condition,
    min_rt = as.numeric(min_rt[as.character(seq_along(subj_ids))]),
    grainsize = 1,
    subj_ids = subj_ids
  )
}

# EZ statistics for the pilot prior derivation
make_ez_classic_stats <- function(dat) {
  cells <- dat %>%
    mutate(k = match(condition, COND_LEVELS)) %>%
    group_by(subj, k) %>%
    summarise(J = n(), C = sum(acc == 1),
              MRT = mean(rt[acc == 1]), VRT = var(rt[acc == 1]), .groups = "drop") %>%
    arrange(subj, k)
  I <- length(unique(cells$subj)); K <- length(COND_LEVELS)
  to_mat <- function(x) matrix(x, nrow = I, ncol = K, byrow = TRUE)
  list(I = I, K = K, J = to_mat(cells$J), C = to_mat(cells$C),
       MRT = to_mat(cells$MRT), VRT = to_mat(cells$VRT))
}


#### priors 

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
  nu_vals        <- as.vector(ez$nu[!is.na(ez$nu)])
  log_alpha_vals <- log(as.vector(ez$alpha[!is.na(ez$alpha) & ez$alpha > 0]))
  tau_vals       <- as.vector(ez$tau[!is.na(ez$tau) & ez$tau > 0.05])
  tau_sd         <- min(max(sd(tau_vals), 0.1),  0.5)     
  tau_scale      <- min(max(sd(tau_vals), 0.03), 0.5)    
  
  nu_sd           <- min(max(sd(nu_vals),        0.3),  1.5)
  log_alpha_sd    <- min(max(sd(log_alpha_vals), 0.1),  0.7)
  log_alpha_scale <- min(max(sd(log_alpha_vals), 0.05), 0.7)
  
  list(
    prior_beta_nu_intercept     = c(mean(nu_vals),        nu_sd),
    prior_beta_nu_c             = c(0,                    nu_sd),
    prior_beta_alpha            = c(mean(log_alpha_vals), 2 * log_alpha_sd),
    prior_beta_tau              = c(mean(tau_vals),       2 * tau_sd),   # seconds
    prior_beta_w                = c(0,                    0.5),          # w near 0.5
    prior_beta_sv               = c(log(0.6),             0.5),          # sv about 0.2 - 1.6
    prior_sigma_nu_intercept    = c(0.3 * nu_sd,          nu_sd),
    prior_sigma_nu_c            = c(0.3 * nu_sd,          nu_sd),
    prior_sigma_alpha_intercept = c(0.3 * log_alpha_scale, log_alpha_scale),
    prior_sigma_tau_intercept   = c(0.3 * tau_scale,      tau_scale),    # seconds 
    prior_sigma_w               = c(0.3 * 0.4,            0.4),          
    prior_sigma_sv              = c(0.3 * 0.3,            0.3)
  )
}

derive_priors_from_pilot <- function(r, I_pilot = 20) {
  pilot_rep <- if (r == n_reps) 1 else r + 1
  dat_p     <- read.csv(file.path(DATA_DIR, sprintf("rep%02d", pilot_rep), "sim_data_full_ddm.csv"))
  set.seed(42 + r)
  pilot_subj <- sample(unique(dat_p$subj), I_pilot)
  sd_p <- make_ez_classic_stats(dat_p[dat_p$subj %in% pilot_subj, ])
  ez   <- ez_point_estimates(sd_p$C, sd_p$J, sd_p$MRT, sd_p$VRT, sd_p$I, sd_p$K)
  pr   <- derive_priors(ez)
  attr(pr, "pilot_rep") <- pilot_rep
  pr
}


#### fit 

reps_to_fit <- if (TEST_MODE) 1 else seq_len(n_reps)
REFIT_REPS  <- NULL       
if (!is.null(REFIT_REPS)) reps_to_fit <- REFIT_REPS

make_inits <- function(priors, I, n_chains, min_rt) {
  t0_init <- 0.7 * min_rt                       # t0 starts below each person's fastest RT
  lapply(seq_len(n_chains), function(ch) list(
    t0                    = t0_init,
    beta_nu_intercept     = priors$prior_beta_nu_intercept[1] + rnorm(1, 0, 0.05),
    beta_nu_c             = 0.3 + rnorm(1, 0, 0.05),
    beta_alpha            = priors$prior_beta_alpha[1] + rnorm(1, 0, 0.05),
    beta_tau              = min(max(mean(t0_init), 0.06), 0.44),
    beta_w                = rnorm(1, 0, 0.05),
    beta_sv               = log(0.6) + rnorm(1, 0, 0.05),
    sigma_nu_intercept    = 0.2,
    sigma_nu_c            = 0.2,
    sigma_alpha_intercept = 0.1,
    sigma_tau_intercept   = 0.03,
    sigma_w               = 0.2,
    sigma_sv              = 0.1,
    z_nu_intercept = rnorm(I, 0, 0.1), z_nu_c = rnorm(I, 0, 0.1),
    z_alpha_intercept = rnorm(I, 0, 0.1),
    z_w = rnorm(I, 0, 0.1), z_sv = rnorm(I, 0, 0.1)
  ))
}

sampler_args <- if (TEST_MODE) {
  list(chains = 2, parallel_chains = 2, iter_warmup = 300, iter_sampling = 300)
} else {
  list(chains = 4, parallel_chains = 4, iter_warmup = 1000, iter_sampling = 1000)
}
THREADS_PER_CHAIN <- max(1, floor((parallel::detectCores() - 1) / sampler_args$parallel_chains))
cat("threads per chain:", THREADS_PER_CHAIN, "\n")

pop_pars    <- c("beta_nu_intercept", "beta_nu_c", "beta_alpha", "beta_w", "beta_sv",
                 "sigma_nu_intercept", "sigma_nu_c", "sigma_alpha_intercept", "sigma_w", "sigma_sv",
                 "beta_tau", "sigma_tau_intercept",    
                 "w_pop_median", "sv_pop_median")
person_pars <- c("nu_intercept", "nu_c", "alpha_intercept", "t0", "w", "sv")

for (r in reps_to_fit) {
  rep_name <- sprintf("rep%02d", r)
  cat("\n ", rep_name, " \n")
  
  dat <- read.csv(file.path(DATA_DIR, rep_name, "sim_data_full_ddm.csv"))
  if (TEST_MODE) dat <- dat[dat$subj %in% sort(unique(dat$subj))[seq_len(TEST_N_PERSON)], ]
  
  stan_input <- make_stan_data(dat)
  priors     <- derive_priors_from_pilot(r)
  cat("priors derived from 20 persons of rep", attr(priors, "pilot_rep"), "\n")
  print(round(do.call(rbind, priors), 3))
  stan_data  <- c(stan_input[setdiff(names(stan_input), "subj_ids")], priors)
  
  set.seed(2026 + r)
  t_start <- Sys.time()
  fit <- do.call(mod$sample, c(list(
    data              = stan_data,
    init              = make_inits(priors, stan_input$I, sampler_args$chains, stan_input$min_rt),
    seed              = 2026 + r,
    threads_per_chain = THREADS_PER_CHAIN,
    adapt_delta       = 0.9,
    max_treedepth     = 10,
    refresh           = 100
  ), sampler_args))
  cat("fit time:", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), "min\n")
  
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
  
  fit$save_object(file.path(out_dir, "fit_ddm5.rds"))
  write.csv(summ,     file.path(out_dir, "summary_ddm5.csv"),     row.names = FALSE)
  write.csv(diag_row, file.path(out_dir, "diagnostics_ddm5.csv"), row.names = FALSE)
  saveRDS(stan_input$subj_ids, file.path(out_dir, "subj_ids.rds"))
  saveRDS(priors,              file.path(out_dir, "priors_ddm5.rds"))
}


#### parameter recovery 

all_reps <- if (TEST_MODE) 1 else seq_len(n_reps)
diag_all <- bind_rows(lapply(all_reps, function(r)
  read.csv(file.path(FIT_DIR, sprintf("rep%02d", r), "diagnostics_ddm5.csv"))))
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
  summ     <- read.csv(file.path(FIT_DIR, rep_name, "summary_ddm5.csv"))
  subj_ids <- readRDS(file.path(FIT_DIR, rep_name, "subj_ids.rds"))
  true_sub <- true_sub[true_sub$subj %in% subj_ids, ]           
  
  true_sd <- c(sigma_nu_intercept = 0.20, sigma_nu_c = 0.30,
               sigma_alpha_intercept = 0.08, sigma_w = 0.25)
  mean_ndt_factor <- 1 + true_pop$st0_raw
  
  # t0: the model has no st0, so its t0 lies somewhere between the true t0 (lower edge of the simulated NDT)
  # and the mean NDT --> "true" = t0, "true_mean_ndt" as second reference
  truth <- data.frame(
    variable = pop_pars,
    true = c(true_pop$beta_nu1, true_pop$nu_c, true_pop$beta_alpha1,
             qlogis(true_pop$w), log(true_pop$sv),
             true_sd[["sigma_nu_intercept"]], true_sd[["sigma_nu_c"]],
             true_sd[["sigma_alpha_intercept"]], true_sd[["sigma_w"]],
             0,                                                  
             true_pop$beta_tau, 0.05,                           
             true_pop$w, true_pop$sv),
    true_mean_ndt = c(rep(NA, 10),
                      true_pop$beta_tau * mean_ndt_factor, 0.05 * mean_ndt_factor,
                      NA, NA)
  )
  
  recovery_pop[[r]] <- summ %>%
    filter(variable %in% pop_pars) %>%
    select(variable, mean, q5, q95, rhat) %>%
    left_join(truth, by = "variable") %>%
    mutate(rep = r, covered_90 = true >= q5 & true <= q95)
  
  person_est <- summ %>%
    filter(grepl("^(nu_intercept|nu_c|alpha_intercept|t0|w|sv)\\[", variable)) %>%
    mutate(par  = sub("\\[.*", "", variable),
           i    = as.integer(sub(".*\\[(\\d+)\\]", "\\1", variable)),
           subj = subj_ids[i]) %>%
    select(par, subj, est = mean, q5, q95)
  
  true_long <- bind_rows(
    data.frame(par = "nu_intercept",    subj = true_sub$subj, true = true_sub$nu_intercept),
    data.frame(par = "nu_c",            subj = true_sub$subj, true = true_sub$nu_c),
    data.frame(par = "alpha_intercept", subj = true_sub$subj, true = true_sub$alpha_intercept),  
    data.frame(par = "t0",              subj = true_sub$subj, true = true_sub$t0),              
    data.frame(par = "w",               subj = true_sub$subj, true = true_sub$w),
    data.frame(par = "sv",              subj = true_sub$subj, true = true_pop$sv)               
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
            true_mean_ndt = first(true_mean_ndt),
            est_mean      = mean(mean),
            est_min       = min(mean),
            est_max       = max(mean),
            coverage_90   = mean(covered_90),
            .groups = "drop") %>%
  print(n = Inf)

# person level
recovery_person_summary <- recovery_person %>%
  group_by(par, rep) %>%
  summarise(r    = if (sd(true) > 0) cor(est, true) else NA_real_,
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
  for (p in c("nu_intercept", "nu_c", "alpha_intercept", "t0", "w")) {
    d <- recovery_person[recovery_person$par == p, ]
    plot(d$true, d$est, pch = 16, cex = 0.5,
         main = sprintf("%s (r = %.2f)", p, cor(d$true, d$est)),
         xlab = switch(p, t0 = "true t0 (s)", w = "true w (0-1)", "true"),
         ylab = "DDM5 estimate (posterior mean)")
    abline(0, 1, lty = 2)
  }
  d <- recovery_person[recovery_person$par == "sv", ]
  
  hist(d$est, breaks = 15, main = "sv per person (dashed = true 0.5)",
       xlab = "DDM5 estimate (posterior mean)",
       xlim = range(c(d$est, true_pop$sv)) + c(-0.02, 0.02))
  abline(v = true_pop$sv, lty = 2, lwd = 2)
  par(mfrow = c(1, 1))
}
draw_scatter()
png(file.path(FIT_DIR, "recovery_person_scatter.png"), width = 2400, height = 1600, res = 200)
draw_scatter()
dev.off()

cat("\nSaved in", FIT_DIR, ":\n",
    " recovery_population_summary.csv, recovery_person_summary.csv,\n",
    " diagnostics_all_reps.csv, recovery_person_scatter.png\n")