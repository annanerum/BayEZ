library(rtdists)
library(dplyr)

n_subj        <- 100
n_raw_trials  <- 101   # one trial lost as warm-up
n_kept_trials <- n_raw_trials - 1   # = 100

true_pop <- list(
  beta_nu1 = 1.8,
  nu_c     = 0.7,
  nu_r     = 0.5,
  nu_cr    = 0.3,
  beta_alpha1 = log(1.8),
  beta_alpha2 = 0.024,
  beta_alpha3 = 0.05,
  beta_tau = 0.30,
  w       = 0.52,
  sw      = 0.15,
  sv      = 0.60,
  st0_raw = 0.15
)

sd_subj <- list(
  sigma_nu_c            = 0.50,
  sigma_nu_r            = 0.30,
  sigma_nu_cr           = 0.30,
  sigma_alpha_intercept = 0.15,
  sigma_tau_intercept   = 0.30
)

## unchanged
simulate_subject <- function(sp) {
  condition <- c(1, sample(rep(c(1, -1), each = n_kept_trials / 2)))
  rt        <- numeric(n_raw_trials)
  acc       <- integer(n_raw_trials)
  resp_type <- integer(n_raw_trials)
  a_trial   <- numeric(n_raw_trials)
  v_trial   <- numeric(n_raw_trials)
  for (i in seq_len(n_raw_trials)) {
    resp_type[i] <- if (i == 1) 1L else acc[i - 1]
    log_a <- sp$alpha_intercept +
      true_pop$beta_alpha2 * resp_type[i] +
      true_pop$beta_alpha3 * condition[i]
    a_trial[i] <- exp(log_a)
    v_trial[i] <- true_pop$beta_nu1 +
      sp$nu_c  * condition[i] +
      sp$nu_r  * resp_type[i] +
      sp$nu_cr * condition[i] * resp_type[i]
    z_abs  <- true_pop$w  * a_trial[i]
    sz_abs <- true_pop$sw * a_trial[i]
    trial_sim <- rdiffusion(
      n = 1,
      a = a_trial[i], v = v_trial[i], t0 = sp$t0,
      z = z_abs, sz = sz_abs, sv = true_pop$sv, st0 = sp$st0
    )
    rt[i]  <- trial_sim$rt
    acc[i] <- ifelse(trial_sim$response == "upper", 1L, -1L)
  }
  data.frame(
    subj = sp$subj, trial = seq_len(n_raw_trials),
    condition = condition, resp_type = resp_type,
    a = a_trial, v = v_trial,
    rt = rt, acc = acc
  )
}

## t0 drawn directly, no fictional min_rt ceiling.
## t0_scale is a fixed constant (not per-subject, not random) that sets the
## overall range t0 can occupy given true_pop$beta_tau and sd_subj$sigma_tau_intercept.
## Picked so mean(t0) stays close to what you had before (~0.20s):
##   plogis(0.30) * t0_scale ~ 0.20  ->  t0_scale ~ 0.35
t0_scale <- 0.35

run_one_rep <- function(rep_id, seed) {
  set.seed(seed)
  z_nu_c  <- rnorm(n_subj); z_nu_r <- rnorm(n_subj); z_nu_cr <- rnorm(n_subj)
  z_alpha <- rnorm(n_subj); z_tau  <- rnorm(n_subj)
  subj_par <- data.frame(
    subj = 1:n_subj,
    nu_c  = true_pop$nu_c  + sd_subj$sigma_nu_c  * z_nu_c,
    nu_r  = true_pop$nu_r  + sd_subj$sigma_nu_r  * z_nu_r,
    nu_cr = true_pop$nu_cr + sd_subj$sigma_nu_cr * z_nu_cr,
    alpha_intercept   = true_pop$beta_alpha1 + sd_subj$sigma_alpha_intercept * z_alpha,
    tau_intercept_raw = true_pop$beta_tau    + sd_subj$sigma_tau_intercept   * z_tau
  )
  
  # no more sp$min_rt ~ N(0.35, 0.03). t0 is drawn directly, on its
  # own generative scale, with no reference to any RT-derived quantity.
  subj_par$t0  <- plogis(subj_par$tau_intercept_raw) * t0_scale
  subj_par$st0 <- true_pop$st0_raw * 2 * subj_par$t0
  
  # NEW: sanity-check t0_scale immediately, per rep, before simulating any trials
  t0_range <- range(subj_par$t0)
  cat(sprintf("  rep%02d t0 check: mean=%.3f, sd=%.3f, range=[%.3f, %.3f]\n",
              rep_id, mean(subj_par$t0), sd(subj_par$t0), t0_range[1], t0_range[2]))
  if (t0_range[1] < 0.05 || t0_range[2] > 0.45) {
    warning(sprintf("rep%02d: t0 range looks off (%.3f-%.3f) -- check t0_scale",
                    rep_id, t0_range[1], t0_range[2]))
  }
  
  sim_list <- lapply(seq_len(n_subj), function(i) simulate_subject(subj_par[i, ]))
  sim_data_raw <- do.call(rbind, sim_list)
  sim_data <- sim_data_raw[sim_data_raw$trial != 1, ]
  sim_data <- sim_data[order(sim_data$subj, sim_data$trial), ]
  rownames(sim_data) <- NULL
  stopifnot(all(table(sim_data$subj) == n_kept_trials))
  
  # min_rt_observed is now purely an IDENTIFICATION device for the fitting
  # model -- it plays no role in how t0 was generated above.
  min_rt_obs <- tapply(sim_data$rt, sim_data$subj, min)
  
  # NEW: confirm t0 stays comfortably below the RTs it's supposed to be a floor for
  # (t0 must physically be less than every observed RT for that subject)
  frac_t0_above_min_rt <- mean(subj_par$t0 > min_rt_obs[as.character(subj_par$subj)])
  if (frac_t0_above_min_rt > 0) {
    warning(sprintf("rep%02d: t0 exceeds min_rt_observed for %.1f%% of subjects -- t0_scale too high",
                    rep_id, 100 * frac_t0_above_min_rt))
  }
  
  out_dir <- file.path("data", sprintf("rep%02d", rep_id))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  write.csv(sim_data, file.path(out_dir, "sim_data_full_ddm.csv"), row.names = FALSE)
  saveRDS(subj_par,   file.path(out_dir, "true_subject_parameters.rds"))
  saveRDS(true_pop,   file.path(out_dir, "true_population_parameters.rds"))
  saveRDS(min_rt_obs, file.path(out_dir, "min_rt_observed.rds"))
  
  cat(sprintf("rep%02d done | seed=%d | trials=%d | mean_rt=%.3f | mean_t0=%.3f | acc=%.3f\n",
              rep_id, seed, nrow(sim_data), mean(sim_data$rt), mean(subj_par$t0), mean(sim_data$acc == 1)))
  invisible(sim_data)
}

n_reps <- 10
base_seed <- 2026
for (r in seq_len(n_reps)) {
  run_one_rep(rep_id = r, seed = base_seed + r)
}

write.csv(subj_par[, c("subj", "t0", "st0")],
          "data/subj_par_export.csv", row.names = FALSE)

# NEW: combined t0_scale sanity check across all 10 reps, run once at the end
all_t0_check <- bind_rows(lapply(1:n_reps, function(r) {
  sp <- readRDS(file.path("data", sprintf("rep%02d", r), "true_subject_parameters.rds"))
  data.frame(rep = r, t0 = sp$t0)
}))
cat("\n== t0_scale sanity check across all reps ==\n")
print(summary(all_t0_check$t0))
cat(sprintf("Any t0 outside [0.05, 0.45]: %s\n",
            any(all_t0_check$t0 < 0.05 | all_t0_check$t0 > 0.45)))

cat("\nAll", n_reps, "saved under data/../\n")