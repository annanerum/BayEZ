recover_subject_level <- function(rep_id, data_dir, fits_dir) {
  rep_dir <- file.path(data_dir, sprintf("rep%02d", rep_id))
  out_dir <- file.path(fits_dir, sprintf("rep%02d_fullscale", rep_id))
  fit <- readRDS(file.path(out_dir, "fit.rds"))
  
  min_rt_obs <- readRDS(file.path(rep_dir, "min_rt_observed.rds"))
  sim_data   <- read.csv(file.path(rep_dir, "sim_data_full_ddm.csv"))
  subj_par   <- readRDS(file.path(rep_dir, "true_subject_parameters.rds"))
  subjects_used <- sort(unique(sim_data$subj))
  min_rt_used   <- as.numeric(min_rt_obs[as.character(subjects_used)])
  
  tau_draws   <- fit$draws(variables = "tau_intercept_raw", format = "matrix")
  t0_hat      <- sweep(plogis(tau_draws), 2, min_rt_used, `*`)
  sv_draws    <- fit$draws(variables = "sv", format = "matrix")
  alpha_draws <- exp(fit$draws(variables = "alpha_intercept", format = "matrix"))
  
  true_sv_single <- readRDS(file.path(rep_dir, "true_population_parameters.rds"))$sv
  
  
  t0_pop_draws    <- rowMeans(t0_hat)
  sv_pop_draws    <- rowMeans(sv_draws)
  alpha_pop_draws <- rowMeans(alpha_draws)
  
  subj_df <- tibble(
    rep = rep_id, subj = subjects_used,
    t0_mean = apply(t0_hat, 2, mean), t0_q05 = apply(t0_hat, 2, quantile, 0.05), t0_q95 = apply(t0_hat, 2, quantile, 0.95),
    sv_mean = apply(sv_draws, 2, mean), sv_q05 = apply(sv_draws, 2, quantile, 0.05), sv_q95 = apply(sv_draws, 2, quantile, 0.95),
    alpha_mean = apply(alpha_draws, 2, mean), alpha_q05 = apply(alpha_draws, 2, quantile, 0.05), alpha_q95 = apply(alpha_draws, 2, quantile, 0.95)
  ) %>%
    left_join(subj_par %>% select(subj, t0_true = t0, alpha_true = alpha_intercept), by = "subj") %>%
    mutate(sv_true = true_sv_single, alpha_true = exp(alpha_true)) %>%
    mutate(
      t0_covered    = t0_true    >= t0_q05    & t0_true    <= t0_q95,
      sv_covered    = sv_true    >= sv_q05    & sv_true    <= sv_q95,
      alpha_covered = alpha_true >= alpha_q05 & alpha_true <= alpha_q95
    )
  
  
  attr(subj_df, "pop_draws") <- list(
    t0    = t0_pop_draws,
    sv    = sv_pop_draws,
    alpha = alpha_pop_draws
  )
  
  subj_df
}