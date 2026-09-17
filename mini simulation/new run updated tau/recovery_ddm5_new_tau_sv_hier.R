library(dplyr)
library(purrr)
library(tidyr)
library(tibble)
library(ggplot2)


DATA_DIR <- "data_full7_globalminrt"
FITS_DIR <- "fits/M2 sv hier globalminrt"


sd_subj_true <- list(
  sigma_nu_c = 0.50, sigma_nu_r = 0.30, sigma_nu_cr = 0.30,
  sigma_alpha_intercept = 0.08, sigma_tau_intercept = 0.05
)


recover_subject_level <- function(rep_id, data_dir, fits_dir) {
  rep_dir <- file.path(data_dir, sprintf("rep%02d", rep_id))
  out_dir <- file.path(fits_dir, sprintf("rep%02d_fullscale", rep_id))
  fit <- readRDS(file.path(out_dir, "fit.rds"))
  
  global_min_rt <- readRDS(file.path(rep_dir, "global_min_rt.rds"))
  sim_data      <- read.csv(file.path(rep_dir, "sim_data_full_ddm.csv"))
  subj_par      <- readRDS(file.path(rep_dir, "true_subject_parameters.rds"))
  subjects_used <- sort(unique(sim_data$subj))
  
  tau_draws   <- fit$draws(variables = "tau_intercept_raw", format = "matrix")
  t0_hat      <- plogis(tau_draws) * global_min_rt   
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
  
  attr(subj_df, "pop_draws") <- list(t0 = t0_pop_draws, sv = sv_pop_draws, alpha = alpha_pop_draws)
  subj_df
}


build_recovery_table <- function(rep_id, data_dir, fits_dir) {
  out_dir  <- file.path(fits_dir, sprintf("rep%02d_fullscale", rep_id))
  summ     <- read.csv(file.path(out_dir, "summary.csv"))
  true_pop <- readRDS(file.path(data_dir, sprintf("rep%02d", rep_id), "true_population_parameters.rds"))
  
  true_vals <- c(
    "beta_nu[1]" = true_pop$beta_nu1, "beta_nu[2]" = true_pop$nu_c,
    "beta_nu[3]" = true_pop$nu_r,     "beta_nu[4]" = true_pop$nu_cr,
    beta_alpha2 = true_pop$beta_alpha2, beta_alpha3 = true_pop$beta_alpha3,
    beta_tau = true_pop$beta_tau,               # raw logit-scale check
    w = true_pop$w,
    sigma_nu_c = sd_subj_true$sigma_nu_c, sigma_nu_r = sd_subj_true$sigma_nu_r,
    sigma_nu_cr = sd_subj_true$sigma_nu_cr,
    sigma_alpha_intercept = sd_subj_true$sigma_alpha_intercept,
    sigma_tau_intercept   = sd_subj_true$sigma_tau_intercept
  )
  
  summ %>%
    filter(variable %in% names(true_vals)) %>%
    mutate(rep = rep_id, true = true_vals[variable], bias = mean - true,
           rel_bias_pct = 100 * bias / true, covered = true >= q5 & true <= q95) %>%
    select(rep, variable, true, mean, sd, q5, q95, rhat, ess_bulk, ess_tail, bias, rel_bias_pct, covered)
}


param_recovery <- map_dfr(1:10, build_recovery_table, data_dir = DATA_DIR, fits_dir = FITS_DIR)
subj_recovery  <- map_dfr(1:10, recover_subject_level, data_dir = DATA_DIR, fits_dir = FITS_DIR)


pop_summaries <- map(1:10, function(rep_id) {
  sr <- recover_subject_level(rep_id, DATA_DIR, FITS_DIR)
  pd <- attr(sr, "pop_draws")
  tibble(
    rep = rep_id,
    t0_mean_pop    = mean(pd$t0),    t0_q05_pop    = quantile(pd$t0, 0.05),    t0_q95_pop    = quantile(pd$t0, 0.95),
    sv_mean_pop    = mean(pd$sv),    sv_q05_pop    = quantile(pd$sv, 0.05),    sv_q95_pop    = quantile(pd$sv, 0.95),
    alpha_mean_pop = mean(pd$alpha), alpha_q05_pop = quantile(pd$alpha, 0.05), alpha_q95_pop = quantile(pd$alpha, 0.95)
  )
})
pop_summaries <- bind_rows(pop_summaries)


alpha_pop <- map_dfr(1:10, function(rep_id) {
  true_pop_r <- readRDS(file.path(DATA_DIR, sprintf("rep%02d", rep_id), "true_population_parameters.rds"))
  pop_summaries %>%
    filter(rep == rep_id) %>%
    transmute(rep, variable = "alpha_mean",
              true = exp(true_pop_r$beta_alpha1),
              mean = alpha_mean_pop, q5 = alpha_q05_pop, q95 = alpha_q95_pop)
})


overall_true_t0 <- mean(subj_recovery$t0_true)

t0_pop <- subj_recovery %>% group_by(rep) %>%
  summarise(true = overall_true_t0) %>%
  left_join(pop_summaries %>% select(rep, mean = t0_mean_pop, q5 = t0_q05_pop, q95 = t0_q95_pop), by = "rep") %>%
  mutate(variable = "t0_mean") %>%
  select(rep, variable, true, mean, q5, q95)


sv_pop <- subj_recovery %>% group_by(rep) %>%
  summarise(true = mean(sv_true)) %>%
  left_join(pop_summaries %>% select(rep, mean = sv_mean_pop, q5 = sv_q05_pop, q95 = sv_q95_pop), by = "rep") %>%
  mutate(variable = "sv_mean") %>%
  select(rep, variable, true, mean, q5, q95)


true_pop_rep1 <- readRDS(file.path(DATA_DIR, "rep01", "true_population_parameters.rds"))

not_estimated <- tibble(
  variable = c("s_beta", "s_tau"),
  true     = c(true_pop_rep1$sw, true_pop_rep1$st0_raw),
  mean = NA_real_, q5 = NA_real_, q95 = NA_real_, covered = NA
)


recovery_by_rep <- bind_rows(
  param_recovery %>% mutate(bias = mean - true, rel_bias_pct = 100 * bias / true),
  t0_pop    %>% mutate(bias = mean - true, rel_bias_pct = 100 * bias / true),
  sv_pop    %>% mutate(bias = mean - true, rel_bias_pct = 100 * bias / true),
  alpha_pop %>% mutate(bias = mean - true, rel_bias_pct = 100 * bias / true)
)

recovery_summary <- recovery_by_rep %>%
  group_by(variable) %>%
  summarise(true = first(true), mean_estimate = mean(mean), mean_bias = mean(bias),
            mean_rel_bias_pct = mean(rel_bias_pct), coverage_90 = mean(covered), n_reps = n())

write.csv(recovery_by_rep,  file.path(FITS_DIR, "m2sv_recovery_table_by_rep.csv"), row.names = FALSE)
write.csv(recovery_summary, file.path(FITS_DIR, "m2sv_recovery_summary.csv"), row.names = FALSE)
write.csv(subj_recovery,    file.path(FITS_DIR, "m2sv_subject_level_recovery.csv"), row.names = FALSE)

print(recovery_summary, n = Inf)


plot_df <- param_recovery %>%
  select(rep, variable, true, mean, q5, q95) %>%
  bind_rows(t0_pop    %>% select(rep, variable, true, mean, q5, q95)) %>%
  bind_rows(sv_pop    %>% select(rep, variable, true, mean, q5, q95)) %>%
  bind_rows(alpha_pop %>% select(rep, variable, true, mean, q5, q95)) %>%
  bind_rows(
    tidyr::crossing(rep = 1:10, not_estimated) %>%
      select(rep, variable, true, mean, q5, q95)
  )


var_labels <- c(
  "beta_nu[1]"            = "nu intercept", # mean drift rate at baseline (condition/resp-type = 0)
  "beta_nu[2]"            = "nu condition slope",# mean effect of condition on drift rate
  "beta_nu[3]"            = "nu resp-type slope", # mean effect of previous-response type on drift rate
  "beta_nu[4]"            = "nu interaction slope",# mean condition × resp-type interaction on drift rate
  "beta_alpha2"           = "alpha resp-type effect (log)", # effect of resp-type on boundary separation
  "beta_alpha3"           = "alpha condition effect (log)",# effect of condition on boundary separation
  "beta_tau"              = "tau intercept (logit)", # mean of the logit-scale parameter that gets transformed into tau 
  "sigma_nu_c"            = "SD nu condition", # between-subject SD of the condition effect on nu
  "sigma_nu_r"            = "SD nu resp-type", # between-subject SD of the resp type effect on nu
  "sigma_nu_cr"           = "SD nu interaction", # between-subject SD of the interaction effect on nu
  "sigma_alpha_intercept" = "SD alpha (log)", # between-subject SD of the boundary separation intercept
  "sigma_tau_intercept"   = "SD tau (logit)", # between-subject SD of the logit scale tau intercept (between-subject variability)
  "w"                     = "beta (starting point)", # z = w * a --> boundary separation relativ to a (0.5) --> 0.52 very tiny bit more to upper bound
  # varies trial by trial so this is important but its where the accumulation starts
  "t0_mean"               = "tau (population mean in seconds)", # mean non-decision time in seconds 
  "sv_mean"               = "drift rate variability (mean)", # mean of the hierarchical trial-to-trial drift rate variability
  # within-trial parameter
  # estimated per subject and hierarchically pooled, its not the between-subject SD of drift
  "alpha_mean"            = "alpha (population mean)", # mean boundary separation
  "s_beta"                = "variability starting point (not estimated)", # 	Trial-to-trial starting-point variability
  "s_tau"                 = "variability NDT (not estimated)" # Trial-to-trial non-decision-time variability
)


fig1 <- ggplot(plot_df, aes(x = factor(rep), y = mean, ymin = q5, ymax = q95)) +
  geom_hline(aes(yintercept = true), color = "red", linetype = "dashed") +
  geom_pointrange(na.rm = TRUE) +
  facet_wrap(~ variable, scales = "free_y", labeller = labeller(variable = var_labels)) +
  labs(x = "Fitting Rep", y = "Posterior mean (with 90% CrI)") +
  theme_minimal()

ggsave(file.path(FITS_DIR, "recovery_M2_sv_hier_gloabalminrt.png"), fig1, width = 14, height = 10, dpi = 150)
print(fig1)