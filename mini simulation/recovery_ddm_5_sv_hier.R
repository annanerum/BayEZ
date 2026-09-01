library(dplyr)
library(purrr)
library(tidyr)
library(tibble)
library(ggplot2)

sd_subj_true <- list(
  sigma_nu_c = 0.50, sigma_nu_r = 0.30, sigma_nu_cr = 0.30,
  sigma_alpha_intercept = 0.15, sigma_tau_intercept = 0.30
)

# t0 and sv are both now subject-level
# need back-transformation before comparison to true values
# beta_tau, sigma_tau_intercept and beta_sv/sigma_sv are all not directly comparable on their
## raw parameter scale 
recover_subject_level <- function(rep_id, data_dir, fits_dir) {
  rep_dir <- file.path(data_dir, sprintf("rep%02d", rep_id))
  out_dir <- file.path(fits_dir, sprintf("rep%02d_fullscale", rep_id))
  fit <- readRDS(file.path(out_dir, "fit.rds"))
  
  min_rt_obs <- readRDS(file.path(rep_dir, "min_rt_observed.rds"))
  sim_data   <- read.csv(file.path(rep_dir, "sim_data_full_ddm.csv"))
  subj_par   <- readRDS(file.path(rep_dir, "true_subject_parameters.rds"))
  subjects_used <- sort(unique(sim_data$subj))
  min_rt_used   <- as.numeric(min_rt_obs[as.character(subjects_used)])
  
  tau_draws <- fit$draws(variables = "tau_intercept_raw", format = "matrix")
  t0_hat    <- sweep(plogis(tau_draws), 2, min_rt_used, `*`)
  
  sv_draws  <- fit$draws(variables = "sv", format = "matrix")   # already on scale, exp() applied in transformed parameters
  
  true_sv_single <- readRDS(file.path(rep_dir, "true_population_parameters.rds"))$sv
  
  tibble(
    rep = rep_id, subj = subjects_used,
    t0_mean = apply(t0_hat, 2, mean), t0_q05 = apply(t0_hat, 2, quantile, 0.05), t0_q95 = apply(t0_hat, 2, quantile, 0.95),
    sv_mean = apply(sv_draws, 2, mean), sv_q05 = apply(sv_draws, 2, quantile, 0.05), sv_q95 = apply(sv_draws, 2, quantile, 0.95)
  ) %>%
    left_join(subj_par %>% select(subj, t0_true = t0), by = "subj") %>%
    mutate(sv_true = true_sv_single) %>%
    mutate(
      t0_covered = t0_true >= t0_q05 & t0_true <= t0_q95,
      sv_covered = sv_true >= sv_q05 & sv_true <= sv_q95
    )
}

build_recovery_table <- function(rep_id, data_dir, fits_dir) {
  out_dir  <- file.path(fits_dir, sprintf("rep%02d_fullscale", rep_id))
  summ     <- read.csv(file.path(out_dir, "summary.csv"))
  true_pop <- readRDS(file.path(data_dir, sprintf("rep%02d", rep_id), "true_population_parameters.rds"))
  
  true_vals <- c(
    "beta_nu[1]" = true_pop$beta_nu1, "beta_nu[2]" = true_pop$nu_c,
    "beta_nu[3]" = true_pop$nu_r,     "beta_nu[4]" = true_pop$nu_cr,
    beta_alpha1 = true_pop$beta_alpha1, beta_alpha2 = true_pop$beta_alpha2,
    beta_alpha3 = true_pop$beta_alpha3, w = true_pop$w,
    sigma_nu_c = sd_subj_true$sigma_nu_c, sigma_nu_r = sd_subj_true$sigma_nu_r,
    sigma_nu_cr = sd_subj_true$sigma_nu_cr,
    sigma_alpha_intercept = sd_subj_true$sigma_alpha_intercept
    # beta_tau, sigma_tau_intercept, beta_sv, sigma_sv excluded --> t0/sv back-transform 
  )
  
  summ %>%
    filter(variable %in% names(true_vals)) %>%
    mutate(rep = rep_id, true = true_vals[variable], bias = mean - true,
           rel_bias_pct = 100 * bias / true, covered = true >= q5 & true <= q95) %>%
    select(rep, variable, true, mean, sd, q5, q95, rhat, ess_bulk, ess_tail, bias, rel_bias_pct, covered)
}


DATA_DIR <- "data_full7"
FITS_DIR <- "fits/M2 sv hier"

param_recovery <- map_dfr(1:10, build_recovery_table, data_dir = DATA_DIR, fits_dir = FITS_DIR)
subj_recovery  <- map_dfr(1:10, recover_subject_level, data_dir = DATA_DIR, fits_dir = FITS_DIR)

t0_pop <- subj_recovery %>% group_by(rep) %>%
  summarise(true = mean(t0_true), mean = mean(t0_mean), covered = mean(t0_covered)) %>%
  mutate(variable = "t0_mean") %>% select(rep, variable, true, mean, covered)

sv_pop <- subj_recovery %>% group_by(rep) %>%
  summarise(true = mean(sv_true), mean = mean(sv_mean), covered = mean(sv_covered),
            sd_recovered = sd(sv_mean)) %>%     
  mutate(variable = "sv_mean") %>% select(rep, variable, true, mean, covered, sd_recovered)

recovery_by_rep <- bind_rows(
  param_recovery %>% mutate(bias = mean - true, rel_bias_pct = 100 * bias / true),
  t0_pop %>% mutate(bias = mean - true, rel_bias_pct = 100 * bias / true),
  sv_pop %>% mutate(bias = mean - true, rel_bias_pct = 100 * bias / true)
)

recovery_summary <- recovery_by_rep %>%
  group_by(variable) %>%
  summarise(true = first(true), mean_estimate = mean(mean), mean_bias = mean(bias),
            mean_rel_bias_pct = mean(rel_bias_pct), coverage_90 = mean(covered), n_reps = n())

write.csv(recovery_by_rep,  file.path(FITS_DIR, "m2sv_recovery_table_by_rep.csv"), row.names = FALSE)
write.csv(recovery_summary, file.path(FITS_DIR, "m2sv_recovery_summary.csv"), row.names = FALSE)
write.csv(subj_recovery,    file.path(FITS_DIR, "m2sv_subject_level_recovery.csv"), row.names = FALSE)

print(recovery_summary, n = Inf)

## does the model invent spurious between-subject sv variabilit where none exists? (true sd_subj$sv is 0)
cat("\n sv between-subject SD (should be small) \n")
print(sv_pop %>% select(rep, sd_recovered))



## plots
plot_df <- param_recovery %>%
  select(rep, variable, true, mean, q5, q95) %>%
  bind_rows(t0_pop %>% rename(mean_val = mean) %>% mutate(q5 = NA, q95 = NA, mean = mean_val) %>% select(rep, variable, true, mean, q5, q95)) %>%
  bind_rows(sv_pop %>% rename(mean_val = mean) %>% mutate(q5 = NA, q95 = NA, mean = mean_val) %>% select(rep, variable, true, mean, q5, q95))

fig1 <- ggplot(plot_df, aes(x = factor(rep), y = mean, ymin = q5, ymax = q95)) +
  geom_hline(aes(yintercept = true), color = "red", linetype = "dashed") +
  geom_pointrange() +
  facet_wrap(~ variable, scales = "free_y") +
  labs(title = "M2 (hierarchical sv) Parameter-Recovery ueber 10 Reps",
       subtitle = "Gestrichelte rote Linie = wahrer Wert",
       x = "Rep", y = "Posterior mean (mit 90%-CrI)") +
  theme_minimal()

ggsave(file.path(FITS_DIR, "m2sv_recovery_by_rep.png"), fig1, width = 12, height = 8, dpi = 150)

print(fig1)



dark_theme <- theme_minimal(base_size = 14) +
  theme(
    plot.background  = element_rect(fill = "#111111", color = NA),
    panel.background = element_rect(fill = "#111111", color = NA),
    panel.grid.major = element_line(color = "#3a3a3a"),
    panel.grid.minor = element_blank(),
    text        = element_text(color = "white"),
    axis.text   = element_text(color = "white"),
    axis.title  = element_text(color = "white"),
    plot.title    = element_text(color = "white"),
    plot.subtitle = element_text(color = "#cccccc"),
    legend.background = element_rect(fill = "#111111"),
    legend.key        = element_rect(fill = "#111111"),
    legend.text  = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    strip.background = element_rect(fill = "#222222"),
    strip.text       = element_text(color = "white")
  )

COL_GREEN <- "#4daf7c"
COL_BLUE  <- "#5b9bd5"
COL_GOLD  <- "#e8a33d"
COL_RED   <- "#d6604d"

caterpillar_plot <- function(subj_mean, subj_q_low, subj_q_high, title,
                             ylab = "Value", color = COL_GREEN,
                             subtitle = "Sorted by posterior mean, gold = population mean") {
  pop_mean <- mean(subj_mean)
  df <- tibble(mean = subj_mean, q_low = subj_q_low, q_high = subj_q_high) %>%
    arrange(mean) %>%
    mutate(rank = row_number())
  
  ggplot(df, aes(x = rank, y = mean, ymin = q_low, ymax = q_high)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "white") +
    geom_hline(yintercept = pop_mean, color = COL_GOLD, linewidth = 1) +
    geom_pointrange(color = color, size = 0.3) +
    labs(title = title, subtitle = subtitle, x = "Participant rank", y = ylab) +
    dark_theme
}

cell_violin_plot <- function(cell_df, value_col, title, ylab,
                             subtitle = "Diamond = cell mean, each dot = one participant") {
  ggplot(cell_df, aes(x = cell, y = .data[[value_col]], color = cell, fill = cell)) +
    geom_violin(alpha = 0.25, linewidth = 0.8) +
    geom_jitter(width = 0.08, alpha = 0.6, size = 1.2) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3,
                 fill = "white", color = "white") +
    labs(title = title, subtitle = subtitle, x = NULL, y = ylab) +
    dark_theme +
    theme(legend.position = "none") +
    scale_x_discrete(labels = function(x) x)
}

forest_plot <- function(df, title, subtitle) {
  df$label <- factor(df$label, levels = rev(df$label))
  ggplot(df, aes(y = label, x = mean, color = family)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "white") +
    geom_errorbarh(aes(xmin = q2.5, xmax = q97.5), height = 0, linewidth = 0.6) +
    geom_errorbarh(aes(xmin = q25,  xmax = q75),  height = 0, linewidth = 2.2) +
    geom_point(size = 3) +
    labs(title = title, subtitle = subtitle, x = "Posterior mean", y = NULL, color = "Parameter") +
    dark_theme
}

make_ddm_diagnostic_plots <- function(fit_path, rep_dir, out_dir, model_label) {
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  fit <- readRDS(fit_path)
  sim_data <- read.csv(file.path(rep_dir, "sim_data_full_ddm.csv"))
  min_rt_obs <- readRDS(file.path(rep_dir, "min_rt_observed.rds"))
  subjects_used <- sort(unique(sim_data$subj))
  min_rt_used   <- as.numeric(min_rt_obs[as.character(subjects_used)])
  
  summ <- fit$summary()
  
  alpha_int_draws <- fit$draws(variables = "alpha_intercept", format = "matrix")
  alpha_nat_draws <- exp(alpha_int_draws)
  
  nu_c_draws  <- fit$draws(variables = "nu_c",  format = "matrix")
  nu_r_draws  <- fit$draws(variables = "nu_r",  format = "matrix")
  nu_cr_draws <- fit$draws(variables = "nu_cr", format = "matrix")
  
  tau_raw_draws <- fit$draws(variables = "tau_intercept_raw", format = "matrix")
  t0_draws <- sweep(plogis(tau_raw_draws), 2, min_rt_used, `*`)
  
  beta_nu1_draws  <- fit$draws(variables = "beta_nu[1]", format = "matrix")[, 1]
  
  ## 1. Caterpillar plots
  fig_alpha <- caterpillar_plot(
    apply(alpha_nat_draws, 2, mean), apply(alpha_nat_draws, 2, quantile, 0.05), apply(alpha_nat_draws, 2, quantile, 0.95),
    title = "Individual intercepts, alpha (natural scale)", ylab = "Value", color = COL_GREEN)
  
  fig_nu_c <- caterpillar_plot(
    apply(nu_c_draws, 2, mean), apply(nu_c_draws, 2, quantile, 0.05), apply(nu_c_draws, 2, quantile, 0.95),
    title = "Condition slope on nu", ylab = "Value", color = COL_BLUE)
  
  fig_nu_cr <- caterpillar_plot(
    apply(nu_cr_draws, 2, mean), apply(nu_cr_draws, 2, quantile, 0.05), apply(nu_cr_draws, 2, quantile, 0.95),
    title = "Interaction slope on nu", ylab = "Value", color = COL_BLUE)
  
  fig_nu_r <- caterpillar_plot(
    apply(nu_r_draws, 2, mean), apply(nu_r_draws, 2, quantile, 0.05), apply(nu_r_draws, 2, quantile, 0.95),
    title = "Resp_type slope on nu", ylab = "Value", color = COL_BLUE)
  
  fig_tau <- caterpillar_plot(
    apply(t0_draws, 2, mean), apply(t0_draws, 2, quantile, 0.05), apply(t0_draws, 2, quantile, 0.95),
    title = "Individual non-decision times, tau (seconds)", ylab = "Value", color = COL_GREEN)
  
  ggsave(file.path(out_dir, "ddm_caterpillar_alpha_intercept.png"), fig_alpha, width = 10, height = 6, dpi = 150, bg = "#111111")
  ggsave(file.path(out_dir, "ddm_caterpillar_nu_c.png"),  fig_nu_c,  width = 10, height = 6, dpi = 150, bg = "#111111")
  ggsave(file.path(out_dir, "ddm_caterpillar_nu_cr.png"), fig_nu_cr, width = 10, height = 6, dpi = 150, bg = "#111111")
  ggsave(file.path(out_dir, "ddm_caterpillar_nu_r.png"),  fig_nu_r,  width = 10, height = 6, dpi = 150, bg = "#111111")
  ggsave(file.path(out_dir, "ddm_caterpillar_tau.png"),   fig_tau,   width = 10, height = 6, dpi = 150, bg = "#111111")
  
  ## 2. Cell-level violins
  beta_alpha2_mean <- summ$mean[summ$variable == "beta_alpha2"]
  beta_alpha3_mean <- summ$mean[summ$variable == "beta_alpha3"]
  beta_nu1_mean  <- mean(beta_nu1_draws)
  alpha_int_mean <- apply(alpha_int_draws, 2, mean)
  nu_c_mean  <- apply(nu_c_draws,  2, mean)
  nu_r_mean  <- apply(nu_r_draws,  2, mean)
  nu_cr_mean <- apply(nu_cr_draws, 2, mean)
  
  cells <- tibble(
    condition = c(1, 1, -1, -1),
    resp_type = c(1, -1, 1, -1),
    cell = c("cong/post-corr", "cong/post-err", "incong/post-corr", "incong/post-err")
  )
  
  person_cell <- map_dfr(1:nrow(cells), function(i) {
    cond <- cells$condition[i]; rt <- cells$resp_type[i]; cell_label <- cells$cell[i]
    tibble(
      subj = subjects_used,
      cell = cell_label,
      alpha = exp(alpha_int_mean + beta_alpha2_mean * rt + beta_alpha3_mean * cond),
      nu    = beta_nu1_mean + nu_c_mean * cond + nu_r_mean * rt + nu_cr_mean * cond * rt
    )
  })
  person_cell$cell <- factor(person_cell$cell, levels = cells$cell)
  
  fig_cell_alpha <- cell_violin_plot(person_cell, "alpha", paste0("Cell-level, alpha posterior means  ", model_label), "alpha")
  fig_cell_nu    <- cell_violin_plot(person_cell, "nu",    paste0("Cell-level, nu posterior means  ", model_label), "nu")
  
  ggsave(file.path(out_dir, "ddm_cells_alpha.png"), fig_cell_alpha, width = 10, height = 8, dpi = 150, bg = "#111111")
  ggsave(file.path(out_dir, "ddm_cells_nu.png"),    fig_cell_nu,    width = 10, height = 8, dpi = 150, bg = "#111111")
  
  ## 3. Forest plots
  q_from_draws <- function(var) {
    d <- fit$draws(variables = var, format = "matrix")[, 1]
    c(mean = mean(d), q2.5 = quantile(d, .025), q25 = quantile(d, .25),
      q75 = quantile(d, .75), q97.5 = quantile(d, .975))
  }
  
  pop_vars <- tribble(
    ~label,                              ~var,          ~family,
    "mu_nu (intercept)",                 "beta_nu[1]",  "nu",
    "b_nu_condition",                    "beta_nu[2]",  "nu",
    "mu_tau (logit-of-min_rt scale)",    "beta_tau",    "tau",
    "b_nu_interaction",                  "beta_nu[4]",  "nu",
    "mu_alpha (log-scale intercept)",    "beta_alpha1", "alpha",
    "b_nu_resp_type",                    "beta_nu[3]",  "nu",
    "b_alpha_condition (log)",           "beta_alpha3", "alpha",
    "b_alpha_resp_type (log)",           "beta_alpha2", "alpha"
  )
  pop_df <- pop_vars %>%
    mutate(qs = map(var, q_from_draws)) %>%
    unnest_wider(qs)
  
  fig_pop <- forest_plot(pop_df,
                         title = paste0("Population effects  ", model_label),
                         subtitle = "Thick bar = 50% CI | thin bar = 95% CI. note: nu/alpha/tau on different scales") +
    scale_color_manual(values = c(alpha = COL_RED, nu = COL_BLUE, tau = COL_GOLD))
  
  sd_vars <- tribble(
    ~label,                    ~var,                    ~family,
    "sigma_nu_c (slope)",      "sigma_nu_c",            "slope",
    "sigma_nu_r (slope)",      "sigma_nu_r",            "slope",
    "sigma_tau (intercept)",   "sigma_tau_intercept",   "intercept",
    "sigma_nu_cr (slope)",     "sigma_nu_cr",           "slope",
    "sigma_alpha (intercept)", "sigma_alpha_intercept", "intercept"
  )
  sd_df <- sd_vars %>%
    mutate(qs = map(var, q_from_draws)) %>%
    unnest_wider(qs)
  
  fig_sd <- forest_plot(sd_df,
                        title = paste0("Between-person SDs  ", model_label),
                        subtitle = "Thick bar = 50% CI, thin bar = 95% CI") +
    scale_color_manual(values = c(intercept = "#3fae5c", slope = COL_BLUE))
  
  ggsave(file.path(out_dir, "ddm_forest_group.png"), fig_pop, width = 12, height = 8, dpi = 150, bg = "#111111")
  ggsave(file.path(out_dir, "ddm_forest_sds.png"),   fig_sd,  width = 12, height = 6, dpi = 150, bg = "#111111")
  
  ## csv export
  ddm_estimates <- bind_rows(
    pop_df %>% mutate(group = "population"),
    sd_df  %>% mutate(group = "between_person_sd")
  ) %>% select(group, label, var, family, mean, q2.5, q25, q75, q97.5)
  write.csv(ddm_estimates, file.path(out_dir, "ddm_estimates.csv"), row.names = FALSE)
  
  write.csv(person_cell, file.path(out_dir, "ddm_person_cell_estimates.csv"), row.names = FALSE)
  
  person_tau <- tibble(
    subj = subjects_used,
    tau_mean = apply(t0_draws, 2, mean),
    tau_q05  = apply(t0_draws, 2, quantile, 0.05),
    tau_q95  = apply(t0_draws, 2, quantile, 0.95)
  )
  write.csv(person_tau, file.path(out_dir, "ddm_person_tau_estimates.csv"), row.names = FALSE)
  
  invisible(list(figs = list(alpha = fig_alpha, nu_c = fig_nu_c, nu_cr = fig_nu_cr, nu_r = fig_nu_r,
                             tau = fig_tau, cell_alpha = fig_cell_alpha, cell_nu = fig_cell_nu,
                             pop = fig_pop, sd = fig_sd),
                 tables = list(estimates = ddm_estimates, person_cell = person_cell, person_tau = person_tau)))
}


make_ddm_diagnostic_plots(
  fit_path    = file.path(FITS_DIR, "rep01_fullscale", "fit.rds"),
  rep_dir     = file.path(DATA_DIR, "rep01"),
  out_dir     = file.path(FITS_DIR, "rep01_fullscale", "diagnostics"),
  model_label = "M2 sv_hier, rep01"
)