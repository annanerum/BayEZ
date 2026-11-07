
library(dplyr)
library(tidyr)
library(ggplot2)

PAPER_DIR   <- "/Users/annegiacobello/Desktop/EZ Paper "    
PROJECT_DIR <- file.path(PAPER_DIR, "sv hier 2 cond")
DATA_DIR    <- file.path(PROJECT_DIR, "data_sim_svhier")
PLOT_DIR    <- file.path(PROJECT_DIR, "plots_model_comparison")
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

reps <- 1:10   



read_summary <- function(file) {
  function(dir, r) {
    f <- file.path(dir, sprintf("rep%02d", r), file)
    if (!file.exists(f)) return(NULL)
    s <- if (grepl("\\.rds$", f)) readRDS(f) else read.csv(f, check.names = FALSE)
    s <- as.data.frame(s)
    if (!"q5"  %in% names(s) && "5%"  %in% names(s)) s$q5  <- s[["5%"]]
    if (!"q95" %in% names(s) && "95%" %in% names(s)) s$q95 <- s[["95%"]]
    s[, c("variable", "mean", "q5", "q95")]
  }
}




MODELS <- list(

  EZ = list(
    label = "EZ", colour = "#E69F00",
    dir = file.path(PROJECT_DIR, "fits_ez_2conditions"),
    reader = read_summary("summary_ez.csv"),
    map = c(nu_intercept = "beta_nu_intercept", nu_c = "beta_nu_c",
            alpha = "beta_alpha", tau = "beta_tau",
            sd_nu_intercept = "sigma_nu_intercept", sd_nu_c = "sigma_nu_c",
            sd_alpha = "sigma_alpha_intercept", sd_tau = "sigma_tau_intercept"),
    transform = c(alpha = "exp"),
    exclude = integer(0)
  ),

  EZ_bias = list(
    label = "EZ + bias", colour = "#56B4E9",
    dir = file.path(PROJECT_DIR, "fits_ez_bias_2conditions"),
    reader = read_summary("summary_ez_bias.csv"),
    map = c(nu_intercept = "beta_nu_intercept", nu_c = "beta_nu_c",
            alpha = "beta_alpha", tau = "beta_tau", w = "w_pop_median",
            sd_nu_intercept = "sigma_nu_intercept", sd_nu_c = "sigma_nu_c",
            sd_alpha = "sigma_alpha_intercept", sd_tau = "sigma_tau_intercept",
            sd_w = "sigma_w"),
    transform = c(alpha = "exp"),
    exclude = integer(0)
  ),

  
  DDM5 = list(
    label = "DDM5", colour = "#009E73",
    dir = file.path(PROJECT_DIR, "fits_ddm5_2conditions"),
    reader = read_summary("summary_ddm5.csv"),
    map = c(nu_intercept = "beta_nu_intercept", nu_c = "beta_nu_c",
            alpha = "beta_alpha", tau = "beta_tau",
            w = "w_pop_median", sv = "sv_pop_median",
            sd_nu_intercept = "sigma_nu_intercept", sd_nu_c = "sigma_nu_c",
            sd_alpha = "sigma_alpha_intercept", sd_tau = "sigma_tau_intercept",
            sd_w = "sigma_w", sd_sv = "sigma_sv"),
    transform = c(alpha = "exp"),     # beta_alpha is log --> exp = population median
    exclude = integer(0)
  ),

  
  DDM6 = list(
    label = "DDM6", colour = "#CC79A7",
    dir = file.path(PROJECT_DIR, "fits_ddm6_st0"),
    reader = read_summary("pop_summary.rds"),
    map = c(nu_intercept = "beta_nu_intercept", nu_c = "beta_nu_c",
            alpha = "alpha_pop_median", tau = "beta_tau",
            w = "w_pop_median", sv = "sv_pop_median", st0 = "st0",
            sd_nu_intercept = "sigma_nu_intercept", sd_nu_c = "sigma_nu_c",
            sd_alpha = "sigma_alpha_intercept", sd_tau = "sigma_tau_intercept",
            sd_w = "sigma_w", sd_sv = "sigma_sv"),
    transform = c(),
    exclude = 8:10     # rep 8: 28% divergences (refit planned), rep 9: borderline, rep 10: running
  ),

  
  DDM7 = list(
    label = "DDM7", colour = "#0072B2",
    dir = file.path(PROJECT_DIR, "fits_ddm7_pop"),
    reader = read_summary("pop_summary.rds"),
    map = c(nu_intercept = "beta_nu_intercept", nu_c = "beta_nu_c",
            alpha = "alpha_pop_median", tau = "beta_tau",
            w = "w_pop_median", sv = "sv_pop_median", sw = "sw", st0 = "st0",
            sd_nu_intercept = "sigma_nu_intercept", sd_nu_c = "sigma_nu_c",
            sd_alpha = "sigma_alpha_intercept", sd_tau = "sigma_tau_intercept",
            sd_w = "sigma_w", sd_sv = "sigma_sv"),
    transform = c(),
    exclude = integer(0)
  )
)


show_model_variables <- function(key, r = 1) {
  m <- MODELS[[key]]
  s <- m$reader(m$dir, r)
  if (is.null(s)) { cat("no summary found for", key, "rep", r, "in", m$dir, "\n"); return(invisible(NULL)) }
  v <- unique(s$variable)
  print(v[!grepl("\\[", v)])
}



read_model <- function(key) {
  m <- MODELS[[key]]
  if (!dir.exists(m$dir)) { message(m$label, ": folder not found --> skipped (", m$dir, ")"); return(NULL) }

  out <- bind_rows(lapply(setdiff(reps, m$exclude), function(r) {
    s <- m$reader(m$dir, r)
    if (is.null(s)) return(NULL)
    s$rep <- r
    s
  }))
  if (nrow(out) == 0) { message(m$label, ": no summary files --> skipped"); return(NULL) }

  not_found <- setdiff(unname(m$map), unique(out$variable))
  if (length(not_found) > 0)
    warning(m$label, ": not in the summaries: ", paste(not_found, collapse = ", "),
            " --> check MODELS$", key, "$map (show_model_variables(\"", key, "\"))")

  out <- out %>%
    inner_join(tibble(param = names(m$map), variable = unname(m$map)), by = "variable") %>%
    as.data.frame()

  for (p in names(m$transform)) {
    f   <- match.fun(m$transform[[p]])
    idx <- out$param == p
    for (cc in c("mean", "q5", "q95")) out[[cc]][idx] <- f(out[[cc]][idx])
  }
  out$model <- m$label
  cat(sprintf("%-10s reps: %s\n", m$label, paste(sort(unique(out$rep)), collapse = " ")))
  out
}

est <- bind_rows(lapply(names(MODELS), read_model))
if (nrow(est) == 0) stop("no model summaries found")




truth_one <- function(r) {
  d  <- file.path(DATA_DIR, sprintf("rep%02d", r))
  tp <- readRDS(file.path(d, "true_population_parameters.rds"))
  ts <- readRDS(file.path(d, "true_sd_parameters.rds"))
  sp <- readRDS(file.path(d, "true_subject_parameters.rds"))
  sw_frac <- tp$sw / (0.99 * 2 * pmin(sp$w, 1 - sp$w))
  tibble(
    rep = r,
    param = c("nu_intercept", "nu_c", "alpha", "tau", "w", "sv", "sw", "st0",
              "sd_nu_intercept", "sd_nu_c", "sd_alpha", "sd_tau", "sd_w", "sd_sv",
              "sd_sw", "sd_st0"),
    truth = c(tp$beta_nu1, tp$nu_c, exp(tp$beta_alpha1), tp$beta_tau, tp$w, tp$sv,
              tp$sw, tp$st0_raw * 2 * tp$beta_tau,
              ts$sigma_nu_intercept, ts$sigma_nu_c, ts$sigma_alpha_intercept,
              ts$sigma_tau_intercept, ts$sigma_w, ts$sigma_sv,
              sd(qlogis(sw_frac)), sd(log(sp$st0)))
  )
}
truth <- bind_rows(lapply(reps, truth_one))


tp1 <- readRDS(file.path(DATA_DIR, "rep01", "true_population_parameters.rds"))
tau_mean_ndt <- tp1$beta_tau * (1 + tp1$st0_raw)




param_labels <- c(
  nu_intercept    = "nu intercept",
  nu_c            = "nu condition slope",
  alpha           = "alpha",
  tau             = "tau (population mean in seconds)",
  w               = "starting point",
  sv              = "sv drift variability",
  sw              = "sw starting point variability",
  st0             = "st0 (NDT variability)",
  sd_nu_intercept = "SD nu intercept",
  sd_nu_c         = "SD nu condition",
  sd_alpha        = "SD alpha (log)",
  sd_tau          = "SD tau (seconds)",
  sd_w            = "SD w (logit)",
  sd_sv           = "SD sv (log)",
  sd_sw           = "SD sw",
  sd_st0          = "SD st0"
)


model_order <- intersect(sapply(MODELS, `[[`, "label"), unique(est$model))
n_mod       <- length(model_order)
offsets     <- setNames((seq_len(n_mod) - (n_mod + 1) / 2) * (0.8 / n_mod), model_order)
colours     <- setNames(sapply(MODELS, `[[`, "colour"), sapply(MODELS, `[[`, "label"))[model_order]

plot_df <- est %>%
  left_join(truth, by = c("rep", "param")) %>%
  filter(param %in% names(param_labels)) %>%
  mutate(model = factor(model, levels = model_order),
         x     = rep + offsets[as.character(model)],
         label = factor(param_labels[param], levels = param_labels)) %>%
  droplevels()

panels <- levels(plot_df$label)

truth_plot <- truth %>%
  filter(param %in% unique(plot_df$param)) %>%
  mutate(label = factor(param_labels[param], levels = panels))
truth_type  <- truth_plot %>% group_by(label) %>%
  summarise(constant = n_distinct(round(truth, 6)) == 1, .groups = "drop")
truth_const <- truth_plot %>% semi_join(filter(truth_type, constant), by = "label") %>%
  distinct(label, truth)
truth_var   <- truth_plot %>% semi_join(filter(truth_type, !constant), by = "label")

ndt_ref <- tibble(label = factor(param_labels["tau"], levels = panels), y = tau_mean_ndt) %>%
  filter(!is.na(label))

p_pop <- ggplot(plot_df, aes(x = x, y = mean, colour = model)) +
  geom_hline(data = truth_const, aes(yintercept = truth),
             colour = "red", linetype = "dashed") +
  geom_segment(data = truth_var,
               aes(x = rep - 0.45, xend = rep + 0.45, y = truth, yend = truth),
               colour = "red", linetype = "dashed", inherit.aes = FALSE) +
  geom_hline(data = ndt_ref, aes(yintercept = y),
             colour = "grey40", linetype = "dotted") +
  geom_linerange(aes(ymin = q5, ymax = q95), linewidth = 0.5) +
  geom_point(size = 1.8) +
  facet_wrap(~ label, scales = "free_y", ncol = 4) +
  scale_x_continuous(breaks = reps) +
  scale_colour_manual(values = colours, name = NULL) +
  labs(x = "Fitting Rep", y = "Posterior mean (with 90% CrI)",
       caption = "Red dashed: true value. Grey dotted (tau): true mean NDT = t0 * (1 + st0_raw).") +
  theme_minimal() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 8),
        panel.grid.minor.x = element_blank())

print(p_pop)
n_rows <- ceiling(length(panels) / 4)
ggsave(file.path(PLOT_DIR, "population_recovery_all_models.png"), p_pop,
       width = 14, height = 2.6 * n_rows + 1, dpi = 300, bg = "white")




rec_table <- plot_df %>%
  group_by(label, model) %>%
  summarise(n_reps     = n(),
            mean_truth = mean(truth),
            mean_est   = mean(mean),
            bias       = mean(mean - truth),
            rel_bias   = mean((mean - truth) / truth),
            coverage90 = mean(truth >= q5 & truth <= q95),
            .groups = "drop") %>%
  arrange(label, model)

print(rec_table, n = Inf)
write.csv(rec_table, file.path(PLOT_DIR, "population_recovery_all_models.csv"), row.names = FALSE)
