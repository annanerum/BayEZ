# Fit the DDM5 (v, a, t0, w, sv free; sw = st0 = 0) to the Krakow Grabowska flanker data
# Same data, filtering and 2x2 effect structure as the grabowska_ez_model_new.stan

# 4 cells, contrast coded +-1
# k=1: cond= 1, resp= 1   k=2: cond= 1, resp=-1
# k=3: cond=-1, resp= 1   k=4: cond=-1, resp=-1
# cond: congruent = 1, incongruent = -1 | resp: post-correct = 1, post-error = -1

library(cmdstanr)
library(dplyr)
library(posterior)

PROJECT_DIR <- "/Users/annegiacobello/Desktop/EZ Paper "
DATA_FILE   <- file.path(PROJECT_DIR, "krakow_data_standardized.csv")
STAN_FILE   <- file.path(PROJECT_DIR, "ddm5_krakow.stan")
FIT_DIR     <- file.path(PROJECT_DIR, "fits_ddm5_krakow")
EZ_FIT_RDS  <- file.path(PROJECT_DIR, "output_krakow_new", "krakow_fit.rds")   

TEST_MODE     <- FALSE
TEST_N_PERSON <- 20

dir.create(FIT_DIR, showWarnings = FALSE, recursive = TRUE)

mod <- cmdstan_model(STAN_FILE, cpp_options = list(stan_threads = TRUE))


#### data 

dat <- read.csv(DATA_FILE)
dat <- dat[dat$is_in_sequence          == "True"  &
           dat$rt_greater_than_1       == "False" &
           dat$log_rt_exceed_threshold == "False", ]

if (TEST_MODE) {
  set.seed(1)
  keep <- sample(unique(dat$participant_index), TEST_N_PERSON)
  dat  <- dat[dat$participant_index %in% keep, ]
}

subj_ids <- sort(unique(dat$participant_index))       
dat$pid  <- match(dat$participant_index, subj_ids)
dat$cond <- as.numeric(dat$condition)                  # congruent = 1, incongruent = -1
dat$resp <- as.numeric(dat$pre_acc)                    # post-correct = 1, post-error = -1
dat$acc_pm <- ifelse(dat$acc == 1, 1L, -1L)            # +1 correct (upper), -1 error (lower)
dat$cell <- with(dat, case_when(
  cond ==  1 & resp ==  1 ~ 1L,
  cond ==  1 & resp == -1 ~ 2L,
  cond == -1 & resp ==  1 ~ 3L,
  cond == -1 & resp == -1 ~ 4L
))
stopifnot(!anyNA(dat$cell), !anyNA(dat$rt))

I      <- length(subj_ids)
min_rt <- as.numeric(tapply(dat$rt, dat$pid, min))     # t0 must stay below this
cat(sprintf("Participants: %d | Trials: %d\n", I, nrow(dat)))
print(table(cell = dat$cell, acc = dat$acc_pm))

stan_input <- list(
  N = nrow(dat), I = I, pid = dat$pid, cell = dat$cell,
  rt = dat$rt, acc = dat$acc_pm, min_rt = min_rt,
  grainsize = 1
)


#### priors


priors <- list(
  # drift
  prior_beta_nu_intercept     = c(0,    2),     # App C: N(0, 2)
  prior_beta_nu_c             = c(0,    2),     # App C: N(0, 2)
  prior_beta_nu_r             = c(0,    0.5),   # App C: N(0, 0.5)
  prior_beta_nu_cr            = c(0,    0.5),   # App C: N(0, 0.5)
  # boundary (natural scale)
  prior_beta_alpha            = c(1,    1),     # App C: N(1, 1)[0, 3] (bounds in Stan)
  prior_beta_alpha_r          = c(0,    0.2),   # App C: N(0, 0.2)
  prior_beta_alpha_c          = c(0,    1),     # App C: N(0, 1)
  # non-decision time
  prior_beta_tau              = c(0.1,  0.2),   # App C: N(0.1, 0.2)[0, 0.4] (bounds in Stan)
  # not in Grabowska
  prior_beta_w                = c(0,    0.5),   # logit: w around 0.5
  prior_beta_sv               = c(log(0.5), 1), # log: median sv 0.5, 95% ~ [0.07, 3.5]

  # between-person SDs
  prior_sigma_nu_intercept    = c(1,    1),     # gamma, App C: gamma(1, 1)
  prior_sigma_nu_c            = c(1,    1),     # gamma, App C: gamma(1, 1)
  prior_sigma_nu_r            = c(0,    0.5),   # half-normal, not in Grabowska (fixed there); as EZ
  prior_sigma_nu_cr           = c(0,    0.5),   # half-normal, not in Grabowska (fixed there); as EZ
  prior_sigma_alpha_intercept = c(1,    1),     # gamma, App C: gamma(1, 1)
  prior_sigma_alpha_r         = c(0,    0.5),   # half-normal, not in Grabowska (fixed there); as EZ
  prior_sigma_tau_intercept   = c(0.3,  1),     # gamma, App C: gamma(0.3, 1)
  prior_sigma_w               = c(0,    0.5),   # half-normal, not in Grabowska (logit scale)
  prior_sigma_sv              = c(0,    0.5)    # half-normal, not in Grabowska (log scale)
)

stan_data <- c(stan_input, priors)


#### fit

make_inits <- function(I, n_chains, min_rt) {
  t0_init <- 0.7 * min_rt                              # t0 starts below each person's fastest RT
  lapply(seq_len(n_chains), function(ch) list(
    t0                    = t0_init,
    beta_nu_intercept     = 2.0   + rnorm(1, 0, 0.05),
    beta_nu_c             = 1.5   + rnorm(1, 0, 0.05),
    beta_nu_r             = 0     + rnorm(1, 0, 0.05),
    beta_nu_cr            = 0     + rnorm(1, 0, 0.05),
    beta_alpha            = 1.0   + rnorm(1, 0, 0.05),       
    beta_alpha_r          = 0     + rnorm(1, 0, 0.02),
    beta_alpha_c          = 0     + rnorm(1, 0, 0.02),
    beta_tau              = min(max(mean(t0_init), 0.06), 0.39),
    beta_w                = rnorm(1, 0, 0.05),
    beta_sv               = log(0.5) + rnorm(1, 0, 0.05),
    sigma_nu_intercept    = 0.3,
    sigma_nu_c            = 0.3,
    sigma_nu_r            = 0.1,
    sigma_nu_cr           = 0.1,
    sigma_alpha_intercept = 0.1,
    sigma_alpha_r         = 0.05,
    sigma_tau_intercept   = 0.03,
    sigma_w               = 0.2,
    sigma_sv              = 0.1,
    z_nu_intercept = rnorm(I, 0, 0.1), z_nu_c  = rnorm(I, 0, 0.1),
    z_nu_r         = rnorm(I, 0, 0.1), z_nu_cr = rnorm(I, 0, 0.1),
    z_alpha_intercept = rnorm(I, 0, 0.1), z_alpha_r = rnorm(I, 0, 0.1),
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

pop_pars <- c("beta_nu_intercept", "beta_nu_c", "beta_nu_r", "beta_nu_cr",
              "beta_alpha", "beta_alpha_r", "beta_alpha_c", "beta_tau", "beta_w", "beta_sv",
              "sigma_nu_intercept", "sigma_nu_c", "sigma_nu_r", "sigma_nu_cr",
              "sigma_alpha_intercept", "sigma_alpha_r", "sigma_tau_intercept",
              "sigma_w", "sigma_sv",
              "w_pop_median", "sv_pop_median", "nu_pop_cell", "alpha_pop_cell")
person_pars <- c("nu_intercept", "nu_c", "nu_r", "nu_cr", "alpha_intercept", "alpha_r",
                 "t0", "w", "sv", "nu", "alpha")


csv_dir <- file.path(FIT_DIR, "csv")
dir.create(csv_dir, showWarnings = FALSE)

set.seed(2026)
t_start <- Sys.time()
fit <- do.call(mod$sample, c(list(
  data              = stan_data,
  init              = make_inits(I, sampler_args$chains, min_rt),
  seed              = 2026,
  threads_per_chain = THREADS_PER_CHAIN,
  adapt_delta       = 0.9,
  max_treedepth     = 10,
  refresh           = 100,
  output_dir        = csv_dir,
  output_basename   = if (TEST_MODE) "ddm5_krakow_test" else "ddm5_krakow"
), sampler_args))
cat("fit time:", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), "min\n")

diag <- fit$diagnostic_summary(quiet = TRUE)
summ <- fit$summary(variables = c(pop_pars, person_pars))
diag_row <- data.frame(
  divergences   = sum(diag$num_divergent),
  treedepth_hit = sum(diag$num_max_treedepth),
  max_rhat      = max(summ$rhat, na.rm = TRUE),
  min_ess_bulk  = min(summ$ess_bulk, na.rm = TRUE),
  min_ess_tail  = min(summ$ess_tail, na.rm = TRUE)
)
print(diag_row)

tag <- if (TEST_MODE) "_test" else ""
fit$save_object(file.path(FIT_DIR, paste0("fit_ddm5_krakow", tag, ".rds")))
write.csv(summ,     file.path(FIT_DIR, paste0("summary_ddm5_krakow", tag, ".csv")),     row.names = FALSE)
write.csv(diag_row, file.path(FIT_DIR, paste0("diagnostics_ddm5_krakow", tag, ".csv")), row.names = FALSE)
saveRDS(subj_ids, file.path(FIT_DIR, paste0("subj_ids", tag, ".rds")))
saveRDS(min_rt,   file.path(FIT_DIR, paste0("min_rt", tag, ".rds")))   # for the t0 plot
saveRDS(priors,   file.path(FIT_DIR, "priors_ddm5_krakow.rds"))

cat("\n population estimates\n")
summ %>% filter(variable %in% pop_pars | grepl("_pop_cell\\[", variable)) %>%
  select(variable, mean, q5, q95, rhat, ess_bulk) %>%
  print(n = Inf)


#### comparison table: EZ vs DDM5 vs Grabowska et al.
grabowska <- data.frame(
  row       = c("nu_intercept", "nu_cond", "nu_resp", "nu_cr",
                "alpha_intercept", "alpha_cond", "alpha_resp",
                "sd_nu_intercept", "sd_nu_cond", "sd_alpha_intercept"),
  grab_mean = c(4.03, 3.07,  0.18, -0.21, 1.48, 0.45, 0.09, 1.13, 1.02, 0.26),
  grab_lo95 = c(3.87, 2.92,  0.12, -0.27, 1.44, 0.42, 0.08, 1.02, 0.92, 0.23),
  grab_hi95 = c(4.19, 3.22,  0.24, -0.15, 1.52, 0.48, 0.10, 1.25, 1.12, 0.29)
)


rows <- data.frame(
  row = c("nu_intercept", "nu_cond", "nu_resp", "nu_cr",
          "alpha_intercept", "alpha_cond", "alpha_resp", "tau",
          "sd_nu_intercept", "sd_nu_cond", "sd_nu_resp", "sd_nu_cr",
          "sd_alpha_intercept", "sd_alpha_resp", "sd_tau", "w", "sv"),
  label = c("Drift: intercept", "Drift: condition", "Drift: response type",
            "Drift: condition x response type",
            "Boundary: intercept", "Boundary: condition", "Boundary: response type",
            "Non-decision time: mean (s)",
            "SD drift: intercept", "SD drift: condition", "SD drift: response type",
            "SD drift: interaction",
            "SD boundary: intercept", "SD boundary: response type", "SD non-decision time",
            "Starting point w", "Drift variability sv"),
  ez = c("beta_nu[1]", "beta_nu[2]", "beta_nu[3]", "beta_nu[4]",
         "beta_alpha[1]", "beta_alpha[3]", "beta_alpha[2]", "beta_tau",
         "sigma_nu", "sigma_nu_c", "sigma_nu_r", "sigma_nu_cr",
         "sigma_alpha_intercept|sigma_alpha", "sigma_alpha_r", "sigma_tau_intercept|sigma_tau",
         "", ""),
  ddm5 = c("beta_nu_intercept", "beta_nu_c", "beta_nu_r", "beta_nu_cr",
           "beta_alpha", "beta_alpha_c", "beta_alpha_r", "beta_tau",
           "sigma_nu_intercept", "sigma_nu_c", "sigma_nu_r", "sigma_nu_cr",
           "sigma_alpha_intercept", "sigma_alpha_r", "sigma_tau_intercept",
           "w_pop_median", "sv_pop_median"),
  note = c("", "", "", "",
           "average over the 4 cells", "EZ/DDM5: fixed, Grabowska: random",
           "EZ/DDM5: random, Grabowska: fixed",
           "Grabowska: not reported",
           "EZ: cell-level SD (no random drift intercept)", "", "Grabowska: fixed effect",
           "Grabowska: fixed effect",
           "EZ: cell-level SD if no random intercept",
           "Grabowska: fixed effect", "EZ: cell-level SD if no random intercept",
           "median person; EZ and Grabowska: not estimated",
           "median person; DDM5 only")
)

# DDM5: draws of the comparison parameters
d5_draws <- fit$draws(unique(rows$ddm5), format = "draws_df")
d5_tab <- summarise_draws(subset_draws(d5_draws, variable = rows$ddm5), "mean",
                          ~ quantile2(.x, probs = c(0.025, 0.975)))
d5_tab <- data.frame(ddm5 = d5_tab$variable, ddm5_mean = d5_tab$mean,
                     ddm5_lo95 = d5_tab$q2.5, ddm5_hi95 = d5_tab$q97.5)

# EZ: rstan summary of the EZ fit (if available)
ez_s <- if (file.exists(EZ_FIT_RDS)) rstan::summary(readRDS(EZ_FIT_RDS))$summary else NULL
if (is.null(ez_s)) cat("\nEZ fit not found at", EZ_FIT_RDS, "--> EZ columns stay NA\n")
pick_ez <- function(cands) {
  if (is.null(ez_s) || cands == "") return(NA_character_)
  hit <- intersect(strsplit(cands, "\\|")[[1]], rownames(ez_s))
  if (length(hit)) hit[1] else NA_character_
}
rows$ez_par <- vapply(rows$ez, pick_ez, character(1))
ez_val <- function(col) vapply(rows$ez_par, function(p)
  if (is.na(p)) NA_real_ else unname(ez_s[p, col]), numeric(1))

comp3 <- rows %>%
  mutate(ez_mean = ez_val("mean"), ez_lo95 = ez_val("2.5%"), ez_hi95 = ez_val("97.5%")) %>%
  left_join(d5_tab,    by = "ddm5") %>%
  left_join(grabowska, by = "row") %>%
  select(row, label, ez_par, ez_mean, ez_lo95, ez_hi95,
         ddm5_par = ddm5, ddm5_mean, ddm5_lo95, ddm5_hi95,
         grab_mean, grab_lo95, grab_hi95, note)

write.csv(comp3, file.path(FIT_DIR, paste0("ez_ddm5_grabowska", tag, ".csv")), row.names = FALSE)

fmt <- function(m, lo, hi) ifelse(is.na(m), "-", sprintf("%.3f [%.3f, %.3f]", m, lo, hi))
comp3_print <- data.frame(
  Parameter = comp3$label,
  EZ        = fmt(comp3$ez_mean,   comp3$ez_lo95,   comp3$ez_hi95),
  DDM5      = fmt(comp3$ddm5_mean, comp3$ddm5_lo95, comp3$ddm5_hi95),
  Grabowska = fmt(comp3$grab_mean, comp3$grab_lo95, comp3$grab_hi95),
  Note      = comp3$note
)
cat("\n EZ vs DDM5 vs Grabowska et al. (2025), population level\n")
print(comp3_print[, c("Parameter", "EZ", "DDM5", "Grabowska")], right = FALSE, row.names = FALSE)
cat("(notes on scales and model differences: column Note in the csv)\n")
write.csv(comp3_print, file.path(FIT_DIR, paste0("ez_ddm5_grabowska_formatted", tag, ".csv")),
          row.names = FALSE)


#### person level: EZ vs DDM5 (same participant_index in both fits)

if (!is.null(ez_s)) {
  d5 <- as.data.frame(summ)
  rownames(d5) <- d5$variable
  person_cmp <- bind_rows(lapply(c("nu_c", "nu_r", "nu_cr", "alpha_r"), function(p) {
    data.frame(par  = p, subj = subj_ids,
               ez   = ez_s[paste0(p, "[", subj_ids, "]"), "mean"],
               ddm5 = d5[paste0(p, "[", seq_along(subj_ids), "]"), "mean"])
  }))
  comp_person <- person_cmp %>% group_by(par) %>%
    summarise(r = cor(ez, ddm5), .groups = "drop")
  cat("\n EZ vs DDM5, person-level correlation of posterior means\n",
      "(all on the natural scale)\n")
  print(comp_person)
  write.csv(person_cmp,  file.path(FIT_DIR, paste0("ez_vs_ddm5_person", tag, ".csv")),     row.names = FALSE)
  write.csv(comp_person, file.path(FIT_DIR, paste0("ez_vs_ddm5_person_cor", tag, ".csv")), row.names = FALSE)
}

cat("\nSaved in", FIT_DIR, "--> plots: run plot_ddm5_krakow.R\n")
