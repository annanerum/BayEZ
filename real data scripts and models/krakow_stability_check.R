library(rstan)
library(dplyr)
library(ggplot2)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores() - 1)

OUT       <- "output_krakow_stability_check"
STAN_FILE <- "grabowska_ez_model_new.stan"
DATA_FILE <- "krakow_data_standardized.csv"
if (!dir.exists(OUT)) dir.create(OUT)

K       <- 4L
P_NU    <- 4L
P_ALPHA <- 3L
P_TAU   <- 1L

CELL_DESIGN <- matrix(
  c( 1,  1,  1,  1,
     1,  1, -1, -1,
     1, -1,  1, -1,
     1, -1, -1,  1),
  nrow = K, byrow = TRUE
)

## priors: 
priors <- list(
  prior_beta_nu = matrix(
    c(0,    0,    0,    0,
      2,    2,    0.5,  0.5),
    nrow = 2, byrow = TRUE
  ),
  prior_beta_alpha = matrix(
    c(1,    0,     0,
      1,    0.2,   0.2),
    nrow = 2, byrow = TRUE
  ),
  prior_beta_tau = matrix(
    c(0.1,
      0.2),
    nrow = 2, byrow = TRUE
  ),
  prior_sigma_nu    = c(0.3,  0.3),
  prior_sigma_alpha = c(0.3,  0.3),
  prior_sigma_tau   = c(0.1,  0.15),
  prior_sigma_nu_c    = c(0, 0.5),
  prior_sigma_nu_r    = c(0, 0.5),
  prior_sigma_nu_cr   = c(0, 0.5),
  prior_sigma_alpha_r = c(0, 0.5)
)

## data prep
cat("Loading data \n")
dat <- read.csv(DATA_FILE)
dat <- dat[dat$is_in_sequence          == "True"  &
           dat$rt_greater_than_1       == "False" &
           dat$log_rt_exceed_threshold == "False", ]

dat$pid     <- as.integer(dat$participant_index)
dat$cond    <- as.numeric(dat$condition)
dat$resp    <- as.numeric(dat$pre_acc)
dat$acc_bin <- ifelse(dat$acc == 1, 1L, 0L)
dat$cell    <- with(dat, dplyr::case_when(
  cond ==  1 & resp ==  1 ~ 1L,
  cond ==  1 & resp == -1 ~ 2L,
  cond == -1 & resp ==  1 ~ 3L,
  cond == -1 & resp == -1 ~ 4L
))

I <- max(dat$pid)
cat(sprintf("Participants: %d , Trials: %d\n", I, nrow(dat)))

J_mat <- matrix(0L, I, K); C_mat <- matrix(0L, I, K)
MRT_mat <- VRT_mat <- matrix(NA_real_, I, K)

for (i in seq_len(I)) for (k in seq_len(K)) {
  sub <- dat[dat$pid == i & dat$cell == k, ]
  n   <- nrow(sub)
  J_mat[i, k] <- max(n, 2L)
  C_mat[i, k] <- sum(sub$acc_bin)
  if (n >= 2) {
    MRT_mat[i, k] <- mean(sub$rt, na.rm = TRUE)
    VRT_mat[i, k] <- var( sub$rt, na.rm = TRUE)
  }
}
for (k in seq_len(K)) {
  MRT_mat[is.na(MRT_mat[, k]), k] <- mean(MRT_mat[, k], na.rm = TRUE)
  VRT_mat[is.na(VRT_mat[, k]), k] <- mean(VRT_mat[, k], na.rm = TRUE)
}
C_mat <- pmax(pmin(C_mat, J_mat - 1L), 1L)

X_nu    <- array(NA_real_, c(I, K, P_NU))
X_alpha <- array(NA_real_, c(I, K, P_ALPHA))
X_tau   <- array(NA_real_, c(I, K, P_TAU))
for (i in seq_len(I)) for (k in seq_len(K)) {
  X_nu[i, k, ]    <- CELL_DESIGN[k, ]
  X_alpha[i, k, ] <- CELL_DESIGN[k, c(1, 3, 2)]
  X_tau[i, k, ]   <- CELL_DESIGN[k, 1]
}

stan_data <- c(
  list(I = I, K = K, P_nu = P_NU, P_alpha = P_ALPHA, P_tau = P_TAU,
       J = J_mat, C = C_mat, MRT = MRT_mat, VRT = VRT_mat,
       X_nu = X_nu, X_alpha = X_alpha, X_tau = X_tau),
  priors
)

GROUP_PARS <- c(
  "beta_nu[1]", "beta_nu[2]", "beta_nu[3]", "beta_nu[4]",
  "beta_alpha[1]", "beta_alpha[2]", "beta_alpha[3]",
  "beta_tau",
  "sigma_nu", "sigma_alpha", "sigma_tau",
  "sigma_nu_c", "sigma_nu_r", "sigma_nu_cr", "sigma_alpha_r"
)
PAR_LABELS <- c(
  "mu_nu", "b_nu_condition", "b_nu_resp_type", "b_nu_interaction",
  "mu_alpha", "b_alpha_resp_type", "b_alpha_condition",
  "mu_tau",
  "sigma_nu", "sigma_alpha", "sigma_tau",
  "sigma_nu_cond", "sigma_nu_resp", "sigma_nu_cr", "sigma_alpha_resp"
)

cat("Compiling \n")
stan_mod <- stan_model(STAN_FILE)

## run one fit with new faster settings
run_one <- function(seed_val, run_id) {
  cat(sprintf("\n=== Run %d (seed = %d) ===\n", run_id, seed_val))
  fit <- sampling(
    stan_mod, data = stan_data,
    chains = 4, iter = 1500, warmup = 750,   
    cores  = getOption("mc.cores"),
    seed   = seed_val,
    refresh = 0
  )
  s <- summary(fit)$summary
  
  # summary matrix is saved but not the full fit object 
  saveRDS(s, file.path(OUT, paste0("summary_run", run_id, ".rds")))
  data.frame(
    run      = run_id,
    param    = PAR_LABELS,
    mean     = s[GROUP_PARS, "mean"],
    se_mean  = s[GROUP_PARS, "se_mean"],   # MCMC standard error of the mean
    sd       = s[GROUP_PARS, "sd"],
    rhat     = s[GROUP_PARS, "Rhat"],
    n_eff    = s[GROUP_PARS, "n_eff"],
    row.names = NULL
  )
}

## run 3 times 
seeds <- c(101, 202, 303)
results <- lapply(seq_along(seeds), function(j) run_one(seeds[j], j))
all_runs <- bind_rows(results)

# stability check, for each parameter:
# across-run SD of the posterior means, vs. the average se_mean (MCMC standard error) within a run
# if the between-run variance is much bigger than what se_mean alone would predict, there is a trade off and not only noise
check_stability <- function(all_runs) {
  all_runs %>%
    group_by(param) %>%
    summarise(
      mean_of_means   = mean(mean),
      between_run_sd  = sd(mean),
      avg_se_mean     = mean(se_mean),
      # ratio > 2-3 = disagreement
      instability_ratio = between_run_sd / avg_se_mean,
      max_rhat        = max(rhat),
      min_n_eff       = min(n_eff),
      .groups = "drop"
    ) %>%
    arrange(desc(instability_ratio))
}

stability_3 <- check_stability(all_runs)
cat("\n Stability check after 3 runs \n")
print(stability_3, n = 20, digits = 3)

write.csv(all_runs,    file.path(OUT, "stability_all_runs_3.csv"), row.names = FALSE)
write.csv(stability_3, file.path(OUT, "stability_summary_3.csv"),  row.names = FALSE)

unstable <- stability_3 %>% filter(instability_ratio > 3 | max_rhat > 1.05)

if (nrow(unstable) > 0) {
  cat("\n Unstable parameter(s) detected (instability_ratio > 3 or Rhat > 1.05) \n")
  print(unstable$param)
  cat("\n Running 2 more \n")

  more_seeds <- c(404, 505)
  more_results <- lapply(seq_along(more_seeds), function(j)
    run_one(more_seeds[j], j + length(seeds)))
  all_runs <- bind_rows(all_runs, bind_rows(more_results))

  stability_5 <- check_stability(all_runs)
  cat("\n Stability check after 5 runs \n")
  print(stability_5, n = 20, digits = 3)

  write.csv(all_runs,    file.path(OUT, "stability_all_runs_5.csv"), row.names = FALSE)
  write.csv(stability_5, file.path(OUT, "stability_summary_5.csv"),  row.names = FALSE)

  still_unstable <- stability_5 %>% filter(instability_ratio > 3 | max_rhat > 1.05)
  if (nrow(still_unstable) > 0) {
    cat("\n Still unstable after 5 runs --> trade-off/\n")
    cat("identifiability problem:\n")
    print(still_unstable$param)
  } else {
    cat("\n Stable across 5 runs --> just MC noise \n")
  }
} else {
  cat("\n Stable across 3 runs --> no trade-offproblem \n")
}

## point estimate +/- se_mean per run, per parameter 
dark <- theme_minimal(base_size = 11) +
  theme(
    plot.background  = element_rect(fill = "#1e1e1e", colour = NA),
    panel.background = element_rect(fill = "#1e1e1e", colour = NA),
    panel.grid.major = element_line(colour = "#333333"),
    panel.grid.minor = element_blank(),
    text  = element_text(colour = "white"),
    axis.text  = element_text(colour = "white"),
    strip.text = element_text(colour = "white")
  )

p_stab <- ggplot(all_runs, aes(x = factor(run), y = mean,
                                ymin = mean - se_mean, ymax = mean + se_mean)) +
  geom_pointrange(colour = "steelblue") +
  facet_wrap(~ param, scales = "free_y") +
  labs(title = "Estimate stability across independent runs",
       subtitle = "Point = posterior mean, whiskers = MCMC se_mean",
       x = "Run", y = "Estimate") +
  dark

ggsave(file.path(OUT, "stability_check.png"), p_stab, width = 12, height = 9, dpi = 300, bg = "#1e1e1e")

cat("\n outputs in:", OUT, "\n")
