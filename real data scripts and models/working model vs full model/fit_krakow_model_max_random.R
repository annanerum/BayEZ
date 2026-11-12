# 4 cells, contrast coded +-1
#   k=1: cond= 1, resp= 1   k=2: cond= 1, resp=-1
#   k=3: cond=-1, resp= 1   k=4: cond=-1, resp=-1


# nu	Intercept	Random (nu_intercept[i])
# nu	Condition (c)	Random (nu_c[i])
# nu	Resp_type (r)	Random (nu_r[i])
# nu	Interaction (cr)	Random (nu_cr[i])
# alpha	Intercept	Random (alpha_intercept[i])
# alpha	Resp_type (r)	Random (alpha_r[i])
# alpha	Condition	Random (alpha_c[i])
# tau	Intercept	Random (tau_intercept[i]
                      
                      
library(rstan)
library(ggplot2)
library(dplyr)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores() - 1)

TEST_RUN <- TRUE

OUT       <- "output_krakow_max_random"
STAN_FILE <- "grabowska_ez_model_max_random.stan"  
DATA_FILE <- "krakow_data_standardized.csv"
if (!dir.exists(OUT)) dir.create(OUT)

K       <- 4L
P_NU    <- 4L   # intercept, condition, response_type, condition*response_type
P_ALPHA <- 3L   # intercept, response_type, condition
P_TAU   <- 1L   # intercept only

# rows = cells
# cols = [intercept, cond, resp, cond*resp])
CELL_DESIGN <- matrix(
  c( 1,  1,  1,  1,
     1,  1, -1, -1,
     1, -1,  1, -1,
     1, -1, -1,  1),
  nrow = K, byrow = TRUE
)


## priors from Appendix C

priors <- list(
  # beta_nu: [intercept, condition, response_type, interaction]
  # intercept ~ N(0,2); condition ~ N(0,2); resp/cr ~ N(0,0.5)
  # Appendix C: mu_theta(delta) ~ N(0,2), mu_theta(c) ~ N(0,2), mu_theta(e), mu_theta(er) ~ N(0,0.5)
  prior_beta_nu = matrix(
    c(0,    0,    0,    0,      # means
      2,    2,    0.5,  0.5),   # SDs
    nrow = 2, byrow = TRUE
  ),

  # beta_alpha: [intercept, response_type, condition]
  # intercept ~ N(1,1); resp ~ N(0,0.2); condition ~ N(0,0.2)
  # Appendix C: mu_theta(alpha) ~ N(1,1)[0,3], mu_theta(e,er,...) ~ N(0,0.2) for boundary
  prior_beta_alpha = matrix(
    c(1,    0,     0,      # means
      1,    0.2,   0.2),   # SDs
    nrow = 2, byrow = TRUE
  ),

  # beta_tau: intercept only
  # ~ N(0.1, 0.2)
  # Appendix C eq. C.4: mu_theta(tau) ~ N(0.1,0.2)[0,0.4]
  prior_beta_tau = matrix(
    c(0.1,           # mean
      0.2),          # SD
    nrow = 2, byrow = TRUE
  ),

  # SDs: half-normal 
  # sigma_nu, sigma_alpha ~ Gamma(1,1) in paper --> approximation: N+(0.3, 0.3)
  # sigma_tau ~ Gamma(0.3,1) --> approximation: N+(0.1, 0.15)
  prior_sigma_nu    = c(0.3,  0.3),
  prior_sigma_alpha = c(0.3,  0.3),
  prior_sigma_tau   = c(0.1,  0.15),

  # random slope SDs: half-normal N+(0, 0.5)
  # Appendix C: sigma_theta(e), sigma_theta(er) ~ N+(0, 0.5)
  prior_sigma_nu_c    = c(0, 0.5),
  prior_sigma_nu_r    = c(0, 0.5),
  prior_sigma_nu_cr   = c(0, 0.5),
  prior_sigma_alpha_r = c(0, 0.5),

  # random intercept SDs (between-person variation in intercepts)
  # own creation, not in Appendix C
  prior_sigma_nu_intercept    = c(0, 1),
  prior_sigma_alpha_intercept = c(0, 0.3),
  prior_sigma_tau_intercept   = c(0, 0.1),

  # alpha's condition effect, random, same scale as alpha_r's prior
  # by symmetry, own creation, not in Appendix C
  prior_sigma_alpha_c = c(0, 0.5)
)


## load and prepare data
cat("Loading data \n")
dat <- read.csv(DATA_FILE)
dat <- dat[dat$is_in_sequence          == "True"  &
           dat$rt_greater_than_1       == "False" &
           dat$log_rt_exceed_threshold == "False", ]

dat$pid     <- as.integer(dat$participant_index)
dat$cond    <- as.numeric(dat$condition)   # congruent=1,  incongruent=-1
dat$resp    <- as.numeric(dat$pre_acc)     # post-correct=1, post-error=-1
dat$acc_bin <- ifelse(dat$acc == 1, 1L, 0L)
dat$cell    <- with(dat, dplyr::case_when(
  cond ==  1 & resp ==  1 ~ 1L,
  cond ==  1 & resp == -1 ~ 2L,
  cond == -1 & resp ==  1 ~ 3L,
  cond == -1 & resp == -1 ~ 4L
))

I <- max(dat$pid)
cat(sprintf("Participants: %d | Trials: %d\n", I, nrow(dat)))


## summary stats
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
# empty cells get cell mean across participants
for (k in seq_len(K)) {
  MRT_mat[is.na(MRT_mat[, k]), k] <- mean(MRT_mat[, k], na.rm = TRUE)
  VRT_mat[is.na(VRT_mat[, k]), k] <- mean(VRT_mat[, k], na.rm = TRUE)
}
C_mat <- pmax(pmin(C_mat, J_mat - 1L), 1L)

cat("Trials per cell (mean):", round(colMeans(J_mat), 1), "\n")


## design matrices
X_nu    <- array(NA_real_, c(I, K, P_NU))
X_alpha <- array(NA_real_, c(I, K, P_ALPHA))
X_tau   <- array(NA_real_, c(I, K, P_TAU))

for (i in seq_len(I)) for (k in seq_len(K)) {
  X_nu[i, k, ]    <- CELL_DESIGN[k, ]          # [1, cond, resp, cond*resp]
  X_alpha[i, k, ] <- CELL_DESIGN[k, c(1, 3, 2)]   # [1, resp, cond]
  X_tau[i, k, ]   <- CELL_DESIGN[k, 1]         # [1]
}


## stan data: data + priors
stan_data <- c(
  list(
    I = I, K = K, P_nu = P_NU, P_alpha = P_ALPHA, P_tau = P_TAU,
    J = J_mat, C = C_mat, MRT = MRT_mat, VRT = VRT_mat,
    X_nu = X_nu, X_alpha = X_alpha, X_tau = X_tau
  ),
  priors
)


# fitting model 
GROUP_PARS <- c(
  "beta_nu[1]", "beta_nu[2]", "beta_nu[3]", "beta_nu[4]",
  "beta_alpha[1]", "beta_alpha[2]", "beta_alpha[3]",
  "beta_tau",
  "sigma_nu", "sigma_alpha", "sigma_tau",
  "sigma_nu_c", "sigma_nu_r", "sigma_nu_cr", "sigma_alpha_r",
  "sigma_nu_intercept",     # NEW: active again
  "sigma_alpha_intercept", "sigma_tau_intercept",
  "sigma_alpha_c"           # NEW
)

cat("Compiling \n")
stan_mod <- stan_model(STAN_FILE)

cat("Fitting \n")
if (TEST_RUN) {
  cat(">>> TEST_RUN = TRUE: short cheap run, just checking for compile/divergence issues\n")
  fit <- sampling(
    stan_mod, data = stan_data,
    chains = 2, iter = 500, warmup = 250,
    cores  = getOption("mc.cores"),
    refresh = 100
  )
} else {
  fit <- sampling(
    stan_mod, data = stan_data,
    chains = 4, iter = 2000, warmup = 1000,
    cores  = getOption("mc.cores"),
    refresh = 200
  )

  s        <- summary(fit)$summary
  rhat_max <- max(s[GROUP_PARS, "Rhat"], na.rm = TRUE)
  cat(sprintf("\nMax Rhat (group): %.3f\n", rhat_max))

  if (rhat_max > 1.05) {
    cat("Refitting with more iterations \n")
    fit <- sampling(
      stan_mod, data = stan_data,
      chains = 6, iter = 4000, warmup = 2000,
      cores  = getOption("mc.cores"),
      control = list(adapt_delta = 0.90, max_treedepth = 12),
      refresh = 500
    )
  }
}

s <- summary(fit)$summary

## HMC diagnostics 
cat("\n=== HMC diagnostics ===\n")
n_div       <- rstan::get_num_divergent(fit)
n_treedepth <- rstan::get_num_max_treedepth(fit)
low_bfmi    <- rstan::get_low_bfmi_chains(fit)
rhat_max    <- max(s[GROUP_PARS, "Rhat"], na.rm = TRUE)
rhat_worst  <- sort(s[GROUP_PARS, "Rhat"], decreasing = TRUE, na.last = TRUE)[1:5]

cat(sprintf("Divergent transitions : %d\n", n_div))
cat(sprintf("Max treedepth hits    : %d\n", n_treedepth))
cat(sprintf("Chains w/ low BFMI    : %s\n",
            if (length(low_bfmi) == 0) "none" else paste(low_bfmi, collapse = ", ")))
cat(sprintf("Max Rhat (group pars) : %.3f\n", rhat_max))
cat("Worst 5 Rhat values:\n")
print(round(rhat_worst, 3))

if (TEST_RUN) {
  cat("\n>>> TEST_RUN complete. Review the diagnostics above before setting TEST_RUN <- FALSE.\n")
  cat(">>> Given the structural expectation (nu saturated, alpha with less slack),\n")
  cat(">>> a high max Rhat / many divergences here would be a documented, expected result,\n")
  cat(">>> not necessarily a bug to chase further.\n")
}

saveRDS(fit, file.path(OUT, "krakow_fit.rds"))

## results 
PAR_LABELS <- c(
  "mu_nu", "b_nu_condition", "b_nu_resp_type", "b_nu_interaction",
  "mu_alpha", "b_alpha_resp_type", "b_alpha_condition",
  "mu_tau",
  "sigma_nu", "sigma_alpha", "sigma_tau",
  "sigma_nu_cond", "sigma_nu_resp", "sigma_nu_cr", "sigma_alpha_resp",
  "sigma_nu_intercept",
  "sigma_alpha_intercept", "sigma_tau_intercept",
  "sigma_alpha_cond"
)

est_df <- data.frame(
  param = PAR_LABELS,
  mean  = s[GROUP_PARS, "mean"],
  sd    = s[GROUP_PARS, "sd"],
  lo95  = s[GROUP_PARS, "2.5%"],
  hi95  = s[GROUP_PARS, "97.5%"],
  rhat  = s[GROUP_PARS, "Rhat"],
  n_eff = s[GROUP_PARS, "n_eff"],
  row.names = NULL
)

cat("\n group estimates\n")
print(est_df[, c("param","mean","sd","lo95","hi95","rhat")],
      digits = 3, row.names = FALSE)
write.csv(est_df, file.path(OUT, "krakow_estimates.csv"), row.names = FALSE)


## plots
dark <- theme_minimal(base_size = 12) +
  theme(
    plot.background   = element_rect(fill = "#1e1e1e", colour = NA),
    panel.background  = element_rect(fill = "#1e1e1e", colour = NA),
    panel.grid.major  = element_line(colour = "#333333"),
    panel.grid.minor  = element_blank(),
    text              = element_text(colour = "white"),
    axis.text         = element_text(colour = "white"),
    strip.text        = element_text(colour = "white")
  )

sv <- function(p, nm, w = 10, h = 7)
  ggsave(file.path(OUT, nm), p, width = w, height = h, dpi = 300, bg = "#1e1e1e")


# forest plot
p_forest <- ggplot(est_df[1:8, ],
                   aes(x = mean, xmin = lo95, xmax = hi95,
                       y = reorder(param, mean))) +
  geom_vline(xintercept = 0, colour = "white", linetype = "dashed") +
  geom_errorbarh(colour = "steelblue", height = 0.3) +
  geom_point(colour = "goldenrod", size = 2.5) +
  labs(title = "Population effects",
       x = "Posterior mean (95% CI)", y = NULL) + dark
sv(p_forest, "krakow_forest.png")

# individual response_type slopes: nu and alpha
for (par_name in c("nu_r", "alpha_r")) {
  label    <- if (par_name == "nu_r") "nu (drift rate)" else "alpha (boundary)"
  colour   <- if (par_name == "nu_r") "steelblue" else "#e07070"
  pop_mean <- if (par_name == "nu_r") est_df$mean[est_df$param == "b_nu_resp_type"] else
                                       est_df$mean[est_df$param == "b_alpha_resp_type"]
  means <- s[paste0(par_name, "[", seq_len(I), "]"), "mean"]
  lo    <- s[paste0(par_name, "[", seq_len(I), "]"), "2.5%"]
  hi    <- s[paste0(par_name, "[", seq_len(I), "]"), "97.5%"]

  df_sl <- data.frame(mean = means, lo = lo, hi = hi) |>
    arrange(mean) |> mutate(rank = seq_len(I))

  p_sl <- ggplot(df_sl, aes(x = rank, y = mean, ymin = lo, ymax = hi)) +
    geom_hline(yintercept = 0,        colour = "white",    linetype = "dashed") +
    geom_hline(yintercept = pop_mean, colour = "goldenrod", linewidth = 0.8) +
    geom_errorbar(colour = colour, width = 0, alpha = 0.5) +
    geom_point(size = 0.8, colour = colour) +
    labs(title = paste("Individual response_type slopes ", label),
         subtitle = "Sorted, gold = population mean",
         x = "Rank", y = "Slope") + dark
  sv(p_sl, paste0("krakow_slopes_", par_name, ".png"))
}

cat("\n outputs saved to:", OUT, "\n")


## more plots (later)
FIT_RDS <- file.path(OUT, "krakow_fit.rds")

cat("Loading fit \n")
fit <- readRDS(FIT_RDS)
s   <- summary(fit)$summary
I   <- 222L
K   <- 4L

CELL_LABELS  <- c("cong/post-corr", "cong/post-err",
                  "incong/post-corr", "incong/post-err")
CELL_COLOURS <- c("goldenrod", "steelblue", "#e07070", "#a070e0")

dark <- theme_minimal(base_size = 12) +
  theme(
    plot.background   = element_rect(fill = "#1e1e1e", colour = NA),
    panel.background  = element_rect(fill = "#1e1e1e", colour = NA),
    panel.grid.major  = element_line(colour = "#333333"),
    panel.grid.minor  = element_blank(),
    text              = element_text(colour = "white"),
    axis.text         = element_text(colour = "white"),
    strip.text        = element_text(colour = "white"),
    legend.background = element_rect(fill = "#1e1e1e"),
    legend.text       = element_text(colour = "white"),
    axis.title        = element_text(colour = "white")
  )

sv <- function(p, nm, w = 10, h = 7)
  ggsave(file.path(OUT, nm), p, width = w, height = h, dpi = 300, bg = "#1e1e1e")

extract_summary <- function(pars, labels) {
  data.frame(
    param = labels,
    mean  = s[pars, "mean"],
    lo50  = s[pars, "25%"],
    hi50  = s[pars, "75%"],
    lo95  = s[pars, "2.5%"],
    hi95  = s[pars, "97.5%"],
    row.names = NULL
  )
}


# forest plot: group-level effects (betas)
group_pars   <- c("beta_nu[1]","beta_nu[2]","beta_nu[3]","beta_nu[4]",
                  "beta_alpha[1]","beta_alpha[2]","beta_alpha[3]","beta_tau")
group_labels <- c("mu_nu (intercept)","b_nu_condition","b_nu_resp_type",
                  "b_nu_interaction","mu_alpha (intercept)",
                  "b_alpha_resp_type","b_alpha_condition","mu_tau")

grp_df <- extract_summary(group_pars, group_labels) |>
  mutate(type = c("nu","nu","nu","nu","alpha","alpha","alpha","tau"))

p_forest_group <- ggplot(grp_df,
                         aes(x = mean, xmin = lo95, xmax = hi95,
                             y = reorder(param, mean), colour = type)) +
  geom_vline(xintercept = 0, colour = "white", linetype = "dashed") +
  geom_errorbarh(aes(xmin = lo95, xmax = hi95), height = 0, linewidth = 1.2, alpha = 0.4) +
  geom_errorbarh(aes(xmin = lo50, xmax = hi50), height = 0, linewidth = 2.5, alpha = 0.8) +
  geom_point(size = 3) +
  scale_colour_manual(values = c(nu="steelblue", alpha="#e07070", tau="goldenrod")) +
  labs(title = "Population effects — Krakow EZ-DDM (maximal random effects)",
       subtitle = "Thick bar = 50% CI | thin bar = 95% CI",
       x = "Posterior mean", y = NULL, colour = "Parameter") +
  dark
sv(p_forest_group, "forest_group.png", h = 6)


# forest plot: SD parameters -- sigma_alpha_c added, sigma_nu_intercept restored
sd_pars   <- c("sigma_nu","sigma_alpha","sigma_tau",
               "sigma_nu_c","sigma_nu_r","sigma_nu_cr","sigma_alpha_r","sigma_alpha_c",
               "sigma_nu_intercept","sigma_alpha_intercept","sigma_tau_intercept")
sd_labels <- c("sigma_nu (residual)","sigma_alpha (residual)","sigma_tau (residual)",
               "sigma_nu_cond (slope)","sigma_nu_resp (slope)",
               "sigma_nu_cr (slope)","sigma_alpha_resp (slope)","sigma_alpha_cond (slope)",
               "sigma_nu (intercept)","sigma_alpha (intercept)","sigma_tau (intercept)")

sd_df <- extract_summary(sd_pars, sd_labels) |>
  mutate(type = c("residual","residual","residual",
                   "slope","slope","slope","slope","slope",
                   "intercept","intercept","intercept"))

p_forest_sds <- ggplot(sd_df,
                       aes(x = mean, xmin = lo95, xmax = hi95,
                           y = reorder(param, mean), colour = type)) +
  geom_vline(xintercept = 0, colour = "white", linetype = "dashed") +
  geom_errorbarh(aes(xmin = lo95, xmax = hi95), height = 0, linewidth = 1.2, alpha = 0.4) +
  geom_errorbarh(aes(xmin = lo50, xmax = hi50), height = 0, linewidth = 2.5, alpha = 0.8) +
  geom_point(size = 3) +
  scale_colour_manual(values = c(residual = "goldenrod", slope = "steelblue",
                                  intercept = "seagreen")) +
  labs(title = "Between-person SDs (maximal random effects)",
       subtitle = "Thick bar = 50% CI, thin bar = 95% CI",
       x = "Posterior mean", y = NULL, colour = NULL) +
  dark
sv(p_forest_sds, "forest_sds.png", h = 6)


# caterpillar plots -- nu_intercept restored, alpha_c added
make_caterpillar <- function(par_stem, pop_mean, title_str, col) {
  pars  <- paste0(par_stem, "[", seq_len(I), "]")
  df <- data.frame(
    mean = s[pars, "mean"],
    lo95 = s[pars, "2.5%"],
    hi95 = s[pars, "97.5%"]
  ) |> arrange(mean) |> mutate(rank = seq_len(I))
  
  ggplot(df, aes(x = rank, y = mean, ymin = lo95, ymax = hi95)) +
    geom_hline(yintercept = 0,        colour = "white",   linetype = "dashed") +
    geom_hline(yintercept = pop_mean, colour = "goldenrod", linewidth = 0.8) +
    geom_errorbar(colour = col, width = 0, alpha = 0.4) +
    geom_point(size = 0.8, colour = col) +
    labs(title = title_str,
         subtitle = "Sorted by posterior mean, gold = population mean",
         x = "Participant rank", y = "Slope") +
    dark
}

pm <- function(par) s[par, "mean"]

sv(make_caterpillar("nu_c",    pm("beta_nu[2]"),    "Condition slope on nu",       "steelblue"),
   "caterpillar_nu_c.png")
sv(make_caterpillar("nu_r",    pm("beta_nu[3]"),    "Resp_type slope on nu",       "steelblue"),
   "caterpillar_nu_r.png")
sv(make_caterpillar("nu_cr",   pm("beta_nu[4]"),    "Interaction slope on nu",     "steelblue"),
   "caterpillar_nu_cr.png")
sv(make_caterpillar("alpha_r", pm("beta_alpha[2]"), "Resp_type slope on alpha",    "#e07070"),
   "caterpillar_alpha_r.png")
sv(make_caterpillar("alpha_c", pm("beta_alpha[3]"), "Condition slope on alpha",    "#e07070"),
   "caterpillar_alpha_c.png")
sv(make_caterpillar("nu_intercept",    pm("beta_nu[1]"),    "Individual intercepts, nu",    "seagreen"),
   "caterpillar_nu_intercept.png")
sv(make_caterpillar("alpha_intercept", pm("beta_alpha[1]"), "Individual intercepts, alpha", "seagreen"),
   "caterpillar_alpha_intercept.png")
sv(make_caterpillar("tau_intercept",   pm("beta_tau"),      "Individual intercepts, tau",   "seagreen"),
   "caterpillar_tau_intercept.png")


# violin + jitter of posterior means
make_cell_df <- function(par_stem) {
  bind_rows(lapply(seq_len(K), function(k) {
    pars <- paste0(par_stem, "[", seq_len(I), ",", k, "]")
    data.frame(
      value      = s[pars, "mean"],
      cell       = k,
      cell_label = factor(CELL_LABELS[k], levels = CELL_LABELS)
    )
  }))
}

plot_cells <- function(df, par_name, col_vec) {
  ggplot(df, aes(x = cell_label, y = value, fill = cell_label, colour = cell_label)) +
    geom_violin(alpha = 0.25, linewidth = 0.4) +
    geom_jitter(width = 0.15, size = 0.5, alpha = 0.4) +
    stat_summary(fun = mean, geom = "point", size = 3,
                 colour = "white", shape = 18) +
    scale_fill_manual(values   = col_vec) +
    scale_colour_manual(values = col_vec) +
    labs(title = paste0("Cell-level, ", par_name ,
                        "posterior means across participants"),
         subtitle = "Diamond = cell mean, each dot = one participant",
         x = NULL, y = par_name) +
    dark + theme(legend.position = "none",
                 axis.text.x = element_text(angle = 20, hjust = 1))
}

sv(plot_cells(make_cell_df("nu"),    "nu",    CELL_COLOURS), "cells_nu.png")
sv(plot_cells(make_cell_df("alpha"), "alpha", CELL_COLOURS), "cells_alpha.png")
sv(plot_cells(make_cell_df("tau"),   "tau",   CELL_COLOURS), "cells_tau.png")

cat("\nAll plots saved to:", OUT, "\n")
cat("Files:\n")
cat(paste0("  ", list.files(OUT, pattern="\\.png$|\\.csv$"), collapse="\n"), "\n")
