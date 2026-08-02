# 3 MODEL VERSIONS 
# in ALL three, w = 0.5 and sw = 0 
#   "quicktest" -> 7% of participants, 25% of each participant's trials, short warmup/sampling 
#   "full"      -> 100% of participants, 100% of trials 

library(cmdstanr)
library(dplyr)
library(tidyr)

set.seed(42)

MODEL_VERSION <- "sv_only"   # "zero", "free", or "sv_only"
SCALE         <- "full"      # "quicktest", "tryout" or "full"


SCALE_CONFIG <- list(
  
  zero = list(
    quicktest = list(frac_participants = 0.07, frac_trials = 0.25, warmup = 50,   sampling = 30),
    tryout    = list(frac_participants = 0.30, frac_trials = 0.50, warmup = 500,  sampling = 300),
    full      = list(frac_participants = 1.00, frac_trials = 1.00, warmup = 3000, sampling = 3000)
  ),
  
  free = list(
    quicktest = list(frac_participants = 0.07, frac_trials = 0.25, warmup = 50,  sampling = 30),
    tryout    = list(frac_participants = 0.30, frac_trials = 0.50, warmup = 100, sampling = 60),
    full      = list(frac_participants = 1.00, frac_trials = 1.00, warmup = 300, sampling = 200)
  ),
  
  sv_only = list(
    quicktest = list(frac_participants = 0.07, frac_trials = 0.25, warmup = 50,  sampling = 30),
    tryout    = list(frac_participants = 0.30, frac_trials = 0.30, warmup = 300, sampling = 200),
    full      = list(frac_participants = 1.00, frac_trials = 1.00, warmup = 2000, sampling = 3000)
  )
  
)

cfg <- SCALE_CONFIG[[MODEL_VERSION]][[SCALE]]
FRAC_PARTICIPANTS <- cfg$frac_participants
FRAC_TRIALS       <- cfg$frac_trials
WARMUP            <- cfg$warmup
SAMPLING          <- cfg$sampling


TREEDEPTH   <- 10    
ADAPT_DELTA <- 0.8   
REFRESH     <- 10    

CHAINS          <- 2
PARALLEL_CHAINS <- 2  
N_THREADS <- max(1, parallel::detectCores() - 1)
THREADS_PER_CHAIN <- max(1, floor(N_THREADS / PARALLEL_CHAINS))

data_path <- "krakow_data_standardized.csv"
stan_file <- switch(MODEL_VERSION,
                    zero    = "hierarchical_fullddm_no_intertrial_variability_nu_alpha.stan",
                    free    = "hierarchical_fullddm_intertrial_variability_nu_alpha.stan",
                    sv_only = "hierarchical_fullddm_sv_only_nu_alpha.stan"
)

out_dir <- sprintf("output_ddm_v2_%s_%s", MODEL_VERSION, SCALE)
if (!dir.exists(out_dir)) dir.create(out_dir)

message(sprintf(
  "Config: MODEL_VERSION=%s, SCALE=%s, treedepth=%d, adapt_delta=%.2f, stan_file=%s, out_dir=%s",
  MODEL_VERSION, SCALE, TREEDEPTH, ADAPT_DELTA, stan_file, out_dir
))
message(sprintf(
  "FRAC_PARTICIPANTS=%.2f, FRAC_TRIALS=%.2f, WARMUP=%d, SAMPLING=%d",
  FRAC_PARTICIPANTS, FRAC_TRIALS, WARMUP, SAMPLING
))


raw_full <- read.csv(data_path) %>%
  transmute(
    participant_orig = participant_index,
    condition  = condition,
    resp_type  = pre_acc,
    acc        = acc,
    rt         = rt
  )

stopifnot(all(raw_full$acc %in% c(-1, 1)))
stopifnot(all(raw_full$condition %in% c(-1, 1)))
stopifnot(all(raw_full$resp_type %in% c(-1, 1)))
stopifnot(all(raw_full$rt > 0))

# randomly select a fraction of participants 
all_ids <- unique(raw_full$participant_orig)
n_keep  <- round(FRAC_PARTICIPANTS * length(all_ids))
keep_ids <- sample(all_ids, n_keep)
raw <- raw_full %>% filter(participant_orig %in% keep_ids)

message(sprintf("Selected %d / %d participants (%.0f%%)",
                n_keep, length(all_ids), 100 * FRAC_PARTICIPANTS))

#  within each kept participant, keep a fraction of their trials
stratified_fraction <- function(df, frac) {
  df %>%
    group_by(participant_orig, condition, resp_type) %>%
    group_modify(~ {
      n_cell <- nrow(.x)
      target <- max(1, round(frac * n_cell))  # keep at least 1 per cell
      if (n_cell <= target) return(.x)
      .x[sample.int(n_cell, target), , drop = FALSE]
    }) %>%
    ungroup()
}

n_before <- nrow(raw)
raw <- stratified_fraction(raw, FRAC_TRIALS)
message(sprintf("Trial fraction %.0f%%: %d -> %d trials", 100 * FRAC_TRIALS, n_before, nrow(raw)))

# re-index participants to a contiguous 1..I 
raw <- raw %>% mutate(pid = as.integer(factor(participant_orig, levels = unique(participant_orig))))
I <- max(raw$pid)
N <- nrow(raw)

# per-participant minimum RT, used to bound non-decision time
# (t0 = min_rt * inv_logit(...) in the Stan model, guaranteeing t0 < min_rt
# <= every trial's RT for that person, which wiener_lpdf requires)
min_rt <- raw %>%
  group_by(pid) %>%
  summarise(min_rt = min(rt), .groups = "drop") %>%
  arrange(pid) %>%
  pull(min_rt)
stopifnot(length(min_rt) == I)

message(sprintf("Final dataset: I=%d participants, N=%d trials", I, N))


## priors from Grabowska et al. (2025) Appendix C
# Grabowska's priors are specified on the natural scale for alpha (raw boundary separation, 
# e.g. Normal(1,1)) and raw seconds for tau
# my model uses a log scale for alpha (for positivity) and a logit of min_rt scale for tau (to keep 0 < t0 < min_rt)
# the numbers below are therefore approximated

# nu (drift rate): intercept fixed, cond/resp/interaction random 
# mu_theta(delta) [intercept] ~ Normal(0, 2)      
# mu_theta(c)     [condition] ~ Normal(0, 2)      
# theta(r)        [resp_type, their fixed effect] ~ Normal(0, 0.5)

prior_beta_nu <- rbind(
  c(0, 0, 0, 0),      # means: intercept, cond, resp, cond:resp
  c(2, 2, 0.5, 0.5)   # sds
)

# sigma_theta(delta), sigma_theta(c) ~ Gamma(1,1) in Grabowska 
prior_sigma_nu_c  <- c(1, 1)  
prior_sigma_nu_r  <- c(1, 1)
prior_sigma_nu_cr <- c(1, 1)  

# alpha (boundary separation, log scale): intercept random, resp/cond fixed 
# Grabowska (natural scale): mu_theta(alpha) ~ Normal(1, 1)[0,3]
#                            mu_theta(c)     ~ Normal(0, 1)
#                            theta(r)        ~ Normal(0, 0.2)
# translated to log scale: natural Normal(1,1) truncated [0,3] has its mass
# concentrated around alpha=1 (log(1)=0), so we center the log-intercept at 0

prior_beta_alpha1 <- c(0, 0.5)    # log-alpha intercept, centered at log(1)=0
prior_beta_alpha2 <- c(0, 0.2)    # resp_type (fixed), translated from Normal(0,0.2)
prior_beta_alpha3 <- c(0, 0.3)    # condition (fixed), translated from Normal(0,1)
prior_sigma_alpha_intercept <- c(1, 1)  # Gamma(1,1), matching Grabowska's sigma_theta(alpha)

# tau (non-decision time, logit of min_rt scale): intercept random only 
# Grabowska: mu_theta(tau) ~ Normal(0.1, 0.2)[0, 0.4]
#           sigma_theta(tau) ~ Gamma(0.3, 1)
prior_beta_tau <- c(0.4, 0.5)
prior_sigma_tau_intercept <- c(0.3, 1)  # Gamma(0.3, 1) 

# sv, st0: sv is used by "free" AND "sv_only"; st0 only by "free"
# sw stays fixed at 0 in all three, Grabowska fixed all three (sv, sw, st0) to 0
prior_sv  <- c(0.5, 0.5)
prior_st0 <- c(0.3, 0.2)   # prior on st0_raw, the (0,1) fraction of the max allowed

# STAN DATA LIST 
stan_data <- list(
  N = N, I = I,
  pid = raw$pid,
  rt = raw$rt,
  acc = as.integer(raw$acc),
  condition = raw$condition,
  resp_type = raw$resp_type,
  min_rt = min_rt,
  grainsize = 1L,
  
  prior_beta_nu = prior_beta_nu,
  prior_beta_alpha1 = prior_beta_alpha1,
  prior_beta_alpha2 = prior_beta_alpha2,
  prior_beta_alpha3 = prior_beta_alpha3,
  prior_beta_tau    = prior_beta_tau,
  
  prior_sigma_nu_c  = prior_sigma_nu_c,
  prior_sigma_nu_r  = prior_sigma_nu_r,
  prior_sigma_nu_cr = prior_sigma_nu_cr,
  prior_sigma_alpha_intercept = prior_sigma_alpha_intercept,
  prior_sigma_tau_intercept   = prior_sigma_tau_intercept
)

#   "zero"    -> no sv/st0 
#   "sv_only" -> prior_sv only (sv is free, st0 hardcoded)
#   "free"    -> both prior_sv and prior_st0 free
if (MODEL_VERSION %in% c("free", "sv_only")) {
  stan_data$prior_sv <- prior_sv
}
if (MODEL_VERSION == "free") {
  stan_data$prior_st0 <- prior_st0
}

# INITIAL VALUES
# centered near the prior means (small jitter per chain) so the sampler
# starts somewhere plausible rather than Stan's default random unconstrained draw
init_fun <- function() {
  base <- list(
    beta_nu     = c(0, 0, 0, 0) + rnorm(4, 0, 0.05),
    beta_alpha1 = 0 + rnorm(1, 0, 0.05),
    beta_alpha2 = 0 + rnorm(1, 0, 0.02),
    beta_alpha3 = 0 + rnorm(1, 0, 0.02),
    beta_tau    = 0.4 + rnorm(1, 0, 0.05),
    
    sigma_nu_c  = 0.5 + runif(1, 0, 0.1),
    sigma_nu_r  = 0.3 + runif(1, 0, 0.05),
    sigma_nu_cr = 0.3 + runif(1, 0, 0.05),
    sigma_alpha_intercept = 0.3 + runif(1, 0, 0.05),
    sigma_tau_intercept   = 0.3 + runif(1, 0, 0.05),
    
    z_nu_c  = rnorm(I, 0, 0.1),
    z_nu_r  = rnorm(I, 0, 0.1),
    z_nu_cr = rnorm(I, 0, 0.1),
    z_alpha_intercept = rnorm(I, 0, 0.1),
    z_tau_intercept   = rnorm(I, 0, 0.1)
    # no w init needed in either version, w=0.5 is hardcoded in the files' likelihood

  )
  if (MODEL_VERSION %in% c("free", "sv_only")) {
    base$sv <- 0.3 + runif(1, 0, 0.1)
  }
  if (MODEL_VERSION == "free") {
    base$st0_raw <- 0.1 + runif(1, 0, 0.05)
  }
  # "zero" has no sv/st0 parameters at all
  # "sv_only" has sv but not st0_raw
  # none of the three versions has an sw parameter (0)
  base
}


## COMPILE AND SAMPLE
mod <- cmdstan_model(stan_file, cpp_options = list(stan_threads = TRUE))

fit <- mod$sample(
  data = stan_data,
  chains = CHAINS,
  parallel_chains = PARALLEL_CHAINS,
  threads_per_chain = THREADS_PER_CHAIN,
  iter_warmup = WARMUP,
  iter_sampling = SAMPLING,
  max_treedepth = TREEDEPTH,
  adapt_delta = ADAPT_DELTA,
  refresh = REFRESH,
  init = init_fun,
  save_warmup = TRUE,   # so warmup progress is recoverable if this crashes
  output_dir = out_dir
)

fit$save_object(file.path(out_dir, "fit.rds"))


## DIAGNOSTICS
# variable list depends on model version 
group_pars <- c("beta_nu[1]", "beta_nu[2]", "beta_nu[3]", "beta_nu[4]",
                "beta_alpha1", "beta_alpha2", "beta_alpha3", "beta_tau",
                "sigma_nu_c", "sigma_nu_r", "sigma_nu_cr",
                "sigma_alpha_intercept", "sigma_tau_intercept")
if (MODEL_VERSION %in% c("free", "sv_only")) {
  group_pars <- c(group_pars, "sv")
}
if (MODEL_VERSION == "free") {
  group_pars <- c(group_pars, "st0")
}

summ <- fit$summary(variables = group_pars)
print(summ, n = Inf)

cat("\n=== Sampler diagnostics ===\n")
diag <- fit$diagnostic_summary()
print(diag)

n_total_transitions <- SAMPLING * CHAINS
# diagnostic_summary() only counts post-warmup sampling transitions
cat(sprintf("\nMax treedepth hit rate: %d / %d (%.0f%%)\n",
            sum(diag$num_max_treedepth), n_total_transitions,
            100 * sum(diag$num_max_treedepth) / n_total_transitions))
cat(sprintf("Divergences: %d\n", sum(diag$num_divergent)))
cat(sprintf("Max Rhat: %.4f\n", max(summ$rhat, na.rm = TRUE)))
cat(sprintf("Min ESS bulk: %.0f\n", min(summ$ess_bulk, na.rm = TRUE)))

cat(sprintf(
  "\nRun complete: MODEL_VERSION=%s, SCALE=%s, I=%d, N=%d, %d+%d iterations.\n",
  MODEL_VERSION, SCALE, I, N, WARMUP, SAMPLING
))


## figures and estimates saving 

library(cmdstanr)
library(posterior)
library(dplyr)
library(tidyr)
library(ggplot2)

MODEL_VERSION <- "sv_only"   # must match whichever run produced the fit below
SCALE         <- "full"      # must match whichever run produced the fit below
FIT_PATH      <- "output_ddm_v2_sv_only_full/fit.rds"
DATA_PATH     <- "krakow_data_standardized.csv"

OUT <- "output_ddm_v2_sv_only_full/figures"
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)

message("Loading fit from disk: ", FIT_PATH)
fit <- readRDS(FIT_PATH)
stopifnot(inherits(fit, "CmdStanFit"))  # checks whether if this is an RStan object instead

# regenerate the same participant/trial subset used to produce this fit
# using the same seed and fraction logic as before
# must match that script's section above or I/min_rt will be wrong
set.seed(42)
FRAC_PARTICIPANTS <- switch(SCALE, quicktest = 0.07, tryout = 0.30, real = 0.50, full = 1.00)
FRAC_TRIALS       <- if (MODEL_VERSION == "sv_only" && SCALE == "tryout") {
  0.30  # your explicit sv_only tryout request
} else {
  switch(SCALE, quicktest = 0.25, tryout = 0.50, real = 0.65, full = 1.00)
}

raw_full <- read.csv(DATA_PATH) %>%
  transmute(participant_orig = participant_index, condition = condition,
            resp_type = pre_acc, acc = acc, rt = rt)

all_ids <- unique(raw_full$participant_orig)
n_keep  <- round(FRAC_PARTICIPANTS * length(all_ids))
keep_ids <- sample(all_ids, n_keep)
raw <- raw_full %>% filter(participant_orig %in% keep_ids)

stratified_fraction <- function(df, frac) {
  df %>%
    group_by(participant_orig, condition, resp_type) %>%
    group_modify(~ {
      n_cell <- nrow(.x)
      target <- max(1, round(frac * n_cell))
      if (n_cell <= target) return(.x)
      .x[sample.int(n_cell, target), , drop = FALSE]
    }) %>%
    ungroup()
}
raw <- stratified_fraction(raw, FRAC_TRIALS)
raw <- raw %>% mutate(pid = as.integer(factor(participant_orig, levels = unique(participant_orig))))

I <- max(raw$pid)
min_rt <- raw %>% group_by(pid) %>% summarise(min_rt = min(rt), .groups = "drop") %>%
  arrange(pid) %>% pull(min_rt)

message(sprintf("Regenerated data subset: I=%d participants", I))

# Sanity check: does this I actually match the fit's own dimensions? 
# if not the config above (MODEL_VERSION/SCALE/seed) doesn't match FIT_PATH and everything would be wrong 
fit_I <- ncol(fit$draws(variables = "alpha_intercept", format = "matrix"))
if (I != fit_I) {
  stop(sprintf(
    "MISMATCH: regenerated data has I=%d participants, but the fit itself has I=%d.\n",
    I, fit_I),
    "MODEL_VERSION/SCALE above don't match FIT_PATH.\n"
  )
}


## POSTERIORS 
# as.matrix()/as.numeric() to get rid of the `posterior` package's draws_matrix class 
# otherwise error 

draws_beta_nu    <- as.matrix(fit$draws(variables = "beta_nu",     format = "matrix"))  # [draws x 4]
draws_nu_c       <- as.matrix(fit$draws(variables = "nu_c",        format = "matrix"))  # [draws x I]
draws_nu_r       <- as.matrix(fit$draws(variables = "nu_r",        format = "matrix"))
draws_nu_cr      <- as.matrix(fit$draws(variables = "nu_cr",       format = "matrix"))
draws_alpha_int  <- as.matrix(fit$draws(variables = "alpha_intercept", format = "matrix"))
draws_alpha2     <- as.numeric(fit$draws(variables = "beta_alpha2", format = "matrix")[, 1])  # scalar per draw
draws_alpha3     <- as.numeric(fit$draws(variables = "beta_alpha3", format = "matrix")[, 1])
draws_t0         <- as.matrix(fit$draws(variables = "t0",          format = "matrix"))  # [draws x I]

n_draws <- nrow(draws_beta_nu)


##   k=1: cond= 1, resp= 1   k=2: cond= 1, resp=-1
##   k=3: cond=-1, resp= 1   k=4: cond=-1, resp=-1

cond_k <- c( 1,  1, -1, -1)
resp_k <- c( 1, -1,  1, -1)
CELL_LABELS <- c("cong/post-corr", "cong/post-err", "incong/post-corr", "incong/post-err")

# nu[i,k] draws: intercept + person's slopes x cell design (natural scale)
nu_draws <- array(NA_real_, dim = c(n_draws, I, 4))
for (k in 1:4) {
  nu_draws[, , k] <- matrix(draws_beta_nu[, 1], nrow = n_draws, ncol = I) +
    draws_nu_c  * cond_k[k] +
    draws_nu_r  * resp_k[k] +
    draws_nu_cr * (cond_k[k] * resp_k[k])
}

# alpha[i,k] draws: exp(person's log-intercept + fixed effects x cell design)
# exponentiated so the output is on the natural scale and comparable to EZ's alpha 
alpha_draws <- array(NA_real_, dim = c(n_draws, I, 4))
for (k in 1:4) {
  log_alpha_k <- draws_alpha_int + draws_alpha2 * resp_k[k] + draws_alpha3 * cond_k[k]
  alpha_draws[, , k] <- exp(log_alpha_k)
}

# tau: t0 is already in natural seconds 
tau_draws <- draws_t0  # [draws x I]


## POPULATION-LEVEL ESTIMATES CSV 
group_pars <- c("beta_nu[1]", "beta_nu[2]", "beta_nu[3]", "beta_nu[4]",
                "beta_alpha1", "beta_alpha2", "beta_alpha3", "beta_tau",
                "sigma_nu_c", "sigma_nu_r", "sigma_nu_cr",
                "sigma_alpha_intercept", "sigma_tau_intercept")
if ("sv" %in% fit$metadata()$stan_variables) group_pars <- c(group_pars, "sv")
if ("st0" %in% fit$metadata()$stan_variables) group_pars <- c(group_pars, "st0")

summ <- fit$summary(
  variables = group_pars,
  mean = mean, sd = sd,
  ~quantile(.x, probs = c(0.025, 0.975)),
  rhat = rhat, ess_bulk = ess_bulk
)
names(summ)[names(summ) == "2.5%"]  <- "lo95"
names(summ)[names(summ) == "97.5%"] <- "hi95"

write.csv(summ, file.path(OUT, "ddm_estimates.csv"), row.names = FALSE)
message("Saved population-level estimates to ddm_estimates.csv")


# PERSON x CELL ESTIMATES CSV 
person_cell_df <- bind_rows(lapply(1:4, function(k) {
  data.frame(
    participant = 1:I,
    cell = k,
    cell_label = CELL_LABELS[k],
    nu_mean    = colMeans(nu_draws[, , k]),
    nu_lo95    = apply(nu_draws[, , k], 2, quantile, 0.025),
    nu_hi95    = apply(nu_draws[, , k], 2, quantile, 0.975),
    alpha_mean = colMeans(alpha_draws[, , k]),
    alpha_lo95 = apply(alpha_draws[, , k], 2, quantile, 0.025),
    alpha_hi95 = apply(alpha_draws[, , k], 2, quantile, 0.975)
  )
}))
person_tau_df <- data.frame(
  participant = 1:I,
  tau_mean = colMeans(tau_draws),
  tau_lo95 = apply(tau_draws, 2, quantile, 0.025),
  tau_hi95 = apply(tau_draws, 2, quantile, 0.975)
)

write.csv(person_cell_df, file.path(OUT, "ddm_person_cell_estimates.csv"), row.names = FALSE)
write.csv(person_tau_df,  file.path(OUT, "ddm_person_tau_estimates.csv"),  row.names = FALSE)
message("Saved person x cell estimates to ddm_person_cell_estimates.csv / ddm_person_tau_estimates.csv")

## PLOTS 
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
sv_plot <- function(p, nm, w = 10, h = 7)
  ggsave(file.path(OUT, nm), p, width = w, height = h, dpi = 300, bg = "#1e1e1e")

# population effects forest plot 
# note: nu here is natural, alpha is log scale, tau (beta_tau) is logit of min_rt scale
# they do not compare across types (nu vs alpha vs tau) so only comparable within a type or against the CI crossing zero
grp_vars   <- c("beta_nu[1]","beta_nu[2]","beta_nu[3]","beta_nu[4]",
                "beta_alpha1","beta_alpha2","beta_alpha3","beta_tau")
grp_labels <- c("mu_nu (intercept)","b_nu_condition","b_nu_resp_type","b_nu_interaction",
                "mu_alpha (log-scale intercept)","b_alpha_resp_type (log)","b_alpha_condition (log)",
                "mu_tau (logit-of-min_rt scale)")
grp_types  <- c("nu","nu","nu","nu","alpha","alpha","alpha","tau")

grp_summ <- fit$summary(variables = grp_vars, mean = mean,
                        ~quantile(.x, probs = c(0.125, 0.875)),
                        ~quantile(.x, probs = c(0.025, 0.975)))
grp_df <- data.frame(
  param = grp_labels, type = grp_types,
  mean = grp_summ$mean,
  lo50 = grp_summ$`12.5%`, hi50 = grp_summ$`87.5%`,
  lo95 = grp_summ$`2.5%`,  hi95 = grp_summ$`97.5%`
)

p_forest_group <- ggplot(grp_df, aes(x = mean, xmin = lo95, xmax = hi95,
                                     y = reorder(param, mean), colour = type)) +
  geom_vline(xintercept = 0, colour = "white", linetype = "dashed") +
  geom_errorbarh(aes(xmin = lo95, xmax = hi95), height = 0, linewidth = 1.2, alpha = 0.4) +
  geom_errorbarh(aes(xmin = lo50, xmax = hi50), height = 0, linewidth = 2.5, alpha = 0.8) +
  geom_point(size = 3) +
  scale_colour_manual(values = c(nu = "steelblue", alpha = "#e07070", tau = "goldenrod")) +
  labs(title = "Population effects full DDM",
       subtitle = "Thick bar = 50% CI | thin bar = 95% CI. note: nu/alpha/tau on different scales",
       x = "Posterior mean", y = NULL, colour = "Parameter") +
  dark
sv_plot(p_forest_group, "ddm_forest_group.png", h = 6)

# random-effect SDs forest plot 
sd_vars   <- c("sigma_nu_c","sigma_nu_r","sigma_nu_cr","sigma_alpha_intercept","sigma_tau_intercept")
sd_labels <- c("sigma_nu_cond (slope)","sigma_nu_resp (slope)","sigma_nu_cr (slope)",
               "sigma_alpha (intercept)","sigma_tau (intercept)")
sd_types  <- c("slope","slope","slope","intercept","intercept")

sd_summ <- fit$summary(variables = sd_vars, mean = mean,
                       ~quantile(.x, probs = c(0.125, 0.875)),
                       ~quantile(.x, probs = c(0.025, 0.975)))
sd_df <- data.frame(
  param = sd_labels, type = sd_types,
  mean = sd_summ$mean, lo50 = sd_summ$`12.5%`, hi50 = sd_summ$`87.5%`,
  lo95 = sd_summ$`2.5%`, hi95 = sd_summ$`97.5%`
)

p_forest_sds <- ggplot(sd_df, aes(x = mean, xmin = lo95, xmax = hi95,
                                  y = reorder(param, mean), colour = type)) +
  geom_vline(xintercept = 0, colour = "white", linetype = "dashed") +
  geom_errorbarh(aes(xmin = lo95, xmax = hi95), height = 0, linewidth = 1.2, alpha = 0.4) +
  geom_errorbarh(aes(xmin = lo50, xmax = hi50), height = 0, linewidth = 2.5, alpha = 0.8) +
  geom_point(size = 3) +
  scale_colour_manual(values = c(slope = "steelblue", intercept = "seagreen")) +
  labs(title = "Between-person SDs -- full DDM",
       subtitle = "Thick bar = 50% CI, thin bar = 95% CI. No EZ-style residual SD exists in this model (see caveats).",
       x = "Posterior mean", y = NULL, colour = NULL) +
  dark
sv_plot(p_forest_sds, "ddm_forest_sds.png", h = 6)

# individual (caterpillar) plots 
make_caterpillar <- function(values_matrix, pop_mean, title_str, col) {
  df <- data.frame(
    mean = colMeans(values_matrix),
    lo95 = apply(values_matrix, 2, quantile, 0.025),
    hi95 = apply(values_matrix, 2, quantile, 0.975)
  ) %>% arrange(mean) %>% mutate(rank = seq_len(I))
  
  ggplot(df, aes(x = rank, y = mean, ymin = lo95, ymax = hi95)) +
    geom_hline(yintercept = 0,        colour = "white",    linetype = "dashed") +
    geom_hline(yintercept = pop_mean, colour = "goldenrod", linewidth = 0.8) +
    geom_errorbar(colour = col, width = 0, alpha = 0.4) +
    geom_point(size = 0.8, colour = col) +
    labs(title = title_str, subtitle = "Sorted by posterior mean, gold = population mean",
         x = "Participant rank", y = "Value") +
    dark
}

pm <- function(vname) fit$summary(variables = vname, mean = mean)$mean

sv_plot(make_caterpillar(draws_nu_c,  pm("beta_nu[2]"), "Condition slope on nu",   "steelblue"),
        "ddm_caterpillar_nu_c.png")
sv_plot(make_caterpillar(draws_nu_r,  pm("beta_nu[3]"), "Resp_type slope on nu",   "steelblue"),
        "ddm_caterpillar_nu_r.png")
sv_plot(make_caterpillar(draws_nu_cr, pm("beta_nu[4]"), "Interaction slope on nu", "steelblue"),
        "ddm_caterpillar_nu_cr.png")

# alpha intercept and tau plotted on natural scale for interpretability
alpha_int_natural <- exp(draws_alpha_int)  # exp() since alpha_intercept is log-scale
sv_plot(make_caterpillar(alpha_int_natural, exp(pm("beta_alpha1")),
                         "Individual intercepts, alpha (natural scale)", "seagreen"),
        "ddm_caterpillar_alpha_intercept.png")
sv_plot(make_caterpillar(tau_draws, mean(colMeans(tau_draws)),
                         "Individual non-decision times, tau (seconds)", "seagreen"),
        "ddm_caterpillar_tau.png")

# cell-level violin plots (nu and alpha only)
make_cell_df <- function(draws_arr, par_name) {
  bind_rows(lapply(1:4, function(k) {
    data.frame(
      value = colMeans(draws_arr[, , k]),
      cell = k,
      cell_label = factor(CELL_LABELS[k], levels = CELL_LABELS)
    )
  }))
}
CELL_COLOURS <- c("goldenrod", "steelblue", "#e07070", "#a070e0")

plot_cells <- function(df, par_name) {
  ggplot(df, aes(x = cell_label, y = value, fill = cell_label, colour = cell_label)) +
    geom_violin(alpha = 0.25, linewidth = 0.4) +
    geom_jitter(width = 0.15, size = 0.5, alpha = 0.4) +
    stat_summary(fun = mean, geom = "point", size = 3, colour = "white", shape = 18) +
    scale_fill_manual(values = CELL_COLOURS) +
    scale_colour_manual(values = CELL_COLOURS) +
    labs(title = paste0("Cell-level, ", par_name, " posterior means across participants -- full DDM"),
         subtitle = "Diamond = cell mean, each dot = one participant",
         x = NULL, y = par_name) +
    dark + theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1))
}

sv_plot(plot_cells(make_cell_df(nu_draws,    "nu"),    "nu"),    "ddm_cells_nu.png")
sv_plot(plot_cells(make_cell_df(alpha_draws, "alpha"), "alpha"), "ddm_cells_alpha.png")

cat(sprintf("\nAll estimates and plots saved to: %s\n", OUT))
