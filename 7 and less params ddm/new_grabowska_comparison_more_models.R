
library(cmdstanr)
library(posterior)
library(dplyr)
library(tidyr)

set.seed(42)

MODEL_VERSION <- "free_wfree_st0fix"   
SCALE         <- "full"  


# free_scalars : welche der vier DDM-Parameter (w, sv, sw, st0) im Modell frei geschätzt werden (Rest ist im Stan-Code fest fixiert)
# st0_percapita: TRUE bei den Modellen, die personen-spezifische st0[p] = st0_raw * 2 * t0[p] Reparametrisierung nutzen
# (kein globales min(t0) mehr, um das Problem mit dem min-RT von allen zu beheben)
# random : TRUE nur bei der Variante in der sv selbst hierarchisch (mit Random Intercept pro Person) statt als einen 
# population Parameter geschätzt wird 

MODEL_REGISTRY <- list(
  zero               = list(stan_file = "hierarchical_fullddm_no_intertrial_variability_nu_alpha.stan",
                             free_scalars = c(),                    random = FALSE),
  sv_only            = list(stan_file = "hierarchical_fullddm_sv_only_nu_alpha.stan",
                             free_scalars = c("sv"),                random = FALSE),
  free               = list(stan_file = "hierarchical_fullddm_intertrial_variability_nu_alpha.stan",
                             free_scalars = c("sv","st0"),          random = FALSE),
  sv_wfree_swfree    = list(stan_file = "hierarchical_fullddm_sv_wfree_swfree_nu_alpha.stan",
                             free_scalars = c("w","sv","sw"),       random = FALSE),
  sv_wfree           = list(stan_file = "hierarchical_fullddm_sv_wfree_nu_alpha.stan",
                             free_scalars = c("w","sv"),            random = FALSE),
  free_wfree_swfree  = list(stan_file = "hierarchical_fullddm_free_wfree_swfree_nu_alpha.stan",
                             free_scalars = c("w","sv","sw","st0"), random = FALSE),
  
  free_wfree_swfree_st0fix = list(stan_file = "hierarchical_fullddm_free_wfree_swfree_st0fix_nu_alpha.stan",
                             free_scalars = c("w","sv","sw"), st0_percapita = TRUE, random = FALSE),
  free_wfree         = list(stan_file = "hierarchical_fullddm_free_wfree_nu_alpha.stan",
                             free_scalars = c("w","sv","st0"),      random = FALSE),

  free_wfree_st0fix  = list(stan_file = "hierarchical_fullddm_free_wfree_st0fix_nu_alpha.stan",
                             free_scalars = c("w","sv"), st0_percapita = TRUE, random = FALSE),
  zero_wfree_swfree  = list(stan_file = "hierarchical_fullddm_zero_wfree_swfree_nu_alpha.stan",
                             free_scalars = c("w","sw"),            random = FALSE),
  zero_wfree         = list(stan_file = "hierarchical_fullddm_zero_wfree_nu_alpha.stan",
                             free_scalars = c("w"),                 random = FALSE),
  sv_random          = list(stan_file = "hierarchical_fullddm_sv_random_intercept_nu_alpha.stan",
                             free_scalars = c(),                    random = TRUE)
)
reg <- MODEL_REGISTRY[[MODEL_VERSION]]
stan_file <- reg$stan_file

# Übersetzung in die stan-Variablennamen/Datennamen/Defaults
# sw und st0 werde im Stan-Code nicht direkt geschätzt, sondern über die raw-Reparametrisierung(sw_raw, st0_raw) 
PARAM_STAN_NAME <- list(w = "w", sv = "sv", sw = "sw_raw", st0 = "st0_raw")
PRIOR_DATA_NAME <- list(w = "prior_w", sv = "prior_sv", sw = "prior_sw", st0 = "prior_st0")
DEFAULT_PRIOR   <- list(w = c(0.5, 0.1), sv = c(0.5, 0.5), sw = c(0.3, 0.2), st0 = c(0.3, 0.2))
DEFAULT_INIT    <- list(w = 0.5, sv = 0.3, sw = 0.1, st0 = 0.1)


SCALE_CONFIG <- list(
  zero    = list(quicktest = list(fp=0.07, ft=0.25, wu=50,  sp=30),
                  tryout    = list(fp=0.30, ft=0.50, wu=500, sp=300),
                  real      = list(fp=0.50, ft=0.65, wu=3000, sp=3000),
                  full      = list(fp=1.00, ft=1.00, wu=3000, sp=3000)),
  sv_only = list(quicktest = list(fp=0.07, ft=0.25, wu=50,  sp=30),
                  tryout    = list(fp=0.30, ft=0.30, wu=300, sp=200),
                  real      = list(fp=0.50, ft=0.65, wu=2000, sp=3000),
                  full      = list(fp=1.00, ft=1.00, wu=2000, sp=3000)),
  free    = list(quicktest = list(fp=0.07, ft=0.25, wu=50,  sp=30),
                  tryout    = list(fp=0.30, ft=0.50, wu=100, sp=60),
                  real      = list(fp=0.50, ft=0.65, wu=300, sp=200),
                  full      = list(fp=1.00, ft=1.00, wu=300, sp=200))
)

CHEAP_GUESS <- list(quicktest = list(fp=0.07, ft=0.25, wu=50,  sp=30),
                     tryout    = list(fp=0.30, ft=0.30, wu=300, sp=200),
                     real      = list(fp=0.50, ft=0.65, wu=2000, sp=2000),
                     full      = list(fp=1.00, ft=1.00, wu=2000, sp=2000))
EXPENSIVE_GUESS <- list(quicktest = list(fp=0.07, ft=0.25, wu=50, sp=30),
                         tryout    = list(fp=0.20, ft=0.20, wu=80, sp=50),
                         real      = list(fp=0.30, ft=0.30, wu=100, sp=60),
                         full      = list(fp=0.30, ft=0.30, wu=100, sp=60))


# eigene kleinschrittigere Stufen für sv_wfree_swfree 
SCALE_CONFIG$sv_wfree_swfree <- list(
  quicktest = list(fp=0.07, ft=0.25, wu=50,  sp=30),
  tryout    = list(fp=0.30, ft=0.30, wu=200, sp=150),   
  real      = list(fp=0.50, ft=0.50, wu=500, sp=400),
  full      = list(fp=1.00, ft=0.65, wu=500, sp=400)
)

# zuordnung der restlichen Modelle zu einer der beiden Vorlagen
# je nachdem wie rechenintensiv sie sind (vermutlich)
SCALE_CONFIG$sv_wfree          <- CHEAP_GUESS      
SCALE_CONFIG$free_wfree_swfree <- EXPENSIVE_GUESS  
SCALE_CONFIG$free_wfree        <- EXPENSIVE_GUESS  
SCALE_CONFIG$free_wfree_st0fix <- EXPENSIVE_GUESS

SCALE_CONFIG$free_wfree_swfree_st0fix <- EXPENSIVE_GUESS
SCALE_CONFIG$zero_wfree_swfree <- EXPENSIVE_GUESS  
SCALE_CONFIG$zero_wfree        <- CHEAP_GUESS      
SCALE_CONFIG$sv_random         <- CHEAP_GUESS      

# Config für den aktiven lauf
cfg <- SCALE_CONFIG[[MODEL_VERSION]][[SCALE]]
FRAC_PARTICIPANTS <- cfg$fp
FRAC_TRIALS       <- cfg$ft
WARMUP            <- cfg$wu
SAMPLING          <- cfg$sp

# sampler settings
TREEDEPTH   <- 10
ADAPT_DELTA <- 0.8
REFRESH     <- 10

CHAINS          <- 2
PARALLEL_CHAINS <- 2

# threads of verfügare cores verteilen 
N_THREADS <- max(1, parallel::detectCores() - 1)
THREADS_PER_CHAIN <- max(1, floor(N_THREADS / PARALLEL_CHAINS))

DATA_PATH <- "krakow_data_standardized.csv"
out_dir <- sprintf("output_ddm_v2_%s_%s", MODEL_VERSION, SCALE)
if (!dir.exists(out_dir)) dir.create(out_dir)

message(sprintf("Config: MODEL_VERSION=%s, SCALE=%s, stan_file=%s, out_dir=%s",
                 MODEL_VERSION, SCALE, stan_file, out_dir))
message(sprintf("FRAC_PARTICIPANTS=%.2f, FRAC_TRIALS=%.2f, WARMUP=%d, SAMPLING=%d",
                 FRAC_PARTICIPANTS, FRAC_TRIALS, WARMUP, SAMPLING))


# load data und alle studen samplen die nicht das full sample sind
raw_full <- read.csv(DATA_PATH) %>%
  transmute(participant_orig = participant_index, condition = condition,
            resp_type = pre_acc, acc = acc, rt = rt)

# Sanity-Checks auf erwartete Kodierung
stopifnot(all(raw_full$acc %in% c(-1, 1)), all(raw_full$condition %in% c(-1, 1)),
          all(raw_full$resp_type %in% c(-1, 1)), all(raw_full$rt > 0))

# stude 1 sub sampling (zufällige personen)
all_ids <- unique(raw_full$participant_orig)
n_keep  <- round(FRAC_PARTICIPANTS * length(all_ids))
keep_ids <- sample(all_ids, n_keep)
raw <- raw_full %>% filter(participant_orig %in% keep_ids)

# stufe 2 subsampling: innerhalb jeder person stratifiziert 
# pro condition x resp_type-Zelle einen Bruchteil (ft) der Trials ziehen
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
# fortlaufende personen IDs für stan
raw <- raw %>% mutate(pid = as.integer(factor(participant_orig, levels = unique(participant_orig))))

I <- max(raw$pid)
N <- nrow(raw)
# min_rt pro Person aus den (ggf subgesampleten) Daten 
min_rt <- raw %>% group_by(pid) %>% summarise(min_rt = min(rt), .groups = "drop") %>%
  arrange(pid) %>% pull(min_rt)

message(sprintf("Final dataset: I=%d participants, N=%d trials", I, N))


# Priors: gemeinsam fuer alle Modelle der Registry --> nu/alpha/tau- Effektstruktur ist ja in jedem Modell gleich
# nur die w/sv/sw/st0 unterscheiden sich je nach Modell

# # prior_beta_nu: Zeile 1 = Means, Zeile 2 = SDs, Spalten = [intercept,
# condition, resp_type, interaction]
prior_beta_nu <- rbind(c(0, 0, 0, 0), c(2, 2, 0.5, 0.5))
prior_beta_alpha1 <- c(0, 0.5)
prior_beta_alpha2 <- c(0, 0.2)
prior_beta_alpha3 <- c(0, 0.3)
prior_sigma_nu_c  <- c(1, 1)
prior_sigma_nu_r  <- c(1, 1)
prior_sigma_nu_cr <- c(1, 1)
prior_sigma_alpha_intercept <- c(1, 1)
prior_beta_tau <- c(0.4, 0.5)
prior_sigma_tau_intercept <- c(0.3, 1)

# Priorsfür die sv_random-Variante
# sv hierarchisch, auf log-Skala, deshalb der negative Mean
# entspricht exp(-2) ~ 0.14 --> sv-Baseline natürliche skala
prior_beta_sv <- c(-2, 1)              
prior_sigma_sv_intercept <- c(1, 1)    


stan_data <- list(
  N = N, I = I, pid = raw$pid, rt = raw$rt, acc = as.integer(raw$acc),
  condition = raw$condition, resp_type = raw$resp_type, min_rt = min_rt,
  grainsize = 1L,
  prior_beta_nu = prior_beta_nu,
  prior_beta_alpha1 = prior_beta_alpha1, prior_beta_alpha2 = prior_beta_alpha2,
  prior_beta_alpha3 = prior_beta_alpha3, prior_beta_tau = prior_beta_tau,
  prior_sigma_nu_c = prior_sigma_nu_c, prior_sigma_nu_r = prior_sigma_nu_r,
  prior_sigma_nu_cr = prior_sigma_nu_cr,
  prior_sigma_alpha_intercept = prior_sigma_alpha_intercept,
  prior_sigma_tau_intercept = prior_sigma_tau_intercept
)

# nur für die in diese Modell freien Parameter die passende Prior-Data-Variable ergänzen
for (nm in reg$free_scalars) {
  stan_data[[PRIOR_DATA_NAME[[nm]]]] <- DEFAULT_PRIOR[[nm]]
}
# st0_percapita-Modelle brauchen den st0-Prior auch wenn "st0"  icht in free_scalars auftaucht 
# (weil dort z.B. nur "w","sv","sw" gelistet sind, aber st0_raw trotzdem als eigener Parameter existiert)
if (isTRUE(reg$st0_percapita)) {
  stan_data$prior_st0 <- DEFAULT_PRIOR$st0  
}
if (reg$random) {
  stan_data$prior_beta_sv <- prior_beta_sv
  stan_data$prior_sigma_sv_intercept <- prior_sigma_sv_intercept
}


# init_fun: eigene Startwerte statt Stan Defaults
# wichtig bei den bounded Parametern w/sw_raw/st0_raw nahe 0 oder 1, wo die Wiener-Likelihood numerisch instabil werden kann
# jede Chain bekommt trotzdem eigenes Jitter um den Default, damit die Chains nicht alle exakt am selben Punkt starten
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
  )
  # nur für die in diesem Modell freien Parameter Startwerte setzen
  
  # w bekommt ein kleineres symmetrisches Jitter um 0.5, weil w auf [0,1] bounded ist und zu
  # grosses Jitter schnell Richtung Rand rutschen kann
  for (nm in reg$free_scalars) {
    jitter <- if (nm %in% c("w")) runif(1, -0.02, 0.02) else runif(1, 0, 0.05)
    base[[PARAM_STAN_NAME[[nm]]]] <- DEFAULT_INIT[[nm]] + jitter
  }
  if (isTRUE(reg$st0_percapita)) {
    base$st0_raw <- DEFAULT_INIT$st0 + runif(1, 0, 0.05)
  }
  if (reg$random) {
    base$beta_sv <- -2 + rnorm(1, 0, 0.1)
    base$sigma_sv_intercept <- 0.3 + runif(1, 0, 0.05)
    base$z_sv <- rnorm(I, 0, 0.1)
  }
  base
}


# modell comp and fit

mod <- cmdstan_model(stan_file, cpp_options = list(stan_threads = TRUE))

fit <- mod$sample(
  data = stan_data, chains = CHAINS, parallel_chains = PARALLEL_CHAINS,
  threads_per_chain = THREADS_PER_CHAIN, iter_warmup = WARMUP, iter_sampling = SAMPLING,
  max_treedepth = TREEDEPTH, adapt_delta = ADAPT_DELTA, refresh = REFRESH,
  init = init_fun, save_warmup = TRUE, output_dir = out_dir
)
fit$save_object(file.path(out_dir, "fit.rds"))


# diagnostics: Populations-Parameter-Summary & sampler check 

group_pars <- c("beta_nu[1]", "beta_nu[2]", "beta_nu[3]", "beta_nu[4]",
                 "beta_alpha1", "beta_alpha2", "beta_alpha3", "beta_tau",
                 "sigma_nu_c", "sigma_nu_r", "sigma_nu_cr",
                 "sigma_alpha_intercept", "sigma_tau_intercept")
for (nm in reg$free_scalars) group_pars <- c(group_pars, nm)  
if (isTRUE(reg$st0_percapita)) group_pars <- c(group_pars, "st0_raw")  
if (reg$random) group_pars <- c(group_pars, "beta_sv", "sigma_sv_intercept")

summ <- fit$summary(variables = group_pars)
print(summ, n = Inf)

cat("\n Sampler diagnostics \n")
diag <- fit$diagnostic_summary()
print(diag)

n_total_transitions <- SAMPLING * CHAINS
cat(sprintf("\nMax treedepth hit rate: %d / %d (%.0f%%)\n",
            sum(diag$num_max_treedepth), n_total_transitions,
            100 * sum(diag$num_max_treedepth) / n_total_transitions))
cat(sprintf("Divergences: %d\n", sum(diag$num_divergent)))
cat(sprintf("Max Rhat: %.4f\n", max(summ$rhat, na.rm = TRUE)))
cat(sprintf("Min ESS bulk: %.0f\n", min(summ$ess_bulk, na.rm = TRUE)))

cat(sprintf("\nRun complete: MODEL_VERSION=%s, SCALE=%s, I=%d, N=%d, %d+%d iterations.\n",
            MODEL_VERSION, SCALE, I, N, WARMUP, SAMPLING))

library(posterior)
library(ggplot2)

OUT <- file.path(out_dir, "figures")
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)

# Draws extrahieren für Personen und Zellen
draws_beta_nu   <- as.matrix(fit$draws(variables = "beta_nu", format = "matrix"))
draws_nu_c      <- as.matrix(fit$draws(variables = "nu_c", format = "matrix"))
draws_nu_r      <- as.matrix(fit$draws(variables = "nu_r", format = "matrix"))
draws_nu_cr     <- as.matrix(fit$draws(variables = "nu_cr", format = "matrix"))
draws_alpha_int <- as.matrix(fit$draws(variables = "alpha_intercept", format = "matrix"))
draws_alpha2    <- as.numeric(fit$draws(variables = "beta_alpha2", format = "matrix")[, 1])
draws_alpha3    <- as.numeric(fit$draws(variables = "beta_alpha3", format = "matrix")[, 1])
draws_t0        <- as.matrix(fit$draws(variables = "t0", format = "matrix"))

n_draws <- nrow(draws_beta_nu)

# die gleiche 4-Zellen-Kodierung wie im EZ-Modell
# k=1 cong/post-corr, k=2 cong/post-err, k=3 incong/post-corr, k=4 incong/post-err
cond_k <- c( 1,  1, -1, -1)
resp_k <- c( 1, -1,  1, -1)
CELL_LABELS <- c("cong/post-corr", "cong/post-err", "incong/post-corr", "incong/post-err")

# nu pro Draw x Person x Zelle: linearer Prädiktor direkt aus den Draws (nicht aus Posterior-Means!)
nu_draws <- array(NA_real_, dim = c(n_draws, I, 4))
for (k in 1:4) {
  nu_draws[, , k] <- matrix(draws_beta_nu[, 1], nrow = n_draws, ncol = I) +
    draws_nu_c  * cond_k[k] + draws_nu_r * resp_k[k] + draws_nu_cr * (cond_k[k] * resp_k[k])
}

# alpha pro Draw x Person x Zelle: exp() des log-linearen Praediktors
alpha_draws <- array(NA_real_, dim = c(n_draws, I, 4))
for (k in 1:4) {
  log_alpha_k <- draws_alpha_int + draws_alpha2 * resp_k[k] + draws_alpha3 * cond_k[k]
  alpha_draws[, , k] <- exp(log_alpha_k)
}

# # tau/t0 variiert nihct nach Zelle (nur nach Person)
tau_draws <- draws_t0


group_pars <- c("beta_nu[1]", "beta_nu[2]", "beta_nu[3]", "beta_nu[4]",
                "beta_alpha1", "beta_alpha2", "beta_alpha3", "beta_tau",
                "sigma_nu_c", "sigma_nu_r", "sigma_nu_cr",
                "sigma_alpha_intercept", "sigma_tau_intercept")
group_pars <- c(group_pars, reg$free_scalars)  
if (isTRUE(reg$st0_percapita)) group_pars <- c(group_pars, "st0_raw")  
if (reg$random) group_pars <- c(group_pars, "beta_sv", "sigma_sv_intercept")

# Populations-Level-Summary mit 95% CI
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


person_cell_df <- bind_rows(lapply(1:4, function(k) {
  data.frame(
    participant = 1:I, cell = k, cell_label = CELL_LABELS[k],
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


draws_sv_person <- NULL
if (reg$random) {
  draws_sv_person <- as.matrix(fit$draws(variables = "sv", format = "matrix"))
  person_sv_df <- data.frame(
    participant = 1:I,
    sv_mean = colMeans(draws_sv_person),
    sv_lo95 = apply(draws_sv_person, 2, quantile, 0.025),
    sv_hi95 = apply(draws_sv_person, 2, quantile, 0.975)
  )
  write.csv(person_sv_df, file.path(OUT, "ddm_person_sv_estimates.csv"), row.names = FALSE)
  message("Saved person-level sv estimates to ddm_person_sv_estimates.csv")
}


draws_st0_person <- NULL
if (isTRUE(reg$st0_percapita)) {
  draws_st0_person <- as.matrix(fit$draws(variables = "st0", format = "matrix"))
  person_st0_df <- data.frame(
    participant = 1:I,
    st0_mean = colMeans(draws_st0_person),
    st0_lo95 = apply(draws_st0_person, 2, quantile, 0.025),
    st0_hi95 = apply(draws_st0_person, 2, quantile, 0.975)
  )
  write.csv(person_st0_df, file.path(OUT, "ddm_person_st0_estimates.csv"), row.names = FALSE)
  message("Saved person-level st0 estimates to ddm_person_st0_estimates.csv")
}


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
  param = grp_labels, type = grp_types, mean = grp_summ$mean,
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
  labs(title = paste0("Population effects -- ", MODEL_VERSION),
       subtitle = "Thick bar = 50% CI | thin bar = 95% CI. note: nu/alpha/tau on different scales",
       x = "Posterior mean", y = NULL, colour = "Parameter") +
  dark
sv_plot(p_forest_group, "ddm_forest_group.png", h = 6)

# random-effect SDs forest 
sd_vars   <- c("sigma_nu_c","sigma_nu_r","sigma_nu_cr","sigma_alpha_intercept","sigma_tau_intercept")
sd_labels <- c("sigma_nu_cond (slope)","sigma_nu_resp (slope)","sigma_nu_cr (slope)",
               "sigma_alpha (intercept)","sigma_tau (intercept)")
sd_types  <- c("slope","slope","slope","intercept","intercept")
sd_summ <- fit$summary(variables = sd_vars, mean = mean,
                       ~quantile(.x, probs = c(0.125, 0.875)),
                       ~quantile(.x, probs = c(0.025, 0.975)))
sd_df <- data.frame(
  param = sd_labels, type = sd_types, mean = sd_summ$mean,
  lo50 = sd_summ$`12.5%`, hi50 = sd_summ$`87.5%`,
  lo95 = sd_summ$`2.5%`, hi95 = sd_summ$`97.5%`
)
p_forest_sds <- ggplot(sd_df, aes(x = mean, xmin = lo95, xmax = hi95,
                                  y = reorder(param, mean), colour = type)) +
  geom_vline(xintercept = 0, colour = "white", linetype = "dashed") +
  geom_errorbarh(aes(xmin = lo95, xmax = hi95), height = 0, linewidth = 1.2, alpha = 0.4) +
  geom_errorbarh(aes(xmin = lo50, xmax = hi50), height = 0, linewidth = 2.5, alpha = 0.8) +
  geom_point(size = 3) +
  scale_colour_manual(values = c(slope = "steelblue", intercept = "seagreen")) +
  labs(title = paste0("Between-person SDs -- ", MODEL_VERSION),
       subtitle = "Thick bar = 50% CI, thin bar = 95% CI. No EZ-style residual SD exists in this model.",
       x = "Posterior mean", y = NULL, colour = NULL) +
  dark
sv_plot(p_forest_sds, "ddm_forest_sds.png", h = 6)


extras_vars <- c(reg$free_scalars, if (isTRUE(reg$st0_percapita)) "st0_raw")
if (length(extras_vars) > 0) {
  extras_summ <- fit$summary(variables = extras_vars, mean = mean,
                             ~quantile(.x, probs = c(0.125, 0.875)),
                             ~quantile(.x, probs = c(0.025, 0.975)))
  extras_df <- data.frame(
    param = extras_vars, mean = extras_summ$mean,
    lo50 = extras_summ$`12.5%`, hi50 = extras_summ$`87.5%`,
    lo95 = extras_summ$`2.5%`, hi95 = extras_summ$`97.5%`
  )
  p_extras <- ggplot(extras_df, aes(x = mean, xmin = lo95, xmax = hi95, y = param)) +
    geom_vline(xintercept = 0, colour = "white", linetype = "dashed") +
    geom_errorbarh(aes(xmin = lo95, xmax = hi95), height = 0, linewidth = 1.2, alpha = 0.4, colour = "orchid") +
    geom_errorbarh(aes(xmin = lo50, xmax = hi50), height = 0, linewidth = 2.5, alpha = 0.8, colour = "orchid") +
    geom_point(size = 3, colour = "orchid") +
    labs(title = paste0("Free DDM-extra parameters ", MODEL_VERSION),
         subtitle = "w, sv, sw, st0(_raw)",
         x = "Posterior mean", y = NULL) +
    dark
  sv_plot(p_extras, "ddm_forest_extras.png", h = 4)
  
}


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

alpha_int_natural <- exp(draws_alpha_int)
sv_plot(make_caterpillar(alpha_int_natural, exp(pm("beta_alpha1")),
                         "Individual intercepts, alpha (natural scale)", "seagreen"),
        "ddm_caterpillar_alpha_intercept.png")
sv_plot(make_caterpillar(tau_draws, mean(colMeans(tau_draws)),
                         "Individual non-decision times, tau (seconds)", "seagreen"),
        "ddm_caterpillar_tau.png")

# sv_random version only: person-level sv caterpillar plot
if (reg$random) {
  sv_plot(make_caterpillar(draws_sv_person, exp(pm("beta_sv")),
                           "Individual trial-to-trial drift noise, sv (natural scale)", "orchid"),
          "ddm_caterpillar_sv.png")
}

# st0_percapita versions only: person-level st0 caterpillar plot
if (isTRUE(reg$st0_percapita)) {
  sv_plot(make_caterpillar(draws_st0_person, mean(colMeans(draws_st0_person)),
                           "Individual inter-trial non-decision variability, st0 (seconds)", "orchid"),
          "ddm_caterpillar_st0.png")
}

# cell-level violin plots (nu and alpha only, tau has no cell variation)
make_cell_df <- function(draws_arr) {
  bind_rows(lapply(1:4, function(k) {
    data.frame(value = colMeans(draws_arr[, , k]), cell = k,
               cell_label = factor(CELL_LABELS[k], levels = CELL_LABELS))
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
    labs(title = paste0("Cell-level, ", par_name, " posterior means  ", MODEL_VERSION),
         subtitle = "Diamond = cell mean, each dot = one participant",
         x = NULL, y = par_name) +
    dark + theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1))
}
sv_plot(plot_cells(make_cell_df(nu_draws),    "nu"),    "ddm_cells_nu.png")
sv_plot(plot_cells(make_cell_df(alpha_draws), "alpha"), "ddm_cells_alpha.png")

cat(sprintf("\nAll estimates and plots saved to: %s\n", OUT))
