# draw population parameters from the priors
# draw one person from that population
# simulate N_TRIALS trials for that person (half in each condition)
# repeat N_PERSONS times
# --> each simulated person is one draw from the prior predictive distribution for a single
# person (population + person-level uncertainty together)


# drift: nu[i, cond] = nu_intercept[i] + nu_c[i] * cond, cond = +1 / -1
# alpha: alpha[i] = exp(alpha_intercept[i]) (log scale, no condition effect)
# tau/t0: beta_tau + sigma_tau_intercept * z (seconds, no condition effect), truncated at 0
# all SD priors are normal truncated at 0
# w = 0.5, no trial-to-trial variability (sv = sw = st0 = 0)
# upper boundary = correct response (acc = 1)
# EZ likelihood: C ~ binomial(J, Pc); MRT ~ normal(MRT_pred, sqrt(VRT_pred / C))
# VRT ~ normal(VRT_pred, sqrt(2 * VRT_pred^2 / (C - 1))) (correct RTs only)

library(dplyr)
library(tidyr)
library(ggplot2)

PROJECT_DIR <- "/Users/annegiacobello/Desktop/EZ Paper "
OUT_DIR     <- file.path(PROJECT_DIR, "model complexity", "ez")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

N_PERSONS   <- 1000     # simulated persons
N_TRIALS    <- 1000     # trials per person, split evenly over the 2 conditions
COND_LEVELS <- c(1, -1)

set.seed(101)

## priors 
priors <- list(
  beta_nu_intercept     = c(2,         1),     # mean drift across conditions
  beta_nu_c             = c(0,         0.5),   # condition effect on drift (+1 / -1 coding)
  beta_alpha            = c(log(1.2),  0.3),   # log boundary separation
  beta_tau              = c(0.3,       0.1),   # non-decision time (s)
  sigma_nu_intercept    = c(0,         0.5),   # between-person SDs (half-normal)
  sigma_nu_c            = c(0,         0.3),
  sigma_alpha_intercept = c(0,         0.2),
  sigma_tau_intercept   = c(0,         0.05)
)

#### draws: population parameters, then one person per population draw

# normal truncated to [lower, upper] (inverse CDF), vectorised
rtnorm <- function(n, mean, sd, lower = -Inf, upper = Inf) {
  p <- runif(n, pnorm(lower, mean, sd), pnorm(upper, mean, sd))
  qnorm(p, mean, sd)
}
nd <- function(nm, lower = -Inf, upper = Inf) {
  rtnorm(N_PERSONS, priors[[nm]][1], priors[[nm]][2], lower, upper)
}

persons <- data.frame(
  person                = seq_len(N_PERSONS),
  beta_nu_intercept     = nd("beta_nu_intercept"),
  beta_nu_c             = nd("beta_nu_c"),
  beta_alpha            = nd("beta_alpha"),
  beta_tau              = nd("beta_tau",              lower = 0),
  sigma_nu_intercept    = nd("sigma_nu_intercept",    lower = 0),
  sigma_nu_c            = nd("sigma_nu_c",            lower = 0),
  sigma_alpha_intercept = nd("sigma_alpha_intercept", lower = 0),
  sigma_tau_intercept   = nd("sigma_tau_intercept",   lower = 0)
) %>%
  mutate(
    nu_intercept    = beta_nu_intercept + sigma_nu_intercept    * rnorm(N_PERSONS),
    nu_c            = beta_nu_c         + sigma_nu_c            * rnorm(N_PERSONS),
    alpha_intercept = beta_alpha        + sigma_alpha_intercept * rnorm(N_PERSONS),
    alpha           = exp(alpha_intercept),
    tau             = rtnorm(N_PERSONS, beta_tau, sigma_tau_intercept, lower = 0),
    nu_cond1        = nu_intercept + nu_c * COND_LEVELS[1],
    nu_cond2        = nu_intercept + nu_c * COND_LEVELS[2]
  )

#### trial-level simulation

# vectorised Euler simulation of the Wiener first-passage time 
# boundaries moved inward by 0.5826 * sqrt(dt) (Broadie-Glasserman continuity correction)
rwiener_euler <- function(v, a, z, dt = 1e-3, t_max = 20) {
  n     <- length(v)
  shift <- 0.5826 * sqrt(dt)
  lo    <- rep(shift, n); hi <- a - shift
  x     <- pmin(pmax(z, lo + 1e-9), hi - 1e-9)
  dec_t <- rep(NA_real_, n); resp <- rep(NA_integer_, n)
  idx <- seq_len(n); xa <- x; va <- v * dt; la <- lo; ha <- hi
  sq <- sqrt(dt); step <- 0L; n_max <- ceiling(t_max / dt)
  while (length(idx) > 0 && step < n_max) {
    step <- step + 1L
    xa   <- xa + va + sq * rnorm(length(xa))
    up   <- xa >= ha; dn <- xa <= la
    done <- up | dn
    if (any(done)) {
      dec_t[idx[done]] <- step * dt
      resp[idx[done]]  <- ifelse(up[done], 1L, -1L)
      keep <- !done
      idx <- idx[keep]; xa <- xa[keep]; va <- va[keep]; la <- la[keep]; ha <- ha[keep]
    }
  }
  list(dt = dec_t, resp = resp)   # NA = not absorbed within t_max seconds
}

sim_person_trials <- function(pp) {
  n_cond <- N_TRIALS / length(COND_LEVELS)
  idx    <- rep(seq_len(nrow(pp)), each = N_TRIALS)
  cc     <- pp[idx, ]
  cond   <- rep(rep(COND_LEVELS, each = n_cond), nrow(pp))
  v      <- cc$nu_intercept + cc$nu_c * cond
  x      <- rwiener_euler(v, cc$alpha, 0.5 * cc$alpha)
  data.frame(person = cc$person, condition = cond,
             rt  = x$dt + cc$tau,
             acc = ifelse(x$resp == 1L, 1, -1))
}

# in batches of 100 persons 
batches <- split(seq_len(N_PERSONS), ceiling(seq_len(N_PERSONS) / 100))
trials  <- bind_rows(lapply(seq_along(batches), function(b) {
  cat("batch", b, "of", length(batches), "\n")
  sim_person_trials(persons[batches[[b]], ])
}))

n_unfinished <- sum(is.na(trials$rt))
cat("trials not finished within 20 s (dropped):", n_unfinished,
    sprintf("(%.4f%%)", 100 * n_unfinished / nrow(trials)), "\n")
trials <- filter(trials, !is.na(rt))

#### summary statistics from the EZ likelihood (one cell per person and condition)

# ez equations
ez_equations <- function(v, a) {
  z  <- a / 2
  av <- pmax(abs(v), 1e-2)
  Pu <- ifelse(abs(v) < 1e-8, 0.5, (1 - exp(-2 * v * z)) / (1 - exp(-2 * v * a)))
  Mu <- (a / tanh(a * av) - z / tanh(z * av)) / av
  Vu <- (a^2 / sinh(a * av)^2 - z^2 / sinh(z * av)^2) / av^2 +
    (a / tanh(a * av) - z / tanh(z * av)) / av^3
  data.frame(Pu = Pu, MDT = Mu, VDT = Vu)
}

J_cell <- N_TRIALS / length(COND_LEVELS)
lik <- persons %>%
  select(person, nu_intercept, nu_c, alpha, tau) %>%
  crossing(condition = COND_LEVELS) %>%
  mutate(nu = nu_intercept + nu_c * condition)
mo  <- ez_equations(lik$nu, lik$alpha)
n   <- nrow(lik)
lik <- lik %>%
  mutate(Pc_pred  = mo$Pu,
         MRT_pred = tau + mo$MDT,
         VRT_pred = mo$VDT,
         C        = rbinom(n, J_cell, Pc_pred),
         acc      = C / J_cell,
         MRT_correct = ifelse(C >= 1, rnorm(n, MRT_pred, sqrt(VRT_pred / pmax(C, 1))), NA),
         VRT_correct = ifelse(C >= 2, rnorm(n, VRT_pred, sqrt(2 * VRT_pred^2 / pmax(C - 1, 1))), NA))

cat("EZ likelihood: share of negative VRT:", mean(lik$VRT_correct < 0, na.rm = TRUE),
    "| share of cells with acc < 0.5:", mean(lik$acc < 0.5), "\n")

#### summaries per person and condition (from the simulated trials)

# acc is computed last: summarise() would otherwise overwrite it before the RT statistics use it
cells <- trials %>%
  mutate(correct = acc == 1) %>%
  group_by(person, condition) %>%
  summarise(MRT_correct  = if (sum(correct)  >= 1) mean(rt[correct])  else NA_real_,
            MRT_error    = if (sum(!correct) >= 1) mean(rt[!correct]) else NA_real_,
            SDRT_correct = if (sum(correct)  >= 2) sd(rt[correct])    else NA_real_,
            q10_correct  = if (any(correct)) quantile(rt[correct], 0.1, names = FALSE) else NA_real_,
            q90_correct  = if (any(correct)) quantile(rt[correct], 0.9, names = FALSE) else NA_real_,
            min_rt       = min(rt),
            acc          = mean(correct),
            .groups = "drop") %>%
  mutate(err_minus_cor = MRT_error - MRT_correct)       # > 0: slow errors, < 0: fast errors

# condition effects per person (cond +1 vs cond -1)
effects <- cells %>%
  group_by(person) %>%
  summarise(acc_effect = acc[condition == COND_LEVELS[1]] - acc[condition == COND_LEVELS[2]],
            MRT_effect = MRT_correct[condition == COND_LEVELS[2]] - MRT_correct[condition == COND_LEVELS[1]],
            .groups = "drop")

# distribution over the simulated persons
stat_table <- cells %>%
  select(condition, acc, MRT_correct, MRT_error, err_minus_cor, SDRT_correct,
         q10_correct, q90_correct, min_rt) %>%
  pivot_longer(-condition, names_to = "stat") %>%
  bind_rows(effects %>% pivot_longer(-person, names_to = "stat") %>%
              mutate(condition = NA) %>% select(-person)) %>%
  group_by(stat, condition) %>%
  summarise(q01    = quantile(value, 0.01, na.rm = TRUE),
            q05    = quantile(value, 0.05, na.rm = TRUE),
            median = median(value, na.rm = TRUE),
            q95    = quantile(value, 0.95, na.rm = TRUE),
            q99    = quantile(value, 0.99, na.rm = TRUE),
            share_NA = mean(is.na(value)),
            .groups = "drop")
print(stat_table, n = Inf)
write.csv(stat_table, file.path(OUT_DIR, "pp_stat_table.csv"), row.names = FALSE)

#### plots

MODEL_LABEL <- "EZ"
theme_set(theme_minimal(base_size = 11))
save_plot <- function(p, name, w = 10, h = 6) {
  print(p); ggsave(file.path(OUT_DIR, name), p, width = w, height = h, dpi = 200, bg = "white")
}
cond_lab <- function(x) factor(x, levels = COND_LEVELS, labels = paste("cond", COND_LEVELS))

# implied person parameters
p1 <- persons %>%
  select(person, nu_cond1, nu_cond2, alpha, tau) %>%
  pivot_longer(-person, names_to = "par") %>%
  ggplot(aes(value)) +
  geom_histogram(bins = 60, fill = "grey60") +
  facet_wrap(~ par, scales = "free") +
  labs(title = paste(MODEL_LABEL, "- person parameters implied by the priors"),
       subtitle = "nu_cond1 / nu_cond2 = drift in cond +1 / cond -1", x = NULL, y = NULL)
save_plot(p1, "pp_person_parameters.png")

# summary statistics per person and condition (upper / lower 0.5% cut for readability)
p2 <- cells %>%
  select(person, condition, acc, MRT_correct, MRT_error, err_minus_cor, SDRT_correct, min_rt) %>%
  pivot_longer(-c(person, condition), names_to = "stat") %>%
  filter(is.finite(value)) %>%
  group_by(stat) %>%
  filter(value >= quantile(value, 0.005), value <= quantile(value, 0.995)) %>%
  ungroup() %>%
  ggplot(aes(value, fill = cond_lab(condition))) +
  geom_histogram(bins = 60, alpha = 0.6, position = "identity") +
  facet_wrap(~ stat, scales = "free") +
  labs(title = paste(MODEL_LABEL, "- one simulated person per prior draw,", N_TRIALS, "trials each"),
       subtitle = "err_minus_cor > 0: errors slower than correct responses", x = NULL, y = NULL, fill = NULL)
save_plot(p2, "pp_person_stats.png", w = 11, h = 7)



# accuracy vs mean correct RT, one point per person and condition
p4 <- ggplot(cells, aes(MRT_correct, acc, colour = cond_lab(condition))) +
  geom_point(size = 0.6, alpha = 0.4) +
  labs(title = paste(MODEL_LABEL, "- accuracy vs mean correct RT"),
       x = "mean correct RT (s)", y = "accuracy", colour = NULL)
save_plot(p4, "pp_acc_vs_rt.png", w = 8, h = 6)


cat("\nSaved in", OUT_DIR, "\n")