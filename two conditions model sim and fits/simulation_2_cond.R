library(rtdists)
library(dplyr)
PROJECT_DIR <- "/Users/annegiacobello/Desktop/EZ Paper "   
DATA_DIR    <- file.path(PROJECT_DIR, "data_sim_fullddm_2condtions")

# condition with two levels (+1 / -1)
# balanced
# condition effect on drift, alpha and tau are intercept-only
# starting point bias at 0.55 
# hierarchical: subject-level random effects on nu_intercept, nu_c, alpha_intercept, tau_intercept_raw and w

n_subj   <- 100
n_trials <- 200  


true_pop <- list(
  beta_nu1    = 1.30,       # drift intercept, average drift when condition is 0
  # condition +1: drift --> 1.30 + 0.5 = 1.80
  # condition -1: drift --> 1.30 − 0.5 = 0.80
  nu_c        = 0.5,        # average condition effect on drift, shift between conditions
  beta_alpha1 = log(1.5),   # 0.405 --> boundary sep intercept on the log scale
  # --> exp(beta_alpha1) = 1.5 = average boundary sep (always positive via exp)
  beta_tau    = 0.20,       # population mean non-decision time t0 (in seconds, natural scale)
  w           = 0.55,       # population mean relative starting point (bias)
  # light biased toward upper boundary
  # must be between 0 and 1
  # proportion of the boundary separation 
  # (0 = starting at the lower boundary, 1 = at the upper boundary)
  sw          = 0.15,       # relative range of trial-to-trial starting point variability
  # (jitter of w within a person)
  sv          = 0.50,       # between-trial SD (variability) of drift rate (in the person)
  st0_raw     = 0.15        # relative non-decision time variability (trial-to-trial t0 jitter)
  # st0 = 2 * st0_raw * t0 = 0.3 * t0 (mean NDT = 1.15 * t0)
) 

sd_subj <- list(
  sigma_nu_intercept    = 0.20, # SD of drift intercept across subjects (to be confirmed)
  sigma_nu_c            = 0.30, # SD of condition effect on drift
  # share of people with negative drift in condition -1:
  # pnorm(0, 1.30 - 0.5, sqrt(0.20^2 + 0.30^2)) --> about 1.3%
  # share of people with reversed condition effect: pnorm(0, 0.5, 0.30) --> about 4.8%
  sigma_alpha_intercept = 0.08, # SD of log boundary sep intercept
  sigma_tau_intercept   = 0.05, # SD of t0 across subjects (in seconds)
  sigma_w               = 0.25  # SD of starting point bias (on the logit scale, keeps it in (0,1))
  # centered on w = 0.55 --> qlogis(0.55) = 0.2007
  # +-1 SD on the logit scale: plogis(0.2007 -/+ 0.25) --> 0.488 - 0.611
  # +-2 SD on the logit scale: plogis(0.2007 -/+ 0.50) --> 0.426 - 0.668
  # plogis(qlogis(p)) == p
)

simulate_subject <- function(sp) {
  # sp = one row of subj_par --> one persons parameters
  
  condition <- sample(rep(c(1, -1), each = n_trials / 2))
  # randomly shuffled and balanced per person 100 trials at +1 and 100 at -1
  
  a_subj  <- exp(sp$alpha_intercept)  # constant per person across trials (varies by person)
  # boundary separation (real scale)
  # intercept-only (no trial-to-trial variation and no condition effect)
  
  
  v_trial <- sp$nu_intercept + sp$nu_c * condition  # person's drift rate per trial
  # person-specific intercept +/- person-specific condition effect
  
  
  z_abs  <- sp$w * a_subj  # person specific starting point
  # converts relative (0-1) to real units 
  
  sz_abs <- true_pop$sw * a_subj
  # trial-to-trial starting point variability (range), also converted 
  
  sim <- rdiffusion(
    n = n_trials, a = a_subj, v = v_trial, t0 = sp$t0,
    z = z_abs, sz = sz_abs, sv = true_pop$sv, st0 = sp$st0
  )
  # drawing 200 random RT and response pairs from the diffusion model 
  # this persons a, v (one per trial), t0 (NDT in seconds), z (bias), sz (var in bias), sv (var in drift), st0 (var in NDT)
  
  data.frame(
    subj      = sp$subj, # ppn id
    trial     = seq_len(n_trials), # trial within subject
    condition = condition, # +1 or -1 for that trial
    a         = a_subj, # this ppns boundary, repeated every row (constant per ppn)
    v         = v_trial, # this trials drift rate (one of two values per ppn depending on condition)
    rt        = sim$rt, # simulated reaction time in seconds
    acc       = ifelse(sim$response == "upper", 1L, -1L) # +1 = correct/upper, -1 = error/lower
  )
}

# DATA_DIR is defined at the top (absolute path) -- do not redefine it here


run_one_rep <- function(rep_id, seed) {
  set.seed(seed)  # each rep uses a different seed (base_seed + rep_id below)
  
  z_nu_c  <- rnorm(n_subj)
  z_alpha <- rnorm(n_subj)
  z_tau   <- rnorm(n_subj)
  z_w     <- rnorm(n_subj)
  z_nu1   <- rnorm(n_subj)
  # standard-normal random draws, one set per ppn
  # random effect for each parameter
  # z_nu1 drawn last so the other draws stay identical to earlier runs with the same seed
  
  subj_par <- data.frame(
    subj              = 1:n_subj,
    nu_c              = true_pop$nu_c        + sd_subj$sigma_nu_c            * z_nu_c,
    # each ppns condition effect on drift
    nu_intercept      = true_pop$beta_nu1    + sd_subj$sigma_nu_intercept    * z_nu1,
    # each ppns drift intercept
    alpha_intercept   = true_pop$beta_alpha1 + sd_subj$sigma_alpha_intercept * z_alpha,
    # each ppns log boundary sep intercept
    tau_intercept_raw = true_pop$beta_tau    + sd_subj$sigma_tau_intercept   * z_tau,
    # each ppns t0 in seconds
    w                 = plogis(qlogis(true_pop$w) + sd_subj$sigma_w * z_w)
    # each ppns starting point --> kept in (0,1) by logit transform
  )
  
  subj_par$t0  <- subj_par$tau_intercept_raw
  # tau_intercept_raw is already on the natural scale (seconds) --> used directly as t0
  stopifnot(all(subj_par$t0 > 0)) 
  
  subj_par$st0 <- true_pop$st0_raw * 2 * subj_par$t0
  # each subject's trial-to-trial NDT variability, scaled relative to their own t0
  
  sim_list <- lapply(seq_len(n_subj), function(i) simulate_subject(subj_par[i, ]))
  
  sim_data <- do.call(rbind, sim_list)
  
  rownames(sim_data) <- NULL
  stopifnot(all(table(sim_data$subj) == n_trials))
  
  min_rt_obs <- tapply(sim_data$rt, sim_data$subj, min)
  gap_min_rt_t0 <- min_rt_obs[as.character(subj_par$subj)] - subj_par$t0
  # how far each person's fastest RT lies above their true t0 (always > 0)
  
  out_dir <- file.path(DATA_DIR, sprintf("rep%02d", rep_id))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  write.csv(sim_data, file.path(out_dir, "sim_data_full_ddm.csv"), row.names = FALSE)
  # trial-level data for fitting
  saveRDS(subj_par,      file.path(out_dir, "true_subject_parameters.rds"))
  # true person-level parameters
  saveRDS(true_pop,      file.path(out_dir, "true_population_parameters.rds"))
  # true population-level parameters
  saveRDS(min_rt_obs,    file.path(out_dir, "min_rt_observed.rds"))
  # per person observed minimum RT 
  
  cat(sprintf("rep%02d done | seed=%d | mean_rt=%.3f | mean_t0=%.3f | acc=%.3f | min_gap=%.3f\n",
              rep_id, seed, mean(sim_data$rt), mean(subj_par$t0),
              mean(sim_data$acc == 1), min(gap_min_rt_t0)))
  invisible(sim_data)
}

n_reps    <- 10   # how many simulated datasets 
base_seed <- 8181 # starting seed
# each rep uses base_seed + rep_id --> every rep is reproducible 
for (r in seq_len(n_reps)) run_one_rep(rep_id = r, seed = base_seed + r)


#### plotting 

# all reps into one data frame
dat_all <- bind_rows(lapply(1:10, function(r)
  read.csv(file.path(DATA_DIR, sprintf("rep%02d", r), "sim_data_full_ddm.csv")) %>%
    mutate(rep = r)))

# condition effect per rep
cond_tab <- dat_all %>%
  group_by(rep, condition) %>%
  summarise(mean_rt = mean(rt), median_rt = median(rt), acc = mean(acc == 1), .groups = "drop")
print(cond_tab, n = 20)

# cond effect summarised over reps (mean, min, max)
cond_tab %>%
  group_by(condition) %>%
  summarise(across(c(mean_rt, median_rt, acc), list(mean = mean, min = min, max = max)))

# per-person summaries 
subj_all <- dat_all %>%
  group_by(rep, subj) %>%
  summarise(rt_effect = mean(rt[condition == -1]) - mean(rt[condition == 1]),
            acc_pos   = mean(acc[condition ==  1] == 1),
            acc_neg   = mean(acc[condition == -1] == 1),
            acc_all   = mean(acc == 1),
            .groups = "drop")

# per rep at or below chance in condition -1
subj_all %>%
  group_by(rep) %>%
  summarise(share_rt_effect_neg = mean(rt_effect < 0),
            share_acc_neg_le_05 = mean(acc_neg <= 0.5),
            min_acc_neg         = min(acc_neg),
            min_acc_all         = min(acc_all))

# pooled plots over all reps
par(mfrow = c(2, 3))
hist(subj_all$rt_effect, breaks = 40, main = "RT condition effect per person", xlab = "RT(-1) - RT(+1) in s")
abline(v = 0, lty = 2)
hist(dat_all$rt, breaks = 100, main = "RT distribution (all trials, all reps)", xlab = "RT in s")
hist(subj_all$acc_all, breaks = 30, main = "Accuracy (overall)", xlab = "proportion correct", xlim = c(0, 1))
hist(subj_all$acc_pos, breaks = 30, main = "Accuracy (cond +1)", xlab = "proportion correct", xlim = c(0, 1))
hist(subj_all$acc_neg, breaks = 30, main = "Accuracy (cond -1)", xlab = "proportion correct", xlim = c(0, 1))
abline(v = 0.5, lty = 2)
par(mfrow = c(1, 1))