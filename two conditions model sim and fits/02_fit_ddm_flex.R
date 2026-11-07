
options(error = NULL)

PAPER_DIR   <- "/Users/annegiacobello/Desktop/EZ Paper "     
PROJECT_DIR <- file.path(PAPER_DIR, "sv hier 2 cond")
source(file.path(PROJECT_DIR, "ddm_flex_helpers.R"))


ESTIMATE_SW   <- FALSE
ESTIMATE_ST0  <- TRUE
REPS          <- 1:10      
N_PERSON      <- 100
ITER_WARMUP   <- 500
ITER_SAMPLING <- 500
PRECISION     <- 1e-3


model_tag <- if (!ESTIMATE_SW && ESTIMATE_ST0) "ddm6_st0" else
             if (ESTIMATE_SW && ESTIMATE_ST0) "ddm7_pop" else
             sprintf("ddm_flex_sw%d_st0%d", ESTIMATE_SW, ESTIMATE_ST0)
FIT_DIR <- file.path(PROJECT_DIR, paste0("fits_", model_tag))
dir.create(FIT_DIR, showWarnings = FALSE, recursive = TRUE)

n_chains          <- 4
n_cores           <- parallel::detectCores()
threads_per_chain <- max(1, floor((n_cores - 1) / n_chains))   # leave one core for the Mac

cat(sprintf("model: %s | reps: %s | %d persons | %d + %d iterations | precision %g\ncores: %d --> %d chains x %d threads\n",
            model_tag, paste(REPS, collapse = " "), N_PERSON, ITER_WARMUP, ITER_SAMPLING,
            PRECISION, n_cores, n_chains, threads_per_chain))

cm  <- compile_flex()
mod <- cm$mod
if (!cm$has_precision && PRECISION != 1e-4) {
  cat("precision argument not available --> Stan default 1e-4 is used\n")
  PRECISION <- 1e-4
}

pop_pars  <- c("beta_nu_intercept", "beta_nu_c", "beta_alpha", "alpha_pop_median",
               "beta_tau", "beta_w", "w_pop_median", "beta_sv", "sv_pop_median",
               if (ESTIMATE_SW) "sw", if (ESTIMATE_ST0) "st0",
               "sigma_nu_intercept", "sigma_nu_c", "sigma_alpha_intercept",
               "sigma_tau_intercept", "sigma_w", "sigma_sv")
subj_pars <- c("nu_intercept", "nu_c", "a", "t0", "w", "sv", "mean_ndt")


fit_one_rep <- function(r) {
  rep_tag <- sprintf("rep%02d", r)
  out_dir <- file.path(FIT_DIR, rep_tag)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  if (file.exists(file.path(out_dir, "pop_summary.rds"))) {
    cat(rep_tag, "already fitted --> skipped\n"); return(invisible(NULL))
  }

  d <- make_flex_data(r, N_PERSON, ESTIMATE_SW, ESTIMATE_ST0, PRECISION)
  set.seed(9100 + r)
  init_fun <- make_flex_init_fun(r, d, ESTIMATE_SW, ESTIMATE_ST0)

  t_start <- Sys.time()
  fit <- mod$sample(
    data = d$stan_data, init = init_fun, seed = 9100 + r,
    chains = n_chains, parallel_chains = n_chains, threads_per_chain = threads_per_chain,
    iter_warmup = ITER_WARMUP, iter_sampling = ITER_SAMPLING,
    adapt_delta = 0.90, max_treedepth = 8, step_size = 0.1,
    refresh = 5, show_exceptions = FALSE,
    output_dir = out_dir,                       # not the macOS temp folder
    output_basename = sprintf("%s_%s", model_tag, rep_tag)
  )
  runtime_h <- as.numeric(difftime(Sys.time(), t_start, units = "hours"))

  pop_summ <- fit$summary(pop_pars, "mean", "median", "sd", "quantile2",
                          "rhat", "ess_bulk", "ess_tail") %>% mutate(rep = r)
  subj_summ <- fit$summary(subj_pars, "mean", "median", "sd", "quantile2",
                           "rhat", "ess_bulk", "ess_tail") %>%
    mutate(param = sub("\\[.*", "", variable),
           pid   = as.integer(gsub(".*\\[|\\]", "", variable)),
           subj  = d$subj_ids[pid], rep = r)

  diag <- fit$diagnostic_summary(quiet = TRUE)
  diag_df <- data.frame(rep = r, model = model_tag,
                        divergences   = sum(diag$num_divergent),
                        max_treedepth = sum(diag$num_max_treedepth),
                        min_ebfmi     = min(diag$ebfmi),
                        max_rhat      = max(c(pop_summ$rhat, subj_summ$rhat), na.rm = TRUE),
                        min_ess_bulk  = min(c(pop_summ$ess_bulk, subj_summ$ess_bulk), na.rm = TRUE),
                        runtime_h     = runtime_h,
                        n_person      = d$stan_data$I,
                        precision     = PRECISION)

  saveRDS(pop_summ,   file.path(out_dir, "pop_summary.rds"))
  saveRDS(subj_summ,  file.path(out_dir, "subj_summary.rds"))
  saveRDS(diag_df,    file.path(out_dir, "diagnostics.rds"))
  saveRDS(d$subj_ids, file.path(out_dir, "subj_ids.rds"))      # which persons were fitted

  cat(sprintf("%s done | %.1f h | div=%d | treedepth=%d | max Rhat=%.3f | min ESS=%.0f\n",
              rep_tag, runtime_h, diag_df$divergences, diag_df$max_treedepth,
              diag_df$max_rhat, diag_df$min_ess_bulk))
  print(pop_summ %>% select(variable, mean, q5, q95, rhat), n = Inf)
  invisible(diag_df)
}

for (r in REPS) {
  tryCatch(fit_one_rep(r),
           error = function(e) cat(sprintf("rep%02d FAILED: %s\n", r, conditionMessage(e))))
}

diag_all <- bind_rows(lapply(REPS, function(r) {
  f <- file.path(FIT_DIR, sprintf("rep%02d", r), "diagnostics.rds")
  if (file.exists(f)) readRDS(f)
}))
print(diag_all)
