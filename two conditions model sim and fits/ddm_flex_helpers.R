
suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
})

DATA_DIR  <- file.path(PROJECT_DIR, "data_sim_svhier")
DDM5_DIR  <- file.path(PROJECT_DIR, "fits_ddm5_2conditions")
STAN_FILE <- file.path(PROJECT_DIR, "ddm_flex.stan")
if (!file.exists(STAN_FILE)) stop("Stan file not found: ", STAN_FILE)



compile_flex <- function() {
  mod <- tryCatch(
    cmdstan_model(STAN_FILE, cpp_options = list(stan_threads = TRUE)),
    error = function(e) { message("compile error: ", conditionMessage(e)); NULL }
  )
  if (!is.null(mod)) {
    cat("\nCOMPILE CHECK: precision argument ACCEPTED by this Stan version\n\n")
    return(list(mod = mod, has_precision = TRUE))
  }
  cat("\nCOMPILE CHECK: precision argument NOT accepted --> compiling a copy without it",
      "(Stan default 1e-4 is always used)\n\n")
  code <- readLines(STAN_FILE)
  code <- gsub(", sw, st0, prec)", ", sw, st0)", code, fixed = TRUE)
  f <- file.path(PROJECT_DIR, "ddm_flex_noprec.stan")
  writeLines(code, f)
  list(mod = cmdstan_model(f, cpp_options = list(stan_threads = TRUE)), has_precision = FALSE)
}



flex_priors <- list(
  prior_beta_nu_intercept     = c(1.5, 1.0),
  prior_beta_nu_c             = c(0.0, 1.0),
  prior_beta_alpha            = c(log(1.5), 0.5),     # log scale
  prior_beta_tau              = c(0.25, 0.10),        # seconds
  prior_beta_w                = c(0.0, 0.5),          # logit scale
  prior_beta_sv               = c(log(0.5), 0.5),     # log scale
  prior_sigma_nu_intercept    = c(0, 0.5),
  prior_sigma_nu_c            = c(0, 0.5),
  prior_sigma_alpha_intercept = c(0, 0.2),
  prior_sigma_tau_intercept   = c(0, 0.1),
  prior_sigma_w               = c(0, 0.5),
  prior_sigma_sv              = c(0, 0.3),
  prior_sw                    = c(0.15, 0.10),        # used only if sw is estimated
  prior_st0                   = c(0.08, 0.05),        # seconds; used only if st0 is estimated
  sw_upper      = 0.4,   # true sw 0.15; true w 0.33-0.75
  guard_a_max   = 5,     # true a 1.1-2.0
  guard_sv_max  = 3,     # true sv 0.3-0.9
  guard_v_max   = 8      # true |drift| < 4.8
)



make_flex_data <- function(r, n_person, estimate_sw, estimate_st0, precision) {
  dat <- read.csv(file.path(DATA_DIR, sprintf("rep%02d", r), "sim_data_full_ddm.csv"))
  all_subj  <- sort(unique(dat$subj))
  keep_subj <- seq_len(min(n_person, length(all_subj)))
  dat <- dat[dat$subj %in% all_subj[keep_subj], ]
  dat$pid <- match(dat$subj, sort(unique(dat$subj)))
  min_rt <- as.numeric(tapply(dat$rt, dat$pid, min))
  stan_data <- c(list(N = nrow(dat), I = max(dat$pid), pid = dat$pid, rt = dat$rt,
                      acc = as.integer(dat$acc), condition = dat$condition,
                      min_rt = min_rt, grainsize = 1),
                 flex_priors,
                 list(precision_derivatives = precision,
                      estimate_sw = as.integer(estimate_sw),
                      estimate_st0 = as.integer(estimate_st0)))
  list(stan_data = stan_data, subj_ids = all_subj[keep_subj],
       keep_subj = keep_subj, min_rt = min_rt)
}



make_flex_init_fun <- function(r, d, estimate_sw, estimate_st0,
                               sw_init = 0.10, st0_init = 0.05, verbose = TRUE) {
  I <- d$stan_data$I; min_rt <- d$min_rt; keep_subj <- d$keep_subj
  f5 <- file.path(DDM5_DIR, sprintf("rep%02d", r), "summary_ddm5.csv")
  if (file.exists(f5)) {
    s5 <- read.csv(f5)
    g  <- function(v) s5$mean[s5$variable == v]
    pv <- function(base) s5$mean[grepl(paste0("^", base, "\\["), s5$variable)][keep_subj]
    zf <- function(x, mu, sig) pmin(pmax((x - mu) / max(sig, 0.05), -3), 3)
    b  <- list(nu1 = g("beta_nu_intercept"), nuc = g("beta_nu_c"), al = g("beta_alpha"),
               tau = g("beta_tau"), w = g("beta_w"), sv = g("beta_sv"),
               s_nu1 = g("sigma_nu_intercept"), s_nuc = g("sigma_nu_c"),
               s_al = g("sigma_alpha_intercept"), s_tau = g("sigma_tau_intercept"),
               s_w = g("sigma_w"), s_sv = g("sigma_sv"))
    z0 <- list(nu1 = zf(pv("nu_intercept"), b$nu1, b$s_nu1),
               nuc = zf(pv("nu_c"),         b$nuc, b$s_nuc),
               al  = zf(pv("alpha_intercept"), b$al, b$s_al),
               w   = zf(qlogis(pv("w")),    b$w,   b$s_w),
               sv  = zf(log(pv("sv")),      b$sv,  b$s_sv))
    t0_0 <- pmax(pmin(pv("t0") - 0.02, 0.9 * min_rt), 0.05)
    if (verbose) cat("inits from DDM5 rep", r, "\n")
  } else {
    b  <- list(nu1 = 1.5, nuc = 0.5, al = log(1.5), tau = 0.2, w = 0, sv = log(0.5),
               s_nu1 = 0.3, s_nuc = 0.2, s_al = 0.1, s_tau = 0.05, s_w = 0.2, s_sv = 0.1)
    z0 <- list(nu1 = rep(0, I), nuc = rep(0, I), al = rep(0, I), w = rep(0, I), sv = rep(0, I))
    t0_0 <- pmin(0.15, 0.8 * min_rt)
    if (verbose) cat("default inits (no DDM5 summary for rep", r, ")\n")
  }
  jit <- function(x, sd = 0.05) x + rnorm(length(x), 0, sd)
  function() {
    init <- list(
      beta_nu_intercept = jit(b$nu1), beta_nu_c = jit(b$nuc), beta_alpha = jit(b$al, 0.02),
      beta_tau = min(max(jit(b$tau, 0.01), 0.06), 0.44),
      beta_w = jit(b$w, 0.02), beta_sv = jit(b$sv),
      sigma_nu_intercept = max(b$s_nu1, 0.05), sigma_nu_c = max(b$s_nuc, 0.05),
      sigma_alpha_intercept = max(b$s_al, 0.03), sigma_tau_intercept = max(b$s_tau, 0.02),
      sigma_w = max(b$s_w, 0.05), sigma_sv = max(b$s_sv, 0.05),
      z_nu_intercept = jit(z0$nu1), z_nu_c = jit(z0$nuc), z_alpha_intercept = jit(z0$al),
      z_w = jit(z0$w), z_sv = jit(z0$sv),
      t0 = pmax(pmin(t0_0 + rnorm(I, 0, 0.005), 0.95 * min_rt), 0.02)
    )
   
    if (estimate_sw)  init$sw_par  <- array(sw_init,  dim = 1)
    if (estimate_st0) init$st0_par <- array(st0_init, dim = 1)
    init
  }
}
