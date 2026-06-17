# Model 2 NCP: nu=cond+slope (non-centred), alpha=cov, tau=intercept

library(rstan)
library(bayesplot)
library(ggplot2)
library(dplyr)
library(tidyr)
library(knitr)
library(parallel)

rstan_options(auto_write = TRUE)
options(mc.cores = 1)

if (!dir.exists("output_model2_v2_ncp_new_params")) dir.create("output_model2_v2_ncp_new_params")

## simulation grid
I_vals <- c(80, 160)
J_vals <- c(80, 160)
N_sim  <- 10
K      <- 2

## true parameters
true_m2 <- list(
  mu_nu    = 0.7, #1.0,
  b1_nu    = 0.2,
  mu_alpha = 1.5,
  b_alpha  = 0.3,
  mu_tau   = 0.3,
  sigma_nu    = 0.18, #0.3
  sigma_alpha = 0.25,
  sigma_tau   = 0.15,
  sigma_v     = 0.1  #0.2
)


## EZ point estimates
ez_point_estimates <- function(C_arr, J_arr, MRT_mat, VRT_mat, I, K) {
  nu_hat <- matrix(NA, I, K); alpha_hat <- matrix(NA, I, K); tau_hat <- matrix(NA, I, K)
  for (i in seq_len(I)) for (k in seq_len(K)) {
    Pc  <- C_arr[i,k] / J_arr[i,k]
    Pc  <- max(min(Pc, 1 - 1/(2*J_arr[i,k])), 1/(2*J_arr[i,k]))
    mrt <- MRT_mat[i,k]; vrt <- VRT_mat[i,k]
    if (is.na(mrt) || is.na(vrt) || vrt <= 0) next
    L   <- log(Pc/(1-Pc))
    nu4 <- L*(Pc^2*L - Pc*L + Pc - 0.5)/vrt
    if (is.na(nu4) || nu4 <= 0) next
    nu_hat[i,k]    <- sign(Pc-0.5)*nu4^(1/4)
    alpha_hat[i,k] <- L/nu_hat[i,k]
    tau_hat[i,k]   <- mrt - (alpha_hat[i,k]/(2*nu_hat[i,k]))*(2*Pc-1)
  }
  list(nu=nu_hat, alpha=alpha_hat, tau=tau_hat)
}


## derive priors
derive_priors <- function(ez) {
  nu_vals    <- as.vector(ez$nu[!is.na(ez$nu)])
  alpha_vals <- as.vector(ez$alpha[!is.na(ez$alpha) & ez$alpha > 0])
  tau_vals   <- as.vector(ez$tau[!is.na(ez$tau) & ez$tau > 0.05])
  
  nu_sd    <- min(max(sd(nu_vals),    0.3), 1.5)
  alpha_sd <- min(max(sd(alpha_vals), 0.2), 1.0)
  tau_sd   <- min(max(sd(tau_vals),   0.1), 0.5)
  
  beta_nu    <- matrix(c(mean(nu_vals), 0,
                         nu_sd,         nu_sd),        nrow=2, byrow=TRUE)
  beta_alpha <- matrix(c(mean(alpha_vals), 0,
                         alpha_sd*2,       alpha_sd*2), nrow=2, byrow=TRUE)
  beta_tau   <- matrix(c(mean(tau_vals),
                         tau_sd*2),                    nrow=2, byrow=TRUE)
  
  list(
    beta_nu     = beta_nu,
    beta_alpha  = beta_alpha,
    beta_tau    = beta_tau,
    sigma_nu    = c(nu_sd,    nu_sd/2),
    sigma_alpha = c(alpha_sd, alpha_sd/2),
    sigma_tau   = c(tau_sd,   tau_sd/2),
    sigma_v     = c(nu_sd*0.3, nu_sd)
  )
}


## pilot prior derivation
derive_priors_from_pilot <- function(J_per, true, I_pilot=20) {
  K <- 2; cond_dummy <- c(0,1)
  J_p   <- matrix(J_per, I_pilot, K)
  C_p   <- matrix(NA_integer_, I_pilot, K)
  MRT_p <- matrix(NA_real_,    I_pilot, K)
  VRT_p <- matrix(NA_real_,    I_pilot, K)
  
  x_p      <- rnorm(I_pilot, 0, 1)
  true_v_p <- rnorm(I_pilot, 0, true$sigma_v)
  
  for (i in seq_len(I_pilot)) for (k in seq_len(K)) {
    mean_nu_ik    <- true$mu_nu + (true$b1_nu + true_v_p[i]) * cond_dummy[k]
    mean_alpha_ik <- true$mu_alpha + true$b_alpha * x_p[i]
    mean_tau_ik   <- true$mu_tau
    
    repeat { a <- rnorm(1, mean_alpha_ik, true$sigma_alpha)
    if (a > 0.1 && a < 3.0) { ap <- a; break } }
    np <- rnorm(1, mean_nu_ik, true$sigma_nu)
    repeat { t <- rnorm(1, mean_tau_ik, true$sigma_tau)
    if (t > 0.05) { tp <- t; break } }
    
    e      <- exp(-ap*np); pi_c <- 1/(1+e)
    mu_rt  <- tp + (ap/(2*np))*((1-e)/(1+e))
    sig2rt <- (ap/(2*np^3))*((1-2*ap*np*e-e^2)/(1+e)^2)
    C_ik   <- max(1L, min(rbinom(1,J_per,pi_c), J_per-1L))
    MRT_p[i,k] <- rnorm(1, mu_rt, sqrt(sig2rt/C_ik))
    VRT_p[i,k] <- max(rnorm(1, sig2rt, sqrt(2*sig2rt^2/max(C_ik-1,1))), 1e-6)
    C_p[i,k]   <- C_ik
  }
  
  ez <- ez_point_estimates(C_p, J_p, MRT_p, VRT_p, I_pilot, K)
  derive_priors(ez)
}


## run one simulation and fit
run_one <- function(I, J_per, true, model, priors,
                    rhat_threshold=1.05, hopeless_threshold=10,
                    iter1=1000, warmup1=300, iter2=3000, warmup2=1000) {
  
  K <- 2; cond_dummy <- c(0,1)
  prs    <- priors
  x      <- rnorm(I, 0, 1)
  true_v <- rnorm(I, 0, true$sigma_v)
  
  P_nu <- 2; P_alpha <- 2; P_tau <- 1
  
  X_nu    <- array(NA_real_, dim=c(I, K, P_nu))
  X_alpha <- array(NA_real_, dim=c(I, K, P_alpha))
  X_tau   <- array(NA_real_, dim=c(I, K, P_tau))
  
  for (i in seq_len(I)) for (k in seq_len(K)) {
    X_nu[i,k,]    <- c(1, cond_dummy[k])
    X_alpha[i,k,] <- c(1, x[i])
    X_tau[i,k,]   <- c(1)
  }
  
  # true individual parameters
  alpha_true <- matrix(NA_real_, I, K)
  nu_true    <- matrix(NA_real_, I, K)
  tau_true   <- matrix(NA_real_, I, K)
  
  for (i in seq_len(I)) for (k in seq_len(K)) {
    mean_nu_ik    <- true$mu_nu + (true$b1_nu + true_v[i]) * cond_dummy[k]
    mean_alpha_ik <- true$mu_alpha + true$b_alpha * x[i]
    mean_tau_ik   <- true$mu_tau
    
    repeat { a <- rnorm(1, mean_alpha_ik, true$sigma_alpha)
    if (a > 0.1 && a < 3.0) { alpha_true[i,k] <- a; break } }
    nu_true[i,k] <- rnorm(1, mean_nu_ik, true$sigma_nu)
    repeat { t <- rnorm(1, mean_tau_ik, true$sigma_tau)
    if (t > 0.05) { tau_true[i,k] <- t; break } }
  }
  
  # observed data
  J_mat   <- matrix(J_per, I, K)
  C_mat   <- matrix(NA_integer_, I, K)
  MRT_mat <- matrix(NA_real_,    I, K)
  VRT_mat <- matrix(NA_real_,    I, K)
  
  for (i in seq_len(I)) for (k in seq_len(K)) {
    a <- alpha_true[i,k]; v <- nu_true[i,k]; tau_ik <- tau_true[i,k]
    e      <- exp(-a*v); pi_c <- 1/(1+e)
    mu_rt  <- tau_ik + (a/(2*v))*((1-e)/(1+e))
    sig2rt <- (a/(2*v^3))*((1-2*a*v*e-e^2)/(1+e)^2)
    C_ik   <- max(1L, min(rbinom(1,J_per,pi_c), J_per-1L))
    MRT_mat[i,k] <- rnorm(1, mu_rt, sqrt(sig2rt/C_ik))
    VRT_mat[i,k] <- max(rnorm(1, sig2rt, sqrt(2*sig2rt^2/max(C_ik-1,1))), 1e-6)
    C_mat[i,k]   <- C_ik
  }
  
  stan_data <- list(
    I=I, K=K, P_nu=P_nu, P_alpha=P_alpha, P_tau=P_tau,
    J=J_mat, C=C_mat, MRT=MRT_mat, VRT=VRT_mat,
    X_nu=X_nu, X_alpha=X_alpha, X_tau=X_tau,
    prior_beta_nu     = prs$beta_nu,
    prior_beta_alpha  = prs$beta_alpha,
    prior_beta_tau    = prs$beta_tau,
    prior_sigma_nu    = prs$sigma_nu,
    prior_sigma_alpha = prs$sigma_alpha,
    prior_sigma_tau   = prs$sigma_tau,
    prior_sigma_v     = prs$sigma_v
  )
  
  group_pars <- c("beta_nu[1]", "beta_nu[2]",
                  "beta_alpha[1]", "beta_alpha[2]",
                  "beta_tau[1]",
                  "sigma_nu", "sigma_alpha", "sigma_tau", "sigma_v")
  
  hopeless <- function(fit) {
    rh <- tryCatch(summary(fit)$summary[group_pars, "Rhat"], error=function(e) NA_real_)
    any(is.na(rh)) || any(is.nan(rh)) || any(is.infinite(rh)) ||
      any(rh > hopeless_threshold, na.rm=TRUE)
  }
  
  fit <- tryCatch(
    sampling(model, data=stan_data,
             iter=iter1, warmup=warmup1, chains=4, cores=1, refresh=0,
             control=list(adapt_delta=0.85, max_treedepth=10)),
    error=function(e) NULL)
  
  if (is.null(fit)) return(list(status="error", I=I, J=J_per))
  if (hopeless(fit)) return(list(status="excluded", reason="degenerate_fit1", I=I, J=J_per))
  
  rhat_max  <- max(summary(fit)$summary[group_pars, "Rhat"], na.rm=TRUE)
  converged <- "first_fit"
  
  if (rhat_max > rhat_threshold) {
    message(sprintf("Refit triggered (I=%d, J=%d) -- Rhat=%.2f", I, J_per, rhat_max))
    fit <- tryCatch(
      sampling(model, data=stan_data,
               iter=iter2, warmup=warmup2, chains=4, cores=1, refresh=0,
               control=list(adapt_delta=0.95, max_treedepth=12)),
      error=function(e) NULL)
    
    if (is.null(fit)) return(list(status="error", I=I, J=J_per))
    if (hopeless(fit)) return(list(status="excluded", reason="degenerate_refit", I=I, J=J_per))
    rhat_max  <- max(summary(fit)$summary[group_pars, "Rhat"], na.rm=TRUE)
    converged <- "refit"
    if (rhat_max > rhat_threshold)
      return(list(status="excluded", reason="non_converged", I=I, J=J_per))
  }
  
  s <- summary(fit)$summary
  
  # extract mean_ik
  mean_nu_est    <- matrix(NA_real_, I, K)
  mean_alpha_est <- matrix(NA_real_, I, K)
  mean_tau_est   <- matrix(NA_real_, I, K)
  for (i in seq_len(I)) for (k in seq_len(K)) {
    mean_nu_est[i,k]    <- s[paste0("mean_nu_ik[",    i, ",", k, "]"), "mean"]
    mean_alpha_est[i,k] <- s[paste0("mean_alpha_ik[", i, ",", k, "]"), "mean"]
    mean_tau_est[i,k]   <- s[paste0("mean_tau_ik[",   i, ",", k, "]"), "mean"]
  }
  
  # group-level results
  group_results <- data.frame(
    I = I, J = J_per,
    param    = c("mu_nu", "b1_nu", "mu_alpha", "b_alpha", "mu_tau",
                 "sigma_nu", "sigma_alpha", "sigma_tau", "sigma_v"),
    true_val = c(true$mu_nu, true$b1_nu, true$mu_alpha, true$b_alpha, true$mu_tau,
                 true$sigma_nu, true$sigma_alpha, true$sigma_tau, true$sigma_v),
    estimate  = s[c("beta_nu[1]", "beta_nu[2]",
                    "beta_alpha[1]", "beta_alpha[2]",
                    "beta_tau[1]",
                    "sigma_nu", "sigma_alpha", "sigma_tau", "sigma_v"), "mean"],
    converged = converged,
    rhat_max  = rhat_max,
    stringsAsFactors = FALSE
  )
  
  # individual results
  alpha_est <- matrix(NA_real_, I, K)
  nu_est    <- matrix(NA_real_, I, K)
  tau_est   <- matrix(NA_real_, I, K)
  for (i in seq_len(I)) for (k in seq_len(K)) {
    alpha_est[i,k] <- s[paste0("alpha[",i,",",k,"]"), "mean"]
    nu_est[i,k]    <- s[paste0("nu[",   i,",",k,"]"), "mean"]
    tau_est[i,k]   <- s[paste0("tau[",  i,",",k,"]"), "mean"]
  }
  
  indiv_results <- data.frame(
    I = I, J = J_per,
    condition  = rep(seq_len(K), each=I),
    true_alpha = as.vector(alpha_true), est_alpha = as.vector(alpha_est),
    true_nu    = as.vector(nu_true),    est_nu    = as.vector(nu_est),
    true_tau   = as.vector(tau_true),   est_tau   = as.vector(tau_est)
  )
  
  # mean_ik verification
  mean_ik_results <- data.frame(
    I = I, J = J_per,
    participant    = rep(seq_len(I), times=K),
    condition      = rep(seq_len(K), each=I),
    est_mean_nu    = as.vector(mean_nu_est),
    est_mean_alpha = as.vector(mean_alpha_est),
    est_mean_tau   = as.vector(mean_tau_est),
    true_mean_nu   = as.vector(matrix(
      sapply(seq_len(K), function(k)
        sapply(seq_len(I), function(i)
          true$mu_nu + (true$b1_nu + true_v[i]) * cond_dummy[k]
        )), nrow=I)),
    true_mean_alpha = as.vector(matrix(
      sapply(seq_len(K), function(k)
        rep(true$mu_alpha + true$b_alpha * x, 1)), nrow=I)),
    true_mean_tau = true$mu_tau
  )
  
  # random slope results
  v_results <- data.frame(
    I = I, J = J_per,
    participant = seq_len(I),
    true_v      = true_v,
    est_v       = s[paste0("v[", seq_len(I), "]"), "mean"]
  )
  
  list(status="ok", group=group_results, indiv=indiv_results,
       mean_ik=mean_ik_results, v=v_results,
       rhat_max=rhat_max, converged=converged, I=I, J=J_per)
}


## compile Stan model
stan_file <- "model2_v2_ncp.stan"
cache <- paste0(tools::file_path_sans_ext(stan_file), ".rds")
if (file.exists(cache)) file.remove(cache)

## derive priors per J
set.seed(42)
priors_by_J <- setNames(
  lapply(J_vals, function(j) derive_priors_from_pilot(j, true_m2, I_pilot=20)),
  as.character(J_vals)
)

## build job list
jobs <- expand.grid(I=I_vals, J=J_vals, sim=seq_len(N_sim))
jobs_list <- split(jobs, seq_len(nrow(jobs)))

## parallel run
cl <- makeCluster(detectCores() - 1)
clusterSetRNGStream(cl, 42)

clusterExport(cl, varlist=c(
  "run_one", "derive_priors_from_pilot", "ez_point_estimates", "derive_priors",
  "true_m2", "priors_by_J", "stan_file", "K"
))

clusterEvalQ(cl, {
  library(rstan)
  library(dplyr)
  rstan_options(auto_write=FALSE)
  options(mc.cores=1)
  stan_model2 <<- stan_model(stan_file)
})

results_list <- parLapply(cl, jobs_list, function(row) {
  run_one(
    I      = row$I,
    J_per  = row$J,
    true   = true_m2,
    model  = stan_model2,
    priors = priors_by_J[[as.character(row$J)]]
  )
})

stopCluster(cl)
saveRDS(results_list, "output_model2_v2_ncp_new_params/simulation_results.rds")


## post-processing
ok_idx     <- sapply(results_list, function(r) !is.null(r) && r$status=="ok")
ok_results <- results_list[ok_idx]

results_group   <- bind_rows(lapply(ok_results, `[[`, "group"))
results_indiv   <- bind_rows(lapply(ok_results, `[[`, "indiv"))
results_mean_ik <- bind_rows(lapply(ok_results, `[[`, "mean_ik"))
results_v       <- bind_rows(lapply(ok_results, `[[`, "v"))

# exclusion report
n_total    <- length(results_list)
n_ok       <- sum(ok_idx)
n_excluded <- sum(sapply(results_list, function(r) !is.null(r) && r$status=="excluded"))
n_error    <- sum(sapply(results_list, function(r) is.null(r) || r$status=="error"))
cat("\n=== OVERALL ===\n")
print(data.frame(Total=n_total, Converged=n_ok, Excluded=n_excluded, Error=n_error,
                 `Excluded (%)`=round(100*n_excluded/n_total,1), check.names=FALSE))

status_df <- bind_rows(lapply(seq_along(results_list), function(i) {
  r <- results_list[[i]]; row <- jobs[i,]
  data.frame(I=row$I, J=row$J, sim=row$sim,
             status = if (is.null(r)) "error" else r$status,
             reason = if (is.null(r) || r$status=="ok") NA else r$reason)
}))

cat("\n=== EXCLUSIONS BY CELL ===\n")
status_df %>%
  group_by(I, J) %>%
  summarise(n=n(), excluded=sum(status=="excluded"),
            pct=round(100*excluded/n,1),
            reasons=paste(table(reason[status=="excluded"]), collapse=", "),
            .groups="drop") %>%
  print(n=Inf)

cat("\n=== DESIGN MATRIX VERIFICATION ===\n")
cat(sprintf("nu: r=%.3f | alpha: r=%.3f | tau: r=%.3f\n",
            cor(results_mean_ik$true_mean_nu,    results_mean_ik$est_mean_nu,    use="complete.obs"),
            cor(results_mean_ik$true_mean_alpha, results_mean_ik$est_mean_alpha, use="complete.obs"),
            cor(results_mean_ik$true_mean_tau,   results_mean_ik$est_mean_tau,   use="complete.obs")))

## plots
I_vals_f <- paste0("I = ", I_vals)
J_vals_f <- paste0("J = ", J_vals)

results_indiv <- results_indiv %>%
  mutate(
    I_label   = factor(paste0("I = ", I), levels=I_vals_f),
    J_label   = factor(paste0("J = ", J), levels=J_vals_f),
    condition = factor(condition, labels=c("Condition 1", "Condition 2"))
  )

results_v <- results_v %>%
  mutate(
    true_v  = true_v + true_m2$b1_nu,
    I_label = factor(paste0("I = ", I), levels=I_vals_f),
    J_label = factor(paste0("J = ", J), levels=J_vals_f)
  )

dark_theme <- theme_minimal(base_size=12) +
  theme(
    plot.background   = element_rect(fill="#1e1e1e", colour=NA),
    panel.background  = element_rect(fill="#1e1e1e", colour=NA),
    panel.grid.major  = element_line(colour="#333333"),
    panel.grid.minor  = element_blank(),
    text              = element_text(colour="white"),
    axis.text         = element_text(colour="white"),
    strip.text        = element_text(colour="white"),
    legend.background = element_rect(fill="#1e1e1e"),
    legend.text       = element_text(colour="white")
  )

save_plot <- function(p, filename, width=10, height=8) {
  ggsave(file.path("output_model2_v2_ncp_new_params", filename), p,
         width=width, height=height, dpi=300, bg="#1e1e1e")
}

p_alpha <- ggplot(results_indiv, aes(x=true_alpha, y=est_alpha, colour=condition)) +
  geom_abline(intercept=0, slope=1, linetype="dashed", colour="white", linewidth=0.5) +
  geom_point(alpha=0.35, size=0.9) +
  scale_colour_manual(values=c("Condition 1"="goldenrod","Condition 2"="steelblue")) +
  facet_grid(I_label ~ J_label) +
  labs(title="M2 NCP | alpha recovery", colour="Condition") +
  dark_theme + theme(legend.position="bottom")

p_nu <- ggplot(results_indiv, aes(x=true_nu, y=est_nu, colour=condition)) +
  geom_abline(intercept=0, slope=1, linetype="dashed", colour="white", linewidth=0.5) +
  geom_point(alpha=0.35, size=0.9) +
  scale_colour_manual(values=c("Condition 1"="goldenrod","Condition 2"="steelblue")) +
  facet_grid(I_label ~ J_label) +
  labs(title="M2 NCP | nu recovery", colour="Condition") +
  dark_theme + theme(legend.position="bottom")

p_tau <- ggplot(results_indiv, aes(x=true_tau, y=est_tau, colour=condition)) +
  geom_abline(intercept=0, slope=1, linetype="dashed", colour="white", linewidth=0.5) +
  geom_point(alpha=0.35, size=0.9) +
  scale_colour_manual(values=c("Condition 1"="goldenrod","Condition 2"="steelblue")) +
  facet_grid(I_label ~ J_label) +
  labs(title="M2 NCP | tau recovery", colour="Condition") +
  dark_theme + theme(legend.position="bottom")

p_v <- ggplot(results_v, aes(x=true_v, y=est_v)) +
  geom_abline(intercept=0, slope=1, linetype="dashed", colour="white", linewidth=0.5) +
  geom_point(colour="goldenrod", alpha=0.4, size=1.2) +
  facet_grid(I_label ~ J_label) +
  labs(title="M2 NCP | random slope v[i] recovery",
       x="True v", y="Estimated v") +
  dark_theme

save_plot(p_alpha, "m2_ncp_recovery_alpha.png")
save_plot(p_nu,    "m2_ncp_recovery_nu.png")
save_plot(p_tau,   "m2_ncp_recovery_tau.png")
save_plot(p_v,     "m2_ncp_recovery_v.png")

cat("\nDone. Results saved to output_model2_v2_ncp_new_params/\n")