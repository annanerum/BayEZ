#===============================================================================
# Helper functions to automate design matrix simulations
# to be used in the simulation study JAGS report
# (inside function run_one_simulation())
#===============================================================================

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
# extract_beta: Extract beta parameters from JAGS output
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
# If P=1 we need a scalar name, otherwise we use array names
# Returns: Vector of posterior means for all beta coefficients
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
extract_beta <- function(sims, param_name, P) {
    if (P == 1) {
        return(mean(sims[, param_name]))
    } else {
        sapply(1:P, function(p) mean(sims[, paste0(param_name, "[", p, "]")]))
    }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
# get_beta_names: Identifies the correct JAGS parameter names based on P
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
get_beta_names <- function(param_name, P) {
    if (P == 1) {
        return(param_name)
    } else {
        paste0(param_name, "[", 1:P, "]")
    }
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
# get_rhat_params: Get all parameter names for convergence diagnostics
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
get_rhat_params <- function(P_nu, P_alpha, P_tau) {
    c(get_beta_names("beta_nu", P_nu),
      get_beta_names("beta_alpha", P_alpha),
      get_beta_names("beta_tau", P_tau),
      "sigma_nu", "sigma_alpha", "sigma_tau", "sigma_v")
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
# build_pop_estimates: Extracts and names all population-level estimates
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
build_pop_estimates <- function(sims, P_nu, P_alpha, P_tau) {
    
    beta_nu_est    <- extract_beta(sims, "beta_nu", P_nu)
    beta_alpha_est <- extract_beta(sims, "beta_alpha", P_alpha)
    beta_tau_est   <- extract_beta(sims, "beta_tau", P_tau)
    
    c(  setNames(beta_nu_est, paste0("beta_nu_", 1:P_nu)),
        setNames(beta_alpha_est, paste0("beta_alpha_", 1:P_alpha)),
        setNames(beta_tau_est, paste0("beta_tau_", 1:P_tau)),
        sigma_nu    = mean(sims[, "sigma_nu"]),
        sigma_alpha = mean(sims[, "sigma_alpha"]),
        sigma_tau   = mean(sims[, "sigma_tau"]),
        sigma_v     = mean(sims[, "sigma_v"]))
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
# build_pop_true: Creates a named vector of true pop pars matching the estimates
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
build_pop_true <- function(true_pars, P_nu, P_alpha, P_tau) {
    
    # Build beta vectors based on P values
    beta_nu_true <- if (P_nu >= 2) {
        c(true_pars$mu_nu, true_pars$b1_nu)
    } else {
        true_pars$mu_nu
    }
    
    beta_alpha_true <- if (P_alpha >= 2) {
        c(true_pars$mu_alpha, true_pars$b_alpha)
    } else {
        true_pars$mu_alpha
    }
    
    beta_tau_true <- if (P_tau >= 2) {
        c(true_pars$mu_tau, true_pars$b1_tau)
    } else {
        true_pars$mu_tau
    }
    
    c(setNames(beta_nu_true, paste0("beta_nu_", 1:P_nu)),
      setNames(beta_alpha_true, paste0("beta_alpha_", 1:P_alpha)),
      setNames(beta_tau_true, paste0("beta_tau_", 1:P_tau)),
      sigma_nu    = true_pars$sigma_nu,
      sigma_alpha = true_pars$sigma_alpha,
      sigma_tau   = true_pars$sigma_tau,
      sigma_v     = true_pars$sigma_v)
}

