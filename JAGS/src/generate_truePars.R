#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
# This file contains functions to generate individual-level true parameters
# from a given simulation design.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

#------------------------------------------------------------------------------#
# generate_truePars_simple: Generate parameters without conditions
#------------------------------------------------------------------------------#
generate_truePars_simple <- function(n_sub,
                                     drift_mean = 0.5, drift_sd = 0.2,
                                     bound_mean = 1.5, bound_sd = 0.25,
                                     nondt_mean = 0.3, nondt_sd = 0.05) {
    
    # Initialize output
    out <- data.frame(sub   = 1:n_sub, drift = numeric(n_sub),
                      bound = numeric(n_sub),nondt = numeric(n_sub))
    
    # Generate parameters for each subject
    for (i in 1:n_sub) {        
        # Drift rate: unconstrained
        out$drift[i] <- rnorm(1, drift_mean, drift_sd)
        
        # Boundary separation: truncated to [0.1, 3.0]
        repeat {
            b <- rnorm(1, bound_mean, bound_sd)
            if (b > 0.1 && b < 3.0) {
                out$bound[i] <- b
                break
            }
        }
        
        # Non-decision time: truncated to [0.05, Inf)
        repeat {
            t <- rnorm(1, nondt_mean, nondt_sd)
            if (t > 0.05) {
                out$nondt[i] <- t
                break
            }
        }
    }
    
    return(out)
}


#------------------------------------------------------------------------------#
# generate_truePars_condition: Generate parameters with condition effects
#------------------------------------------------------------------------------#
generate_truePars_condition <- function(n_sub, n_cond = 2,
                                        drift_mean = 0.5, drift_cond = 0.2, drift_sd = 0.18,
                                        slope_sd = 0.2,
                                        bound_mean = 1.5, bound_cov = 0.3, bound_sd = 0.25,
                                        nondt_mean = 0.3, nondt_cond = 0, nondt_sd = 0.15) {
    
    # Condition dummy coding: 0 for condition 1, 1 for condition 2, etc.
    cond_dummy <- 0:(n_cond - 1)
    
    # Generate subject-level random effects
    v <- rnorm(n_sub, 0, slope_sd)  # Random slopes for condition effect
    x <- rnorm(n_sub, 0, 1)         # Standardized covariate
    
    # Total rows: n_sub x n_cond
    n_rows <- n_sub * n_cond
    
    # Initialize output
    out <- data.frame(sub   = integer(n_rows),
                    cond  = integer(n_rows),
                    drift = numeric(n_rows),
                    bound = numeric(n_rows),
                    nondt = numeric(n_rows))
    
    # Generate parameters for each subject x condition
    row_idx <- 1
    for (i in 1:n_sub) {
        for (k in 1:n_cond) {
            
            out$sub[row_idx]  <- i
            out$cond[row_idx] <- cond_dummy[k]
            
            # Mean for drift rate: intercept + (fixed + random slope) * condition
            mean_drift_ik <- drift_mean + (drift_cond + v[i]) * cond_dummy[k]
            out$drift[row_idx] <- rnorm(1, mean_drift_ik, drift_sd)
            
            # Mean for boundary: intercept + covariate effect
            mean_bound_ik <- bound_mean + bound_cov * x[i]
            repeat {
                b <- rnorm(1, mean_bound_ik, bound_sd)
                if (b > 0.1 && b < 3.0) {
                    out$bound[row_idx] <- b
                    break
                }
            }
            
            # Mean for non-decision time: intercept + condition effect
            mean_nondt_ik <- nondt_mean + nondt_cond * cond_dummy[k]
            repeat {
                t <- rnorm(1, mean_nondt_ik, nondt_sd)
                if (t > 0.05) {
                    out$nondt[row_idx] <- t
                    break
                }
            }
            
            row_idx <- row_idx + 1
        }
    }
    
    # Store true values and random effects for later comparison
    true_values <- list(
        drift_mean = drift_mean,
        drift_cond = drift_cond,
        drift_sd   = drift_sd,
        slope_sd   = slope_sd,
        bound_mean = bound_mean,
        bound_cov  = bound_cov,
        bound_sd   = bound_sd,
        nondt_mean = nondt_mean,
        nondt_cond = nondt_cond,
        nondt_sd   = nondt_sd)
    
    random_effects <- data.frame(
        sub = 1:n_sub,
        v   = v,      # Random slope
        x   = x       # Covariate value
    )
    
    return(list(pars   = out,
                true_values    = true_values,
                random_effects = random_effects))
}