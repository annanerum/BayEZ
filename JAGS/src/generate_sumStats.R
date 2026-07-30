#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
# This file contains functions to generate summary statistics directly from
# true DDM parameters using the EZ-DDM forward equations.
#
# Two main functions:
# 1. EZ_forward(): Computes EXPECTED summary statistics (deterministic)
#    - Returns the theoretical values of Pc, MRT, VRT given parameters
#
# 2. generate_sumStats(): Samples summary statistics from sampling distributions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

#------------------------------------------------------------------------------#
# EZ_forward: Compute expected summary statistics from DDM parameters
#------------------------------------------------------------------------------#
EZ_forward <- function(bound, drift, nondt) {
    
    # Intermediate calculation: exponential term
    ey  <- exp(-bound * drift)
    
    # Proportion correct (accuracy)
    Pc  <- 1 / (1 + ey)
    
    # RT precision (inverse variance)
    # From EZ-DDM equations
    PRT <- 2 * drift^3 / bound * (ey + 1)^2 / 
           (2 * -bound * drift * ey - ey^2 + 1)
    
    # Mean decision time
    MDT <- (bound / (2 * drift)) * (1 - ey) / (1 + ey)
    
    # Mean reaction time = mean decision time + non-decision time
    MRT <- MDT + nondt
    
    # Variance of RT (inverse of precision)
    VRT <- 1 / PRT
    
    return(list(Pc = Pc, MRT = MRT, VRT = VRT, PRT = PRT))
}


#------------------------------------------------------------------------------#
# generate_sumStats: Sample summary statistics from sampling distributions
#------------------------------------------------------------------------------#
generate_sumStats <- function(true_pars, nTrials) {

    # Handle nTrials: expand to vector if single value
    n_rows <- nrow(true_pars)
    if (length(nTrials) == 1) {
        nTrials <- rep(nTrials, n_rows)
    }
    
    # Check if condition column exists
    has_condition <- "cond" %in% colnames(true_pars)
    
    # Initialize output matrix
    if (has_condition) {
        out <- matrix(NA, nrow = n_rows, ncol = 7)
        colnames(out) <- c("sub", "cond", "nTrials", "sum_correct", 
                           "meanAccuracy", "meanRT", "varRT")
        out[, "cond"] <- true_pars$cond
    } else {
        out <- matrix(NA, nrow = n_rows, ncol = 6)
        colnames(out) <- c("sub", "nTrials", "sum_correct", 
                           "meanAccuracy", "meanRT", "varRT")
    }

    out[, "sub"]  <- true_pars$sub
    out[, "nTrials"] <- nTrials
    
    # Generate summary statistics for each row
    for (i in 1:n_rows) {
        # Extract true parameters for this row
        drift_i <- true_pars$drift[i]
        bound_i <- true_pars$bound[i]
        nondt_i <- true_pars$nondt[i]
        n_i     <- nTrials[i]
        
        # Compute expected summary statistics
        expected <- EZ_forward(bound_i, drift_i, nondt_i)
        
        #----------------------------------------------------------------------#
        # Sample from sampling distributions
        #----------------------------------------------------------------------#
        
        # Number correct: Binomial(nTrials, Pc)
        # Clamp to [1, n-1] to avoid edge cases (following Anne's code)
        sum_correct_i <- rbinom(1, n_i, expected$Pc)
        sum_correct_i <- max(1, min(sum_correct_i, n_i - 1))
        
        # Mean RT: Normal(MRT, sqrt(VRT / n_correct))
        # Standard error of the mean
        se_mrt <- sqrt(expected$VRT / sum_correct_i)
        meanRT_i <- rnorm(1, expected$MRT, se_mrt)
        
        # RT Variance: sampled from its sampling distribution
        # Var(s^2) = 2 * sigma^4 / (n - 1) for normal data
        # So s^2 ~ Normal(sigma^2, sqrt(2 * sigma^4 / (n-1)))
        df_var <- max(sum_correct_i - 1, 1)
        se_var <- sqrt(2 * expected$VRT^2 / df_var)
        varRT_i <- rnorm(1, expected$VRT, se_var)
        varRT_i <- max(varRT_i, 1e-6)  # Ensure positive
        
        #----------------------------------------------------------------------#
        # Store results
        #----------------------------------------------------------------------#
        out[i, "sum_correct"]   <- sum_correct_i
        out[i, "meanAccuracy"]  <- sum_correct_i / n_i
        out[i, "meanRT"]        <- meanRT_i
        out[i, "varRT"]         <- varRT_i
    }
    
    return(out)
}

