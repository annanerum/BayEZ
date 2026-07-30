#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
#                     TRIAL-LEVEL DATA GENERATION FOR DDM
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
# This file contains functions to generate trial-level data by emulating the
# random walk (Wiener) process implied by the Drift Diffusion Model.
#
# FUNCTIONS IN THIS FILE:
# -----------------------
# 1. sample_trial(a, v, dt, max_steps)
#    Generates a SINGLE trial. Returns RT and final choice (C).
#
# 2. sample_many_trials(a, v, t, n)
#    Generates N trials for one participant/condition using sample_trial().
#    Adds non-decision time to RTs and returns data frame with RT and accuracy.
#
# 3. sample_data(nPart, nTrials, parameter_set, nTrialsPerCondition, 
#                prevent_zero_accuracy)
#    Generates complete dataset for simulation studies.
#    Handles both simple (no conditions) and within-subject (2 conditions) designs.
#
# AUXILIARY FUNCTIONS:
# --------------------
# 4. identify_design(nPart, nTrials, nTrialsPerCondition, parameter_set)
#    Determines design (with/without conditions) and prepares
#    data matrix and parameter indexing.
#
# 5. get_data(nPart, N, params, prevent_zero_accuracy)
#    Wrapper that generates data for a single parameter set.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
# sample_trial: Generate a SINGLE TRIAL by emulating the DDM random walk
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
sample_trial <- function(a, v, dt, max_steps){
    x <- 0  # Initialize the evidence accumulator at 0 (starting point)

    # Generate random step-size noise in advance - Noise around the drift rate
    random_dev <- rnorm(max_steps)

    # Scale the random noise and drift rate to approximate a continuous Wiener process
    noise <- random_dev * sqrt(dt)
    drift <- v * dt
    
    # Compute the evidence sampled on each step
    steps <- drift + noise
    
    # Start the random walk
    for(i in 2:max_steps){
          # Update the evidence accumulator on each step
          x = x + steps[i]
          
          # Check if a decision boundary has been reached
          # The boundaries are at +a/2 and -a/2
          if(abs(x)>=(a/2)){
            break  # Stop the simulation if a boundary is reached
          }
    }

    # Create the output:
    # - RT: Reaction time (number of steps (i) - initial step * time step size (dt))
    # - C: Final evidence accumulated (x)  
    output <- list("RT" = (i-1)*dt, "C"  = x)
return(output)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
# sample_many_trials: Generate multiple trials for one participant/condition
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
sample_many_trials <- function(a,v,t,n){
    # Set random walk emulation parameters
    dt = 0.001  # Time step size (in seconds)
    max_steps = 10 / dt  # Maximum number of steps (~10 seconds)
    
    # Initialize vectors to store results
    rt = rep(NA,n)        # Reaction times
    accuracy = rep(NA,n)  # Accuracy (1 = correct, 0 = incorrect)
    
    # Generate n trials
    for(i in 1:n){
          # Generate a single DDM trial
          X <- sample_trial(a, v, dt, max_steps)
          
          # Store the reaction time (without non-decision time for now)
          rt[i] <- X$RT
          
          # Determine accuracy based on the final evidence value
          if(X$C > 0){
              # Hitting the upper boundary is considered a "correct" response
              accuracy[i] <- 1
          } else {
              # Hitting the lower boundary is considered an "incorrect" response
              accuracy[i] <- 0
          }
    }

    RT <- rt + t  # Add the non-decision time (t) to all reaction times

    # Create output data frame with reaction times and accuracy
    output <- data.frame("RT" = RT, "accuracy" = accuracy)
return(output)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
# sample_data: MAIN FUNCTION - Generate complete dataset for simulation studies
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
sample_data <- function(nPart, nTrials = NA, parameter_set, 
                        nTrialsPerCondition = NA, prevent_zero_accuracy = TRUE){
  
  design <- identify_design(nPart, nTrials, nTrialsPerCondition, parameter_set)
            data <- design$data
            nObs <- design$nObs
            nCols <- design$nCols
            n_subIndex <- design$n_subIndex
            col_accuracy <- design$col_accuracy
            col_rt <- design$col_rt
            N <- design$N
            colNames <- design$colNames
            adjusted_parameter_set <- design$adjusted_parameter_set
            cell_index <- design$cell_index
            n_par_sets <- length(adjusted_parameter_set$bound)

  for(i in 1:n_par_sets){
        # Identify rows for current participant
        this.cell <- which(cell_index==i)

        params <- list(bound = parameter_set$bound[i],
                      drift = parameter_set$drift[i],
                      nondt = parameter_set$nondt[i])
        # Generate dataset for this participant using their specific parameters
        # First generate the dataset once
        temp <- get_data(nPart, N, params, prevent_zero_accuracy)
        data[this.cell,col_accuracy] <- temp$accuracy
        data[this.cell,col_rt] <- temp$RT
  }
  # Convert to matrix and add column names
  data <- as.matrix(data)
  colnames(data) <- colNames
  
return(data)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
#                          AUXILIARY FUNCTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
# identify_design: Determine design structure and prepare data matrix
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
identify_design <- function(nPart, nTrials, nTrialsPerCondition, parameter_set){
      withinSubject <- !is.na(nTrialsPerCondition)
      if(withinSubject){
          nObs <- nPart*nTrialsPerCondition*2
          nCols <- 4
          colNames <- c("sub","cond","rt", "accuracy")
          n_subIndex <- nTrialsPerCondition*2
          col_accuracy <- 4
          col_rt <- 3
          N <- nTrialsPerCondition
          adjusted_parameter_set <- list(
                  bound = rep(parameter_set$bound, each=2),
                  drift = parameter_set$drift,
                  nondt = rep(parameter_set$nondt, each=2))
      } else {
          nObs <- nPart*nTrials
          nCols <- 3
          colNames <- c("sub","rt", "accuracy")
          n_subIndex <- nTrials
          col_accuracy <- 3
          col_rt <- 2
          N <- nTrials
          adjusted_parameter_set <- parameter_set
      }

      data <- matrix(NA, ncol=nCols, nrow=nObs)  
      data[,1] <- rep(1:nPart, each=n_subIndex)
      cell_index <- data[,1]

      if(withinSubject){
          data[,2] <- rep(rep(c(1,0), each=nTrialsPerCondition), nPart)
          cell_index <- rep(1:nPart*2, each=nTrialsPerCondition)
      }

  return(list(data = data, nObs = nObs,
              nCols = nCols, N = N,
              n_subIndex = n_subIndex,
              col_accuracy = col_accuracy,
              col_rt = col_rt,
              colNames = colNames,
              adjusted_parameter_set = adjusted_parameter_set,
              cell_index = cell_index))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
# get_data: Generate data for a single parameter set with edge case handling
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
get_data <- function(nPart, N, params, prevent_zero_accuracy){
      # Generate dataset for this participant using their specific parameters
      # First generate the dataset once
      temp <- sample_many_trials(a = params$bound, v = params$drift, t = params$nondt, n = N)
      accuracy <- temp$accuracy

      # If prevent_zero_accuracy is TRUE and we got all zeros, keep trying
      while(sum(accuracy)==0 && prevent_zero_accuracy){
        temp <- sample_many_trials(a = params$bound, v = params$drift, t = params$nondt, n = N)
        accuracy <- temp$accuracy
      }

return(temp)
}