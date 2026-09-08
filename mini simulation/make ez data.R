build_ez_data <- function(sim_data) {
  ez_long <- sim_data %>%
    group_by(subj, condition, resp_type) %>%
    summarise(
      J = n(),
      C = sum(acc == 1),
      MRT = mean(rt[acc == 1]),
      VRT = var(rt[acc == 1]),
      .groups = "drop"
    ) %>%
    filter(C >= 5)  # need enough correct trials to estimate MRT/VRT reliably
  
  # keep only subjects with all 4 cells present
  complete_subjects <- ez_long %>%
    count(subj) %>%
    filter(n == 4) %>%
    pull(subj)
  
  n_dropped <- length(unique(sim_data$subj)) - length(complete_subjects)
  if (n_dropped > 0) {
    warning(sprintf("Dropping %d subjects with incomplete cells (C < 5 in at least one cell)", n_dropped))
  }
  
  ez_long %>% filter(subj %in% complete_subjects)
}

build_stan_data_ez <- function(rep_id, data_dir) {
  rep_dir <- file.path(data_dir, sprintf("rep%02d", rep_id))
  sim_data <- read.csv(file.path(rep_dir, "sim_data_full_ddm.csv"))
  
  ez_data <- build_ez_data(sim_data)
  
  subjects_used <- sort(unique(ez_data$subj))
  I <- length(subjects_used)
  
  # define the 4 cells in a fixed order 
  cells <- tibble(
    k = 1:4,
    condition = c(1, 1, -1, -1),
    resp_type = c(1, -1, 1, -1)
  )
  K <- nrow(cells)
  
  ez_data$pid <- as.integer(factor(ez_data$subj, levels = subjects_used))
  ez_data <- ez_data %>%
    left_join(cells, by = c("condition", "resp_type"))
  
  # initialize I x K arrays
  J_arr   <- matrix(NA_integer_, I, K)
  C_arr   <- matrix(NA_integer_, I, K)
  MRT_arr <- matrix(NA_real_,    I, K)
  VRT_arr <- matrix(NA_real_,    I, K)
  
  for (r in seq_len(nrow(ez_data))) {
    i <- ez_data$pid[r]; k <- ez_data$k[r]
    J_arr[i, k]   <- ez_data$J[r]
    C_arr[i, k]   <- ez_data$C[r]
    MRT_arr[i, k] <- ez_data$MRT[r]
    VRT_arr[i, k] <- ez_data$VRT[r]
  }
  
  stopifnot(!anyNA(J_arr), !anyNA(C_arr), !anyNA(MRT_arr), !anyNA(VRT_arr))
  
  # design matrices, one row per (i,k), stored as array[I,K] row_vector[P]
  X_nu <- array(0, dim = c(I, K, 4))     # intercept, cond, resp, cond*resp
  X_alpha <- array(0, dim = c(I, K, 3))  # intercept, resp, cond
  X_tau <- array(1, dim = c(I, K, 1))    # intercept only
  
  for (k in seq_len(K)) {
    cond <- cells$condition[k]; resp <- cells$resp_type[k]
    X_nu[, k, ]    <- matrix(c(1, cond, resp, cond*resp), I, 4, byrow = TRUE)
    X_alpha[, k, ] <- matrix(c(1, resp, cond),             I, 3, byrow = TRUE)
  }
  
  list(
    I = I, K = K, P_nu = 4, P_alpha = 3, P_tau = 1,
    J = J_arr, C = C_arr, MRT = MRT_arr, VRT = VRT_arr,
    X_nu = X_nu, X_alpha = X_alpha, X_tau = X_tau,
    subjects_used = subjects_used
  )
}


sd <- build_stan_data_ez(rep_id = 1, data_dir = "data_full7")


str(sd, max.level = 1)
cat("I =", sd$I, "| K =", sd$K, "\n")


sapply(list(J=sd$J, C=sd$C, MRT=sd$MRT, VRT=sd$VRT), function(x) sum(is.na(x)))
# should all be 0


range(sd$J)     
range(sd$C)    
all(sd$C <= sd$J)   # should be TRUE
range(sd$MRT)   # should look like plausible RTs (seconds)
range(sd$VRT)   # should be positive
any(sd$VRT <= 0)  # should be FALSE


sd$X_nu[1, , ]     # print participant 1's 4 rows: should be
#      intercept  cond  resp  cond*resp
# k=1:     1        1     1        1
# k=2:     1        1    -1       -1
# k=3:     1       -1     1       -1
# k=4:     1       -1    -1        1

sd$X_alpha[1, , ]  # should be
#      intercept  resp  cond
# k=1:     1        1     1
# k=2:     1       -1     1
# k=3:     1        1    -1
# k=4:     1       -1    -1

sd$X_tau[1, , ]    # should just be a column of 1s, 4 rows


length(sd$subjects_used)   
# if far below 100, check the warning message from build_ez_data() for how many were dropped and why