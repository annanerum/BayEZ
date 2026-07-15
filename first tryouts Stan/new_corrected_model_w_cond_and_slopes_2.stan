data {
  int I;        // number of participants
  int K;        // number of conditions
  int P_nu;     // numer of predictors for nu --> 3 intercept + covariate + condition
  int P_alpha;  // numer of predictors for alpha --> 2 intercept + condition 
  int P_tau;    // number of predictors for tau --> 2 intercept + condition 
  
  // two dimensional arrays, I rows, K columns
  array[I, K] int  J;   // trials per participant per condition
  // J[i,k] = how many trials ppn I did in cond K
  // because its integers it must be an array 
  
  array[I, K] int  C;   // correct responses per participant per condition
  // C[i,k] = how many of those trials are correct 
  
  matrix[I, K] MRT;     // mean RT of correct responses
  // observed mean RT of correct responses for for ppn i in cond k 
  // could also be array [i,k] real but stan prefers matrix 
  
  matrix[I, K] VRT;     // variance of RT of correct responses
  // observed variance of RT for ppn i iin cond k 


  // design matrices 
  // each cell [i,k] holds a row vector = a rowof predictors for that ppn-cond combi 
  array[I, K] row_vector[P_nu]    X_nu;    // design matrix for drift rate (log scale)
  // X_nu[i,k] = row vector with 3 predictor: 1, x_i and cond_dummy_k
  array[I, K] row_vector[P_alpha] X_alpha; // design matrix for boundary separation
  // X_alpha[i,k] = row vector with 2: 1, cond_dummy_k
  array[I, K] row_vector[P_tau]   X_tau;   // design matrix for non-decision time
  // multiplying this row vector by the column vector of the coefficient later (like X_nu[i,k]*beta_nu) w
  // gives us a dot product = the linear predictor for that cell 

  // priors passed from R, row 1 = prior means for each regression coeff, row 2 = prior SDs
  // 2xP matrices of priors hyperparameters from pilot data
  matrix[2, P_nu]    prior_beta_nu;    // beta_nu priors on log scale
  matrix[2, P_alpha] prior_beta_alpha;
  matrix[2, P_tau]   prior_beta_tau;
  // prior_beta_nu[1, 1] = prior mean for the nu intercept
  // prior_beta_nu[2, 1] = prior SD
  // prior_beta_nu[1, 2] = prior mean for covariate slope
  // prior_beta_nu[2, 2] = prior SD for slope
  //                intercept   covariate   condition
  // row 1 (means)  [  0.0        0.0         0.0   ]
  // row 2 (SDs)    [  0.58       0.58        0.58  ]

  
  // length 2 vectors --> SD priors parameters
  //  element 1 = prior mean, 2 = prior SD
  vector[2] prior_sigma_nu;
  // prior_sigma_nu[1] --> prior mean (center), prior_sigma_nu[2] --> prior SD (how wide)
  vector[2] prior_sigma_alpha;
  vector[2] prior_sigma_tau;
  vector[2] prior_sigma_v;
  // v = every persons deviation from the mean condition effect on boundary sep (b_alpha = 0.2)
  // b_alpha + v[i] (v[i] centered at 0 and some become + and some -)
  // --> ind differences in how much ppl adapt their caution
  // sigma_v = SD of all the v[i] values across ppns --> how much ppl differ from each other in how (much) their adjust their caution
}



// sampled from posterior
parameters {
  
  // column vectors of regression coeff --> one per predictor for each parameter 
  // can be any real number bc regression coeff can be - or + 
  vector[P_nu]    beta_nu;    // raw-scale drift rate regression coefficients
  // beta_nu --> 3: raw-scale intercept for nu, covariate effect on nu, condition effect on nu
  vector[P_alpha] beta_alpha; // boundary separation regression coefficients
  // beta_alpha --> 2: boundary sep intercept, condition effect on boundary sep 
  vector[P_tau]   beta_tau;   // non-decision time regression coefficients
  // also intercept and boundary sep

  
  // constraint
  real<lower=0.1, upper=0.8>  sigma_alpha; // between-participant SD for alpha
  real<lower=0>               sigma_nu;    // between-participant SD for nu (raw scale)
  real<lower=0.05, upper=0.6> sigma_tau;   // between-participant SD for tau
  real<lower=0>               sigma_v;     // SD of participant-level random slopes


  vector[I] v; // by-participant random slopes on condition effect for alpha
  // one value per ppn (length-I vector)
  // personal deviation from average cond effect on boundary sep 
  
  
  // individual parameters
  // IxK matrix of individual parameters with constraints where needed
  matrix<lower=0.1, upper=3.0>[I, K] alpha; // individual boundary separations
  matrix[I, K]                        nu;    // individual drift rates -- unconstrained, negative = drift toward wrong boundary
  matrix<lower=0.05>[I, K]            tau;   // individual non-decision times
}



model {
  // priors on regression intercepts
  beta_nu[1]    ~ normal(prior_beta_nu[1, 1],    prior_beta_nu[2, 1]);
  beta_alpha[1] ~ normal(prior_beta_alpha[1, 1], prior_beta_alpha[2, 1]);
  beta_tau[1]   ~ normal(prior_beta_tau[1, 1],   prior_beta_tau[2, 1]);
  // each intercept (element [1] of its coefficient-vector) gets a normal prior 
  // mean and SD come from the pilot-derived prior matrices
  // row 1 column 1 = mean, row 2 column 1 = SD

  // priors on regression slopes
  if (P_nu > 1)
    beta_nu[2:P_nu] ~ normal(to_vector(prior_beta_nu[1, 2:P_nu]),
                             to_vector(prior_beta_nu[2, 2:P_nu]));
  if (P_alpha > 1)
    beta_alpha[2:P_alpha] ~ normal(to_vector(prior_beta_alpha[1, 2:P_alpha]),
                                   to_vector(prior_beta_alpha[2, 2:P_alpha]));
  if (P_tau > 1)
    beta_tau[2:P_tau] ~ normal(to_vector(prior_beta_tau[1, 2:P_tau]),
                               to_vector(prior_beta_tau[2, 2:P_tau]));


  // priors on group-level SDs
  sigma_nu    ~ normal(prior_sigma_nu[1],    prior_sigma_nu[2]);
  sigma_alpha ~ normal(prior_sigma_alpha[1], prior_sigma_alpha[2]);
  sigma_tau   ~ normal(prior_sigma_tau[1],   prior_sigma_tau[2]);
  // half-normal because sigma parameters are constrained to > 0 so already truncated

  // random slope model
  sigma_v ~ normal(prior_sigma_v[1], prior_sigma_v[2]); // hyperprior on slope variability
  v       ~ normal(0, sigma_v);                         // participant slopes centred at zero

  // participant × condition loop
  for (i in 1:I) {
    for (k in 1:K) {
      // linear predictors (raw scale for nu now -- not log scale)
      real mean_nu_ik    = X_nu[i, k]    * beta_nu;
      // dot product: [1, x_i, cond_dummy_k] · [beta_nu[1], beta_nu[2], beta_nu[3]]
      // negative mean_nu_ik allowed --> drift toward wrong boundary

      real mean_alpha_ik = X_alpha[i, k] * beta_alpha + v[i] * X_alpha[i, k][2];
      // fixed effect + ppn random slope (only active in cond 2 where dummy = 1)

      real mean_tau_ik   = X_tau[i, k]   * beta_tau;


      // individual parameter priors
      nu[i, k]    ~ normal(mean_nu_ik,    sigma_nu);    // normal again
      alpha[i, k] ~ normal(mean_alpha_ik, sigma_alpha);
      tau[i, k]   ~ normal(mean_tau_ik,   sigma_tau);


      // EZ formulas
      real e         = exp(-alpha[i, k] * nu[i, k]);
      // when nu < 0: -alpha*nu > 0, so e > 1 --> pi_c < 0.5 (below-chance accuracy)

      real pi_c      = 1.0 / (1.0 + e);

      real mu_rt     = tau[i, k] + (alpha[i, k] / (2.0 * nu[i, k])) *
                       ((1.0 - e) / (1.0 + e));

      real sigma2_rt = (alpha[i, k] / (2.0 * pow(nu[i, k], 3))) *
                       ((1.0 - 2.0 * alpha[i, k] * nu[i, k] * e - e^2) /
                        square(1.0 + e));
      // sigma2_rt can go negative when nu is negative and large 
      // fmax prevents sqrt() from receiving a negative value
      real sigma2_rt_safe = fmax(sigma2_rt, 1e-10);


      // likelihood
      C[i, k]   ~ binomial(J[i, k], pi_c);
      MRT[i, k] ~ normal(mu_rt, sqrt(sigma2_rt_safe / C[i, k]));
      VRT[i, k] ~ normal(sigma2_rt_safe,
                         sqrt(2.0 * square(sigma2_rt_safe) / fmax(C[i, k] - 1, 1)));
    }
  }
}
