data {
  int I;        // number of participants
  int K;        // number of conditions
  int P_nu;     // number of predictors for nu --> 3: intercept + covariate + condition
  int P_alpha;  // number of predictors for alpha --> 2: intercept + condition
  int P_tau;    // number of predictors for tau --> 2: intercept + condition

  // two dimensional arrays, I rows, K columns
  array[I, K] int  J;   // trials per participant per condition
  // J[i,k] = how many trials ppn i did in cond k
  // because its integers it must be an array

  array[I, K] int  C;   // correct responses per participant per condition
  // C[i,k] = how many of those trials are correct

  matrix[I, K] MRT;     // mean RT of correct responses
  // observed mean RT of correct responses for ppn i in cond k

  matrix[I, K] VRT;     // variance of RT of correct responses
  // observed variance of RT for ppn i in cond k


  // design matrices
  // each cell [i,k] holds a row vector = a row of predictors for that ppn-cond combination
  array[I, K] row_vector[P_nu]    X_nu;    // design matrix for drift rate (raw scale)
  // X_nu[i,k] = row vector with 3 predictors: 1, x_i, cond_dummy_k
  array[I, K] row_vector[P_alpha] X_alpha; // design matrix for boundary separation
  // X_alpha[i,k] = row vector with 2: 1, cond_dummy_k
  array[I, K] row_vector[P_tau]   X_tau;   // design matrix for non-decision time
  // multiplying this row vector by the column vector of coefficients gives the linear predictor


  // priors passed from R: row 1 = prior means, row 2 = prior SDs
  matrix[2, P_nu]    prior_beta_nu;    // beta_nu priors (raw scale)
  matrix[2, P_alpha] prior_beta_alpha;
  matrix[2, P_tau]   prior_beta_tau;

  // length-2 vectors: element 1 = prior mean, element 2 = prior SD
  vector[2] prior_sigma_nu;
  vector[2] prior_sigma_alpha;
  vector[2] prior_sigma_tau;
  vector[2] prior_sigma_v;
}


parameters {

  // regression coefficients (unconstrained: can be negative)
  vector[P_nu]    beta_nu;    // raw-scale drift rate regression coefficients
  vector[P_alpha] beta_alpha; // boundary separation regression coefficients
  vector[P_tau]   beta_tau;   // non-decision time regression coefficients

  // group-level SDs
  real<lower=0.1, upper=0.8>  sigma_alpha;
  real<lower=0>               sigma_nu;    // raw-scale SD for nu
  real<lower=0.05, upper=0.6> sigma_tau;
  real<lower=0>               sigma_v;

  vector[I] v; // by-participant random slopes on condition effect for alpha

  // individual parameters
  matrix<lower=0.1, upper=3.0>[I, K] alpha; // individual boundary separations
  matrix[I, K]                        nu;    // individual drift rates -- unconstrained
  // negative nu = drift toward wrong boundary (below-chance performance)
  // near-zero nu handled via nu_safe below -- never fed raw nu to EZ formulas
  matrix<lower=0.05>[I, K]            tau;   // individual non-decision times
}


model {

  // safety threshold: nu is never fed to EZ formulas closer to zero than this
  // prevents mu_rt diverging to +/-Inf as nu --> 0
  // value of 0.05 matches the repeat-loop threshold in R data generation
  real eps = 0.05;

  // priors on regression intercepts
  beta_nu[1]    ~ normal(prior_beta_nu[1, 1],    prior_beta_nu[2, 1]);
  beta_alpha[1] ~ normal(prior_beta_alpha[1, 1], prior_beta_alpha[2, 1]);
  beta_tau[1]   ~ normal(prior_beta_tau[1, 1],   prior_beta_tau[2, 1]);

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
  // half-normal because sigma parameters are constrained to > 0
  sigma_nu    ~ normal(prior_sigma_nu[1],    prior_sigma_nu[2]);
  sigma_alpha ~ normal(prior_sigma_alpha[1], prior_sigma_alpha[2]);
  sigma_tau   ~ normal(prior_sigma_tau[1],   prior_sigma_tau[2]);

  // random slope model
  sigma_v ~ normal(prior_sigma_v[1], prior_sigma_v[2]);
  v       ~ normal(0, sigma_v);  // participant slopes centred at zero

  // participant x condition loop
  for (i in 1:I) {
    for (k in 1:K) {

      // linear predictors
      real mean_nu_ik    = X_nu[i, k]    * beta_nu;
      // dot product: [1, x_i, cond_dummy_k] . [beta_nu[1], beta_nu[2], beta_nu[3]]
      // negative mean_nu_ik allowed --> drift toward wrong boundary

      real mean_alpha_ik = X_alpha[i, k] * beta_alpha + v[i] * X_alpha[i, k][2];
      // fixed effect + ppn random slope (only active in cond 2 where dummy = 1)

      real mean_tau_ik   = X_tau[i, k]   * beta_tau;


      // individual parameter priors
      nu[i, k]    ~ normal(mean_nu_ik,    sigma_nu);
      alpha[i, k] ~ normal(mean_alpha_ik, sigma_alpha);
      tau[i, k]   ~ normal(mean_tau_ik,   sigma_tau);


      // nu_safe: keep nu away from zero before feeding into EZ formulas
      // maps nu to the same sign but with |nu| >= eps
      // sign(nu): +1 if nu > 0, -1 if nu < 0 -- but sign(0) = 0 in Stan so we guard with fmax
      // fmax(fabs(nu[i,k]), eps) ensures the magnitude is at least eps
      // together: same direction as nu, magnitude at least eps
      // this is the only place nu=0 would cause problems (mu_rt --> Inf)
      real nu_safe = (nu[i, k] >= 0 ? 1.0 : -1.0) * fmax(fabs(nu[i, k]), eps);


      // EZ formulas -- all use nu_safe, never raw nu[i,k]
      real e         = exp(-alpha[i, k] * nu_safe);
      // when nu_safe > 0 and large: e near 0 --> high accuracy, fast decision
      // when nu_safe < 0: e > 1 --> below-chance accuracy

      real pi_c      = 1.0 / (1.0 + e);
      // proportion correct: > 0.5 when nu_safe > 0, < 0.5 when nu_safe < 0

      real mu_rt     = tau[i, k] + (alpha[i, k] / (2.0 * nu_safe)) *
                       ((1.0 - e) / (1.0 + e));
      // predicted mean RT -- finite for all nu_safe != 0
      // with nu_safe < 0: both (alpha/2*nu_safe) and (1-e)/(1+e) are negative
      // negative * negative = positive decision time, so mu_rt > tau

      real sigma2_rt = (alpha[i, k] / (2.0 * pow(nu_safe, 3))) *
                       ((1.0 - 2.0 * alpha[i, k] * nu_safe * e - e^2) /
                        square(1.0 + e));
      // always positive for nu_safe != 0 (two negatives cancel when nu_safe < 0)
      // no fmax needed here -- sigma2_rt cannot be negative or zero when nu_safe != 0


      // likelihood
      C[i, k]   ~ binomial(J[i, k], pi_c);
      MRT[i, k] ~ normal(mu_rt, sqrt(sigma2_rt / C[i, k]));
      VRT[i, k] ~ normal(sigma2_rt,
                         sqrt(2.0 * square(sigma2_rt) / fmax(C[i, k] - 1, 1)));
    }
  }
}
