data {
  int<lower=1> I;
  int<lower=1> K;
  int<lower=1> P_nu;              // columns in X_nu (here: intercept, covariate, condition)
  int<lower=1> P_alpha;           // columns in X_alpha (here: intercept, condition)
  int<lower=1> P_tau;             // columns in X_tau (here: intercept, condition)
  array[I, K] int<lower=1> J;
  array[I, K] int<lower=0> C;
  matrix[I, K] MRT;
  matrix[I, K] VRT;
  array[I, K] row_vector[P_nu]    X_nu;
  array[I, K] row_vector[P_alpha] X_alpha;
  array[I, K] row_vector[P_tau]   X_tau;
}

parameters {

  // regression coefficients: [intercept, predictors] per parameter
  vector[P_nu]    beta_nu;      // [mu_nu, beta_nu_1, beta_nu_2]
  vector[P_alpha] beta_alpha;   // [mu_alpha, beta_alpha_2]
  vector[P_tau]   beta_tau;     // [mu_tau, beta_tau_2]

  // group-level SDs
  real<lower=0.1, upper=0.8>  sigma_alpha;
  real<lower=0>               sigma_nu;
  real<lower=0.05, upper=0.6> sigma_tau;

  // by-participant random slopes for condition effect on alpha
  vector[I]     v;              // participant deviations from beta_alpha_2
  real<lower=0> sigma_v;        // SD of random slopes

  // individual-level parameters (I x K)
  matrix<lower=0.1, upper=3.0>[I, K] alpha;
  matrix[I, K]                       nu;
  matrix<lower=0.05>[I, K]           tau;
}

model {
  // prior on random-slope SD
  sigma_v ~ normal(0, 0.5);
  v       ~ normal(0, sigma_v);

  for (i in 1:I) {
    for (k in 1:K) {
      // design-matrix-predicted means
      real mean_nu_ik    = X_nu[i, k]    * beta_nu;
      // fixed condition effect + participant-specific deviation
      real mean_alpha_ik = X_alpha[i, k] * beta_alpha + v[i] * X_alpha[i, k][2];
      real mean_tau_ik   = X_tau[i, k]   * beta_tau;

      // hierarchical distributions
      nu[i, k]    ~ normal(mean_nu_ik,    sigma_nu);
      alpha[i, k] ~ normal(mean_alpha_ik, sigma_alpha);
      tau[i, k]   ~ normal(mean_tau_ik,   sigma_tau);

      // EZ likelihood
      real u         = exp(-alpha[i, k] * nu[i, k]);
      real pi_c      = 1.0 / (1.0 + u);
      real mu_rt     = tau[i, k] + (alpha[i, k] / (2.0 * nu[i, k])) *
                       ((1.0 - u) / (1.0 + u));
      real sigma2_rt = (alpha[i, k] / (2.0 * pow(nu[i, k], 3))) *
                       ((1.0 - 2.0 * alpha[i, k] * nu[i, k] * u - u^2) /
                        square(1.0 + u));

      C[i, k]   ~ binomial(J[i, k], pi_c);
      MRT[i, k] ~ normal(mu_rt,     sqrt(sigma2_rt / C[i, k]));
      VRT[i, k] ~ normal(sigma2_rt, sqrt(2.0 * square(sigma2_rt) / (C[i, k] - 1)));
    }
  }
}