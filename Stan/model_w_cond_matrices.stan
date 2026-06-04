data {
  int<lower=1> I;
  int<lower=1> K;
  int<lower=1> P_nu;
  int<lower=1> P_alpha;
  int<lower=1> P_tau;
  array[I, K] int<lower=1> J;
  array[I, K] int<lower=0> C;
  matrix[I, K] MRT;
  matrix[I, K] VRT;
  array[I, K] row_vector[P_nu]    X_nu;
  array[I, K] row_vector[P_alpha] X_alpha;
  array[I, K] row_vector[P_tau]   X_tau;
}

parameters {
  vector[P_nu]    beta_nu;
  vector[P_alpha] beta_alpha;
  vector[P_tau]   beta_tau;

  real<lower=0.1, upper=0.8>  sigma_alpha;
  real<lower=0>               sigma_nu;
  real<lower=0.05, upper=0.6> sigma_tau;

  matrix<lower=0.1, upper=3.0>[I, K] alpha;
  matrix[I, K]                       nu;
  matrix<lower=0.05>[I, K]           tau;
}

model {
  // priors on intercepts
  beta_nu[1]    ~ normal(0.0, 0.5);
  beta_alpha[1] ~ normal(1.5, 0.2);
  beta_tau[1]   ~ normal(0.3, 0.6);

  // priors on covariate and condition effect
  beta_nu[2]    ~ normal(0, 1.0);    // beta_1^(nu): covariate effect
  beta_nu[3]    ~ normal(0, 1.0);    // beta_2^(nu): condition effect
  beta_alpha[2] ~ normal(0, 0.5);    // beta^(alpha): condition effect
  beta_tau[2]   ~ normal(0, 0.05);   // beta^(tau): condition effect

  // priors on group-level SDs
  sigma_alpha ~ normal(0.25, 0.1);
  sigma_nu    ~ normal(0.30, 0.1);
  sigma_tau   ~ normal(0.15, 0.1);

  for (i in 1:I) {
    for (k in 1:K) {
      real mean_nu_ik    = X_nu[i, k]    * beta_nu;
      real mean_alpha_ik = X_alpha[i, k] * beta_alpha;
      real mean_tau_ik   = X_tau[i, k]   * beta_tau;

      nu[i, k]    ~ normal(mean_nu_ik,    sigma_nu);
      alpha[i, k] ~ normal(mean_alpha_ik, sigma_alpha);
      tau[i, k]   ~ normal(mean_tau_ik,   sigma_tau);

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
