data {
  int I;
  int K;
  int P_nu;
  int P_alpha;
  int P_tau;

  array[I, K] int  J;
  array[I, K] int  C;
  matrix[I, K] MRT;
  matrix[I, K] VRT;

  array[I, K] row_vector[P_nu]    X_nu;
  array[I, K] row_vector[P_alpha] X_alpha;
  array[I, K] row_vector[P_tau]   X_tau;

  matrix[2, P_nu]    prior_beta_nu;
  matrix[2, P_alpha] prior_beta_alpha;
  matrix[2, P_tau]   prior_beta_tau;

  vector[2] prior_sigma_nu;
  vector[2] prior_sigma_alpha;
  vector[2] prior_sigma_tau;
  vector[2] prior_sigma_v;
}

parameters {
  // group-level regression coefficients 
  vector[P_nu]    beta_nu;
  vector[P_alpha] beta_alpha;
  vector[P_tau]   beta_tau;

  // group-level SDs 
  // bounds removed from alpha and tau: hard bounds break non-centered
  // the informative pilot priors keep values in sensible range instead
  real<lower=0> sigma_nu;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_tau;
  real<lower=0> sigma_v;

  // raw offsets: these are all standard normal, no structure
  // sampler sees simple geometry regardless of sigma values
  matrix[I, K] nu_raw;
  matrix[I, K] alpha_raw;
  matrix[I, K] tau_raw;
  vector[I]    v_raw;   // random slopes also non-centered
}

transformed parameters {
  // reconstruct individual parameters from raw offsets + group-level mean/SD
  // the actual nu/alpha/tau/v values 
  matrix[I, K] nu;
  matrix[I, K] alpha;
  matrix[I, K] tau;
  vector[I]    v;

  // v first because alpha depends on it
  v = sigma_v * v_raw;
  // v[i] = 0 + sigma_v * v_raw[i]
  // equivalent to v[i] ~ normal(0, sigma_v) but decoupled from sigma_v

  for (i in 1:I) {
    for (k in 1:K) {
      real mean_nu_ik    = X_nu[i, k]    * beta_nu;
      real mean_alpha_ik = X_alpha[i, k] * beta_alpha + v[i] * X_alpha[i, k][2];
      real mean_tau_ik   = X_tau[i, k]   * beta_tau;

      // reconstruct: individual value = group mean + (group SD * raw offset)
      nu[i, k]    = mean_nu_ik    + sigma_nu    * nu_raw[i, k];
      alpha[i, k] = mean_alpha_ik + sigma_alpha * alpha_raw[i, k];
      tau[i, k]   = mean_tau_ik   + sigma_tau   * tau_raw[i, k];
    }
  }
}

model {
  // priors on regression coefficients 
  beta_nu[1]    ~ normal(prior_beta_nu[1, 1],    prior_beta_nu[2, 1]);
  beta_alpha[1] ~ normal(prior_beta_alpha[1, 1], prior_beta_alpha[2, 1]);
  beta_tau[1]   ~ normal(prior_beta_tau[1, 1],   prior_beta_tau[2, 1]);

  if (P_nu > 1)
    beta_nu[2:P_nu] ~ normal(to_vector(prior_beta_nu[1, 2:P_nu]),
                             to_vector(prior_beta_nu[2, 2:P_nu]));
  if (P_alpha > 1)
    beta_alpha[2:P_alpha] ~ normal(to_vector(prior_beta_alpha[1, 2:P_alpha]),
                                   to_vector(prior_beta_alpha[2, 2:P_alpha]));
  if (P_tau > 1)
    beta_tau[2:P_tau] ~ normal(to_vector(prior_beta_tau[1, 2:P_tau]),
                               to_vector(prior_beta_tau[2, 2:P_tau]));

  // priors on SDs
  sigma_nu    ~ normal(prior_sigma_nu[1],    prior_sigma_nu[2]);
  sigma_alpha ~ normal(prior_sigma_alpha[1], prior_sigma_alpha[2]);
  sigma_tau   ~ normal(prior_sigma_tau[1],   prior_sigma_tau[2]);
  sigma_v     ~ normal(prior_sigma_v[1],     prior_sigma_v[2]);

  // raw offset priors — replaces the old individual parameter priors
  // no means, no sigmas: sampler always sees standard normal geometry
  to_vector(nu_raw)    ~ normal(0, 1);
  to_vector(alpha_raw) ~ normal(0, 1);
  to_vector(tau_raw)   ~ normal(0, 1);
  v_raw                ~ normal(0, 1);

  // EZ likelihood, nu/alpha/tau come from transformed parameters
  for (i in 1:I) {
    for (k in 1:K) {
      real e         = exp(-alpha[i, k] * nu[i, k]);
      real pi_c      = 1.0 / (1.0 + e);
      real mu_rt     = tau[i, k] + (alpha[i, k] / (2.0 * nu[i, k])) *
                       ((1.0 - e) / (1.0 + e));
      real sigma2_rt = (alpha[i, k] / (2.0 * pow(nu[i, k], 3))) *
                       ((1.0 - 2.0 * alpha[i, k] * nu[i, k] * e - e^2) /
                        square(1.0 + e));

      C[i, k]   ~ binomial(J[i, k], pi_c);
      MRT[i, k] ~ normal(mu_rt, sqrt(sigma2_rt / C[i, k]));
      VRT[i, k] ~ normal(sigma2_rt,
                         sqrt(2.0 * square(sigma2_rt) / fmax(C[i, k] - 1, 1)));
    }
  }
}