data {
  int I;
  int K;
  int P_nu;     // 2: intercept + condition
  int P_alpha;  // 2: intercept + condition
  int P_tau;    // 1: intercept only

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
  // random slope on nu condition effect
}

parameters {
  vector[P_nu]    beta_nu;    // [intercept, condition]
  vector[P_alpha] beta_alpha; // [intercept, condition]
  vector[P_tau]   beta_tau;   // [intercept]

  real<lower=0.1, upper=0.8>  sigma_alpha;
  real<lower=0>               sigma_nu;
  real<lower=0.05, upper=0.6> sigma_tau;
  real<lower=0>               sigma_v;

  vector[I] z_v;  // non-centred: raw N(0,1) slopes; v[i] = sigma_v * z_v[i]

  matrix<lower=0.1, upper=3.0>[I, K] alpha;
  matrix[I, K]                        nu;
  matrix<lower=0.05>[I, K]            tau;
}

transformed parameters {
  vector[I] v = sigma_v * z_v;  // actual random slopes
}

model {
  beta_nu[1]    ~ normal(prior_beta_nu[1, 1],    prior_beta_nu[2, 1]);
  beta_alpha[1] ~ normal(prior_beta_alpha[1, 1], prior_beta_alpha[2, 1]);
  beta_tau[1]   ~ normal(prior_beta_tau[1, 1],   prior_beta_tau[2, 1]);

  if (P_nu > 1)
    beta_nu[2:P_nu] ~ normal(to_vector(prior_beta_nu[1, 2:P_nu]),
                             to_vector(prior_beta_nu[2, 2:P_nu]));
  if (P_alpha > 1)
    beta_alpha[2:P_alpha] ~ normal(to_vector(prior_beta_alpha[1, 2:P_alpha]),
                                   to_vector(prior_beta_alpha[2, 2:P_alpha]));

  sigma_nu    ~ normal(prior_sigma_nu[1],    prior_sigma_nu[2]);
  sigma_alpha ~ normal(prior_sigma_alpha[1], prior_sigma_alpha[2]);
  sigma_tau   ~ normal(prior_sigma_tau[1],   prior_sigma_tau[2]);
  sigma_v     ~ normal(prior_sigma_v[1],     prior_sigma_v[2]);
  z_v         ~ normal(0, 1);  // non-centred prior

  for (i in 1:I) {
    for (k in 1:K) {
      real mu_nu_ik    = X_nu[i, k] * beta_nu + v[i] * X_nu[i, k][2];
      real mu_alpha_ik = X_alpha[i, k] * beta_alpha;
      real mu_tau_ik   = X_tau[i, k]   * beta_tau;

      nu[i, k]    ~ normal(mu_nu_ik,    sigma_nu);
      alpha[i, k] ~ normal(mu_alpha_ik, sigma_alpha);
      tau[i, k]   ~ normal(mu_tau_ik,   sigma_tau);

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

generated quantities {
  matrix[I, K] mean_nu_ik;
  matrix[I, K] mean_alpha_ik;
  matrix[I, K] mean_tau_ik;

  for (i in 1:I) {
    for (k in 1:K) {
      mean_nu_ik[i, k]    = X_nu[i, k] * beta_nu + v[i] * X_nu[i, k][2];
      mean_alpha_ik[i, k] = X_alpha[i, k] * beta_alpha;
      mean_tau_ik[i, k]   = X_tau[i, k]   * beta_tau;
    }
  }
}