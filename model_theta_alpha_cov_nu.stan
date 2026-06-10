// Theta reparameterization for random slope on alpha
// alpha[i,k] = mean_alpha[i] + cond[k] * theta_alpha[i]   (random intercept + random slope)
// nu[i,k]    = mean_nu[i]    + x[i]    * beta_nu[2]        (random intercept + fixed covariate)
// tau[i,k]   = mean_tau[i]                                  (random intercept only)
// No within-cell residual variance -- all variation is participant-level

data {
  int I;
  int K;
  int P_nu;     // 2: intercept + covariate
  int P_alpha;  // 2: [pop mean intercept, pop mean slope]
  int P_tau;    // 1: intercept only

  array[I, K] int  J;
  array[I, K] int  C;
  matrix[I, K] MRT;
  matrix[I, K] VRT;

  array[I, K] row_vector[P_nu]    X_nu;    // [1, x[i]] -- X_nu[i,k][2] = covariate
  array[I, K] row_vector[P_alpha] X_alpha; // [1, cond] -- X_alpha[i,k][2] = cond dummy
  array[I, K] row_vector[P_tau]   X_tau;   // [1]

  matrix[2, P_nu]    prior_beta_nu;
  matrix[2, P_alpha] prior_beta_alpha;
  matrix[2, P_tau]   prior_beta_tau;

  vector[2] prior_sigma_nu;
  vector[2] prior_sigma_alpha;
  vector[2] prior_sigma_tau;
  vector[2] prior_sigma_v;
}

parameters {
  vector[P_nu]    beta_nu;    // [1]=pop mean nu intercept,    [2]=fixed covariate slope
  vector[P_alpha] beta_alpha; // [1]=pop mean alpha intercept, [2]=pop mean alpha slope
  vector[P_tau]   beta_tau;   // [1]=pop mean tau intercept

  real<lower=0>               sigma_nu;    // SD of participant nu intercepts
  real<lower=0.1, upper=0.8>  sigma_alpha; // SD of participant alpha intercepts
  real<lower=0.05, upper=0.6> sigma_tau;   // SD of participant tau intercepts
  real<lower=0>               sigma_v;     // SD of participant alpha slopes

  vector[I]             mean_nu;    // participant nu intercepts
  vector<lower=0.1>[I]  mean_alpha; // participant alpha intercepts (= alpha in cond 1)
  vector<lower=0.05>[I] mean_tau;   // participant tau values

  vector[I] theta_alpha; // participant alpha slopes (total condition effect, fixed + random)
}

transformed parameters {
  matrix[I, K] nu;
  matrix[I, K] alpha;
  matrix[I, K] tau;

  for (i in 1:I) {
    for (k in 1:K) {
      nu[i, k]    = mean_nu[i]    + X_nu[i, k][2]    * beta_nu[2];
      // covariate slope is fixed (beta_nu[2]), intercept varies per participant
      // X_nu[i,k][2] = x[i] -- same across conditions for a given participant
      alpha[i, k] = mean_alpha[i] + X_alpha[i, k][2] * theta_alpha[i];
      // condition slope is random per participant (theta_alpha[i])
      // cond 1 (k=1): alpha[i,1] = mean_alpha[i]
      // cond 2 (k=2): alpha[i,2] = mean_alpha[i] + theta_alpha[i]
      tau[i, k]   = mean_tau[i];
      // no condition effect on tau -- same value in both conditions
    }
  }
}

model {
  // fixed effects priors
  beta_nu[1]    ~ normal(prior_beta_nu[1, 1],    prior_beta_nu[2, 1]);
  beta_alpha[1] ~ normal(prior_beta_alpha[1, 1], prior_beta_alpha[2, 1]);
  beta_tau[1]   ~ normal(prior_beta_tau[1, 1],   prior_beta_tau[2, 1]);

  if (P_nu > 1)
    beta_nu[2:P_nu] ~ normal(to_vector(prior_beta_nu[1, 2:P_nu]),
                             to_vector(prior_beta_nu[2, 2:P_nu]));
  if (P_alpha > 1)
    beta_alpha[2:P_alpha] ~ normal(to_vector(prior_beta_alpha[1, 2:P_alpha]),
                                   to_vector(prior_beta_alpha[2, 2:P_alpha]));

  // variance priors
  sigma_nu    ~ normal(prior_sigma_nu[1],    prior_sigma_nu[2]);
  sigma_alpha ~ normal(prior_sigma_alpha[1], prior_sigma_alpha[2]);
  sigma_tau   ~ normal(prior_sigma_tau[1],   prior_sigma_tau[2]);
  sigma_v     ~ normal(prior_sigma_v[1],     prior_sigma_v[2]);

  // participant-level random effects
  mean_nu     ~ normal(beta_nu[1],    sigma_nu);
  mean_alpha  ~ normal(beta_alpha[1], sigma_alpha);
  mean_tau    ~ normal(beta_tau[1],   sigma_tau);
  theta_alpha ~ normal(beta_alpha[2], sigma_v);
  // theta_alpha[i] is the participant's total condition effect on alpha
  // centred at beta_alpha[2] (population mean), with SD sigma_v

  // EZ likelihood
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

// alpha, nu, tau are already in transformed parameters 
generated quantities {
  matrix[I, K] mean_nu_ik;
  matrix[I, K] mean_alpha_ik;
  matrix[I, K] mean_tau_ik;

  for (i in 1:I) {
    for (k in 1:K) {
      mean_nu_ik[i, k]    = nu[i, k];
      mean_alpha_ik[i, k] = alpha[i, k];
      mean_tau_ik[i, k]   = tau[i, k];
    }
  }
}