//   nu : intercept + condition + response_type + condition*response_type
//   alpha : intercept + response_type
//   tau : intercept only

// 4 cells 
// k=1: cond= 1, resp= 1   k=2: cond= 1, resp=-1
// k=3: cond=-1, resp= 1   k=4: cond=-1, resp=-1

// resp=+1 = post-correct, resp=-1 = post-error

// I expect
// post-error = higher alpha (more caution)
// alpha_r slope should be negative (alpha goes up when resp goes down)

// beta_alpha[2] must be unconstraint (can be negative)

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

  // Priors: row 1 = means, row 2 = SDs
  matrix[2, P_nu]    prior_beta_nu;
  matrix[2, P_alpha] prior_beta_alpha;
  matrix[2, P_tau]   prior_beta_tau;

  vector[2] prior_sigma_nu;
  vector[2] prior_sigma_alpha;
  vector[2] prior_sigma_tau;

  vector[2] prior_sigma_nu_c;
  vector[2] prior_sigma_nu_r;
  vector[2] prior_sigma_nu_cr;
  vector[2] prior_sigma_alpha_r;
}

parameters {
  // population means
  // beta_nu: all four elements can be any sign
  vector[P_nu]    beta_nu;

  // beta_alpha: unconstraint, the slope (element 2) can be negative
  // only alpha[i,k] is constrained to [0.1, 3]
  vector[P_alpha] beta_alpha;

  // beta_tau: intercept only, bounded to plausible range
  real<lower=0.05, upper=0.45> beta_tau;

  // SDs (within-person across cell variation)
  real<lower=0>              sigma_nu;
  real<lower=0, upper=0.8>   sigma_alpha;  // sonst crash
  real<lower=0>              sigma_tau;    

  // random slope SDs (between person variation in effects)
  real<lower=0> sigma_nu_c;
  real<lower=0> sigma_nu_r;
  real<lower=0> sigma_nu_cr;
  real<lower=0> sigma_alpha_r;

  // NCP 
  vector[I] z_nu_c;
  vector[I] z_nu_r;
  vector[I] z_nu_cr;
  vector[I] z_alpha_r;

  // individual 
  matrix[I, K]                        nu;    // drift rates 
  matrix<lower=0.1, upper=3.0>[I, K]  alpha; // boundary sep
  matrix<lower=0.05>[I, K]            tau;   // non decision times
}

transformed parameters {
  // participant specific slopes: population mean + sigma * z_i  (NCP)
  vector[I] nu_c    = beta_nu[2]    + sigma_nu_c    * z_nu_c;
  vector[I] nu_r    = beta_nu[3]    + sigma_nu_r    * z_nu_r;
  vector[I] nu_cr   = beta_nu[4]    + sigma_nu_cr   * z_nu_cr;
  vector[I] alpha_r = beta_alpha[2] + sigma_alpha_r * z_alpha_r;
}

model {
  // population means 
  for (p in 1:P_nu)
    beta_nu[p] ~ normal(prior_beta_nu[1, p], prior_beta_nu[2, p]);
  for (p in 1:P_alpha)
    beta_alpha[p] ~ normal(prior_beta_alpha[1, p], prior_beta_alpha[2, p]);
  beta_tau ~ normal(prior_beta_tau[1, 1], prior_beta_tau[2, 1]);

  // residual SDs 
  sigma_nu    ~ normal(prior_sigma_nu[1],    prior_sigma_nu[2]);
  sigma_alpha ~ normal(prior_sigma_alpha[1], prior_sigma_alpha[2]);
  sigma_tau   ~ normal(prior_sigma_tau[1],   prior_sigma_tau[2]);

  // random slope SDs 
  sigma_nu_c    ~ normal(prior_sigma_nu_c[1],    prior_sigma_nu_c[2]);
  sigma_nu_r    ~ normal(prior_sigma_nu_r[1],    prior_sigma_nu_r[2]);
  sigma_nu_cr   ~ normal(prior_sigma_nu_cr[1],   prior_sigma_nu_cr[2]);
  sigma_alpha_r ~ normal(prior_sigma_alpha_r[1], prior_sigma_alpha_r[2]);

  // NCP  
  z_nu_c    ~ normal(0, 1);
  z_nu_r    ~ normal(0, 1);
  z_nu_cr   ~ normal(0, 1);
  z_alpha_r ~ normal(0, 1);

  // parameters in every call and EZ likelihood 
  for (i in 1:I) {
    for (k in 1:K) {
      // linear predictor for nu:
      // intercept (beta_nu[1] * X[1]=1) + person's condition slope * X[2]
      //   + person's resp_type slope * X[3] + person's cr slope * X[4]
      real mu_nu_ik = beta_nu[1]  * X_nu[i, k][1]
                    + nu_c[i]     * X_nu[i, k][2]
                    + nu_r[i]     * X_nu[i, k][3]
                    + nu_cr[i]    * X_nu[i, k][4];

      // linear predictor for alpha:
      // intercept + person's resp_type slope * X[2]
      // slope alpha_r[i] can be negative --> post-error cells get higher alpha
      real mu_alpha_ik = beta_alpha[1] * X_alpha[i, k][1]
                       + alpha_r[i]    * X_alpha[i, k][2];

      // tau: same across cells per person
      real mu_tau_ik = beta_tau;

      // cell parameters drawn from person means
      nu[i, k]    ~ normal(mu_nu_ik,    sigma_nu);
      alpha[i, k] ~ normal(mu_alpha_ik, sigma_alpha);
      tau[i, k]   ~ normal(mu_tau_ik,   sigma_tau) T[0.05,];

      // EZ equations 
      // fmax prevents division by zero if nu drifts near zero during warmup
      real a = alpha[i, k];
      real v = fmax(nu[i, k], 0.01);  // EZ needs positive drift toward correct

      real e      = exp(-a * v);
      real pi_c   = 1.0 / (1.0 + e);
      real mu_rt  = tau[i, k] + (a / (2.0 * v)) * ((1.0 - e) / (1.0 + e));
      real sig2rt = fmax(
          (a / (2.0 * pow(v, 3)))
            * ((1.0 - 2.0*a*v*e - square(e)) / square(1.0 + e)),
          1e-8);

      C[i, k]   ~ binomial(J[i, k], pi_c);
      MRT[i, k] ~ normal(mu_rt, sqrt(sig2rt / fmax(C[i, k], 1)));
      VRT[i, k] ~ normal(sig2rt,
                          sqrt(2.0 * square(sig2rt) / fmax(C[i, k] - 1, 1)));
    }
  }
}



