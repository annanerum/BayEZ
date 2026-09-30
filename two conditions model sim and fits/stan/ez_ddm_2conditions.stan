//   nu    : intercept + condition         RANDOM (nu_intercept[i], nu_c[i])   
//   alpha : intercept only                RANDOM (alpha_intercept[i]), LOG scale  
//   tau   : intercept only                RANDOM (tau_intercept[i]), seconds

// 2 cells per person                                                              
// k=1: cond = +1     k=2: cond = -1

// EZ uses the mean (MRT) and variance (VRT) of CORRECT RTs, and s = 1
// (same scaling as rtdists::rdiffusion)

data {
  int<lower=1> I;                       // persons
  int<lower=1> K;                       // cells (= 2)
  vector[K] cond;                       // condition code per cell (+1 / -1)

  array[I, K] int<lower=1> J;           // trials per cell
  array[I, K] int<lower=0> C;           // correct responses per cell
  matrix[I, K] MRT;                     // mean correct RT per cell (s)
  matrix[I, K] VRT;                     // variance of correct RTs per cell (s^2)

  // priors: [1] = mean, [2] = SD 
  vector[2] prior_beta_nu_intercept;    // normal
  vector[2] prior_beta_nu_c;            // normal
  vector[2] prior_beta_alpha;           // normal, log scale     
  vector[2] prior_beta_tau;             // normal, seconds

  // normal (half-normal via lower=0) instead of gamma,
  // same form as in the first simulation (model2_v2_ncp)
  vector[2] prior_sigma_nu_intercept;   // normal
  vector[2] prior_sigma_nu_c;           // normal (was prior_sigma_v)
  vector[2] prior_sigma_alpha_intercept;// normal (log scale)
  vector[2] prior_sigma_tau_intercept;  // normal
}

parameters {
  // population means
  real beta_nu_intercept;                       // was beta_nu[1] before (fixed intercept)
  real beta_nu_c;                               // was beta_nu[2]
  real beta_alpha;                              // log boundary sep intercept
  real<lower=0.05, upper=0.45> beta_tau;      

  // between-person SDs
  real<lower=0> sigma_nu_intercept;             
  real<lower=0> sigma_nu_c;
  real<lower=0> sigma_alpha_intercept;          
  real<lower=0> sigma_tau_intercept;

  
  vector[I] z_nu_intercept;                     
  vector[I] z_nu_c;
  vector[I] z_alpha_intercept;
  vector[I] z_tau_intercept;
}

transformed parameters {
  vector[I] nu_intercept    = beta_nu_intercept + sigma_nu_intercept    * z_nu_intercept;
  vector[I] nu_c            = beta_nu_c         + sigma_nu_c            * z_nu_c;
  vector[I] alpha_intercept = beta_alpha        + sigma_alpha_intercept * z_alpha_intercept; // log scale
  vector[I] tau_intercept   = beta_tau          + sigma_tau_intercept   * z_tau_intercept;
}

model {
  // population means
  beta_nu_intercept ~ normal(prior_beta_nu_intercept[1], prior_beta_nu_intercept[2]);
  beta_nu_c         ~ normal(prior_beta_nu_c[1],         prior_beta_nu_c[2]);
  beta_alpha        ~ normal(prior_beta_alpha[1],        prior_beta_alpha[2]);
  beta_tau          ~ normal(prior_beta_tau[1],          prior_beta_tau[2]);

  // between-person SDs 
  // normal as in the first simulation
  sigma_nu_intercept    ~ normal(prior_sigma_nu_intercept[1],    prior_sigma_nu_intercept[2]);
  sigma_nu_c            ~ normal(prior_sigma_nu_c[1],            prior_sigma_nu_c[2]);
  sigma_alpha_intercept ~ normal(prior_sigma_alpha_intercept[1], prior_sigma_alpha_intercept[2]);
  sigma_tau_intercept   ~ normal(prior_sigma_tau_intercept[1],   prior_sigma_tau_intercept[2]);

  // NCP
  z_nu_intercept    ~ std_normal();
  z_nu_c            ~ std_normal();
  z_alpha_intercept ~ std_normal();
  z_tau_intercept   ~ std_normal();

  // EZ likelihood per person and cell
  for (i in 1:I) {
    real a   = exp(alpha_intercept[i]);   
    real tau = tau_intercept[i];          

    for (k in 1:K) {
      // fmax prevents division by zero if drift gets near zero
      real v = fmax(nu_intercept[i] + nu_c[i] * cond[k], 0.01);

      real e      = exp(-a * v);
      real pi_c   = 1.0 / (1.0 + e);
      real mu_rt  = tau + (a / (2.0 * v)) * ((1.0 - e) / (1.0 + e));
      real sig2rt = fmax(
          (a / (2.0 * pow(v, 3)))
            * ((1.0 - 2.0 * a * v * e - square(e)) / square(1.0 + e)),
          1e-8);

      C[i, k]   ~ binomial(J[i, k], pi_c);
      MRT[i, k] ~ normal(mu_rt, sqrt(sig2rt / fmax(C[i, k], 1)));
      VRT[i, k] ~ normal(sig2rt,
                         sqrt(2.0 * square(sig2rt) / fmax(C[i, k] - 1, 1)));
    }
  }
}

generated quantities {
  vector[I] alpha_person = exp(alpha_intercept);
  real alpha_pop_median  = exp(beta_alpha);   
}
