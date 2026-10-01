
// v, a, t0, w, sv free
// sw = 0 and st0 = 0 
// nu: intercept + condition RANDOM (nu_intercept[i], nu_c[i])
// alpha : intercept only RANDOM (alpha_intercept[i])
// t0: intercept only RANDOM (seconds)
// t0[i] ~ normal(beta_tau, sigma_tau), centered and set below the person's fastest RT
// because the Wiener density needs t0 < every RT
// w: intercept only RANDOM (w_logit[i])
// sv: intercept only RANDOM (log_sv[i]), 

// acc = +1 --> upper boundary (correct)
// acc = -1 --> lower boundary (error)

functions {
  real partial_sum(array[] int idx_slice, int start, int end,
                   vector rt, array[] int acc, array[] int pid, vector condition,
                   vector a, vector t0, vector w, vector sv,
                   vector nu_intercept, vector nu_c) {
    real lp = 0;
    for (n in 1:size(idx_slice)) {
      int i = idx_slice[n];
      int p = pid[i];
      real v = nu_intercept[p] + nu_c[p] * condition[i];
      if (acc[i] == 1) {
        lp += wiener_lpdf(rt[i] | a[p], t0[p], w[p], v, sv[p], 0.0, 0.0);
      } else {
        lp += wiener_lpdf(rt[i] | a[p], t0[p], 1 - w[p], -v, sv[p], 0.0, 0.0);
      }
    }
    return lp;
  }
}

data {
  int<lower=1> N;                           // trials
  int<lower=1> I;                           // persons
  array[N] int<lower=1, upper=I> pid;
  vector<lower=0>[N] rt;
  array[N] int<lower=-1, upper=1> acc;      // +1 upper (correct), -1 lower (error)
  vector[N] condition;                      // +1 / -1
  vector<lower=0>[I] min_rt;                // fastest RT per person
  int<lower=1> grainsize;

  // priors: [1] = mean, [2] = SD
  vector[2] prior_beta_nu_intercept;        // normal
  vector[2] prior_beta_nu_c;                // normal
  vector[2] prior_beta_alpha;               // normal, log scale
  vector[2] prior_beta_tau;                 // normal, seconds
  vector[2] prior_beta_w;                   // normal, logit scale
  vector[2] prior_beta_sv;                  // normal, log scale

  vector[2] prior_sigma_nu_intercept;       // normal (half-normal bc of lower=0)
  vector[2] prior_sigma_nu_c;
  vector[2] prior_sigma_alpha_intercept;    // log scale
  vector[2] prior_sigma_tau_intercept;      // seconds
  vector[2] prior_sigma_w;                  // logit scale
  vector[2] prior_sigma_sv;                 // log scale
}

transformed data {
  array[N] int idx;
  for (n in 1:N) idx[n] = n;
}

parameters {
  real beta_nu_intercept;
  real beta_nu_c;
  real beta_alpha;
  real<lower=0.05, upper=0.45> beta_tau;    // population mean t0 (s)
  real beta_w;
  real beta_sv;

  real<lower=0> sigma_nu_intercept;
  real<lower=0> sigma_nu_c;
  real<lower=0> sigma_alpha_intercept;
  real<lower=0> sigma_tau_intercept;
  real<lower=0> sigma_w;
  real<lower=0> sigma_sv;

  vector[I] z_nu_intercept;
  vector[I] z_nu_c;
  vector[I] z_alpha_intercept;
  vector<lower=0, upper=min_rt>[I] t0;      // person t0 in seconds, below each fastest RT
  vector[I] z_w;
  vector[I] z_sv;
}

transformed parameters {
  vector[I] nu_intercept    = beta_nu_intercept + sigma_nu_intercept    * z_nu_intercept;
  vector[I] nu_c            = beta_nu_c         + sigma_nu_c            * z_nu_c;
  vector[I] alpha_intercept = beta_alpha        + sigma_alpha_intercept * z_alpha_intercept;  // log scale
  vector[I] w_logit         = beta_w            + sigma_w               * z_w;
  vector[I] log_sv          = beta_sv           + sigma_sv              * z_sv;

  vector[I] a  = exp(alpha_intercept);
  vector[I] w  = inv_logit(w_logit);
  vector[I] sv = exp(log_sv);
}

model {
  beta_nu_intercept ~ normal(prior_beta_nu_intercept[1], prior_beta_nu_intercept[2]);
  beta_nu_c         ~ normal(prior_beta_nu_c[1],         prior_beta_nu_c[2]);
  beta_alpha        ~ normal(prior_beta_alpha[1],        prior_beta_alpha[2]);
  beta_tau          ~ normal(prior_beta_tau[1],          prior_beta_tau[2]);
  beta_w            ~ normal(prior_beta_w[1],            prior_beta_w[2]);
  beta_sv           ~ normal(prior_beta_sv[1],           prior_beta_sv[2]);

  sigma_nu_intercept    ~ normal(prior_sigma_nu_intercept[1],    prior_sigma_nu_intercept[2]);
  sigma_nu_c            ~ normal(prior_sigma_nu_c[1],            prior_sigma_nu_c[2]);
  sigma_alpha_intercept ~ normal(prior_sigma_alpha_intercept[1], prior_sigma_alpha_intercept[2]);
  sigma_tau_intercept   ~ normal(prior_sigma_tau_intercept[1],   prior_sigma_tau_intercept[2]);
  sigma_w               ~ normal(prior_sigma_w[1],               prior_sigma_w[2]);
  sigma_sv              ~ normal(prior_sigma_sv[1],              prior_sigma_sv[2]);

  z_nu_intercept    ~ std_normal();
  z_nu_c            ~ std_normal();
  z_alpha_intercept ~ std_normal();
  // t0: normal population distribution in seconds
  // t0 < min_rt keeps sampler where Wiener density is defined
  t0 ~ normal(beta_tau, sigma_tau_intercept);
  z_w               ~ std_normal();
  z_sv              ~ std_normal();

  target += reduce_sum(partial_sum, idx, grainsize,
                       rt, acc, pid, condition,
                       a, t0, w, sv, nu_intercept, nu_c);
}

generated quantities {
  // population summaries on the natural scale, for comparison with the simulation
  real w_pop_median     = inv_logit(beta_w);
  real sv_pop_median    = exp(beta_sv);
  real alpha_pop_median = exp(beta_alpha);
}
