// DDM5 for the Krakow (Grabowska) flanker data
// v, a, t0, w, sv free; sw = 0 and st0 = 0

// same effect structure as grabowska_ez_model_new.stan:
// nu: intercept + condition + response_type + condition*response_type (all random)
// alpha: intercept (random) + response_type (random) + condition (fixed)
// t0: intercept only random (seconds), centered, below the person's fastest RT
// w: intercept only random (logit scale)
// sv: intercept only random (log scale)

// Priors: Grabowska et al. (2025) Appendix C when parameter exists in theirs
// gamma for the SDs Grabowska gives Gamma priors, half-normal for the rest

// 4 cells, contrast coded +-1 (same as EZ script):
// k=1: cond= 1, resp= 1   k=2: cond= 1, resp=-1
// k=3: cond=-1, resp= 1   k=4: cond=-1, resp=-1
// cond: congruent = 1, incongruent = -1
// resp: post-correct = 1, post-error = -1

// acc = +1 --> upper boundary (correct)
// acc = -1 --> lower boundary (error)

functions {
  real partial_sum(array[] int idx_slice, int start, int end,
                   vector rt, array[] int acc, array[] int pid, array[] int cell,
                   matrix alpha, matrix nu, vector t0, vector w, vector sv) {
    real lp = 0;
    for (n in 1:size(idx_slice)) {
      int i = idx_slice[n];
      int p = pid[i];
      int k = cell[i];
      if (acc[i] == 1) {
        lp += wiener_lpdf(rt[i] | alpha[p, k], t0[p], w[p], nu[p, k], sv[p], 0.0, 0.0);
      } else {
        lp += wiener_lpdf(rt[i] | alpha[p, k], t0[p], 1 - w[p], -nu[p, k], sv[p], 0.0, 0.0);
      }
    }
    return lp;
  }
}

data {
  int<lower=1> N;                           // trials
  int<lower=1> I;                           // persons
  array[N] int<lower=1, upper=I> pid;
  array[N] int<lower=1, upper=4> cell;      // cell index k (see header)
  vector<lower=0>[N] rt;                    // seconds
  array[N] int<lower=-1, upper=1> acc;      // +1 upper (correct), -1 lower (error)
  vector<lower=0>[I] min_rt;                // fastest RT per person
  int<lower=1> grainsize;

  // priors: [1] = mean, [2] = SD
  vector[2] prior_beta_nu_intercept;        // normal
  vector[2] prior_beta_nu_c;
  vector[2] prior_beta_nu_r;
  vector[2] prior_beta_nu_cr;
  vector[2] prior_beta_alpha;               // normal, natural scale, truncated [0, 3]
  vector[2] prior_beta_alpha_r;             // normal, natural scale
  vector[2] prior_beta_alpha_c;             // normal, natural scale
  vector[2] prior_beta_tau;                 // normal, seconds, truncated [0, 0.4]
  vector[2] prior_beta_w;                   // normal, logit scale
  vector[2] prior_beta_sv;                  // normal, log scale

  vector[2] prior_sigma_nu_intercept;       // gamma
  vector[2] prior_sigma_nu_c;               // gamma
  vector[2] prior_sigma_nu_r;               // normal (half-normal bc of lower=0)
  vector[2] prior_sigma_nu_cr;              // normal (half-normal)
  vector[2] prior_sigma_alpha_intercept;    // gamma, natural scale
  vector[2] prior_sigma_alpha_r;            // normal (half-normal), natural scale
  vector[2] prior_sigma_tau_intercept;      // gamma, seconds
  vector[2] prior_sigma_w;                  // logit scale
  vector[2] prior_sigma_sv;                 // log scale
}

transformed data {
  int K = 4;
  vector[4] cell_cond = [ 1,  1, -1, -1]';
  vector[4] cell_resp = [ 1, -1,  1, -1]';

  array[N] int idx;
  for (n in 1:N) idx[n] = n;
}

parameters {
  // population means
  real beta_nu_intercept;
  real beta_nu_c;
  real beta_nu_r;
  real beta_nu_cr;
  real<lower=0, upper=3> beta_alpha;        // natural scale, App C: N(1,1)[0,3]
  real beta_alpha_r;                        // can be negative (post-error caution?)
  real beta_alpha_c;                        // fixed effect
  real<lower=0, upper=0.4> beta_tau;        // population mean t0 (s), App C: N(0.1,0.2)[0,0.4]
  real beta_w;
  real beta_sv;

  // between-person SDs
  real<lower=0> sigma_nu_intercept;
  real<lower=0> sigma_nu_c;
  real<lower=0> sigma_nu_r;
  real<lower=0> sigma_nu_cr;
  real<lower=0> sigma_alpha_intercept;
  real<lower=0> sigma_alpha_r;
  real<lower=0> sigma_tau_intercept;
  real<lower=0> sigma_w;
  real<lower=0> sigma_sv;

  // NCP
  vector[I] z_nu_intercept;
  vector[I] z_nu_c;
  vector[I] z_nu_r;
  vector[I] z_nu_cr;
  vector[I] z_alpha_intercept;
  vector[I] z_alpha_r;
  vector<lower=0, upper=min_rt>[I] t0;      // person t0 in seconds, below each fastest RT
  vector[I] z_w;
  vector[I] z_sv;
}

transformed parameters {
  // person-level coefficients (names match the EZ model where possible)
  vector[I] nu_intercept    = beta_nu_intercept + sigma_nu_intercept    * z_nu_intercept;
  vector[I] nu_c            = beta_nu_c         + sigma_nu_c            * z_nu_c;
  vector[I] nu_r            = beta_nu_r         + sigma_nu_r            * z_nu_r;
  vector[I] nu_cr           = beta_nu_cr        + sigma_nu_cr           * z_nu_cr;
  vector[I] alpha_intercept = beta_alpha        + sigma_alpha_intercept * z_alpha_intercept;
  vector[I] alpha_r         = beta_alpha_r      + sigma_alpha_r         * z_alpha_r;
  vector[I] w_logit         = beta_w            + sigma_w               * z_w;
  vector[I] log_sv          = beta_sv           + sigma_sv              * z_sv;

  vector[I] w  = inv_logit(w_logit);
  vector[I] sv = exp(log_sv);

  // person x cell drift and boundary 
  matrix[I, K] nu;
  matrix[I, K] alpha;
  for (k in 1:K) {
    nu[, k]    = nu_intercept + nu_c * cell_cond[k] + nu_r * cell_resp[k]
                 + nu_cr * (cell_cond[k] * cell_resp[k]);
    alpha[, k] = alpha_intercept + alpha_r * cell_resp[k] + beta_alpha_c * cell_cond[k];
  }
}

model {
  beta_nu_intercept ~ normal(prior_beta_nu_intercept[1], prior_beta_nu_intercept[2]);
  beta_nu_c         ~ normal(prior_beta_nu_c[1],         prior_beta_nu_c[2]);
  beta_nu_r         ~ normal(prior_beta_nu_r[1],         prior_beta_nu_r[2]);
  beta_nu_cr        ~ normal(prior_beta_nu_cr[1],        prior_beta_nu_cr[2]);
  beta_alpha        ~ normal(prior_beta_alpha[1],        prior_beta_alpha[2]);
  beta_alpha_r      ~ normal(prior_beta_alpha_r[1],      prior_beta_alpha_r[2]);
  beta_alpha_c      ~ normal(prior_beta_alpha_c[1],      prior_beta_alpha_c[2]);
  beta_tau          ~ normal(prior_beta_tau[1],          prior_beta_tau[2]);
  beta_w            ~ normal(prior_beta_w[1],            prior_beta_w[2]);
  beta_sv           ~ normal(prior_beta_sv[1],           prior_beta_sv[2]);

  sigma_nu_intercept    ~ gamma(prior_sigma_nu_intercept[1],     prior_sigma_nu_intercept[2]);
  sigma_nu_c            ~ gamma(prior_sigma_nu_c[1],             prior_sigma_nu_c[2]);
  sigma_nu_r            ~ normal(prior_sigma_nu_r[1],            prior_sigma_nu_r[2]);
  sigma_nu_cr           ~ normal(prior_sigma_nu_cr[1],           prior_sigma_nu_cr[2]);
  sigma_alpha_intercept ~ gamma(prior_sigma_alpha_intercept[1],  prior_sigma_alpha_intercept[2]);
  sigma_alpha_r         ~ normal(prior_sigma_alpha_r[1],         prior_sigma_alpha_r[2]);
  sigma_tau_intercept   ~ gamma(prior_sigma_tau_intercept[1],    prior_sigma_tau_intercept[2]);
  sigma_w               ~ normal(prior_sigma_w[1],               prior_sigma_w[2]);
  sigma_sv              ~ normal(prior_sigma_sv[1],              prior_sigma_sv[2]);

  z_nu_intercept    ~ std_normal();
  z_nu_c            ~ std_normal();
  z_nu_r            ~ std_normal();
  z_nu_cr           ~ std_normal();
  z_alpha_intercept ~ std_normal();
  z_alpha_r         ~ std_normal();
  // t0: normal population distribution in seconds 
  // need t0 < min_rt bc it keeps sampler where density is still defined
  t0 ~ normal(beta_tau, sigma_tau_intercept);
  z_w               ~ std_normal();
  z_sv              ~ std_normal();

  target += reduce_sum(partial_sum, idx, grainsize,
                       rt, acc, pid, cell,
                       alpha, nu, t0, w, sv);
}

generated quantities {
  real w_pop_median  = inv_logit(beta_w);
  real sv_pop_median = exp(beta_sv);

  // population-level cell values (average person)
  vector[4] nu_pop_cell;
  vector[4] alpha_pop_cell;
  for (k in 1:4) {
    nu_pop_cell[k]    = beta_nu_intercept + beta_nu_c * cell_cond[k] + beta_nu_r * cell_resp[k]
                        + beta_nu_cr * cell_cond[k] * cell_resp[k];
    alpha_pop_cell[k] = beta_alpha + beta_alpha_r * cell_resp[k] + beta_alpha_c * cell_cond[k];
  }
}
