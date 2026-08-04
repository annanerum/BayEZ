//   w  : fixed at 0.5 (like Grabowska)
//   sw : fixed at 0
//   sv : random per participant (log scale for positivity)
//   each person has their own level of trial-to-trial drift noise: sv[i] = exp(mu_sv + sigma_sv*z_sv[i])
//   st0: fixed at 0


functions {
  real partial_sum(array[] int idx_slice, int start, int end,
                    vector rt, array[] int acc, array[] int pid,
                    vector condition, vector resp_type,
                    vector alpha_intercept, real beta_alpha2, real beta_alpha3,
                    real beta_nu1, vector nu_c, vector nu_r, vector nu_cr,
                    vector t0, vector sv) {
    real lp = 0;
    for (n in 1:size(idx_slice)) {
      int i = idx_slice[n];
      int p = pid[i];
      real log_a = alpha_intercept[p] + beta_alpha2 * resp_type[i] + beta_alpha3 * condition[i];
      real a     = exp(log_a);
      real v = beta_nu1
             + nu_c[p]  * condition[i]
             + nu_r[p]  * resp_type[i]
             + nu_cr[p] * condition[i] * resp_type[i];
      if (acc[i] == 1) {
        lp += wiener_lpdf(rt[i] | a, t0[p], 0.5, v, sv[p], 0.0, 0.0);
      } else {
        lp += wiener_lpdf(rt[i] | a, t0[p], 0.5, -v, sv[p], 0.0, 0.0);
      }
    }
    return lp;
  }
}

data {
  int<lower=1> N;
  int<lower=1> I;
  array[N] int<lower=1, upper=I> pid;
  vector<lower=0>[N] rt;
  array[N] int<lower=-1, upper=1> acc;
  vector[N] condition;
  vector[N] resp_type;
  vector<lower=0>[I] min_rt;
  int<lower=1> grainsize;

  matrix[2, 4] prior_beta_nu;
  vector[2] prior_beta_alpha1;
  vector[2] prior_beta_alpha2;
  vector[2] prior_beta_alpha3;
  vector[2] prior_beta_tau;
  vector[2] prior_sigma_nu_c;
  vector[2] prior_sigma_nu_r;
  vector[2] prior_sigma_nu_cr;
  vector[2] prior_sigma_alpha_intercept;
  vector[2] prior_sigma_tau_intercept;

  vector[2] prior_beta_sv;             // prior on mu_sv (log scale)
  vector[2] prior_sigma_sv_intercept;  // prior on sigma_sv_intercept
}

transformed data {
  array[N] int idx;
  for (n in 1:N) idx[n] = n;
}

parameters {
  vector[4] beta_nu;
  real beta_alpha1;
  real beta_alpha2;
  real beta_alpha3;
  real beta_tau;
  real<lower=0> sigma_nu_c;
  real<lower=0> sigma_nu_r;
  real<lower=0> sigma_nu_cr;
  real<lower=0> sigma_alpha_intercept;
  real<lower=0> sigma_tau_intercept;
  vector[I] z_nu_c;
  vector[I] z_nu_r;
  vector[I] z_nu_cr;
  vector[I] z_alpha_intercept;
  vector[I] z_tau_intercept;

  real beta_sv;                     // population mean, log scale
  real<lower=0> sigma_sv_intercept; // between-person SD, log scale
  vector[I] z_sv;
}

transformed parameters {
  vector[I] nu_c  = beta_nu[2] + sigma_nu_c  * z_nu_c;
  vector[I] nu_r  = beta_nu[3] + sigma_nu_r  * z_nu_r;
  vector[I] nu_cr = beta_nu[4] + sigma_nu_cr * z_nu_cr;
  vector[I] alpha_intercept   = beta_alpha1 + sigma_alpha_intercept * z_alpha_intercept;
  vector[I] tau_intercept_raw = beta_tau    + sigma_tau_intercept   * z_tau_intercept;
  vector[I] t0 = min_rt .* inv_logit(tau_intercept_raw);

  vector[I] sv_intercept_log = beta_sv + sigma_sv_intercept * z_sv;
  vector[I] sv = exp(sv_intercept_log);  // positivity 
}

model {
  beta_nu ~ normal(to_vector(prior_beta_nu[1]), to_vector(prior_beta_nu[2]));
  beta_alpha1 ~ normal(prior_beta_alpha1[1], prior_beta_alpha1[2]);
  beta_alpha2 ~ normal(prior_beta_alpha2[1], prior_beta_alpha2[2]);
  beta_alpha3 ~ normal(prior_beta_alpha3[1], prior_beta_alpha3[2]);
  beta_tau    ~ normal(prior_beta_tau[1],    prior_beta_tau[2]);
  sigma_nu_c  ~ gamma(prior_sigma_nu_c[1],  prior_sigma_nu_c[2]);
  sigma_nu_r  ~ gamma(prior_sigma_nu_r[1],  prior_sigma_nu_r[2]);
  sigma_nu_cr ~ gamma(prior_sigma_nu_cr[1], prior_sigma_nu_cr[2]);
  sigma_alpha_intercept ~ gamma(prior_sigma_alpha_intercept[1], prior_sigma_alpha_intercept[2]);
  sigma_tau_intercept   ~ gamma(prior_sigma_tau_intercept[1],   prior_sigma_tau_intercept[2]);
  z_nu_c  ~ normal(0, 1);
  z_nu_r  ~ normal(0, 1);
  z_nu_cr ~ normal(0, 1);
  z_alpha_intercept ~ normal(0, 1);
  z_tau_intercept   ~ normal(0, 1);

  beta_sv            ~ normal(prior_beta_sv[1], prior_beta_sv[2]);
  sigma_sv_intercept  ~ gamma(prior_sigma_sv_intercept[1], prior_sigma_sv_intercept[2]);
  z_sv ~ normal(0, 1);

  target += reduce_sum(partial_sum, idx, grainsize,
                        rt, acc, pid, condition, resp_type,
                        alpha_intercept, beta_alpha2, beta_alpha3,
                        beta_nu[1], nu_c, nu_r, nu_cr,
                        t0, sv);
}
