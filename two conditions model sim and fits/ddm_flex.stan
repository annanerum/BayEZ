functions {
  real partial_sum(array[] int idx_slice, int start, int end,
                   vector rt, array[] int acc, array[] int pid, vector condition,
                   vector a, vector t0, vector w, vector sv, real sw, real st0,
                   vector nu_intercept, vector nu_c, data real prec) {
    real lp = 0;
    for (n in 1:size(idx_slice)) {
      int i = idx_slice[n];
      int p = pid[i];
      real v = nu_intercept[p] + nu_c[p] * condition[i];
      if (acc[i] == 1) {
        lp += wiener_lpdf(rt[i] | a[p], t0[p], w[p], v, sv[p], sw, st0, prec);
      } else {
        // lower boundary: mirror the process (start 1 - w, drift -v)
        // sw and st0 are symmetric
        lp += wiener_lpdf(rt[i] | a[p], t0[p], 1 - w[p], -v, sv[p], sw, st0, prec);
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
  array[N] int<lower=-1, upper=1> acc;      // +1 upper (correct), -1 lower (error)
  vector[N] condition;                      // +1 / -1
  vector<lower=0>[I] min_rt;                // fastest RT per person
  int<lower=1> grainsize;

  // priors: [1] = mean, [2] = SD
  vector[2] prior_beta_nu_intercept;
  vector[2] prior_beta_nu_c;
  vector[2] prior_beta_alpha;               // log scale
  vector[2] prior_beta_tau;                 // seconds
  vector[2] prior_beta_w;                   // logit scale
  vector[2] prior_beta_sv;                  // log scale
  vector[2] prior_sw;                       // normal on the natural scale (half-normal via lower=0)
  vector[2] prior_st0;                      // normal on the natural scale, seconds

  vector[2] prior_sigma_nu_intercept;
  vector[2] prior_sigma_nu_c;
  vector[2] prior_sigma_alpha_intercept;
  vector[2] prior_sigma_tau_intercept;
  vector[2] prior_sigma_w;
  vector[2] prior_sigma_sv;

  real<lower=0, upper=1> sw_upper;   // upper bound for sw 

  // guard limits
  real guard_a_max;      // max boundary separation
  real guard_sv_max;     // max drift variability
  real guard_v_max;      // max |drift| in any condition

  real<lower=0> precision_derivatives;   // 1e-4 = Stan default; 1e-3 = faster

  int<lower=0, upper=1> estimate_sw;     // 0 = sw fixed at 0
  int<lower=0, upper=1> estimate_st0;    // 0 = st0 fixed at 0
}

transformed data {
  array[N] int idx;
  for (n in 1:N) idx[n] = n;
}

parameters {
  real beta_nu_intercept;
  real beta_nu_c;
  real beta_alpha;
  real<lower=0.05, upper=0.45> beta_tau;
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
  vector<lower=0, upper=min_rt>[I] t0;      // lower edge of the NDT, below each fastest RT
  vector[I] z_w;
  vector[I] z_sv;

  array[estimate_sw]  real<lower=0, upper=sw_upper> sw_par;   // empty if not estimated
  array[estimate_st0] real<lower=0, upper=0.5>      st0_par;  // 0.5 s far above plausible
}

transformed parameters {
  vector[I] nu_intercept    = beta_nu_intercept + sigma_nu_intercept    * z_nu_intercept;
  vector[I] nu_c            = beta_nu_c         + sigma_nu_c            * z_nu_c;
  vector[I] alpha_intercept = beta_alpha        + sigma_alpha_intercept * z_alpha_intercept;
  vector[I] w_logit         = beta_w            + sigma_w               * z_w;
  vector[I] log_sv          = beta_sv           + sigma_sv              * z_sv;

  real sw  = 0;
  real st0 = 0;
  if (estimate_sw)  sw  = sw_par[1];
  if (estimate_st0) st0 = st0_par[1];

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
  if (estimate_sw)  sw_par[1]  ~ normal(prior_sw[1],  prior_sw[2]);
  if (estimate_st0) st0_par[1] ~ normal(prior_st0[1], prior_st0[2]);

  sigma_nu_intercept    ~ normal(prior_sigma_nu_intercept[1],    prior_sigma_nu_intercept[2]);
  sigma_nu_c            ~ normal(prior_sigma_nu_c[1],            prior_sigma_nu_c[2]);
  sigma_alpha_intercept ~ normal(prior_sigma_alpha_intercept[1], prior_sigma_alpha_intercept[2]);
  sigma_tau_intercept   ~ normal(prior_sigma_tau_intercept[1],   prior_sigma_tau_intercept[2]);
  sigma_w               ~ normal(prior_sigma_w[1],               prior_sigma_w[2]);
  sigma_sv              ~ normal(prior_sigma_sv[1],              prior_sigma_sv[2]);

  z_nu_intercept    ~ std_normal();
  z_nu_c            ~ std_normal();
  z_alpha_intercept ~ std_normal();
  t0 ~ normal(beta_tau, sigma_tau_intercept);   // like in DDM5 no truncation normalisation
  z_w               ~ std_normal();
  z_sv              ~ std_normal();

  // reject implausible proposals and any w for which the start-point
  // range (w - sw/2, w + sw/2) would leave (0, 1), before evaluating the likelihood
  if (max(a) > guard_a_max || max(sv) > guard_sv_max ||
      max(abs(nu_intercept) + abs(nu_c)) > guard_v_max ||
      min(w) <= 0.5 * sw + 0.005 || max(w) >= 1 - 0.5 * sw - 0.005) {
    target += negative_infinity();
  } else {
    target += reduce_sum(partial_sum, idx, grainsize,
                         rt, acc, pid, condition,
                         a, t0, w, sv, sw, st0, nu_intercept, nu_c,
                         precision_derivatives);
  }
}

generated quantities {
  real w_pop_median     = inv_logit(beta_w);
  real sv_pop_median    = exp(beta_sv);
  real alpha_pop_median = exp(beta_alpha);
  vector[I] mean_ndt    = t0 + st0 / 2;     // mean non-decision time per person (this is comparable to EZ tau)
}
