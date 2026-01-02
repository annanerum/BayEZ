functions {
  
  real ddm_decision_time_mean(real a, real v) {
    if (fabs(v) < 1e-6) {
      return 0.25 * square(a);  // smooth limit at drift ~ 0
    } else {
      real half_z = 0.5 * a * v;
      return (a / (2 * v)) * tanh(half_z);
    }
  }
  
  real ddm_decision_time_var(real a, real v) {
    if (fabs(v) < 1e-6) {
      return pow(a, 4) / 24;  // smooth limit at drift ~ 0
    } else {
      real z = a * v;
      real u = exp(-z);
      real numer = 1 - 2 * z * u - u * u;
      real denom = square(1 + u);
      real var_dec = (a / (2 * pow(v, 3))) * (numer / denom);
      return fmax(var_dec, 1e-9);
    }
  }

  real ddm_mu_rt(real a, real v, real t0) {
    return t0 + ddm_decision_time_mean(a, v);
  }

  real ddm_sigma2_rt(real a, real v, real t0) {
    return ddm_decision_time_var(a, v);
  }
}

data {
  int<lower=1> I;              // number of subjects
  int<lower=1> K;              // number of within-subject conditions (usually 2)
  int<lower=1> J[I, K];        // total trials per subject i, condition k
  int<lower=0> C[I, K];        // correct counts per subject-condition
  real MRT[I, K];              // mean RT of correct responses
  real VRT[I, K];              // variance of RTs of correct responses
  vector[I] x;                 // subject-level predictor 
}

parameters {
  // population-level means for condition 1 
  real<lower=0.1> mu_alpha;   
  real mu_nu;      
  real<lower=0.05> mu_tau;     

  // population-level SDs (between-subject variability) -----
  real<lower=0.1,  upper=0.8>   sigma_alpha;
  real<lower=0>sigma_nu;
  real<lower=0.05, upper=0.6>  sigma_tau;

  // beta_cond_ = shift (cond2 - cond1) for each parameter
  // delta = effect of subject predictor x[i] on drift nu
  real beta_cond_alpha;
  real beta_cond_nu;
  real beta_cond_tau;
  
  real delta;

  // individual-level parameters per subject and condition
  matrix<lower=0.1,  upper=3.0>[I, K]   alpha;   
  matrix [I, K]   nu;      
  matrix<lower=0.05>[I, K]   tau;     
}

model {
  // priors

  // population means for condition 1
  mu_alpha ~ normal(1.5, 0.2);
  mu_nu    ~ normal(0.0, 0.5);
  mu_tau   ~ normal(0.3, 0.0.6);

  // between-subject variability 
  sigma_alpha ~ normal(0.25, 0.1);
  sigma_nu    ~ normal(0.3, 0.1);
  sigma_tau   ~ normal(0.15, 0.1);

  // how much condition 2 shifts alpha, nu, tau on average
  beta_cond_alpha ~ normal(0, 0.20);
  beta_cond_nu    ~ normal(0, 1.00);
  beta_cond_tau   ~ normal(0, 0.05);

  // subject-level predictor effect on drift (same across conditions)
  delta ~ normal(0, 1);


  // likelihood 
  for (i in 1:I) {
    for (k in 1:K) {

      // hierarchical priors 
      alpha[i, k] ~ normal(
        mu_alpha + beta_cond_alpha * (k - 1),
        sigma_alpha
      );

      nu[i, k] ~ normal(
        mu_nu + delta * x[i] + beta_cond_nu * (k - 1),
        sigma_nu
      );

      tau[i, k] ~ normal(
        mu_tau + beta_cond_tau * (k - 1),
        sigma_tau
      );

      // predict summary stats from DDM parameters 
      {
        real a  = alpha[i, k];
        real v  = nu[i, k];
        real t0 = tau[i, k];

        // Predicted accuracy (Pc) 
        real pi      = inv_logit(a * v);
        // avoid log(0) in binomial
        real pi_safe = fmin(fmax(pi, 1e-6), 1.0 - 1e-6);

        // accuracy likelihood
        C[i, k] ~ binomial(J[i, k], pi_safe);

  
        real mu_rt     = ddm_mu_rt(a, v, t0);        // mean correct RT
        real sigma2_rt = ddm_sigma2_rt(a, v, t0);    // var of correct RT

        
        // sampling distribution of the sample mean RT
        if (C[i, k] >= 1) {
          real se_mean = sqrt(sigma2_rt / C[i, k]);
          MRT[i, k] ~ normal(mu_rt, se_mean);
        }

        // sampling distribution of the sample variance of RT
        if (C[i, k] > 1) {
          real se_var = sqrt(2 * square(sigma2_rt) / (C[i, k] - 1));
          VRT[i, k] ~ normal(sigma2_rt, se_var);
        }
      }
    }
  }
}

generated quantities {
  
  // baseline group means (condition 1)
  real mu_alpha_cond1 = mu_alpha;
  real mu_nu_cond1    = mu_nu;
  real mu_tau_cond1   = mu_tau;

  // group means in condition 2 (k=2)
  real mu_alpha_cond2 = mu_alpha + beta_cond_alpha;
  real mu_nu_cond2    = mu_nu    + beta_cond_nu;
  real mu_tau_cond2   = mu_tau   + beta_cond_tau;

  // condition differences (cond2 - cond1)
  real delta_alpha = beta_cond_alpha;
  real delta_nu    = beta_cond_nu;
  real delta_tau   = beta_cond_tau;
}
