functions {
  // mean decision time at the upper boundary for start z 
  // coth(x) = 1 / tanh(x)
  real mdt_upper(real v, real a, real z) {
    return (a / tanh(a * v) - z / tanh(z * v)) / v;
  }
  // decision-time variance at the upper boundary for start z 
  // csch(x) = 1 / sinh(x)
  real vdt_upper(real v, real a, real z) {
    return (square(a / sinh(a * v)) - square(z / sinh(z * v))) / square(v)
           + (a / tanh(a * v) - z / tanh(z * v)) / pow(v, 3);
  }
}

data {
  int<lower=1> I;                       // persons
  int<lower=1> K;                       // cells (= 2)
  vector[K] cond;                       // condition code per cell (+1 / -1)

  array[I, K] int<lower=2> J;           // trials per cell
  array[I, K] int<lower=0> U;           // upper (= correct) responses per cell (was C)
  matrix[I, K] MRT_U;                   // mean RT of upper responses (0 if U = 0)   
  matrix[I, K] MRT_L;                   // mean RT of lower responses (0 if L = 0)   
  matrix[I, K] VRT;                     // variance of ALL RTs in the cell          

  // priors: [1] = mean, [2] = SD
  vector[2] prior_beta_nu_intercept;    // normal
  vector[2] prior_beta_nu_c;            // normal
  vector[2] prior_beta_alpha;           // normal, log scale
  vector[2] prior_beta_tau;             // normal, seconds
  vector[2] prior_beta_w;               // normal, logit scale     

  vector[2] prior_sigma_nu_intercept;   // normal (half-normal via lower=0)
  vector[2] prior_sigma_nu_c;           // normal
  vector[2] prior_sigma_alpha_intercept;// normal (log scale)
  vector[2] prior_sigma_tau_intercept;  // normal
  vector[2] prior_sigma_w;              // normal (logit scale)    
}

transformed data {
  array[I, K] int L;                    // lower (= error) responses per cell
  for (i in 1:I) for (k in 1:K) L[i, k] = J[i, k] - U[i, k];
}

parameters {
  // population means
  real beta_nu_intercept;
  real beta_nu_c;
  real beta_alpha;                              // log boundary sep intercept
  real<lower=0.05, upper=0.45> beta_tau;
  real beta_w;                                  // logit of the relative starting point

  // between-person SDs
  real<lower=0> sigma_nu_intercept;
  real<lower=0> sigma_nu_c;
  real<lower=0> sigma_alpha_intercept;
  real<lower=0> sigma_tau_intercept;
  real<lower=0> sigma_w;                        // (logit scale)

  // NCP
  vector[I] z_nu_intercept;
  vector[I] z_nu_c;
  vector[I] z_alpha_intercept;
  vector[I] z_tau_intercept;
  vector[I] z_w;                                
}

transformed parameters {
  vector[I] nu_intercept    = beta_nu_intercept + sigma_nu_intercept    * z_nu_intercept;
  vector[I] nu_c            = beta_nu_c         + sigma_nu_c            * z_nu_c;
  vector[I] alpha_intercept = beta_alpha        + sigma_alpha_intercept * z_alpha_intercept; // log scale
  vector[I] tau_intercept   = beta_tau          + sigma_tau_intercept   * z_tau_intercept;
  vector[I] w_logit         = beta_w            + sigma_w               * z_w;               
}

model {
  // population means
  beta_nu_intercept ~ normal(prior_beta_nu_intercept[1], prior_beta_nu_intercept[2]);
  beta_nu_c         ~ normal(prior_beta_nu_c[1],         prior_beta_nu_c[2]);
  beta_alpha        ~ normal(prior_beta_alpha[1],        prior_beta_alpha[2]);
  beta_tau          ~ normal(prior_beta_tau[1],          prior_beta_tau[2]);
  beta_w            ~ normal(prior_beta_w[1],            prior_beta_w[2]);         

  // between-person SDs
  sigma_nu_intercept    ~ normal(prior_sigma_nu_intercept[1],    prior_sigma_nu_intercept[2]);
  sigma_nu_c            ~ normal(prior_sigma_nu_c[1],            prior_sigma_nu_c[2]);
  sigma_alpha_intercept ~ normal(prior_sigma_alpha_intercept[1], prior_sigma_alpha_intercept[2]);
  sigma_tau_intercept   ~ normal(prior_sigma_tau_intercept[1],   prior_sigma_tau_intercept[2]);
  sigma_w               ~ normal(prior_sigma_w[1],               prior_sigma_w[2]);  

  // NCP
  z_nu_intercept    ~ std_normal();
  z_nu_c            ~ std_normal();
  z_alpha_intercept ~ std_normal();
  z_tau_intercept   ~ std_normal();
  z_w               ~ std_normal();                                               

  // EZ-freeZ likelihood per person and cell
  for (i in 1:I) {
    real a   = exp(alpha_intercept[i]);
    real tau = tau_intercept[i];
    real w   = inv_logit(w_logit[i]);
    real z   = a * w;          // absolute starting point
    real zl  = a - z;          // distance from the start to the upper boundary (used for the lower boundary)

    for (k in 1:K) {
      real v = nu_intercept[i] + nu_c[i] * cond[k];
      // the equations are undefined at exactly v = 0 --> keep v away from 0 (sign kept, so negative drift is fine)
      if (abs(v) < 1e-3) v = (v < 0) ? -1e-3 : 1e-3;

      real p_u = expm1(-2 * v * z) / expm1(-2 * v * a);             
      real M_u = tau + mdt_upper(v, a, z);                           
      real M_l = tau + mdt_upper(v, a, zl);                          
      real V_u = vdt_upper(v, a, z);                                 
      real V_l = vdt_upper(v, a, zl);                               
      real V   = p_u * V_u + (1 - p_u) * V_l
                 + p_u * (1 - p_u) * square(M_u - M_l);              

      p_u = fmin(fmax(p_u, 1e-9), 1 - 1e-9);                         // numerical safety 

      U[i, k] ~ binomial(J[i, k], p_u);                              
      if (U[i, k] > 0)
        MRT_U[i, k] ~ normal(M_u, sqrt(fmax(V_u, 1e-8) / U[i, k])); 
      if (L[i, k] > 0)
        MRT_L[i, k] ~ normal(M_l, sqrt(fmax(V_l, 1e-8) / L[i, k])); 
      VRT[i, k] ~ normal(V, sqrt(2 * square(V) / (J[i, k] - 1)));    
    }
  }
}

generated quantities {
  vector[I] alpha_person = exp(alpha_intercept);
  real alpha_pop_median  = exp(beta_alpha);
  vector[I] w_person     = inv_logit(w_logit);    // relative starting point per person (0-1)
  real w_pop_median      = inv_logit(beta_w);     
}
