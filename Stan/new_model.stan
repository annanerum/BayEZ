data {
  int<lower=1> I;               // number of individuals
  int<lower=1> J[I];            // number of trials per individual
  int<lower=0> C[I];            // correct responses
  vector[I] MRT;                // mean RTs
  vector[I] VRT;                // variance of RTs
  vector[I] x;                  // predictor
}

parameters {
  // group-level means
  real<lower=0.1> mu_alpha; 
  real mu_nu; 
  real<lower=0.05> mu_tau;

  // group-level SDs
  real<lower=0.1, upper=0.8> sigma_alpha; 
  real<lower=0>sigma_nu; 
  real<lower=0.05, upper=0.6> sigma_tau; 

  // regression weight on x
  real<lower=0, upper=1> delta;

  // individual-level parameters
  vector<lower=0.1, upper=3.0>[I] alpha;
  vector [I] nu; 
  vector<lower=0.05>[I] tau;
}

model {
  // Priors
  mu_alpha ~ normal(1.5, 0.2);
  mu_nu    ~ normal(0.0, 0.5);
  mu_tau   ~ normal(0.3, 0.06);

  sigma_alpha ~ normal(0.25, 0.1);   
  sigma_nu    ~ normal(0.3, 0.1);   
  sigma_tau   ~ normal(0.15, 0.1); 

  delta ~ normal(0, 1); // changed 
  

  // hierarchical individual-level parameters
  for (i in 1:I) {
    alpha[i] ~ normal(mu_alpha, sigma_alpha);
    nu[i]    ~ normal(mu_nu + delta * x[i], sigma_nu);
    tau[i]   ~ normal(mu_tau, sigma_tau);
  }

  // likelihood
  for (i in 1:I) {
    real u = exp(-alpha[i] * nu[i]);
    real pi = 1 / (1 + u); // predicted accuracy
    real mu_rt = tau[i] + (alpha[i] / (2 * nu[i])) * ((1 - u) / (1 + u));
    real sigma2_rt = (alpha[i] / (2 * pow(nu[i], 3))) * ((1 - 2 * alpha[i] * nu[i] * u - u^2) / square(1 + u));

    C[i] ~ binomial(J[i], pi);
    MRT[i] ~ normal(mu_rt, sqrt(sigma2_rt / C[i]));
    VRT[i] ~ normal(sigma2_rt, sqrt(2 * square(sigma2_rt) / (C[i] - 1)));
  }
}
