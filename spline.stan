data{
  int<lower=0> nb;
  int<lower=0> ns;
  int<lower=0> ny;
  
  matrix[ny, ns] Y;
  matrix[ny, nb] B;
}

parameters{
  matrix[nb, ns] alpha;
  real<lower=0> sigma[ns];
  real log_sigma_alpha[ns];
  real mu_sigma;
  real<lower=0> tau;
}

transformed parameters{
  matrix[ny, ns] mu;
  mu = B*alpha;
}

model{
  //likelihood
  for(i in 1:ns) {
    for (j in 1:ny) {
      Y[j, i] ~ normal(mu[j, i], sigma[i]);
    }
  }
  
  //priors
  for (i in 1:ns) {
    alpha[1, i] ~ normal(0, exp(log_sigma_alpha[i]));
    alpha[2, i] ~ normal(0, exp(log_sigma_alpha[i]));
    alpha[3:nb, i] ~ normal(2*alpha[2:(nb-1), i] - alpha[1:(nb-2), i], exp(log_sigma_alpha[i]));
  }
  
  log_sigma_alpha ~ normal(mu_sigma, tau);
  mu_sigma ~ normal(0, 1);
  tau ~ normal(0, 1);
  sigma ~ normal(0, 1);
}

generated quantities{
  matrix[ny, ns] Y_rep;
  
  for (i in 1:ns) {
    for(j in 1:ny) {
      Y_rep[j, i]= normal_rng(mu[j, i], sigma[i]);
    }
  }
}


