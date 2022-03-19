
data{
int<lower=1> N;
vector[N] offset;
int<lower=0> deaths[N];
}

parameters{
real log_theta;
real mu;
real<lower=0> sigma_mu;
}


model{
sigma_mu ~ normal(0,1);
mu ~ normal(0,1);
log_theta ~ normal(mu, sigma_mu);

deaths ~ poisson_log(log_theta+log(offset));
}

generated quantities{
  real theta;
  vector[N] death_rep;
  vector[N] log_lik;
  theta=exp(log_theta);
  for (i in 1:N){
  log_lik[i]=poisson_log_lpmf(deaths[i]|log_theta+log(offset[i]));
  death_rep[i]=poisson_log_rng(log_theta+log(offset[i]));
  } 
}

