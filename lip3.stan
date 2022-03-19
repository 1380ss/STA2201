
data{
int<lower=1> N;
vector[N] offset;
int<lower=0> deaths[N];
}

parameters{
vector[N] alpha;
real beta;
real mu;
real<lower=0> sigma_mu;
}

transformed parameters {
vector[N] log_theta;
for (i in 1:N){
 log_theta[i] = alpha[i] + beta;
 }
}

model{


alpha ~ normal(mu, sigma_mu);
mu ~ normal(0,1);
beta ~ normal(0,1);
sigma_mu ~ normal(0,1);

deaths ~ poisson_log(log_theta+log(offset));
}

generated quantities{
  real theta[N];
  vector[N] death_rep;
  vector[N] log_lik;
  for (i in 1:N){
  theta[i]=exp(log_theta[i]);
  log_lik[i]=poisson_log_lpmf(deaths[i]|log_theta[i]+log(offset[i]));
  death_rep[i]=poisson_log_rng(log_theta[i]+log(offset[i]));
  } 
}
