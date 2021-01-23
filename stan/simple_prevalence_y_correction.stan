data {
  int<lower=0> x_positive;
  int<lower=0> n_sample;
  real<lower=0> alpha_prev;
  real<lower=0> beta_prev;
  real<lower=0> alpha_spec;
  real<lower=0> beta_spec;
  real<lower=0> alpha_sens;
  real<lower=0> beta_sens;
  real<lower=0> pop_size;
}
parameters {
  real<lower=0, upper=1> theta;
  real<lower=0, upper=1> spec;
  real<lower=0, upper=1> sens;
}
transformed parameters{
  real<lower=0, upper=1> p_sample = theta * sens + (1 - theta) * (1 - spec);
}
model {
  x_positive ~ binomial(n_sample, p_sample);
  theta ~ beta(alpha_prev, beta_prev);
  spec ~ beta(alpha_spec, beta_spec);
  sens ~ beta(alpha_sens, beta_sens);
}
generated quantities{
  real N = theta*pop_size;
}
