functions{
  real cont_binomial_lpdf(real x, real N, real p){
    real ans = lchoose(N, x) + x*log(p) + (N-x)*log1m(p);
    return(ans);
  }
}
data {
  int<lower=0> x_positive;
  int<lower=0> n_sample;
  real<lower=0> alpha_detect;
  real<lower=0> beta_detect;
  real<lower=0> alpha_prev;
  real<lower=0> beta_prev;
  real<lower=0> alpha_spec;
  real<lower=0> beta_spec;
  real<lower=0> alpha_sens;
  real<lower=0> beta_sens;
  real<lower=0> pop_size;
  real<lower=0, upper=pop_size> Y_obs;
}
parameters {
  real<lower=0, upper=1> theta;
  real<lower=0, upper=1> spec;
  real<lower=0, upper=1> sens;
  real<lower=0, upper=1> p_detect;
}
transformed parameters{
  real<lower=0, upper=1> p_sample = theta * sens + (1 - theta) * (1 - spec);
  real<lower=Y_obs, upper=pop_size> N = theta*pop_size;
}
model {
  x_positive ~ binomial(n_sample, p_sample);
  theta ~ beta(alpha_prev, beta_prev);
  spec ~ beta(alpha_spec, beta_spec);
  sens ~ beta(alpha_sens, beta_sens);
  p_detect ~ beta(alpha_detect, beta_detect);
  Y_obs ~ cont_binomial(N, p_detect);
}
generated quantities{
  real F = 1/p_detect;
}