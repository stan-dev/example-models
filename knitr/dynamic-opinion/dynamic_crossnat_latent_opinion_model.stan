// bayesian latent model for estimating dynamic cross-national opinion

data{
  int<lower=1> N; // number of national survey opinions
  int<lower=1> J; // number of countries
  int<lower=1> K; // number of items
  int<lower=1> P; // number of items-country combinations  
  int<lower=1> R; // number of national opinion estimates
  array[N] int<lower=1,upper=J> jj; // country j for opinion n
  array[N] int<lower=1,upper=K> kk; // item k for opinion n
  array[N] int<lower=1,upper=P> pp; // item-country p for opinion n
  array[N] int<lower=1,upper=R> rr; // estimate r for opinion n
  array[N] int<lower=1> x; // vector of survey responses, count
  array[N] int<lower=1> samp; // vector of sample sizes
  array[K] int<lower=1> it_len; // number of countries for each item
  array[J] int<lower=1> est_pos; // indicator showing cntry start for estimate vector	
  array[J] int<lower=1> len_theta_ts; // indicator showing length of each country's estimates
  real mn_resp_log; // observed response mean proportion on logit scale
}

parameters{
  real<lower=0> sigma_theta; // opinion evolution error SD
  real<lower=0> sigma_delta; // item-country intercept error SD
  row_vector[J] theta_init; // initial latent traits for year 0
  real<lower=0> phi; // dispersion parameter
  corr_matrix[2] Omega; // correlation matrix for item pars
  vector<lower=0>[2] tau; // cor -> cov conversion
  real mu_lambda; // item intercept expectation
  vector[P] delta_ncp; // redundant parameters for item-country effects
  vector[R] theta_ncp; // redundant parameters for latent opinion
  matrix[K,2] Gamma_ncp; // non-centered parameters for item parameters
}

transformed parameters{
  // dynamic model with non-centered parameters
  vector[R] theta; // R-vector of theta values	
  for (j in 1:J) {                  
    theta[est_pos[j]] = theta_init[j] 
	+ sigma_theta * theta_ncp[est_pos[j]];
	for (i in 1:(len_theta_ts[j]-1)) {
	  theta[(est_pos[j]+i)] = theta[(est_pos[j]+i-1)] 
	  + sigma_theta * theta_ncp[(est_pos[j]+i)];
	}
  }
  // variance-covariance matrix for item ints and slopes
  matrix[2,2] Sigma = quad_form_diag(Omega, tau);  
  // non-centered item-country parameters
  vector[P] delta = sigma_delta * delta_ncp; 
  // item parameter models with non-centered parameters
  real mu_gamm = 1; // item slope expectation
  matrix[K,2] Gamma; // matrix of item intercepts and slopes 
  for (k in 1:K) 
    Gamma[k] = [ mu_lambda , mu_gamm ] + Gamma_ncp[k] * Sigma;
  vector[K] lambda = Gamma[,1]; // K estimated item intercepts
  vector[K] gamm = Gamma[,2]; // K estimated item slopes
  // fitted values model
  vector<lower=0,upper=1>[N]eta = inv_logit(lambda[kk] + gamm[kk] .* theta[rr] + delta[pp]);  
  // reparamaterise beta distr parameters
  vector<lower=0>[N] beta_par1 = phi * eta; 
  vector<lower=0>[N] beta_par2 = phi * (1 - eta); 
}

model{
  // response model
  x ~ beta_binomial(samp, beta_par1, beta_par2); 
  // priors
  phi ~ gamma(3, 0.04); 				
  sigma_theta ~ normal(0, 1); 
  sigma_delta ~ normal(0, 1); 
  tau ~ normal(0, 1);
  Omega ~ lkj_corr(2);
  mu_lambda ~ normal(mn_resp_log, 0.5);
  theta_init ~ normal(0, 1);
  theta_ncp ~ normal(0, 1);
  to_vector(Gamma_ncp) ~ normal(0, 1);
  // standard normal prior for item-country effects, centered within items
  int pos; // local variable indicating which item to evaluate	
  pos = 1;
  for (k in 1:K) { 
    segment(delta_ncp, pos, it_len[k]) ~ normal(0, 1);
    pos = pos + it_len[k];
  }
}

generated quantities {
  vector[N] x_pred; // fitted data to check model
  vector[N] log_lik; // log lik for WAIC calc
  for (i in 1:N) {
    x_pred[i] = beta_binomial_rng(samp[i], beta_par1[i], beta_par2[i]);
    log_lik[i] = beta_binomial_lpmf(x[i] | samp[i], beta_par1[i], beta_par2[i]); 
  }
}
