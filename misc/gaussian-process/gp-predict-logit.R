library(rstan)

stan_dat <- read_rdump('gp-predict-logit.data.R')
fit_logit_predict <- stan(file="gp-predict-logit.stan",   
                          data=stan_dat,
                          iter=200, chains=3);
print(fit_logit_predict, c('rho','alpha'));