library(rstan)

stan_dat <- read_rdump('gp-predict.data.R')

fit_predict <- stan(file="gp-predict.stan",   
                    data=stan_dat,
                    iter=200, chains=3);
print(fit_predict, pars = c('rho','alpha','sigma'))