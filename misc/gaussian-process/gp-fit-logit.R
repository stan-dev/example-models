library(rstan)

stan_dat <- read_rdump('gp-fit-logit.data.R')
fit_logit <- stan(file="gp-fit-logit.stan",   
                  data=stan_dat,
                  iter=200, chains=3)
print(fit_logit, c('rho','alpha','a'))
