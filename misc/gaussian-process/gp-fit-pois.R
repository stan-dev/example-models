library(rstan)

stan_dat <- read_rdump('gp-fit-pois.data.R')
fit_pois <- stan(file="gp-fit-pois.stan",   
                 data=stan_dat,
                 iter=200, chains=3)
print(fit_pois, c('rho','alpha','a'))
