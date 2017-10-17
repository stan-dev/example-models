library(rstan)

stan_dat <- read_rdump('gp-fit-ARD.data.R')
fit_ARD <- stan(file="gp-fit-ARD.stan",   
                data=stan_dat,
                iter=200, chains=3)
print(fit_ARD, c('rho','alpha','sigma'))
