library(rstan)
stan_dat <- read_rdump('gp-fit.data.R')

fit_gp <- stan(file="gp-fit.stan", data=stan_dat,
               iter=200, chains=3);

print(fit_gp, pars = c('rho','alpha','sigma'))

fit_lat_gp <- stan(file="gp-fit-latent.stan", data=stan_dat,
                   iter=200, chains=3);

print(fit_lat_gp, pars = c('rho','alpha','sigma'))
