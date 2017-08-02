
# FIRST RUN SIMULATION to get a value to fit:
# source("gp-sim.R");  # first run simulation 
library(rstan)
stan_dat <- read_rdump('gp-fit.data.R')

fit_fit <- stan(file="gp-fit.stan", data=stan_dat,
                 iter=200, chains=3);

print(fit_fit, pars = c('rho','alpha','sigma'))

fit_fit_lat <- stan(file="gp-fit-latent.stan", data=stan_dat,
                 iter=200, chains=3);

print(fit_fit_lat, pars = c('rho','alpha','sigma'))
