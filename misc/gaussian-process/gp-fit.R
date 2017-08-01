
# FIRST RUN SIMULATION to get a value to fit:
# source("gp-sim.R");  # first run simulation 
library(rstan)
stan_dat <- read_rdump('gp-fit.data.R')

fit_fit <- stan(file="gp-fit.stan", data=stan_dat,
                 iter=200, chains=3);

fit_fit_ss <- extract(fit_fit, permuted=TRUE);

traceplot(fit_fit);

# histograms of fits nice here
