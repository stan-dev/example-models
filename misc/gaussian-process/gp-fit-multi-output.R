library(rstan)

stan_dat <- read_rdump('gp-fit-multi-output.data.R')
fit_mult_out <- stan(file="gp-fit-multi-output.stan",   
                     data=stan_dat,
                     iter=200, chains=3,
                     control = list(max_treedepth = 15))
print(fit_mult_out, c('rho','alpha','sigma','Omega'))

print(apply(rstan::extract(fit_mult_out)$Omega,c(2,3), mean))
print(stan_dat$Omega)
