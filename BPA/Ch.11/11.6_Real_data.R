## 11. Integrated population models
# 11.6. Real data example: hoopoe population dynamics 

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Read data
## The data generation code is in bpa-code.txt, available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
stan_data <- read_rdump("ipm_hoopoe.data.R")

## Parameters monitored
params <- c("phij", "phia", "f", "omega", "p", "lambda",
            "mphij", "mphia", "mfec", "mim", "mlam",
            "sig_phij", "sig_phia", "sig_fec", "sig_im",
            "N1", "NadSurv", "Nadimm", "Ntot")

## MCMC settings
ni <- 12000
nt <- 6
nb <- 6000
nc <- 4

## Initial values
inits <- lapply(1:nc, function(i)
    list(l_mphij = rnorm(1, 0.2, 0.5),
         l_mphia = rnorm(1, 0.2, 0.5),
         l_mfec = rnorm(1, 0.2, 0.5),
         l_mim = rnorm(1, 0.2, 0.5),
         l_p = rnorm(1, 0.2, 1),
         sig_phij = runif(1, 0.1, 10),
         sig_phia = runif(1, 0.1, 10),
         sig_fec = runif(1, 0.1, 10),
         sig_im = runif(1, 0.1, 10),
         N1 = round(runif(stan_data$nyears, 1, 50), 0),
         NadSurv = round(runif(stan_data$nyears, 5, 50), 0),
         Nadimm = round(runif(stan_data$nyears, 1, 50), 0)))

## Call Stan from R
ipm_hoopoe <- stan("ipm_hoopoe.stan",
                   data = stan_data, init = inits, pars = params,
                   chains = nc, iter = ni, warmup = nb, thin = nt,
                   seed = 2,
                   control = list(adapt_delta = 0.95),
                   open_progress = FALSE)
## Divergent transitions after warmup may occur.
print(ipm_hoopoe, digits = 3)
