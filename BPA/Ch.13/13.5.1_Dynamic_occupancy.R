## 13. Metapopulation modeling of species distributions using
##     hierarchical logistic regression: Site-occupancy models
## 13.5. Dynamic (multi-season) site-occupancy models
## 13.5.1. Generation and analysis of simulated data

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Read data
## The data generation code is in bpa-code.txt, available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
stan_data <- read_rdump("Dynocc.data.R")

## True values for this data
## psi: 0.6000000 0.5251555 0.5420522 0.6103358 0.5279381
##      0.5787685 0.5348172 0.5196947 0.5649945 0.5470323
## phi: 0.7785394 0.8627761 0.7752497 0.7761624 0.7529837
##      0.7878669 0.7915214 0.8081415 0.8331360
## gamma: 0.1450796 0.1873469 0.4151345 0.1391412 0.3839322
##        0.1871292 0.2071774 0.3019072 0.1754348
## p: 0.1554887 0.7542202 0.8540974 0.3155055 0.2354785
##    0.1271165 0.2430280 0.6133323 0.1183022 0.1066599

## Parameters monitored
params <- c("psi", "phi", "gamma", "p", "n_occ", "growthr", "turnover")

## MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 4

## Initial values
inits <- lapply(1:nc, function(i)
    list(psi1 = runif(1, 0, 1),
         p = runif(stan_data$nyear, 0, 1)))

## Call Stan from R
out <- stan("Dynocc.stan",
            data = stan_data,
            init = inits, pars = params,
            chains = nc, iter = ni, warmup = nb, thin = nt,
            seed = 1,
            control = list(adapt_delta = 0.8),
            open_progress = FALSE)
print(out, digits = 2)
