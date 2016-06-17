## 5. State-space models
## 5.4. Real example: House martin population counts in the village of
## Magden

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Data generation code is transplanted from original bpa-code.txt
## House martin population data from Magden
pyears <- 6 # Number of future years with predictions
hm <- c(271, 261, 309, 318, 231, 216, 208, 226, 195, 226, 233,
        209, 226, 192, 191, 225, 245, 205, 191, 174)
year <- 1990:2009

## Bundle data
stan_data <- list(y = log(hm), T = length(year), pyears = pyears)

## Parameters monitored
params <- c("r", "mean_r", "sigma2_obs", "sigma2_proc",
            "N_est")

## MCMC settings
ni <- 30000
nt <- 10
nb <- 20000
nc <- 4

## Initial values
inits <- lapply(1:nc, function(i) {
    list(sigma_proc = runif(1, 0, 1),
         mean_r = rnorm(1),
         sigma_obs = runif(1, 0, 1))})

## Call Stan from R
hm_ssm <- stan("ssm2.stan",
               data = stan_data, init = inits, pars = params,
               chains = nc, thin = nt, iter = ni, warmup = nb,
               seed = 1,
               control = list(adapt_delta = 0.999),
               open_progress = FALSE)

## Note: there may be divergent transitions after warmup.

## Summarize posteriors
print(hm_ssm, digits = 3)
