## 5. State-space models
## 5.2. A simple model

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Data generation code is transplanted from bpa-code.txt
n.years <- 25           # Number of years
N1 <- 30                # Initial population size
mean.lambda <- 1.02     # Mean annual population growth rate
sigma2.lambda <- 0.02   # Process (temporal) variation of the growth rate
sigma2.y <- 20          # Variance of the observation error

y <- N <- numeric(n.years)
N[1] <- N1
lambda <- rnorm(n.years - 1, mean.lambda, sqrt(sigma2.lambda))
for (t in 1:(n.years-1)) {
  N[t + 1] <- N[t] * lambda[t]
}

for (t in 1:n.years) {
  y[t] <- rnorm(1, N[t], sqrt(sigma2.y))
}

## Bundle data
stan_data <- list(y = y, T = n.years)

## Parameters monitored
params <- c("lambda", "mean_lambda", "sigma2_obs", "sigma2_proc",
            "N_est")

## MCMC settings
ni <- 10000
nt <- 5
nb <- 5000
nc <- 4

## Initial values
inits <- lapply(1:nc, function(i) {
    list(sigma_proc = runif(1, 0, 5),
         mean_lambda = runif(1, 0.1, 2),
         sigma_obs = runif(1, 0, 10),
         N_est1 = runif(1, 20, 40))})

## Call Stan from R
ssm <- stan("ssm.stan",
            data = stan_data, init = inits, pars = params,
            chains = nc, iter = ni, warmup = nb, thin = nt,
            seed = 1,
            control = list(adapt_delta = 0.999),
            open_progress = FALSE)

## Note: there may be divergent transitions after warmup.

## Summarize posteriors
print(ssm)
