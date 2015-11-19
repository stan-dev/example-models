# BUGS rats example (Vol 1, Example 1)
# http://www.mrc-bsu.cam.ac.uk/bugs/winbugs/Vol1.pdf
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores(logical = FALSE))
devAskNewPage(ask = FALSE)
chains = 4
iter = 1000

sourceToList = function(file){
  source(file, local = TRUE)
  d = mget(ls())
  d$file = NULL
  d
}
# Data are the same for all models
data = sourceToList("rats.data.R")
init = rep(list(sourceToList("rats.init.R")), chains)

# Indexed data for use with rats_stanified.stan
dataIndexed = with(data, list(
  N = N,
  Npts = length(y),
  rat = rep(1:nrow(y), ncol(y)),
  x = rep(x, each = nrow(y)),
  y = as.numeric(y),
  xbar = xbar
))


# With BUGS init
rats = stan(file = "rats.stan", data = data, chains = chains,
            init = init, iter = iter)
plot(rats, ask = FALSE)

# With random init
rats = stan(file = "rats.stan", data = data, chains = chains,  iter = iter)
plot(rats, ask = FALSE)
# traceplot(rats, inc_warmup = FALSE)

# rats_vec.stan
rats_vec = stan(file = "rats_vec.stan", data = data, chains = chains,
                init = init, iter = iter)
plot(rats_vec, ask = FALSE)

# rats_vec.stan_unit

# TODO: This gives an increment_log_prob() warning.
rats_vec_unit = stan(file = "rats_vec_unit.stan", data = data, chains = chains,
                init = init, iter = iter)
plot(rats_vec_unit, ask = FALSE)


# Indexed and simplified rats example, using random initialization
rats_stanified = stan(file = "rats_stanified.stan", data = dataIndexed, chains = chains,
                      iter = iter)
plot(rats_stanified, ask = FALSE)

# Compare with BUGS and Gelfand
c(
  rats_mu_beta = mean(extract(rats, "mu_beta")[[1]]),
  rats_vec_mu_beta = mean(extract(rats_vec, "mu_beta")[[1]]),
  rats_vec_unit_mu_beta = mean(extract(rats_vec_unit, "mu_beta")[[1]])*100,
  rats_stanified_mu_beta = mean(extract(rats_stanified, "mu_beta")[[1]]),
  BUGS_mu_beta = 6.183,
  Gelfand_mu_beta =  6.19)

c(
  rats_sigma_beta = mean(extract(rats, "sigma_beta")[[1]]),
  rats_vec_sigma_beta = mean(extract(rats_vec, "sigma_beta")[[1]]),
  rats_vec_unit_sigma_beta = mean(extract(rats_vec_unit, "sigma_beta")[[1]]),
  rats_stanified_sigma_beta = mean(extract(rats_stanified, "sigma_beta")[[1]])
)

