# BUGS pump example (Vol 1, Example 2)
# http://www.mrc-bsu.cam.ac.uk/bugs/winbugs/Vol1.pdf
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores(logical = FALSE))

chains = 4
iter = 1500

sourceToList = function(file){
  source(file, local = TRUE)
  d = mget(ls())
  d$file = NULL
  d
}
# Data are the same for all models
data = sourceToList("pump.data.R")
init = rep(list(sourceToList("pump.init.R")), chains)

# Original init
pump = stan(file = "pump.stan", data = data, chains = chains,
            warmup = 500, init = init, iter = iter, seed = 4711)
plot(pump, ask = FALSE)
pump

# Random (default) init
pump_random = stan(file = "pump.stan", data = data, chains = chains,
                   warmup = 500, iter = iter, seed = 4711)
plot(pump_random, ask = FALSE)
pump_random

# Get parameters as listed in manual
pars = c("alpha","beta","theta[1]","theta[2]","theta[9]", "theta[10]")
bugs_results = as.matrix(data.frame(
  bugs_mean = c(0.73, 0.98, 0.06, 0.10, 1.58, 1.97),
  bugs_2_5 = c(0.28, 0.24, 0.02, 0.01, 0.47, 1.24),
  bugs_97_5 = c(1.38, 2.36, 0.12, 0.30, 3.39, 2.93)
))

pump_random_summary =  summary(pump_random, probs = c(.025,0.975),
  pars = pars)$summary[,c(1,3,4,5,6)]

cbind(pump_random_summary, bugs_results)[,c(1,6,3,7,4,8,5)]

#             mean bugs_mean    2.5% bugs_2_5 97.5% bugs_97_5 n_eff
# alpha     0.6891      0.73 0.28885     0.28 1.314      1.38  2372
# beta      0.9250      0.98 0.18769     0.24 2.295      2.36  2316
# theta[1]  0.0596      0.06 0.02045     0.02 0.117      0.12  3537
# theta[2]  0.1012      0.10 0.00754     0.01 0.316      0.30  4000
# theta[9]  1.5872      1.58 0.50311     0.47 3.487      3.39  4000
# theta[10] 1.9848      1.97 1.22120     1.24 2.883      2.93  4000
