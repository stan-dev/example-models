# BUGS seeds example (Vol 1, Example 3)
# http://www.mrc-bsu.cam.ac.uk/bugs/winbugs/Vol1.pdf
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores(logical = FALSE))
chains = 4
iter = 1000
set.seed(4711)

sourceToList = function(file){
  source(file, local = TRUE)
  d = mget(ls())
  d$file = NULL
  d
}
# Data are the same for all models
data = sourceToList("seeds.data.R")
init = rep(list(sourceToList("seeds.init.R")), chains)

# Original init
seeds = stan(file = "seeds.stan", data = data, chains = chains,
             init = init, iter = iter)
plot(seeds, ask = FALSE)
pars = c("alpha0", "alpha1","alpha2", "alpha12","sigma")
seeds_summary = summary(seeds,probs = c(.025,0.975),
              pars = pars)$summary[,c(1,3,4,5,6)]

# Default random init
seeds_random = stan(file = "seeds.stan", data = data, chains = chains,
             iter = iter)
plot(seeds_random, ask = FALSE)
seeds_random_summary =  summary(seeds_random,probs = c(.025,0.975),
          pars = pars)$summary[,c(1,3,4,5,6)]

# Removed tau, using sigma direct
seeds_stanified = stan(file = "seeds_stanified.stan", data = data, chains = chains,
                    iter = iter)
plot(seeds_stanified, ask = FALSE)
seeds_stanified_summary =  summary(seeds_stanified, probs = c(.025,0.975),
                                pars = pars)$summary[,c(1,3,4,5,6)]


(seeds_summary)
(seeds_random_summary)
(seeds_stanified_summary)

# alpha1: seed type
# alpha2: root extract
#
# (seeds_summary)
#            mean    sd    2.5%  97.5% n_eff
# alpha0  -0.5643 0.207 -0.9630 -0.156   470
# alpha1   0.0917 0.335 -0.5626  0.727   435
# alpha2   1.3642 0.285  0.7954  1.957   500
# alpha12 -0.8303 0.450 -1.6577  0.088   552
# sigma    0.2982 0.137  0.0869  0.621   163

#  (seeds_random_summary)
#            mean    sd    2.5%   97.5% n_eff
# alpha0  -0.5511 0.201 -0.9639 -0.1554   665
# alpha1   0.0785 0.320 -0.5598  0.6876   762
# alpha2   1.3459 0.291  0.7932  1.9299   689
# alpha12 -0.8084 0.447 -1.7063  0.0684   767
# sigma    0.3066 0.136  0.0859  0.6213   118

# (seeds_stanified_summary)
#            mean    sd   2.5%   97.5% n_eff
# alpha0  -0.5399 0.230 -1.004 -0.0897   563
# alpha1   0.0597 0.353 -0.683  0.7364   627
# alpha2   1.3545 0.316  0.766  2.0405   618
# alpha12 -0.8432 0.488 -1.833  0.1024   740
# sigma    0.3638 0.140  0.141  0.6681   183

bugs_results = as.matrix(data.frame(
        bugs_mean = c(-0.542, 0.028, 1.37, -0.79, 0.29),
        bugs_se   = c(0.18, .34,.25, .43, .15)))

cbind(seeds_random_summary[,1:2], bugs_results)
#            mean    sd bugs_mean bugs_se
# alpha0  -0.5511 0.201    -0.542    0.18
# alpha1   0.0785 0.320     0.028    0.34
# alpha2   1.3459 0.291     1.370    0.25
# alpha12 -0.8084 0.447    -0.790    0.43
# sigma    0.3066 0.136     0.290    0.15

# 3.1 Constraining random effects to sum to zero
# There is an example in the bugs manual how to constrain parameters to sum
# up to 0.
# This does not improve convergence and n_eff
seeds_centered = stan(file = "seeds_centered.stan", data = data, chains = chains,
                       iter = iter)
plot(seeds_centered, ask = FALSE)
seeds_centered_summary =  summary(seeds_centered, probs = c(.025,0.975),
                                   pars = pars)$summary[,c(1,3,4,5,6)]
(seeds_centered_summary)
par(mfcol = c(3,1))
b_centered = as.numeric(extract(seeds_centered, pars = "b")$b)
hist(b_centered, main = paste("Mean of centered b=", signif(mean(b_centered),2)),
     xlim = c(-1.5,1.5))
b_stan = as.numeric(extract(seeds, pars = "b")$b)
hist(b_stan, main = paste("Mean of default model b=", signif(mean(b_stan),2)),
     xlim =  c(-1.5,1.5))
b_stanified = as.numeric(extract(seeds_stanified, pars = "b")$b)
hist(b_stanified, main = paste("Mean of stanified b=", signif(mean(b_stanified),2)),
     xlim = c(-1.5,1.5))


