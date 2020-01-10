library(rstan)

mod <- stan_model("bernoulli-bandits.stan")

# Simulate bandits
K <- 2
theta <- c(0.5, 0.4)

MAX_N <- 200

theta_hat <- matrix(0, MAX_N, 2)
p_best <- matrix(0, MAX_N, 2)

y <- array(0.0, 0)
z <- array(0.0, 0)
prefix <- function(y, n) array(y, dim = n - 1)
for (n in 1:MAX_N) {
  data <- list(K = K, N = n - 1, y = prefix(y, n), z = prefix(z, n))
  fit <- sampling(mod, data = data, chains = 1, refresh=0,
                  init = 0, control = list(stepsize=0.1, adapt_delta=0.99))
  p_best[n, ] <-
    summary(fit, pars="is_best", probs = c())$summary[ , "mean"]
  z[n] <- sample(K, 1, replace = TRUE, p_best[n, ])
  y[n] <- rbinom(1, 1, theta[z[n]])
}

library(ggplot)
ggplot(data.frame(trial = 1:MAX_N, prob_best = p_best[1:MAX_N, 1])) +
  geom_line(aes(trial, prob_best))
