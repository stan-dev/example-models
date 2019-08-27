library(rstan);
library(MCMCpack);

rcat <- function(n,theta)
          sample(length(theta), n, replace = TRUE, prob = theta);
K <- 3;
I <- 500;
J <- 5;
y <- array(0, c(I, J));
alpha <- rep(3, K);       # modest regularization
beta <- matrix(1, K, K);
for (k in 1:K)
  beta[k, k] <- 2.5 * K;  # 75%-ish accuracy
pi <- rdirichlet(1, alpha)
theta <- array(0, c(J, K, K))
for (j in 1:J)
  for (k in 1:K)
    theta[j, k, ] <- rdirichlet(1, beta[k, ]);
z <- rcat(I, pi);
for (i in 1:I)
  for (j in 1:J)
    y[i ,j] <- rcat(1,theta[j, z[i], ]);

theta_init <- array(.2 / (K - 1), c(J, K, K))
for (j in 1:J)
  for (k in 1:K)
      theta_init[j, k, k] <- 0.8
pi_init <- rep(1 / K, K)

fit <- stan("dawid-skene.stan",
            data = c("K", "I", "J", "y", "alpha", "beta"),
            init = function(n) list(theta = theta_init, pi = pi_init),
            chains = 4, iter = 2000);

fit_ss <- extract(fit);
for (j in 1:J)
  for (k_ref in 1:K)
    for (k_resp in 1:K)
      print(sprintf("(%d, %d, %d) theta=%5.3f hat-theta=%5.3f",
                  j, k_ref, k_resp, theta[j, k_ref, k_resp],
                  mean(fit_ss$theta[, j, k_ref, k_resp])), quote=FALSE);
