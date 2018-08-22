inv_logit <- function(u) 1 / (1 + exp(-u))

N <- 500
R <- 10
D <- 20

x <- matrix(rnorm(N * D), N, D)
w <- rnorm(D)
w0 <- -2

alpha <- runif(R, 0.65, 0.95)
beta <- runif(R, 0.65, 0.95)

z <- rep(0, N);
for (n in 1:N)
  z[n] <- rbinom(1, 1, inv_logit(w0 + (x[n, ] %*% w)))

y <- matrix(0, N, R);
for (n in 1:N)
  for (r in 1:R)
    y[n, r] <- rbinom(1, 1, ifelse(z[n], alpha[r], 1 - beta[r]))
