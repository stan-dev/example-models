library("rstan")

N <- 100
K <- 4
x <- array(rnorm(N*K), c(N, K))
alpha <- -3
beta <- 1:K
sigma <- 5
y <- c(alpha + x %*% beta + rnorm(N, 0, sigma))  # Needed to do c() to convert y from a 1-column matrix into a vector
stan_data <- list(N=N, K=K, x=x, y=y)

model <- stan_model("multiple_regression_qr.stan")
fit <- sampling(model, data=stan_data)
print(fit)

cat(
    "alpha =", alpha, "\n",
    "beta =", beta, "\n",
    "sigma =", sigma, "\n"
)

