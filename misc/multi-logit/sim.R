K <- 4;  # categories
D <- 8;  # predictors including intercept
N <- 1024; # observations 

x <- matrix(rnorm(N * D, 0, 1) ,N, D);
x[1:N,1] <- rep(1,N);  # intercepts
beta <- matrix(rnorm((K - 1) * D, 0, 1), K-1, D);
lambda <- rep(0, K);
y <- rep(NA, N);
for (n in 1:N) {
  for (k in 1:(K-1)) {
    lambda[k] <- beta[k,1:D] %*% x[n,1:D];
  }
  y_multi <- rmultinom(1, 1, exp(lambda));  # rmultinom normalizes
  y[n] <- sort((1:K) * y_multi)[K]; # pull out result
}
stan_rdump( c("K","D","N","x","y"), "multi_logit.data.R");
stan_rdump( c("beta"), "multi_logit.parameters.R");
