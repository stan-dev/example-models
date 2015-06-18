library(MASS);

K <- 8;
D <- 4;
N <- 500;

# stick breaking construction to generate a random correlation matrix
L_Omega <- matrix(0,D,D);
L_Omega[1,1] <- 1;
for (i in 2:D) {
  bound <- 1;
  for (j in 1:(i-1)) {
    L_Omega[i,j] <- runif(1, -sqrt(bound), sqrt(bound));
    bound <- bound - L_Omega[i,j]^2;
  }
  L_Omega[i,i] <- sqrt(bound);
}

Omega <- L_Omega %*% t(L_Omega);

x <- matrix(rnorm(N * K, 0, 1), N, K);

beta <- matrix(rnorm(D * K, 0, 1), D, K);

z <- matrix(NA, N, D);
for (n in 1:N) {
  z[n,] <- mvrnorm(1, x[n,] %*% t(beta), Omega);
}

y <- matrix(0, N, D);
for (n in 1:N) {
  for (d in 1:D) {
    y[n,d] <- ifelse((z[n,d] > 0), 1, 0);
  }
}

library(rstan);
fit <- stan("probit-multi.stan", 
            data = c("K", "D", "N", "y", "x"),
            control = list(stepsize=0.01, adapt_delta=0.99));
