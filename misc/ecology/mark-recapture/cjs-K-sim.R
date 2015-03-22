K <- 10;                      # captures
p <- runif(K, 0.7, 0.8);      # capture at k
phi <- runif(K-1, 0.8, 0.9);  # survive k to k+1

I <- 1000;                    # individual
X <- matrix(0,I,K);           # capture of individual
for (i in 1:I) {
  alive <- 1;
  for (k in 1:K) {
    if (alive)
      X[i,k] <- rbinom(1,1,p[k]);
    if (alive && (k < K))
      alive <- rbinom(1,1,phi[k]);
  }
}

chi <- rep(0,K);              # uncaptured after given alive
chi[K] <- 1;
k <- K - 1;
while (k > 0) {
  chi[k] <- (1 - phi[k]) + phi[k] * (1 - p[k+1]) * chi[k+1];
  k <- k - 1;
}

beta <- phi[K-1] * p[K];      # not individually identifiable

# to fit:
# library('rstan');
# fit <- stan('cjs-K.stan',data=c("K","I","X"), chains=4,iter=1000)

pop <- rep(0,K);
for (k in 1:K)
  pop[k] <- sum(X[,k]) / p[k];


