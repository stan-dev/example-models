## ---- sim-1pl ----
inv_logit <- function(u) {
  return(1 / (1 + exp(-u)));
}

I <- 10;
J <- 20;
theta <- rnorm(J, 0, 1);
b <- rnorm(I, -1, 2);

y <- matrix(I * J, I, J);
for (i in 1:I) {
  for (j in 1:J) {
    y[i,j] <- rbinom(1, 1, inv_logit(theta[j] - b[i]));
  }
}
