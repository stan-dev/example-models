library(rstan);
library(R2jags);

inv_logit <- function(u) { return(1 / (1 + exp(-u))); };

irt_stan <- stan_model("irt.stan");

I <- 10;
J <- 20;
alpha <- rnorm(J, 0, 1);
delta <- exp(rnorm(I, 0, 0.5));
beta <- rnorm(I, -1, 2);
y <- matrix(I * J, I, J);
for (i in 1:I) {
  for (j in 1:J) {
    y[i,j] <- rbinom(1,1,inv_logit(delta[i] * (alpha[j] - beta[i])));
  }
}

init_fun = function(chain_id) {
  return(list(alpha=runif(J, -2, 2), beta=runif(I, -2, 2), delta=exp(runif(I, -2, 2))));
}

fit_stan <- sampling(irt_stan, data = c("I", "J", "y"), init=init_fun,
                     chains=4, iter=2000);

start_time <- proc.time();
fit_jags <- jags(data=list(I=I, J=J, y=y), , parameters.to.save=c("alpha","beta","delta"),
                 model.file="irt.jags", n.chains=4, n.iter=10000, init=init_fun);
elapsed_time = proc.time() - start_time;
print(elapsed_time);

fit_jags_mcmc <- as.mcmc(fit_jags);
arr <- as.array(fit_jags_mcmc);
arr <- aperm(arr, c(1,3,2));
monitor(arr, warmup = 0);
