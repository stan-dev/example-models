library(rstan);
library(R2jags);

get_n_eff_jags <- function(fit_jags, num_params) {
  return(n_eff_jags[1:num_params, ]);
}


inv_logit <- function(u) { return(1 / (1 + exp(-u))); };

I <- 20;
J <- 200;
alpha <- rnorm(J, 0, 1);
delta <- exp(rnorm(I, 0, 0.5));
beta <- rnorm(I, -1, 2);
num_params <- J + I + I
y <- matrix(I * J, I, J);
for (i in 1:I) {
  for (j in 1:J) {
    y[i,j] <- rbinom(1, 1, inv_logit(delta[i] * (alpha[j] - beta[i])));
  }
}


## this provides an init_fun matching Stan's controlled by an R seed
num_chains <- 4;
init_seed <- 1234;
init_list <- list();
for (chain_id in 1:num_chains) {
  set.seed(chain_id * init_seed);
  list_chain_id <-  list(alpha = runif(J, -2, 2),
                         beta = runif(I, -2, 2),
                         delta = exp(runif(I, -2, 2)));
  init_list <- append(init_list, list(list_chain_id));
}


## STAN
# order: alpha[1:J], beta[1:I], delta[1:I], lp__

start_time_stan <- proc.time();

irt_stan <- stan_model("irt.stan");
fit_stan <- sampling(irt_stan,
                     data = c("I", "J", "y"),
                     init = init_list,
                     chains = num_chains,
                     iter = 2000,
                     seed = init_seed);

elapsed_time_stan = proc.time() - start_time_stan;
print(c("stan", elapsed_time_stan));

n_eff_stan <- summary(fit_stan)$summary[1:num_params, "n_eff"];
sec_stan <- elapsed_time_stan[3];
n_eff_per_sec_stan <- n_eff_stan / sec_stan;

## JAGS
# order: alpha[1:J] , beta[1:I], delta[1:I], deviance
start_time_jags <- proc.time();

fit_jags <- jags(model.file = "irt.jags",
                 data = list(I=I, J=J, y=y),
                 init = init_list,
                 parameters.to.save = c("alpha","beta","delta"),
                 n.chains = num_chains,
                 n.iter = 10000,
                 jags.seed = init_seed);

elapsed_time_jags = proc.time() - start_time_jags;
print(c("jags", elapsed_time_jags));

fit_jags_mcmc <- as.mcmc(fit_jags);
arr <- as.array(fit_jags_mcmc);
arr <- aperm(arr, c(1,3,2));
mon <- monitor(arr, inc_warmup = FALSE);

n_eff_jags <- mon[1:num_params, "n_eff"];
sec_jags <- elapsed_time_jags[3];
n_eff_per_sec_jags <- n_eff_jags / sec_jags;

# translate back to traditional names for display (should fix models)
param_names <- c(rep("theta", J), rep("b", I), rep("a", I));

df_stan <- data.frame(n_eff_per_sec = n_eff_per_sec_stan,
                      parameter = param_names,
                      mcmc = rep("stan", num_params));

df_jags <- data.frame(n_eff_per_sec = n_eff_per_sec_jags,
                      parameter = param_names,
                      mcmc = rep("jags", num_params));

df <- rbind(df_stan, df_jags);

library(ggplot2);
plot <- ggplot(df, aes(x = n_eff_per_sec, fill = mcmc)) +
        geom_histogram() +
        scale_x_log10() +
        facet_wrap(~ parameter) +
        ggtitle("Stan vs. JAGS: IRT 2PL, 20 questions (I), 200 students (J)") +
        xlab("n_eff / sec") +
        ylab("# parameters");
