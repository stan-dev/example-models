library(devtools)
if(!require(cmdstanr)){
  devtools::install_github("stan-dev/cmdstanr", dependencies=c("Depends", "Imports"))
}
library(cmdstanr)   
options(digits=3)

source("mungeCARdata4stan.R")  
source("scotland_data.R")
y = data$y;
x = 0.1 * data$x;
E = data$E;

nbs = mungeCARdata4stan(data$adj, data$num);
N = nbs$N;
node1 = nbs$node1;
node2 = nbs$node2;
N_edges = nbs$N_edges;

data = list(N=N,
            N_edges=N_edges,
            node1=node1,
            node2=node2,
            y=y,
            x=x,
            E=E);

bym_model = cmdstan_model("bym_predictor_plus_offset.stan");

bym_scot_stanfit = bym_model$sample(
         data = data,
         parallel_chains = 4,
         refresh=0);
         
bym_scot_stanfit$summary(variables = c("lp__", "beta0", "beta1",
                                       "sigma_phi", "tau_phi",
                                       "sigma_theta", "tau_theta",
                                       "mu[5]","phi[5]","theta[5]"));

bym_scot_stanfit$summary(variables = c("lp__", "beta0", "beta1",
                                       "sigma_phi", "tau_phi",
                                       "sigma_theta", "tau_theta",
                                       "mu[5]","phi[5]","theta[5]"),
                         ~quantile(.x, probs = c(0.025, 0.5, 0.975)));




