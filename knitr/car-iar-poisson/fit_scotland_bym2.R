#setwd("~/Documents/example-models/knitr/car-iar-poisson/")
library(rstan)   
options(mc.cores = parallel::detectCores())  

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

scot_stanfit = stan("bym2_predictor_plus_offset.stan",
         data=list(N,N_edges,node1,node2,y,x,E),
         iter=10000,diagnostic_file = "diag.txt");

print(scot_stanfit,
      pars=c("lp__", "beta0", "beta1", "mixing_parameter", "sigma_random", "mu[5]"),
      probs=c(0.025, 0.5, 0.975));

