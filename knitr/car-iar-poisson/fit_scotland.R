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

scot_stanfit = stan("bym_predictor_plus_offset.stan",
         data=list(N,N_edges,node1,node2,y,x,E),
         iter=10000);

print(scot_stanfit,
      pars=c("lp__", "beta0", "beta1", "sigma_phi", "tau_phi", "sigma_theta", "tau_theta","mu[5]","phi[5]","theta[5]"),
      probs=c(0.025, 0.5, 0.975),digits=3);

