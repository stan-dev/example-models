setwd("~/github/stan-dev/example-models/knitr/car-iar-poisson")
library(rstan)   

source("mungeCARdata4stan.R")  
source("carlin_data.R")

nbs = mungeCARdata4stan(data$adj, data$num);

y = data$O;
x = data$aff;
E = data$E;
N = nbs$N;
node1 = nbs$node1;
node2 = nbs$node2;
N_edges = nbs$N_edges;

bym_stanfit = stan("bym_carlin_scotland.stan",
        data=list(N,N_edges,node1,node2,y,x,E));

#options(max.print = 5000);
#print(bym_stanfit, digits=4);

print(bym_stanfit, digits=4,
pars=c("lp__", "beta0", "beta1", "tau_theta", "tau_phi", "sigma_theta", "sigma_phi", "psi", "eta[1]", "eta[2]"))
