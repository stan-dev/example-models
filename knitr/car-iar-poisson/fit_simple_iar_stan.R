library(rstan)   
options(mc.cores = parallel::detectCores())  

source("mungeCARdata4stan.R")  
source("scotland_data.R")

nbs = mungeCARdata4stan(data$adj, data$num);
N = data$N;
node1 = nbs$node1;
node2 = nbs$node2;
N_edges = nbs$N_edges;
y = data$y
x = data$x
E = data$E

#icar_stan = stan_model("icar_prior.stan");
#fit_icar_stan = sampling(icar_stan, data=list(N,N_edges,node1,node2), chains=3, warmup=4000, iter=5000, save_warmup=FALSE);

#ofile = "icar_prior_stan_5K_iters.txt";
#capture.output(print(fit_stan, digits=3, probs=c(0.025, 0.975)),file=ofile);

pois_icar_stan = stan_model("pois_icar.stan");
fit_pois_icar = sampling(pois_icar_stan, data=list(N,N_edges,node1,node2,y,x,E), chains=3, warmup=4000, iter=5000, save_warmup=FALSE);


ofile = "icar_pois_stan_5K_iters.txt";
capture.output(print(fit_pois_icar, digits=3, probs=c(0.025, 0.975)),file=ofile);
