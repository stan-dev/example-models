setwd("~/github/stan-dev/example-models/knitr/car-iar-poisson")
library(rstan)   
options(mc.cores = parallel::detectCores())  

library(maptools);   
library(spdep);

source("nyc_ped_subset.data.R")  
source("nbdata4stan.R")  

nyc_all_tracts.shp<-readShapePoly("nycTracts10/nycTracts10");
bklyn_tracts <- nyc_all_tracts.shp$GEOID10 %in% bklyn_tractIDs;
bklyn_tracts.shp <- nyc_all_tracts.shp[bklyn_tracts,]
bklyn_tracts.shp <- bklyn_tracts.shp[order(bklyn_tracts.shp$GEOID10),]

nb_bk = poly2nb(bklyn_tracts.shp);
nbs=nbdata4stan(nb_bk);

y = events_all_2001[all_tractIDs %in% bklyn_tractIDs];
x = pop_adj_2001[all_tractIDs %in% bklyn_tractIDs];

N = nbs$N;
node1 = nbs$node1
node2 = nbs$node2
N_edges = nbs$N_edges

gc()
start = proc.time()
bk_stanfit = stan("bym_predictor_only.stan",
             data=list(N,N_edges,node1,node2,y,x),
             iter=5000);
end = proc.time()
# proc time reports *seconds*
end - start;

print(bk_stanfit,
      pars=c("lp__", "beta0", "beta1", "tau_theta", "tau_phi", "sigma_theta", "sigma_phi", "psi", "eta[1]", "eta[2]", "eta[3]", "eta[4]", "eta[709]"),
      probs=c(0.025, 0.5, 0.975));

# compare to simple poisson regression
r_glm_poisson = glm(y~x,family="poisson")
print(r_glm_poisson)
