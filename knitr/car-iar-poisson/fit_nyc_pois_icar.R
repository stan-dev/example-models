library(maptools);
library(spdep);
library(rgdal)
library(rstan);
options(mc.cores = 3);

load("nyc_subset.data.R");

nyc_shp<-readOGR("nycTracts10", layer="nycTracts10");
geoids <- nyc_shp$GEOID10 %in% nyc_tractIDs;
nyc_subset_shp <- nyc_shp[geoids,];
nyc_subset_shp <- nyc_subset_shp[order(nyc_subset_shp$GEOID10),];
nb_nyc_subset = poly2nb(nyc_subset_shp);

y = events_2001
E = pop_2001;
## set pop > 0 so we can use log(pop) as offset
E[E < 10] = 10;

source("nb_data_funs.R");
nbs=nb2graph(nb_nyc_subset);
N = nbs$N;
node1 = nbs$node1;
node2 = nbs$node2;
N_edges = nbs$N_edges;
scaling_factor = scale_nb_components(nb_nyc_subset)[1];

pois_icar_stan = stan_model("pois_icar.stan");
pois_icar_fit = sampling(pois_icar_stan, data=list(N,N_edges,node1,node2,y,E,scaling_factor), chains=3, save_warmup=FALSE);
print(pois_icar_fit, digits=3, pars=c("beta0", "mu[1]", "mu[2]", "mu[3]", "mu[500]", "mu[1000]", "mu[1500]", "mu[1900]", "phi[1]", "phi[2]", "phi[3]", "phi[500]", "phi[1000]", "phi[1500]", "phi[1900]", "sigma"), probs=c(0.025, 0.5, 0.975));

save(pois_icar_fit, file="nyc_pois_icar_fit.data.R");
