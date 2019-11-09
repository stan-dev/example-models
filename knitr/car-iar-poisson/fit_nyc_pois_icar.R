library(maptools);
library(spdep);
library(rgdal)
library(rstan);
options(mc.cores = parallel::detectCores());

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

pois_icar_stan = stan_model("pois_icar.stan");
pois_icar_fit = sampling(pois_icar_stan, data=list(N,N_edges,node1,node2,y,E), chains=3, save_warmup=FALSE);
