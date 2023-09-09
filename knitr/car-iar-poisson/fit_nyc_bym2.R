library(maptools);
library(spdep);
library(rgdal)
library(cmdstanr);

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

data = list(N=N,
            N_edges=N_edges,
            node1=node1,
            node2=node2,
            y=y,
            E=E,
            scaling_factor=scaling_factor);

bym2_model = cmdstan_model("bym2_offset_only.stan");
bym2_fit = bym2_model$sample(data=data, parallel_chains=4, refresh=0);

bym2_fit$summary(
             variables = c(
                 "beta0", "rho", "sigma",
                 "mu[1]", "mu[2]", "mu[3]", "mu[500]", "mu[1000]", "mu[1500]", "mu[1900]",
                 "phi[1]", "phi[2]", "phi[3]", "phi[500]", "phi[1000]", "phi[1500]", "phi[1900]",
                 "theta[1]", "theta[2]", "theta[3]", "theta[500]", "theta[1000]", "theta[1500]", "theta[1900]"));

save(bym2_fit, file="nyc_bym2_fit.data.R");
