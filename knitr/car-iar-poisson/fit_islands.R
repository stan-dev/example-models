## template of an R script to fit BYM2 model
## for spatial data where there are disconnected components
## and/or islands 

library(rstan);
options(mc.cores = 2);

library(maptools);   
library(spdep);
library(rgdal)

# load R script which contains functions
# for munging nb objects
source("scale_nb_components.R");

# load the data
# need observed counts y
# exposure E
# x - matrix of covariates

# load your geodata - shp files or similar
nyc_all_tracts_shp<-readOGR("datasets/nyct2010", layer="nyct2010");

# use spdep "poly2nb" to get nb object 
nb_nyc = poly2nb(nyc_poptracts_shp, row.names=nyc_poptracts_shp$CT2010full)

## check spatial structure - strictly optional
#coords<-coordinates(nyc_poptracts_shp);
#plot(nb_nyc, coords, pch = ".", col = "blue");
#components <- n.comp.nb(nb_nyc)
#table(components$comp.id)


# calculate scaling factor for all components
scaling_factor = scale_nb_components(nb_nyc);

# identify islands
singletons = which(scaling_factor == -1);
N_singletons = length(singletons);

# entire graph represented as node pairs
nbs=nb2graph(nb_nyc);
N = nbs$N;
node1 = nbs$node1;
node2 = nbs$node2;
N_edges = nbs$N_edges;

# fit model
bym2_nyc_stanfit = stan("models/bym2_islands.stan", data=list(N,N_edges,node1,node2,y,E,scaling_factor,N_singletons,singletons), iter=5000, warmup=4000, chains= 4, control = list(adapt_delta = 0.95), save_warmup=FALSE);
