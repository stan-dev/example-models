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
source("nb_data_funs.R");

# load data 
# source("my.data.R")
# y = mydata$count;
# E = mydata$expected;

# load your geodata - shp files or similar
# using nyc census tract data version 2010
nyc_tracts_shp<-readOGR("nyct2010", layer="nyct2010");

# use spdep "poly2nb" to get nb object 
nb_nyc = poly2nb(nyc_tracts_shp);

## bym2 model requires symmetric neighborhood structure
validateNb(nb_nyc);

## only need bym2_islands model if neighborhood is disconnected
## or contains nodes with no neighbords (singletons)
isDisconnected(nb_nyc);

## order neighborhood graph
new_order = orderByComponent(nb_nyc);
nb_new = new_order[1][[1]];
rMap = new_order[2][[1]];

## get component sizes for model
comps_new = as.matrix(table(n.comp.nb(nb_new)[[2]]));
nodes_per_component = as.vector(comps_new[,1]);
N_components = length(nodes_per_component);

## check reordering is correct
DEBUG = TRUE;
if (DEBUG) {
  comps_old = as.matrix(table(n.comp.nb(nb_nyc)[[2]]));
  sum(comps_old[,1]) == sum(comps_new[,1])
  nrow(comps_old) == nrow(comps_new);
  coords = coordinates(nyc_poptracts_shp);
  coords_new = reorderMatrix(coords, rMap);

  ## the two plots should look the same
  plot(nb_nyc, coords, pch = ".", col = "blue");
  plot(nb_new, coords_new, pch = ".", col = "red");
}

# calculate scaling factor for all components
scales = scale_nb_components(nb_new);

# identify islands
singletons = which(scales == 0);
N_singletons = length(singletons);

# entire graph represented as node pairs
nbs=nb2graph(nb_new);
N = nbs$N;
node1 = nbs$node1;
node2 = nbs$node2;
N_edges = nbs$N_edges;

# use rMap to reorder the data
# need observed counts y
# exposure E

# we don't really have data
# y = reorderVector(y, rMap);
# E = reorderVector(E, rMap);

# fit model
#bym2_nyc_stanfit = stan("bym2_islands.stan",
#  data=list(N,N_edges,node1,node2,y,E,N_components,N_singletons,nodes_per_component,scales),
#  iter=5000, warmup=4000, chains= 4, control = list(adapt_delta = 0.95), save_warmup=FALSE);
