## template of an R script to fit BYM2 model
## for spatial data where there are disconnected components
## and/or islands 

library(rstan);
options(mc.cores = 2);

library(maptools);   
library(spdep);
library(rgdal)

# using nb_object helper functions
source("nb_data_funs.R");

# load data - need to restrict to populated tracts
source("nyc_data/nyckidped.data.R")
nyc_poptractIDs<-nyckidped_tract$ct2010full
y = nyckidped_tract$count;
E = nyckidped_tract$expected;

# if co-variates, set up K and x
#K = 0;
#x = matrix(nrow=0,ncol=0);

# load the fully connected nb object
# based on shapefiles in "nyct2010"
source("nyc_data/nb_nyc_mm.R");

## bym2 model requires symmetric neighborhood structure
validateNb(nb_nyc_mm);

## only need bym2_islands model if neighborhood is disconnected
## or contains nodes with no neighbords (singletons)
!(isDisconnected(nb_nyc_mm))

# calculate scaling factor for all components
scaling_factor = scale_nb_components(nb_nyc_mm)[1];

# entire graph represented as node pairs
nbs=nb2graph(nb_nyc_mm);
N = nbs$N;
node1 = nbs$node1;
node2 = nbs$node2;
N_edges = nbs$N_edges;

# fit model
bym2_nyc_stanfit = stan("bym2.stan",
  data=list(N,N_edges,node1,node2,y,E,scaling_factor),
  iter=5000, warmup=4000, chains= 4, control = list(adapt_delta = 0.95), save_warmup=FALSE);
