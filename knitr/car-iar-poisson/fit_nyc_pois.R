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

N = length(events_2001);
y = events_2001;
E = pop_2001;
## set pop > 0 so we can use log(pop) as offset
E[E < 10] = 10;

pois_stan = stan_model("pois.stan");
pois_fit = sampling(pois_stan, data=list(N,y,E), chains=3, save_warmup=FALSE);
print(pois_fit, digits=3, pars=c("beta0", "mu[1]", "mu[2]", "mu[3]", "mu[500]", "mu[1000]", "mu[1500]", "mu[1900]"), probs=c(0.025, 0.5, 0.975));

save(pois_fit, file="nyc_pois_fit.data.R");
