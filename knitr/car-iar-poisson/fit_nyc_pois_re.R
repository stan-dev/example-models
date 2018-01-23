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
y = events_2001
E = pop_2001;
## set pop > 0 so we can use log(pop) as offset
E[E < 10] = 10;

pois_re_stan = stan_model("pois_re.stan");
pois_re_fit = sampling(pois_re_stan, data=list(N,y,E), chains=3, warmup=4000, iter=5000, save_warmup=FALSE);
print(pois_re_fit, digits=3, pars=c("beta0", "mu[1]", "mu[2]", "mu[3]", "mu[500]", "mu[1000]", "mu[1500]", "mu[1900]", "theta[1]", "theta[2]", "theta[3]", "theta[500]", "theta[1000]", "theta[1500]", "theta[1900]"), probs=c(0.025, 0.5, 0.975));

save(pois_re_fit, file="nyc_pois_re_fit.data.R");
