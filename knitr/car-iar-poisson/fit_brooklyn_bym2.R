library(rstan);
options(mc.cores = parallel::detectCores());
library(maptools);   
library(spdep);

load("nyc_subset.data.R");

bklny_tractIDs <- NYC_tractIDs[substr(NYC_tractIDs,1,5)=="36047"];
y = events_2001[all_tractIDs %in% bklyn_tractIDs];
E = pop_2001[all_tractIDs %in% bklyn_tractIDs];

nyc_shp<-readShapePoly("nycTracts10/nycTracts10");
bklyn_tracts <- nyc_shp$GEOID10 %in% bklyn_tractIDs;
bklyn_shp <- nyc_shp[bklyn_tracts,];
bklyn_shp <- bklyn_shp[order(bklyn_shp$GEOID10),];
nb_bklyn = poly2nb(bklyn_shp);

coords<-coordinates(bklyn_shp);
plot(nb_bklyn, coords, pch = ".", col = "blue");

components <- n.comp.nb(nb_bklyn)
table(components$comp.id)


pois_stan = stan_model("pois.stan");
pois_fit = sampling(pois_stan, data=list(N=length(y), y=y, E=E));

print(pois_fit, digits=3, pars=c("beta0", "mu[1]", "mu[2]", "mu[3]", "mu[500]", "mu[600]", "mu[700]"), probs=c(0.025, 0.5, 0.975));

re_stan = stan_model("pois_re.stan");
re_fit = sampling(re_stan, data=list(N=length(y), y=y, E=E));
print(re_fit, digits=3, pars=c("beta0", "mu[1]", "mu[2]", "mu[3]", "mu[500]", "mu[600]", "mu[700]"), probs=c(0.025, 0.5, 0.975));

source("nb_data_funs.R");
nbs=nb2graph(nb_bklyn);
N = nbs$N;
node1 = nbs$node1;
node2 = nbs$node2;
N_edges = nbs$N_edges;

icar_stan = stan_model("pois_icar.stan");
icar_fit = sampling(icar_stan, data=list(N,N_edges,node1,node2,y,E=pop), iter=5000, warmup=4000, chains= 4, control = list(adapt_delta = 0.95), save_warmup=FALSE);

print(icar_fit, digits=3, pars=c("beta0", "mu[1]", "mu[2]", "mu[3]", "mu[500]", "mu[600]", "mu[700]"), probs=c(0.025, 0.5, 0.975));

scaling_factor = scale_nb_components(nb_bklyn)[1];

bym2_fit = stan("bym2_offset_only.stan",
data=list(N,N_edges,node1,node2,y,E,scaling_factor), control = list(adapt_delta = 0.95), warmup=5000, iter=6000);

print(bym2_fit, pars=c("beta0", "rho", "sigma", "mu[1]", "mu[2]", "mu[3]", "mu[500]", "mu[600]", "mu[700]"), probs=c(0.025, 0.5, 0.975));


pois_samples = as.data.frame(pois_fit, pars=c("mu"))
pois_mu = apply(pois_samples, 2, mean)

re_samples = as.data.frame(re_fit, pars=c("mu"))
re_mu = apply(re_samples, 2, mean)

icar_samples = as.data.frame(icar_fit, pars=c("mu"))
icar_mu = apply(icar_samples, 2, mean)

bym2_samples = as.data.frame(bym2_fit, pars=c("mu"))
bym2_mu = apply(bym2_samples, 2, mean)

df_model_fits = data.frame(id=bklyn_tractIDs, y, pois_mu, icar_mu, bym2_mu);

df_bklyn_shp = tidy(bklyn_shp, region="GEOID10");
df_wide = left_join(df_bklyn_shp, df_model_fits);
#df_long = melt(df_wide,id.vars=c("long","lat","order","hole","piece","group","id"), measure.vars=c("y", "pois_mu", "icar_mu", "bym2_mu"), variable.name="model");
df_long = melt(df_wide,id.vars=c("long","lat","order","hole","piece","group","id"), measure.vars=c("icar_mu", "bym2_mu"), variable.name="model");

p1 = ggplot() + geom_polygon(data=df_long, aes(x=long, y=lat, group=group, fill=value)) + facet_wrap(~model) + coord_map() + coord_fixed()
p1 = p1 + scale_fill_gradientn(limits=c(1,20), colors=blues9, oob=scales::squish);
p1
