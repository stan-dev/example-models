library(rstan);
source("schools.data.R");

fit <- stan("schools-3.stan",
            data=c("N","M","LRT",
                   "school","School_denom","School_gender",
                   "VR","Y","Gender","R"),
            # control=list(adapt_delta=0.9),
            chains=1, iter=400, 
            init=0);
