library(rstan)

### Data

source("ARM/Ch.3/kidiq.data.R", echo = TRUE)

### Model (kidiq_multi_preds.stan): kid_score ~ mom_hs + mom_iq

data.list <- c("N", "kid_score", "mom_hs", "mom_iq")
kidiq_multi_preds <- stan(file='ARM/Ch.3/kidiq_multi_preds.stan', data=data.list,
                          iter=500, chains=4)
kidiq_multi_preds
pairs(kidiq_multi_preds)
