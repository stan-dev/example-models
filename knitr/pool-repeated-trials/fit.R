# DATA FROM:
# http://www.swarthmore.edu/NatSci/peverso1/Sports%20Data/JamesSteinData/Efron-Morris%20Baseball/EfronMorrisBB.txt
# ORIGINAL PAPER:
# http://www.stat.cmu.edu/~acthomas/724/Efron-Morris.pdf

library(rstan);

df <- read.csv("data.tsv", sep="\t");
N <- dim(df)[1];
K1 <- df[[3]][1];
y1 <- df[[4]];
K2 <- df[[6]];
y <- df[[9]];
fit <- stan("hier.stan", data=c("N", "K1", "y1", "K2"));
fit_ss <- extract(fit)
L <- 0.1;
H <- 0.9;
print(sprintf("%12s  (%3.2f, %3.2f)  %3s", "name", L, H, "y"), 
      quote=FALSE);
for (n in 1:N) {
  quantiles <- quantile(fit_ss$y[,n], prob=c(L, H))
  print(sprintf("%12s  (%3.0f, %3.0f)  %3.0f",
                df[[2]][n], quantiles[1], quantiles[2], y[n]),
        quote=FALSE);
}
fit_logistic <-
  stan("hier-logit.stan", data=c("N", "K1", "y1", "K2"),
       iter=5000,
       control=list(stepsize=0.01, adapt_delta=0.99));


# STEP 1:  Hierarchical model as in BDA ch. 5
# STEP 2:  Logistic version --- centered
# STEP 3:  logistic --- non-centered
# STEP 4:  logistic --- non-centered + control

# CF:   no pooling
#       complete pooling
