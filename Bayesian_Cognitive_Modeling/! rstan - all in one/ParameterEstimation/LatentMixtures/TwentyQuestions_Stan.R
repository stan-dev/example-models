# clears workspace: 
rm(list=ls()) 

library(rstan)

#### Notes to Stan model #######################################################
## 1) There are two models in this script. The first one is not able to 
##    incorporate missing values as is required for exercise 6.3.2. The second
##    model is able to do that. Actually you can use the second model on both 
##    datasets - with or without NAs. 
## 2) Models with missing values can be tricky. It's important to know that
##    Stan treats variables either as know or unknown, therefore you have to
##    separate both parts first and then make it work somehow. Chapter on
##    Missing Data in Stan manual is useful.
################################################################################
model1 <- "
// Twenty Questions
data { 
  int<lower=0> np;  # rows
  int<lower=0> nq;  # columns
  int k[np,nq];
}
parameters {
  real<lower=0,upper=1> p[np];
  real<lower=0,upper=1> q[nq];
} 
transformed parameters {
  real<lower=0,upper=1> theta[np, nq];

  // Probability Correct Is Product Of Question By Person Rates
  for (i in 1:np)
    for (j in 1:nq)
      theta[i,j] <- p[i] * q[j];
}
model {
  // Priors For People and Questions
  p ~ beta(1, 1);
  q ~ beta(1, 1);
    
  // Correctness Of Each Answer Is Bernoulli Trial
  for (i in 1:np)
    for (j in 1:nq)
      k[i, j] ~ bernoulli(theta[i,j]);
}"
model2 <- "
# Twenty Questions
data { 
  int<lower=0> np;  // rows
  int<lower=0> nq;  // columns
  int k[np,nq];
  int k_na[np,nq];  // locating NAs in k
  int n_na;         // number of NAs
}
parameters {
  real<lower=0,upper=1> p[np];
  real<lower=0,upper=1> q[nq];
} 
transformed parameters {
  real<lower=0,upper=1> theta[np,nq];

  // Probability Correct Is Product Of Question By Person Rates
  for (i in 1:np)
    for (j in 1:nq)
      theta[i,j] <- p[i] * q[j];
}
model {
  // Priors For People and Questions
  p ~ beta(1, 1);
  q ~ beta(1, 1);
    
  // Correctness Of Each Answer Is Bernoulli Trial
  for (i in 1:np)
    for (j in 1:nq)
      if (k_na[i,j] == 0)     // If k[i,j] is not missing
        k[i,j] ~ bernoulli(theta[i,j]);
}
generated quantities {
  int na_array[n_na];
  int index;
  
  index <- 1;
  for (i in 1:np) {
    for (j in 1:nq) {   
      if (k_na[i,j] == 1) {   // If k[i,j] is missing
        na_array[index] <- bernoulli_rng(theta[i,j]);
        index <- index + 1;
      }
    }
  }
}"

dset <- 1  # Chooes dataset/model

if (dset == 1) {
  k <- c(1,1,1,1,0,0,1,1,0,1,0,0,1,0,0,1,0,1,0,0,
        0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,1,0,0,0,1,1,0,0,0,0,1,0,0,0,0,0,0,0,
        0,0,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,
        1,0,1,1,0,1,1,1,0,1,0,0,1,0,0,0,0,1,0,0,
        1,1,0,1,0,0,0,1,0,1,0,1,1,0,0,1,0,1,0,0,
        0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,1,1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,1,0,1,
        1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,0)
}
if (dset == 2) {
  k <- c(1,1,1,1,0,0,1,1,0,1,0,0,NA,0,0,1,0,1,0,0,
        0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,1,0,0,0,1,1,0,0,0,0,1,0,0,0,0,0,0,0,
        0,0,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,
        1,0,1,1,0,1,1,1,0,1,0,0,1,0,0,0,0,1,0,0,
        1,1,0,1,0,0,0,1,0,1,0,1,1,0,0,1,0,1,0,0,
        0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
        0,0,0,0,NA,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,1,1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,1,0,1,
        1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,0,0,NA,0,0)
}

k <- matrix(k, nrow=10, byrow=T)
np <- nrow(k)
nq <- ncol(k)

if (dset == 1) {
  model <- model1
  
  data <- list(k=k, np=np, nq=nq)  # to be passed on to Stan
  parameters <- c("p", "q")  # parameters to be monitored
  
} else if (dset == 2) {
  model <- model2
  k_na <- is.na(k) + 0  # Matrix locating missing values: 1 = missing
  n_na <- sum(k_na)  # number of missing values
  k[is.na(k)] <- 99  # some numeric value, since Stan doesn't eat NAs
  
  # to be passed on to Stan:
  data <- list(k=k, np=np, nq=nq, k_na=k_na, n_na=n_na)  
  parameters <- c("p", "q", "na_array")  # parameters to be monitored   
}

myinits <- list(
  list(q=runif(nq), p=runif(np)))

# The following command calls Stan with specific options.
# For a detailed description type "?rstan".
samples <- stan(model_code=model,   
                data=data, 
                init=myinits,  # If not specified, gives random inits
                pars=parameters,
                iter=10000, 
                chains=1, 
                thin=1,
                # warmup = 100,  # Stands for burn-in; Default = iter/2
                # seed = 123  # Setting seed; Default is random seed
)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

print(samples, digits=3)
