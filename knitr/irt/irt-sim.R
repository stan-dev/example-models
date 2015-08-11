inv_logit <- 
function(x) { 
  return(1 / (1 + exp(-x))); 
}

logit <-
function(y) {
  return(log(y/(1 - y)));
}

I <- 20;        # questions
K <- 100;       # schools
L <- 20         # students / school
J <- K * L;     # students

N <- I * J;

ii <- rep(1:I,J);
jj <- rep(1:J,I);
kk <- rep(NA,N);
n <- 1;
for (k in 1:K) {
  for (l in 1:L) {
    for (i in 1:I) {
      kk[n] <- k;
      n <- n + 1;
    }
  }
}
# ll <- rep(NA,J); 
# for (n in 1:N)
#  ll[jj[n]] <- kk[n];

alpha <- rnorm(J,0,1);
beta <- rnorm(I,0,0.8);
sigma_gamma <- 0.62;
gamma <- rnorm(K,0,sigma_gamma);
delta <- runif(I,0.65,1.35);
y <- rep(NA,N);
for (n in 1:N) {
  y[n] <- rbinom(1,1,inv_logit(delta[ii[n]] * (alpha[jj[n]] + gamma[kk[n]] - beta[ii[n]])));
}
dump(c("I","J","K","N","ii","jj","kk","y"),
     file="irt.sim.data.R");

# Inits matching --init=0 in (R)Stan 1.3
sigma_gamma <- 1;
alpha <- rep(0,J);
beta <- rep(0,I);
gamma <- rep(0,K);
delta <- rep(1,I);

dump(c("sigma_gamma","alpha","beta","gamma","delta","zeta"),  
     file="irt.sim.inits.R");
