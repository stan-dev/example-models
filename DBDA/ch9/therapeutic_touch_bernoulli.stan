data {
    int<lower=0> P;                      // # practitioners
    int<lower=0> N;                      // # trials
    int<lower=0, upper=1> y[N];          // outcomes (must be sliced to obtain observations for a particular practitioner)
    int<lower=1, upper=N+1> start[P+1];  // starting indexes for slicing vector 'outcomes' 
    int<lower=0> A;                      // alpha parameter for beta distribution (omega)
    int<lower=0> B;                      // beta parameter for beta distribution (omega)
    real<lower=0> S;                     // shape parameter for gamma distribution (kappa)
    real<lower=0> R;                     // rate parameter for gamma distribution (kappa)
}
parameters {
    real<lower=0, upper=1> theta[P];  // coin bias
    real<lower=0, upper=1> omega;     // mint bias
    real<lower=0> kappaMinusTwo;      // needed because value of kappa - 2 must be non-negative
}
transformed parameters {
    real<lower=0> alpha;  // alpha parameter for prior on coin bias
    real<lower=0> beta;   // beta parameter for prior on coin bias
    real<lower=2> kappa;  // concentration; dependence between omega and theta
    
    kappa = kappaMinusTwo + 2;
    alpha = omega * (kappa - 2) + 1;
    beta = (1 - omega) * (kappa - 2) + 1;
    
}
model {
    kappaMinusTwo ~ gamma(S, R);
    omega ~ beta(A, B);
    theta ~ beta(alpha, beta);
    
    for (p in 1:P) {
        y[start[p]:start[p+1]-1] ~ bernoulli(theta[p]);
    }
}
