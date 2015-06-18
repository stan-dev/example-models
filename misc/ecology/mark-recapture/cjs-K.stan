/**
 * Cormack-Jolly-Seber Model
 * 
 * following section 1.2.1 of:
 * http://www.maths.otago.ac.nz/home/resources/theses/PhD_Matthew_Schofield.pdf
 *
 */
data {
  int<lower=2> K;                      // capture events
  int<lower=0> I;                      // number of individuals
  int<lower=0,upper=1> X[I,K];         // X[i,k]: individual i captured at k
}
transformed data {
  int<lower=0,upper=K+1> first[I];     // first[i]: ind i first capture
  int<lower=0,upper=K+1> last[I];      // last[i]:  ind i last capture
  int<lower=0,upper=I> n_captured[K];  // n_capt[k]: num aptured at k

  first <- rep_array(K+1,I);
  last <- rep_array(0,I);
  for (i in 1:I) {
    for (k in 1:K) {
      if (X[i,k] == 1) {
        if (k < first[i]) 
          first[i] <- k;
        if (k > last[i]) 
          last[i] <- k;
      }
    }
  }

  n_captured <- rep_array(0,K);
  for (i in 1:I)
    for (k in 1:K)
      n_captured[k] <- n_captured[k] + X[i,k];
}
parameters {
  vector<lower=0,upper=1>[K-1] phi;  // phi[k]: Pr[alive at k + 1 | alive at k]
  vector<lower=0,upper=1>[K] p;      // p[k]: Pr[capture at k]

  // note:  p[1] not used in model and hence not identified
}
transformed parameters {
  vector<lower=0,upper=1>[K] chi;   // chi[k]: Pr[no capture >  k | alive at k]
  {
    int k;
    chi[K] <- 1.0;              
    k <- K - 1;
    while (k > 0) {
      chi[k] <- (1 - phi[k]) + phi[k] * (1 - p[k+1]) * chi[k+1]; 
      k <- k - 1;
    }
  }
}
model {
  for (i in 1:I) {
    if (last[i] > 0) {
      for (k in (first[i]+1):last[i]) {
        increment_log_prob(log(phi[k-1]));     // i survived from k-1 to k
        if (X[i,k] == 1)
          increment_log_prob(log(p[k]));       // i captured at k
        else
          increment_log_prob(log1m(p[k]));     // i not captured at k
      }
      increment_log_prob(log(chi[last[i]]));   // i not seen after last[i]
    }
  }
}
generated quantities {
  // phi[K-1] and p(K) not identified, but product is
  real beta;
  vector<lower=0>[K] pop_hat;  // population

  beta <- phi[K-1] * p[K];

  for (k in 1:K)
    pop_hat[k] <- n_captured[k] / p[k];  
}
