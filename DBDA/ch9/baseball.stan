data {
    int<lower=0> P;                                // number of players
    int<lower=0> N_positions;                      // number of positions
    int<lower=1, upper=N_positions> positions[P];  // position for player

    int<lower=1, upper=P+1> start[N_positions+1];  // start index for each group of players
    int<lower=0> N[P];                             // number of attemps per player
    int<lower=0> y[P];                             // number of hits per player
    
    real<lower=0> A;
    real<lower=0> B;
    real<lower=0> S;
    real<lower=0> R;
}
parameters {
    real<lower=0, upper=1> omega0;
    real<lower=0, upper=1> omega[N_positions];  // group level hit rate
    real<lower=0> kappaMinusTwo0;
    real<lower=0> kappaMinusTwo[N_positions];
    real<lower=0, upper=1> theta[P];            // player level hit rate
} 
transformed parameters {
    real<lower=0> kappa0;
    real<lower=0> alpha0;
    real<lower=0> beta0;
    real<lower=0> kappa[N_positions];
    real<lower=0> alpha[N_positions];
    real<lower=0> beta[N_positions];

    kappa0 = kappaMinusTwo0 + 2;
    alpha0 = omega0 * (kappa0 - 2) + 1;
    beta0 = (1 - omega0) * (kappa0 - 2) + 1;

    for (c in 1:N_positions) {
        kappa[c] = kappaMinusTwo[c] + 2;
        alpha[c] = omega[c] * (kappa[c] - 2) + 1;
        beta[c] = (1 - omega[c]) * (kappa[c] - 2) + 1;
    }
}
model {
    omega0 ~ beta(A, B);
    kappaMinusTwo0 ~ gamma(S, R);
    
    for (c in 1:N_positions) {
        omega[c] ~ beta(alpha0, beta0);
        kappaMinusTwo[c] ~ gamma(S, R);
    }
    
    for (i in 1:P) {
        theta[i] ~ beta(alpha[positions[i]], beta[positions[i]]);
        y[i] ~ binomial(N[i], theta[i]);
    }
}
