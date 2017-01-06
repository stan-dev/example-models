data {
  int n_age;
  int marr_j[n_age + 1];
  int marr_a[n_age + 1];
}

parameters {
  real<lower=0,upper=1> sjuv;   // Juv. survival
  real<lower=0,upper=1> ssub;   // Subad. survival
  real<lower=0,upper=1> sad;    // Ad. survival
  real<lower=0,upper=1> rjuv;   // Juv. recovery
  real<lower=0,upper=1> rad;    // Ad. recovery
}

transformed parameters {
  simplex[n_age + 1] pr_a;
  simplex[n_age + 1] pr_j;

  // Define the cell probabilities of the juvenile m-array
  // First element
  pr_j[1] = (1 - sjuv) * rjuv;

  // Second element
  pr_j[2] = sjuv * (1 - ssub) * rad;

  // Third and further elements
  for (t in 3:n_age)
    pr_j[t] = sjuv * ssub * pow(sad, t - 3) * (1 - sad) * rad;

  // Probability of non-recovery
  pr_j[n_age + 1] = 1 - sum(pr_j[:n_age]);

  // Define the cell probabilities of the adult m-array
  // All elements
  for (t in 1:n_age)
    pr_a[t] = pow(sad, t - 1) * (1 - sad) * rad;

  // Probability of non-recovery
  pr_a[n_age + 1] = 1 - sum(pr_a[:n_age]);
}

model {
  // Priors
  sjuv ~ beta(4.2, 2.8); // Informative prior for juv. survival: Analysis A
  //sjuv ~ uniform(0, 1);  // Non-informative for juv. survival prior: Analysis B
  // Uniform priors are implicitly defined.
  //  ssub ~ uniform(0, 1);
  //  sad ~ uniform(0, 1);
  //  rjuv ~ uniform(0, 1);
  //  rad ~ uniform(0, 1);

  // Define the multinomial likelihood
  marr_j ~ multinomial(pr_j);
  marr_a ~ multinomial(pr_a);
}
