data {
  int<lower=1> N;  // total number of tests
  int<lower=1> N_age;
  int<lower=1> N_eth;
  int<lower=1> N_edu;
  real baseline;
  real<lower=0, upper=1> sens;
  real<lower=0, upper=1> spec;
}
transformed data {
  int strata = 2 * N_age * N_eth * N_edu;
}
generated quantities {
  real beta_0 = baseline;
  // some difference by sex, unequal observations
  real beta_sex = normal_rng(0, 0.5);
  vector[2] pct_sex = [0.4, 0.6]';

  vector[N_age] pct_age = dirichlet_rng(rep_vector(2, N_age));
  vector[N_age] beta_age;
  for (n in 1:N_age) {
    beta_age[n] = std_normal_rng();
  }
  // increased risk with age, some variation in observations
  // {
  // vector[N_age] beta_age_tmp;
  // for (n in 1:N_age) {
  // beta_age_tmp[n] = normal_rng(0, 1);
  // }
  // beta_age = sort_asc(beta_age_tmp);
  // }

  vector[N_eth] pct_eth = dirichlet_rng(rep_vector(1, N_eth));
  vector[N_eth] beta_eth;
  for (n in 1:N_eth) {
    beta_eth[n] = std_normal_rng();
  }

  vector[N_edu] pct_edu = dirichlet_rng(rep_vector(1.5, N_edu));
  vector[N_edu] beta_edu;
  for (n in 1:N_edu) {
    beta_edu[n] = std_normal_rng();
  }
  // decreased risk with edu, some variation in observations
  // {
  // vector[N_edu] beta_edu_tmp;
  // for (n in 1:N_edu) {
  // beta_edu_tmp[n] = normal_rng(0, 1);
  // }
  // beta_edu = sort_desc(beta_edu_tmp);
  // }

  array[strata] int sex, age, eth, edu, pos_tests, tests;
  array[strata] real p;
  array[strata] real p_sample;

  int idx = 1;
  for (i_sex in 1:2) {
    for (i_age in 1:N_age) {
      for (i_eth in 1:N_eth) {
        for (i_edu in 1:N_edu) {

	  // corresponds to unmodeled data inputs
          sex[idx] = i_sex; age[idx] = i_age; eth[idx] = i_eth; edu[idx] = i_edu;
          tests[idx] = to_int(pct_sex[i_sex] * pct_age[i_age] * pct_eth[i_eth] * pct_edu[i_edu] * N);

	  // corresponds to transformed parameters
          p[idx] = inv_logit(beta_0 + beta_sex * (i_sex)
                    + beta_age[i_age] + beta_eth[i_eth] +  beta_edu[i_edu]);
          p_sample[idx] = p[idx] * sens + (1 - p[idx]) * (1 - spec);

	  // corresponds to likelihood
          pos_tests[idx] = binomial_rng(tests[idx], p_sample[idx]);
          idx += 1;
        }
      }
    }
  }
}
