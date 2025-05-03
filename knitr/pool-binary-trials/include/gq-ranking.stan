/** include in generated quantities block of model with parameter
    vector<lower=0, upper=1>[N] theta; // chance of success
    and data N (number of observations)
*/
  array[N] int<lower=1, upper=N> rnk; // rank of player n
  {
    array[N] int dsc;
    dsc = sort_indices_desc(theta);
    for (n in 1 : N) {
      rnk[dsc[n]] = n;
    }
  }
  array[N] int<lower=0, upper=1> is_best; // Pr[player n highest chance of success]
  for (n in 1 : N) {
    is_best[n] = rnk[n] == 1;
  }
