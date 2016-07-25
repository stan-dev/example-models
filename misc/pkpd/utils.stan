/*
 * Author: S. Weber
 *
 * Various utlities
 */

// helper function for rle_int
int rle_elem_count(int[] set) {
  int U;
    U <- 1;
    for(i in 2:num_elements(set)) {
      if(set[i-1] != set[i])
	U <- U + 1;
    }
    return(U);
}

// repeated length encoding, see rle in R
int[] rle_int(int[] set) {
  int res[rle_elem_count(set)];
  int c;
  res[1] <- 1;
  c <- 1;
  for(i in 2:num_elements(set)) {
    if(set[i-1] == set[i]) {
      res[c] <- res[c] + 1;
    } else {
      c <- c + 1;
      res[c] <- 1;
    }
  }
  return(res);
}

/* calculate the absolute value of a - b in log-space with log(a)
   and log(b) given. Does so by squaring and taking the root, i.e.
   
   la = log(a)
   lb = log(b)
   
   sqrt( (a - b)^2 ) = sqrt( a^2 - 2 * a * b + b^2 )
   
   <=> 0.5 * log_diff_exp(log_sum_exp(2*la, 2*lb), log(2) + la + lb)
*/
real log_diff_exp_abs(real la, real lb) {
  return(0.5 * log_diff_exp(log_sum_exp(2*la, 2*lb), log(2) + la + lb));
}


/* find_interval, see findInterval from R
 * i.e. find the ranks of x in sorted; sorted is assumed to be weakly-sorted.
 */
int[] find_interval_slow(vector x, vector sorted) {
  int res[num_elements(x)];
  // very brutal and ineffcient way of doing this, but well, it's
  // C++ speed...
  for(i in 1:num_elements(x)) {
    res[i] <- rank(append_row(rep_vector(x[i], 1), sorted), 1);
  }
  return(res);
}

/* faster version which uses bisectioning search
 */
int find_interval_elem(real x, vector sorted, int start_ind) {
  int res;
  int N;
  int max_iter;
  real left;
  real right;
  int left_ind;
  int right_ind;
  int iter;
    
  N <- num_elements(sorted);
  
  if(N == 0) return(0);
  
  left_ind  <- start_ind;
  right_ind <- N;
  
  max_iter <- 100 * N;
  left  <- sorted[left_ind ] - x;
  right <- sorted[right_ind] - x;
  
  if(0 <= left)  return(left_ind-1);
  if(0 == right) return(N-1);
  if(0 >  right) return(N);
  
  iter <- 1;
  while((right_ind - left_ind) > 1  && iter != max_iter) {
    int mid_ind;
    real mid;
    // is there a controlled way without being yelled at with a
    // warning?
    mid_ind <- (left_ind + right_ind) / 2;
    mid <- sorted[mid_ind] - x;
    if (mid == 0) return(mid_ind-1);
    if (left  * mid < 0) { right <- mid; right_ind <- mid_ind; }
    if (right * mid < 0) { left  <- mid; left_ind  <- mid_ind; }
    iter <- iter + 1;
  }
  if(iter == max_iter)
    print("Maximum number of iterations reached.");
  return(left_ind);
}

int[] find_interval(vector x, vector sorted) {
  int res[num_elements(x)];
  for(i in 1:num_elements(x)) {
    res[i] <- find_interval_elem(x[i], sorted, 1);
  }
  return(res);
}

// takes as input x an ascending sorted vector x which allows to
// move the left starting index to be moved
int[] find_interval_asc(vector x, vector sorted) {
  int res[num_elements(x)];
  int last;
  last <- 1;
  for(i in 1:num_elements(x)) {
    res[i] <- find_interval_elem(x[i], sorted, last);
    if(res[i] > 0) last <- res[i];
  }
  return(res);
}

int[] find_interval_blocked(int[] vals_M, vector vals, int[] sorted_M, vector sorted) {
  int res[num_elements(vals)];
  int M;
  int v;
  int s;
  M <- num_elements(vals_M);
  v <- 1;
  s <- 1;
  for(m in 1:M) {
    int temp[vals_M[m]];
    temp <- find_interval(segment(vals, v, vals_M[m]), segment(sorted, s, sorted_M[m]));
    for(n in 1:vals_M[m])
      res[v + n - 1] <- temp[n];
    v <- v + vals_M[m];
    s <- s + sorted_M[m];
  }
  return(res);
}

// count number times elem appears in test set
int count_elem(int[] test, int elem) {
  int count;
  count <- 0;
  for(i in 1:num_elements(test))
    if(test[i] == elem)
      count <- count + 1;
  return(count);
}

// count number times elems appears in test set
int[] count_elems(int[] test, int[] elems) {
  int counts[num_elements(elems)];
  for(i in 1:num_elements(elems))
    counts[i] <- count_elem(test, elems[i]);
  return(counts);
}

// find elements in test which are equal to elem
int[] which_elem(int[] test, int elem) {
  int res[count_elem(test, elem)];
  int ci;
  ci <- 1;
  for(i in 1:num_elements(test))
    if(test[i] == elem) {
      res[ci] <- i;
      ci <- ci + 1;
    }
  return(res);
}

// divide fac by div and return the rounded down integer
int floor_div_int(real fac, real div) {
  int count;
  if(fac < 0)
    reject("floor_div_int only works for positive values.");
  count <- 1;
  while(count * div <= fac) { count <- count + 1; }
  count <- count - 1;
  return count;
}

int[] count_obs_event_free(int[] obs_timeRank, int ndose) {
  int dose_next_obs[ndose];
  int o;
  int O;
  dose_next_obs <- rep_array(0, ndose);
  o <- 0;
  O <- size(obs_timeRank);
  while (o < O && obs_timeRank[o+1] == 0) { o <- o + 1; }
  for (i in 1:ndose) {
    int count;
    count <- 0;
    while(o < O && obs_timeRank[o+1] == i) {
      o <- o + 1;
      count <- count + 1;
    }
    dose_next_obs[i] <- count;
  }
  return(dose_next_obs);
}

int[] count_obs_event_free_blocked(int[] M, int[] obs_timeRank, int[] ndose) {
  int dose_next_obs[sum(ndose)];
  int l;
  int ld;
  dose_next_obs <- rep_array(0, sum(ndose));
  l <- 1;
  ld <- 1;
  for (i in 1:size(M)) {
    int u;
    int ud;
    u <- l + M[i] - 1;
    ud <- ld + ndose[i] - 1;
    dose_next_obs[ld:ud] <- count_obs_event_free(obs_timeRank[l:u], ndose[i]);
    l <- u + 1;
    ld <- ud + 1;
  }
  return(dose_next_obs);
}

