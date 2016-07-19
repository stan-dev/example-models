functions {
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
int[] count_dose_given(vector time, vector dose_time, vector dose_tau, int[] dose_addl) {
  int dose_count[num_elements(time)];
  int time_rank[num_elements(time)];
  int o;
  int O;
  o <- 1;
  O <- num_elements(time);
  time_rank <- find_interval(time, dose_time);
  dose_count <- rep_array(0, O);
  //print("time_rank = ", time_rank);
  // first skip all dosings before the first time
  while(o < O && time_rank[o] == 0) { o <- o + 1; }
  //print("o = ", o);
  for(i in o:O) {
    int d;
    d <- time_rank[i];
    if(dose_tau[d] > 0)
      dose_count[i] <- min(floor_div_int(time[i] - dose_time[d], dose_tau[d]), dose_addl[d]);
  }
  return dose_count;
}

int[] count_dose_given_blocked(int[] M, vector time, int[] M_dose, vector dose_time, vector dose_tau, int[] dose_addl) {
  int dose_count[num_elements(time)];
  int B;
  int tl;
  int dl;
  B <- num_elements(M);
  tl <- 1;
  dl <- 1;
  for(b in 1:B) {
    int tu;
    int du;
    tu <- tl + M[b] - 1;
    du <- dl + M_dose[b] - 1;
    dose_count[tl:tu] <- count_dose_given(time[tl:tu], dose_time[dl:du], dose_tau[dl:du], dose_addl[dl:du]);
    tl <- tu + 1;
    dl <- du + 1;
  }
  return dose_count;
}

/*
 * Need to define an operator which develops a state forward in time
 * by an amount Dt. Input is the state at t, the log-coefficients of
 * the exponentials and the exponential exponents and Dt. Output is
 * the new state after Dt time units for one of the cmts (depends on
 * given coefs).
 *
 * lstate_ref is the initial
 *
 * coefsN is a 2-row matrix with the first row being the
 * log-coefficients and the second row the exponents.
 */
matrix evolve_lsystem(int S, vector Dt, matrix coefs, int[] coefs_map) {
  matrix[num_elements(Dt),S] lsystem;
  int T;
  T <- num_elements(Dt);
  
  // initialize to zero
  lsystem <- rep_matrix(-500, T, S);
  
  for(o in 1:cols(coefs)) {
    int s;
    vector[T] term;
    s    <- coefs_map[o];
    term <- coefs[1,o] + Dt * coefs[2,o];
    for(t in 1:T)
      lsystem[t,s] <- log_sum_exp(lsystem[t,s], term[t]);
  }
  return(lsystem);
}

real lgeometric_series(real la, real n) {
  return( log1m_exp(la * n) - log1m_exp(la) );
}

/*
 * models 2 cmts 1 and 2 which each have an elimination rate k1 / k2
 * and we have a flow from 1 to 2, k12. This models a 1cmt oral dosing
 * and/or a building block of a metabolite. Finally all mass exiting
 * cmt 2 is recorded in cmt 3.
 */
matrix pk_1cmt_metabolite_depot(vector lref, vector Dt, real lk1, real lk12, real lk20, real tau, real n) {
  matrix[2,8] coefs1;
  int coefs1_map[8];
  int coefs1_zero[8];
  matrix[2,6] coefs2;
  int coefs2_map[6];
  int coefs2_zero[6];
  // positive terms
  matrix[num_elements(Dt),3] lsystem1;
  // negative terms (due to Bateman function)
  matrix[num_elements(Dt),2] lsystem2;
  real ldeltaSq;
  real nk1;
  real nk20;

  coefs1_zero <- rep_array(0, 8);
  coefs2_zero <- rep_array(0, 6);
  
  ldeltaSq <- 2*log_diff_exp_abs(lk1, lk20);

  nk1  <- -exp(lk1);
  nk20 <- -exp(lk20);

  // setup coefficient matrix and coef index vectors
  coefs1_map[1] <- 1;
  coefs1[1,1] <- lref[1];
  coefs1[2,1] <- nk1;
  
  coefs1_map[2] <- 2;
  coefs1[1,2] <- lref[1] + lk12 - ldeltaSq + lk1;
  coefs1[2,2] <- nk20;
  coefs1_map[3] <- 2;
  coefs1[1,3] <- lref[1] + lk12 - ldeltaSq + lk20;
  coefs1[2,3] <- nk1;
  
  coefs1_map[4] <- 2;
  coefs1[1,4] <- lref[2];
  coefs1[2,4] <- nk20;

  // whatever is in the depot cmt doesn´t go away
  coefs1_map[5] <- 3;
  coefs1[1,5] <- log_sum_exp(lref[3], lref[2]);
  coefs1[2,5] <- 0;
  coefs1_zero[5] <- 1;
  
  coefs1_map[6] <- 3;
  coefs1[1,6] <- lref[1] + lk12 + lk20 - ldeltaSq + log_sum_exp(lk1 - lk20, lk20 - lk1);
  coefs1[2,6] <- 0;
  coefs1_zero[6] <- 1;
  
  coefs1_map[7] <- 3;
  coefs1[1,7] <- lref[1] + lk12 + lk20 - ldeltaSq;
  coefs1[2,7] <- nk1;
  
  coefs1_map[8] <- 3;
  coefs1[1,8] <- lref[1] + lk12 + lk20 - ldeltaSq;
  coefs1[2,8] <- nk20;
  
  // for the negative terms we only use a two cmts; hence 2 is
  // relabeled to 1, and 3 to 2
  coefs2_map[1] <- 1;
  coefs2[1,1] <- lref[1] + lk12 - ldeltaSq + lk1;
  coefs2[2,1] <- nk1;
  coefs2_map[2] <- 1;
  coefs2[1,2] <- lref[1] + lk12 - ldeltaSq + lk20;
  coefs2[2,2] <- nk20;

  coefs2_map[3] <- 2;
  coefs2[1,3] <- lref[2];
  coefs2[2,3] <- nk20;
  
  coefs2_map[4] <- 2;
  coefs2[1,4] <- lref[1] + lk12 - ldeltaSq + lk20 + log(2);
  coefs2[2,4] <- 0;
  coefs2_zero[4] <- 1;

  coefs2_map[5] <- 2;
  coefs2[1,5] <- lref[1] + lk12 - ldeltaSq + lk20 + lk1 - lk20;
  coefs2[2,5] <- nk20;
  
  coefs2_map[6] <- 2;
  coefs2[1,6] <- lref[1] + lk12 - ldeltaSq + lk20 + lk20 - lk1;
  coefs2[2,6] <- nk1;
  
  // in case the initial state is dosed in a regular pattern, we can
  // take advantage of the geometric series here by modifing the
  // coefficients
  if(n>1) {
    real logn;
    logn <- log(n);
    for(i in 1:8) {
      if(coefs1_zero[i]) {
        coefs1[1,i] <- coefs1[1,i] + logn;
      } else {
        coefs1[1,i] <- coefs1[1,i] + lgeometric_series(coefs1[2,i] * tau, n);
      }
    }
    for(i in 1:6) {
      if(coefs2_zero[i]) {
        coefs2[1,i] <- coefs2[1,i] + logn;
      } else {
        coefs2[1,i] <- coefs2[1,i] + lgeometric_series(coefs2[2,i] * tau, n);
      }
    }
  }
    
  //print("AFTER: coefs1 = ", coefs1);
  //print("AFTER: coefs2 = ", coefs2);

  lsystem1 <- evolve_lsystem(3, Dt, coefs1, coefs1_map);
  lsystem2 <- evolve_lsystem(2, Dt, coefs2, coefs2_map);

  //print("lsystem1 = ", lsystem1);
  //print("lsystem2 = ", lsystem2);

  // final system is the difference of the two solutions
  for(t in 1:num_elements(Dt)) {
    lsystem1[t,2] <- log_diff_exp(lsystem1[t,2], lsystem2[t,1]);
    lsystem1[t,3] <- log_diff_exp(lsystem1[t,3], lsystem2[t,2]);
  }

  return(lsystem1);
}

matrix pk_1cmt_metabolite(vector lref, vector Dt, real lk1, real lk12, real lk20, real tau, real n) {
  matrix[2,4] coefs1;
  int coefs1_map[4];
  matrix[2,2] coefs2;
  int coefs2_map[2];
  // positive terms
  matrix[num_elements(Dt),2] lsystem1;
  // negative terms (due to Bateman function)
  matrix[num_elements(Dt),1] lsystem2;
  real ldeltaSq;
  real nk1;
  real nk20;

  ldeltaSq <- 2*log_diff_exp_abs(lk1, lk20);

  nk1  <- -exp(lk1);
  nk20 <- -exp(lk20);

  // setup coefficient matrix and coef index vectors
  coefs1_map[1] <- 1;
  coefs1[1,1] <- lref[1];
  coefs1[2,1] <- nk1;
  
  coefs1_map[2] <- 2;
  coefs1[1,2] <- lref[1] + lk12 - ldeltaSq + lk1;
  coefs1[2,2] <- nk20;
  coefs1_map[3] <- 2;
  coefs1[1,3] <- lref[1] + lk12 - ldeltaSq + lk20;
  coefs1[2,3] <- nk1;
  
  coefs1_map[4] <- 2;
  coefs1[1,4] <- lref[2];
  coefs1[2,4] <- nk20;

  // for the negative terms we only use a two cmts; hence 2 is
  // relabeled to 1, and 3 to 2
  coefs2_map[1] <- 1;
  coefs2[1,1] <- lref[1] + lk12 - ldeltaSq + lk1;
  coefs2[2,1] <- nk1;
  coefs2_map[2] <- 1;
  coefs2[1,2] <- lref[1] + lk12 - ldeltaSq + lk20;
  coefs2[2,2] <- nk20;

  // in case the initial state is dosed in a regular pattern, we can
  // take advantage of the geometric series here by modifing the
  // coefficients
  if(n>1) {
    for(i in 1:4) {
      coefs1[1,i] <- coefs1[1,i] + lgeometric_series(coefs1[2,i] * tau, n);
    }
    for(i in 1:2) {
      coefs2[1,i] <- coefs2[1,i] + lgeometric_series(coefs2[2,i] * tau, n);
    }
  }
    
  //print("AFTER: coefs1 = ", coefs1);
  //print("AFTER: coefs2 = ", coefs2);

  lsystem1 <- evolve_lsystem(2, Dt, coefs1, coefs1_map);
  lsystem2 <- evolve_lsystem(1, Dt, coefs2, coefs2_map);

  //print("lsystem1 = ", lsystem1);
  //print("lsystem2 = ", lsystem2);

  // final system is the difference of the two solutions
  for(t in 1:num_elements(Dt)) {
    lsystem1[t,2] <- log_diff_exp(lsystem1[t,2], lsystem2[t,1]);
  }

  return(lsystem1);
}


// forward declare pk system functions
matrix pk_system(vector lref, vector Dt, vector theta);


matrix pk_system_addl(vector lref, vector Dt, int cmt, real lamt, real tau, int n, vector theta);

// model evaluation function taking dosing (and addl dosing) into
// account for a single patient
matrix pk_model_fast(vector dose_lamt, int[] dose_cmt, vector dose_time, vector dose_tau, int[] dose_addl,
                     vector init_lstate, real init_time,
                     vector obs_time, int[] obs_timeRank, int[] obs_dose_given,
                     vector theta,
                     vector lscale)
{
  int D;
  int O;
  int d;
  int o;
  int active_addl;
  row_vector[num_elements(init_lstate)] lstate_ref;
  row_vector[num_elements(init_lstate)] row_lscale;
  real time_ref;
  matrix[num_elements(obs_time),num_elements(init_lstate)] lstate;
  
  D <- num_elements(dose_lamt);
  O <- num_elements(obs_time);

  row_lscale <- to_row_vector(lscale);
    
  o <- 1;
  d <- 0;
  active_addl <- 0; // 0 = FALSE
  lstate_ref  <- to_row_vector(init_lstate);
  time_ref    <- init_time;
  // first skip all dosings before init_time
  while(d < D && dose_time[d+1] < init_time) { d <- d + 1; }
  // next, process all elements which are past active dosing
  while(o <= O) {
    // first update reference state to be just after the last
    // dosing
    while(d != obs_timeRank[o]) {
      int nd;
      //print("Advancing from dose ", d, " to ", obs_timeRank[o]);
      // the dose of the reference state is not yet the very
      // last dose given before this observation, add it
      nd <- d + 1;
      if(active_addl) {
        // in case of an active addl record, we have to
        // super-impose the developed reference state with the
        // emerged dosing
        lstate_ref <- pk_system_addl(to_vector(lstate_ref), rep_vector(dose_time[nd] - time_ref, 1), dose_cmt[d], dose_lamt[d], dose_tau[d], dose_addl[d], theta)[1];
      } else {
        lstate_ref <- pk_system(to_vector(lstate_ref), rep_vector(dose_time[nd] - time_ref, 1), theta)[1];
      }
      time_ref <- dose_time[nd];
      // add in the dosing, but only if we have a simple dosing
      // event, i.e. no additional dosings
      active_addl <- dose_addl[nd] > 0;
      if(!active_addl) {
        lstate_ref[dose_cmt[nd]] <- log_sum_exp(lstate_ref[dose_cmt[nd]], dose_lamt[nd]);
      }
      d <- nd;
    }
    // ok, evolve from last dose to current observation...
    if(active_addl) {
      int ndose;
      // ...in case of addl dosing, the effect of the multiple
      // dosing has not yet been added
      // note: I would prefer to keep ndose as int, but there is no
      // int floor(int) function available; TODO: pre-compute this counting vector!
      //ndose <- floor_div_int((obs_time[o] - dose_time[d]), dose_tau[d]);
      ndose <- obs_dose_given[o];
      if(ndose >= dose_addl[d]) {
        //print("Adding ", dose_addl[d] + 1 , " doses at once...");
        //print("merging multiple dosing stuff");
        // in this case all doses have been used and we can merge
        // and update the reference state
        lstate_ref <- pk_system_addl(to_vector(lstate_ref), rep_vector(obs_time[o] - time_ref, 1), dose_cmt[d], dose_lamt[d], dose_tau[d], dose_addl[d], theta)[1];
        time_ref <- obs_time[o];
        lstate[o] <- to_row_vector(lstate_ref) - row_lscale;
        active_addl <- 0;
      } else {
        //print("Adding ", ndose+1 , " doses at once...");
        lstate[o] <- pk_system_addl(to_vector(lstate_ref), rep_vector(obs_time[o] - time_ref, 1), dose_cmt[d], dose_lamt[d], dose_tau[d], ndose, theta)[1] - row_lscale;
      }
    } else {
      // ... which is simple for non-addl dosing as dose is
      // already merged
      lstate[o] <- pk_system(to_vector(lstate_ref), rep_vector(obs_time[o] - time_ref, 1), theta)[1] - row_lscale;
    }
    o <- o + 1;
  }
  return(lstate);
}

matrix pk_model(vector dose_lamt, int[] dose_cmt, vector dose_time, vector dose_tau, int[] dose_addl,
                vector init_lstate, real init_time,
                vector obs_time,
                vector theta,
                vector lscale) {
  return(pk_model_fast(dose_lamt, dose_cmt, dose_time, dose_tau, dose_addl,
                       init_lstate, init_time,
                       obs_time, find_interval(obs_time, dose_time), count_dose_given(obs_time, dose_time, dose_tau, dose_addl),
                       theta,
                       lscale));
}

matrix evaluate_model_fast(int[] dose_M, vector dose_lamt, int[] dose_cmt, vector dose_time, vector dose_tau, int[] dose_addl,
                           matrix init_lstate, vector init_time,
                           int[] obs_M, vector obs_time, int[] obs_timeRank, int[] obs_dose_given,
                           matrix theta,
                           matrix lscale) {
  matrix[num_elements(obs_time), cols(init_lstate)] lstate;
  int J;
  int d;
  int o;
  int S;
  
  //print("dose_M = ", dose_M);
  //print("obs_M = ", obs_M);
  
  J <- num_elements(dose_M);
  S <- cols(lscale);
  d <- 1;
  o <- 1;
  for(j in 1:J) {
    matrix[obs_M[j],S] lstate_j;
    int d_m;
    int o_m;
    //print("Processing patient ", j);
    d_m <- dose_M[j];
    o_m <- obs_M[j];
    lstate_j <- pk_model_fast(segment(dose_lamt, d, d_m), segment(dose_cmt, d, d_m), segment(dose_time, d, d_m), segment(dose_tau, d, d_m), segment(dose_addl, d, d_m)
                              ,to_vector(init_lstate[j]), init_time[j]
                              ,segment(obs_time, o, o_m), segment(obs_timeRank, o, o_m), segment(obs_dose_given, o, o_m)
                              ,to_vector(theta[j])
                              ,to_vector(lscale[j]));
    
    for(i in 1:o_m)
      lstate[i + o - 1] <- lstate_j[i];
    
    d <- d + d_m;
    o <- o + o_m;
  }
  return(lstate);
}
  
matrix evaluate_model(int[] dose_M, vector dose_lamt, int[] dose_cmt, vector dose_time, vector dose_tau, int[] dose_addl,
                      matrix init_lstate, vector init_time,
                      int[] obs_M, vector obs_time,
                      matrix theta,
                      matrix lscale) {
  return(evaluate_model_fast(dose_M, dose_lamt, dose_cmt, dose_time, dose_tau, dose_addl,
                             init_lstate, init_time,
                             obs_M, obs_time, find_interval_blocked(obs_M, obs_time, dose_M, dose_time), count_dose_given_blocked(obs_M, obs_time, dose_M, dose_time, dose_tau, dose_addl),
                             theta,
                             lscale));
}

matrix evaluate_model_nm(int[] id, vector time, int[] cmt, int[] evid, vector amt, vector tau, int[] addl, int[] mdv,
                         matrix init_lstate, vector init_time,
                         matrix theta, matrix lscale) {
  int dose_ind[count_elem(evid, 1)];
  dose_ind <- which_elem(evid, 1);
  
  return(evaluate_model(rle_int(id[dose_ind]), log(amt[dose_ind]), cmt[dose_ind], time[dose_ind], tau[dose_ind], addl[dose_ind],
                        init_lstate, init_time,
                        rle_int(id), time,
                        theta,
                        lscale));
}
  
  // we fit a 1-cmt oral dosing situation
  matrix pk_system(vector lref, vector Dt, vector theta) {
    // as we fitting a 1-cmt oral dosing situation such that k1=k12 (all
    // mass out of 1 goes to 2)
    return(pk_1cmt_metabolite(lref, Dt, theta[1], theta[1], theta[2], 0, 0));
  }

  // note that n are the additional doses to be added such that in total
  // n+1 are added
  matrix pk_system_addl(vector lref, vector Dt, int cmt, real lamt, real tau, int n, vector theta) {
    matrix[num_elements(Dt), num_elements(lref)] lstate;
    matrix[num_elements(Dt), num_elements(lref)] lstate_mdose;
    vector[num_elements(lref)] lref_mdose;
    int S;
    
    // evolve reference state freely...
    lstate <- pk_system(lref, Dt, theta);
    
    // ... and add the extra doses correctly time-shifted
    S <- num_elements(lref);
    lref_mdose <- rep_vector(-35, S);
    lref_mdose[cmt] <- lamt;
    //if(prod(Dt - tau * n) < 0) reject("All requested times must be past the last dosing addl event.");
    /*
      for(i in 1:num_elements(Dt))
      if((Dt[i] - tau * n) < 0)
    reject();
    ("All requested times must be past the last dosing addl event.");
    */
    lstate_mdose <- pk_1cmt_metabolite(lref_mdose, Dt - tau * n, theta[1], theta[1], theta[2], tau, n+1);
    for(s in 1:S)
      for(t in 1:num_elements(Dt))
        lstate[t,s] <- log_sum_exp(lstate_mdose[t,s], lstate[t,s]);
    return(lstate);
  }
  
}
data {
  int<lower = 1> N; // number of lines of nm data set
  vector<lower=0>[N] time;
  vector<lower=0>[N] amt;
  int cmt[N];
  int<lower=0, upper=1> mdv[N];
  int<lower=0, upper=2> evid[N];
  int<lower=1, upper=N> id[N];
  int<lower=0> addl[N];
  vector<lower=0>[N] tau;

  vector<lower=0>[N] dv; // observations

  vector[3] prior_theta_mean;
  vector<lower=0>[3] prior_theta_sd;
}
transformed data {
  int dose_ind[count_elem(evid, 1)];
  int obs_ind[count_elem(mdv, 0)];
  int obs_M[rle_elem_count(id)];
  int dose_M[rle_elem_count(id)];
  int obs_time_rank[count_elem(mdv, 0)];
  vector[count_elem(mdv, 0)] obs_time;
  vector[count_elem(mdv, 0)] obs_ldv;
  int obs_dose_given[count_elem(mdv, 0)];
  int obs_cmt[count_elem(mdv, 0)];
  vector[count_elem(evid, 1)] dose_time;
  vector[count_elem(evid, 1)] dose_tau;
  int dose_addl[count_elem(evid, 1)];
  vector[count_elem(evid, 1)] dose_lamt;
  int dose_cmt[count_elem(evid, 1)];
  int J;
  int O;
  row_vector[rle_elem_count(id)] zero;
  matrix[rle_elem_count(id),2] Init_lstate;
  vector[rle_elem_count(id)] init_time;

  dose_ind <- which_elem(evid, 1);
  obs_ind  <- which_elem(mdv , 0);
  
  // note: We implicitly assume here that every patient has at least
  // one dose. If not, this line breaks, but in any PK problem this
  // will be given.
  dose_M <- rle_int(id[dose_ind]);
  obs_M  <- rle_int(id[obs_ind] );

  obs_time <- time[obs_ind];
  obs_ldv  <- log(dv[obs_ind]);
  obs_cmt  <- cmt[obs_ind];

  dose_time <- time[dose_ind];
  dose_tau <- tau[dose_ind];
  dose_lamt <- log(amt[dose_ind]);
  dose_addl <- addl[dose_ind];
  dose_cmt <- cmt[dose_ind];
  
  obs_time_rank <- find_interval_blocked(obs_M, obs_time, dose_M, dose_time);
  obs_dose_given <- count_dose_given_blocked(obs_M, obs_time, dose_M, dose_time, dose_tau, dose_addl);
  
  J <- rle_elem_count(id);
  O <- count_elem(mdv, 0);

  zero <- rep_row_vector(0, J);

  Init_lstate <- rep_matrix(-25, J, 2);
  init_time   <- rep_vector(0, J);
}
parameters {
  ordered[2] theta_lelim;
  real theta_lV;
  vector<lower=0>[2] omega;
  matrix[2,J] xi;
  real<lower=0> sigma_y;
}
transformed parameters {
  vector[3] theta;
  matrix[3,J] Theta;
  row_vector<lower=0>[J] kDelta;

  // log(ka)
  theta[1] <- theta_lelim[2]; // absorbtion must be faster (hence larger) than elimination
  // log(ke)
  theta[2] <- theta_lelim[1];
  // log(V)
  theta[3] <- theta_lV;
  
  Theta[1] <- rep_row_vector(theta[1], J);
  // ncp parametrization
  //Theta[2:3] <- rep_matrix(theta[2:3], J) + diag_pre_multiply(omega, xi);
  // cp parametrization
  Theta[2:3] <- xi;

  // kDelta is only defined to ensure that we have a faster absorption
  // than elimination (avoid "flip-flop") for each patient
  kDelta <- Theta[1] - Theta[2];
}
model {
  vector[O] ipred;

  theta_lelim[2] ~ normal(prior_theta_mean[1], prior_theta_sd[1]);
  theta_lelim[1] ~ normal(prior_theta_mean[2], prior_theta_sd[2]);
  theta_lV       ~ normal(prior_theta_mean[3], prior_theta_sd[3]);

  // ncp parametrization
  //to_vector(xi) ~ normal(0, 1);
  // cp parametrization
  xi[1] ~ normal(theta[2], omega[1]);
  xi[2] ~ normal(theta[3], omega[2]);
  
  omega ~ normal(0, 1);
  sigma_y ~ normal(0, 1);

  {
    matrix[2,J] Lscale;
    matrix[O,2] ly;

    Lscale[1] <- zero;
    Lscale[2] <- Theta[3];
    
    ly <- evaluate_model_fast(dose_M, dose_lamt, dose_cmt, dose_time, dose_tau, dose_addl,
                              Init_lstate, init_time,
                              obs_M, obs_time, obs_time_rank, obs_dose_given,
                              Theta[1:2]',
                              Lscale');
    for (i in 1:O)
      ipred[i] <- ly[i, obs_cmt[i]];
  }

  obs_ldv ~ normal(ipred, sigma_y);
}
generated quantities {
}