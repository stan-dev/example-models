functions { 
  real quantile(vector x, real p){
    int n;            // length of vector x
    real index;       // integer index of p
    int lo;           // lower integer cap of the index
    int hi;           // higher integer cap of the index
    real h;           // generated weight between the lo and hi
    real qs;          // weighted average of x[lo] and x[hi]
    n = num_elements(x);           
    index = 1 + (n - 1)*p;             
    lo = 1; 
    while ((lo + 1) < index) 
      lo = lo + 1; 
    hi = lo + 1; 
    h = index - lo; 
    qs = (1 - h)*sort_asc(x)[lo] + h*sort_asc(x)[hi];     
    return qs;                
  }
}
data {
  int<lower=0> N;                  // sample size
  int<lower=0> N_cov;              // number of covariates
  real y[N];                       // continuous outcome
  int<lower=0,upper=1> y_pos[N];   // indicator for y > 0
  int<lower=0,upper=1> z[N];       // treatment assigned
  vector[N_cov] x[N];              // covariates
  vector[N_cov] xz_inter[N];       // interaction terms
  real<lower=-1,upper=1> rho;       // assumed correlation between the potential outcomes
}

parameters {
  // Binary Sub-Model
  real alpha_bin;                    // intercept
  row_vector[N_cov] beta_bin;        // coefficients for x[N]
  row_vector[N_cov] beta_inter_bin;  // coefficients for x_inter[N] 
  real tau_bin;                      // treatment effect
  
  // Continuous Sub-Model
  real alpha_cont;                   // intercept       
  row_vector[N_cov] beta_cont;       // coefficients for x[N]
  row_vector[N_cov] beta_inter_cont; // coefficients for x_inter[N]
  real tau_cont;                     // treatment effect for continuous part
  
  real<lower=0> sigma_t;             // residual SD for the treated
  real<lower=0> sigma_c;             // residual SD for the control
}

model {
   // Priors for the Binary Sub-Model  
   alpha_bin ~ student_t(5, 0, 2.5);  // student t with high df
   beta_bin ~ student_t(5, 0, 2.5);
   beta_inter_bin ~ student_t(5, 0, 2.5);
   tau_bin ~ student_t(5, 0, 2.5);
   
   // Priors for the Continuous Sub-Model       
   alpha_cont ~ normal(0, 100);       // unscaled normal priors for coefficients
   beta_cont ~ normal(0, 100);  
   beta_inter_cont ~ normal(0, 100); 
   tau_cont ~ normal(0, 100);
   sigma_c ~ normal(0, 100);          
   sigma_t ~ normal(0, 100);
   // sigma_c ~ inv_gamma(1, 0.01);          
   // sigma_t ~ inv_gamma(1, 0.01);

   // LIKELIHOOD
   for(n in 1:N) {
     // Binary Sub-Model
     y_pos[n] ~ bernoulli_logit(alpha_bin + beta_bin * x[n] + beta_inter_bin * xz_inter[n] + tau_bin * z[n]);
     
     // Continuous Sub-Model (if y > 0)
     if(y_pos[n] == 1)
        log(y[n]) ~ normal(alpha_cont + beta_cont * x[n] + beta_inter_cont * xz_inter[n] + tau_cont * z[n], sigma_t*z[n] + sigma_c*(1-z[n]));
   }
}
generated quantities{
  // finite sample average treatment effect (ATE) & mean values
  real y0_bar;
  real y1_bar;
  real tau_samp;   
  
  // finite sample quantile treatment effects (QTEs)
  real tau_qte25;  
  real tau_qte50;
  real tau_qte75;
  { // to create temporary variables
    // Generate science table
    real y0[N];
    real y1[N];
    real tau_ind[N];  // individual treatment effects (ITE)
    
    // temporary variable for positive y_pred
    int y_pred_pos;
    for(n in 1:N){
      // Predicted chance of success: Binary part
      real theta_c = inv_logit(alpha_bin + beta_bin * x[n]); 
      real theta_t = inv_logit(alpha_bin + beta_bin * x[n] + beta_inter_bin * x[n] + tau_bin);
      
      // Predicted log mean: Continuous part
      real mu_c = alpha_cont + beta_cont * x[n];  
      real mu_t = alpha_cont + beta_cont * x[n] + beta_inter_cont * x[n] + tau_cont;   
            
      if(z[n] == 1){  // (1) Impute Y_mis(0)    
        // Binary Sub-model: predict positive y_pred(0)
        y_pred_pos = bernoulli_rng(theta_c);
        if(y_pred_pos == 0){ 
          y0[n] = 0;
        }else{
          y0[n] = exp(normal_rng(mu_c + rho*(sigma_c/sigma_t)*(y[n] - mu_t), sigma_c*sqrt(1 - rho^2)));
        }
        // Fill in Y_obs(1)
        y1[n] = y[n];
      }else{  // z[n] == 0,  (2) Impute Y_mis(1)    
        // Fill in Y_obs(0)
        y0[n] = y[n];
        
        // Binary Sub-model: predict positive y_pred(1)
        y_pred_pos = bernoulli_rng(theta_t);
        if(y_pred_pos == 0){
          y1[n] = 0;
        }else{
          y1[n] = exp(normal_rng(mu_t + rho*(sigma_t/sigma_c)*(y[n] - mu_c), sigma_t*sqrt(1 - rho^2)));
        }
      }
      tau_ind[n] = y1[n] - y0[n];  // generate ITE   
    } // end for loop

    // store mean values and ATE
    y0_bar = mean(y0);
    y1_bar = mean(y1);
    tau_samp = mean(tau_ind);
    
    // store QTEs
    tau_qte25 = quantile(to_vector(y1), 0.25) - quantile(to_vector(y0), 0.25); 
    tau_qte50 = quantile(to_vector(y1), 0.50) - quantile(to_vector(y0), 0.50); 
    tau_qte75 = quantile(to_vector(y1), 0.75) - quantile(to_vector(y0), 0.75); 
  } // end temporary variables
}

