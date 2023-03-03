
functions {
  
  // Create vector with elements 1, ..., N
  vector STAN_seq_len(data int N) {
    vector[N] v = rep_vector(1.0, N);
    for(n in 2:N) v[n] = n;
    return(v);
  }
  
  // Basis function matrix (EQ kernel)
  // - x: covariate vector
  // - B: number of basis functions
  // - L: domain size
  matrix STAN_basisfun_eq(vector x, data int B, real L) {
    int N = num_elements(x);
    vector[B] seq_B = STAN_seq_len(B);
    matrix[N, B] mat_B = transpose(rep_matrix(seq_B, N));
    matrix[N, B] mat_X = rep_matrix(x+L, B);
    matrix[N, B] PHI = 1.0/sqrt(L)*sin(0.5*pi()/L * mat_X .* mat_B);
    return(PHI);
  }
  
  // Create all PHI matrices, stacked horizontally in one matrix
  matrix STAN_create_basisfun_mats(data vector[] X, data int B, vector L) 
  {
    int N = num_elements(X[1]);
    int J = size(X);
    matrix[N, J*B] PHI_hstack;
    for(j in 1:J) {
      int is = 1+(j-1)*B;
      int ie = j*B;
      PHI_hstack[:, is:ie] = STAN_basisfun_eq(X[j], B, L[j]);
    }
    return(PHI_hstack);
  }
  
  // Compute log spectral density of EQ kernel
  //vector STAN_log_spd_eq(real alpha, real ell, vector omega){
  //  return 2*log(alpha)+log(ell)+0.5*log(2*pi())-0.5*ell^2*omega .* omega;
  //}
  
  // Compute the multipliers s_b
  vector STAN_basisfun_eq_multipliers(real alpha, real ell, 
      data int B, real L) {
    vector[B] seq_B = STAN_seq_len(B);
    return alpha^2*ell*sqrt(2*pi())*exp(-(ell^2*pi()^2)/(8.0*L^2) * seq_B .* seq_B);
  }

  // Build the regression coefficients for linear regression
  vector STAN_build_glm_b(data int B, vector L, 
      real[] alpha, real[] ell, vector xi)
  {
    int num_xi = num_elements(xi);
    int J = size(alpha);
    vector[num_xi] weights;

    // Build the weights
    for(j in 1:J) {
      int is = 1 + (j-1)*B;
      int ie = j*B;
      weights[is:ie] = STAN_basisfun_eq_multipliers(alpha[j], ell[j], B, L[j]);
    }
    return(weights .* xi);
  }

  // Same as rep(x, times=N) in R
  vector STAN_rep_vector_times(data vector x, data int N) {
    return to_vector(rep_matrix(x, N));
  }
  
  // Build the GP components
  vector[] STAN_build_gp_components(
    data vector[] X_train, 
    data vector[] X_pred,
    data int B, 
    data vector L, 
    real[] alpha, 
    real[] ell, 
    vector xi
  ){
    int J = size(X_train);
    int num_xi = num_elements(xi);
    int N_pred = num_elements(X_pred[1]);
    vector[num_xi] glm_b = STAN_build_glm_b(B, L, alpha, ell, xi);
    vector[N_pred] f_comps[J];
    for(j in 1:J) {
      int is = 1 + (j-1)*B;
      int ie = j*B;
      matrix[N_pred, B] PHI_pred = STAN_basisfun_eq(X_pred[j], B, L[j]);
      f_comps[j] = PHI_pred*glm_b[is:ie];
    }
    return(f_comps);
  }
  
  // Log likelihood (Gaussian)
  vector log_likelihood_gaussian(vector f, vector y, real sigma) {
    int N = num_elements(y);
    vector[N] ll = rep_vector(0.0, N);
    for(n in 1:N){
      ll[n] = normal_lpdf(y[n] | f[n], sigma);
    } 
    return(ll);
  }

}

data {
  
  // Dimensions
  int<lower=1> N_train; // number of observations
  int<lower=0> N_test;  // number of observations
  int<lower=1> J;       // number of additive components

  // Covariates
  vector[N_train] X_train[J];
  vector[N_test] X_test[J];
  
  // Approximation domain and basis function amount
  int<lower=1> B;         // number of basis functions
  vector<lower=0>[J] L;
  
  // Observations
  vector[N_train] y_train;
  vector[N_test] y_test;
}
