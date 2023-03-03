#include funs_and_data.stan
#include tdata.stan
parameters {
  real<lower=1e-12> alpha[J]; // component magnitudes
  real<lower=1e-12> ell[J]; // lengthscales
  real<lower=1e-12> sigma;
  vector[num_xi] xi;
}

transformed parameters {
  vector[num_xi] glm_b = STAN_build_glm_b(B, L, alpha, ell, xi);
}

model {
  xi ~ normal(0, 1);
  alpha ~ student_t(20, 0, 1);
  ell ~ lognormal(0, 1);
  sigma ~ normal(0, 1); //inv_gamma(2, 1);
  y_train ~ normal_id_glm(PHI_train, glm_a, glm_b, sigma);
}

generated quantities {
  vector[N_train] f_comps_train[J] = STAN_build_gp_components(
    X_train, X_train, B, L, alpha, ell, xi
  );
  vector[N_test] f_comps_test[J] = STAN_build_gp_components(
    X_train, X_test, B, L, alpha, ell, xi
  );
  vector[N_train] f_sum_train = rep_vector(0.0, N_train);
  vector[N_test] f_sum_test = rep_vector(0.0, N_test);
  vector[N_train] loglik_train;
  vector[N_test] loglik_test;
  for(j in 1:J){
    f_sum_train = f_sum_train + f_comps_train[j];
    f_sum_test = f_sum_test + f_comps_test[j];
  } 
  loglik_train = log_likelihood_gaussian(f_sum_train, y_train, sigma);
  loglik_test = log_likelihood_gaussian(f_sum_test, y_test, sigma);
}
