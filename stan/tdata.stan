transformed data {
  matrix[N_train, J*B] PHI_train = STAN_create_basisfun_mats(X_train, B, L);
  matrix[N_test, J*B] PHI_test = STAN_create_basisfun_mats(X_test, B, L);
  vector[N_train] glm_a = rep_vector(0.0, N_train);
  int num_xi = J * B;
}
