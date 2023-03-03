#include funs_and_data.stan
#include tdata.stan

generated quantities {
  matrix[N_train, J*B] PHI_train_gq = PHI_train;
  matrix[N_test, J*B] PHI_test_gq = PHI_test;
}
