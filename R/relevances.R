# Get alpha means and medians
get_alpha_estimates <- function(fit, median = FALSE) {
  aa <- get_alpha_draws(fit)
  fun <- mean
  if (median) {
    fun <- stats::median
  }
  as.numeric(apply(aa, 2, fun))
}

# Get all component variances
comp_vars <- function(fit, num_comps) {
  fv <- rep(0, num_comps)
  for (j in 1:num_comps) {
    cd <- get_f_comp_draws(fit, j, test = FALSE)
    fv[j] <- mean(apply(cd, 1, var))
  }
  fv
}
