# Project the beta
project_weights_lm <- function(X, mu) {
  XtX <- t(X) %*% X
  Xtm <- t(X) %*% mu
  n <- nrow(XtX)
  solve(XtX + 1e-9 * diag(nrow = n, ncol = n), Xtm)
}

# Euclidean norm squared
squared_euc_norm <- function(x) {
  v <- as.vector(x)
  sum(v * v)
}

# Project sigma
project_sigma <- function(sigma, mu, mu_proj) {
  v <- mu_proj - mu
  n <- length(v)
  sigma_proj2 <- sigma^2 + 1 / n * squared_euc_norm(v)
  sqrt(sigma_proj2)
}

# KL divergence of multivariate normals with S = s^2 I
kl_divergence <- function(mu, mu_proj, sigma, sigma_proj) {
  term1 <- log(sigma_proj / sigma) + (sigma^2) / (2 * sigma_proj^2) - 0.5
  term2 <- squared_euc_norm(mu - mu_proj) / (2 * sigma_proj^2)
  length(mu) * term1 + term2
}

# Compute submodel weights for one component
compute_beta_j <- function(xi_j, alpha_j, ell_j, seq_b, L) {
  bbb <- (pi^2 * seq_b^2) / (8 * L^2)
  sqrt_s_j <- alpha_j * sqrt(ell_j * sqrt(2 * pi) * exp(-ell_j^2 * bbb))
  xi_j * sqrt_s_j
}

# Compute all submodel weights
compute_beta <- function(xi_sub, pars, J, B, L) {
  beta <- c()
  seq_B <- seq_len(B)
  for (j in seq_len(J)) {
    xi_j <- xi_sub[inds_submatrix(j, B)]
    beta_j <- compute_beta_j(xi_j, 1.0, pars$ell[j], seq_B, L)
    beta <- c(beta, beta_j)
  }
  beta
}

# Parameter vector on log scale to list of params
logpar_vector_to_params <- function(p, J) {
  list(
    ell = exp(p[1:J]),
    sigma = exp(p[J + 1])
  )
}
# Project by fixing xi
project_xi_fixed <- function(PHI_sub, xi_sub, mu_ref, sigma_ref, J, B) {
  L <- 3 # todo: fix
  target_fun <- function(p) {
    pars <- logpar_vector_to_params(p, J)
    beta <- compute_beta(xi_sub, pars, J, B, L)
    mu <- PHI_sub %*% beta
    kl_divergence(mu_ref, mu, sigma_ref, pars$sigma)
  }
  p0 <- rep(0, J + 1) # init
  res <- optim(p0, target_fun) # run optimizer
  pars_optim <- logpar_vector_to_params(res$par, J)
  beta_optim <- compute_beta(xi_sub, pars_optim, J, B, L)

  # return
  c(pars_optim, list(beta = beta_optim))
}
