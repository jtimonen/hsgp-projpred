# Project to submodel
project_model_kopt <- function(inds, fit, sd) {
  str <- paste(inds, collapse = " ")
  cat("Projecting to submodel [", str, "] \n", sep = "")

  # Project draws
  pd <- project_draws_kopt(inds, fit, sd)

  # Compute mean log predictive densities
  mlpd_test <- mlpd_gaussian(sd$y_test, pd$mu_proj_test, pd$sigma_proj)
  mlpd_train <- mlpd_gaussian(sd$y_train, pd$mu_proj_train, pd$sigma_proj)
  kl <- mean(pd$kl_divs)
  metrics <- data.frame(kl, mlpd_train, mlpd_test)
  all_lpd_test <- log_lik_gaussian(sd$y_test, pd$mu_proj_test, pd$sigma_proj)

  # Output
  list(
    pd = pd,
    model = inds,
    metrics = metrics,
    all_lpd_test = all_lpd_test
  )
}


# Project all draws
project_draws_kopt <- function(inds, fit, sd) {
  X_train <- t(sd$X_train)
  X_test <- t(sd$X_test)
  colnames(X_train) <- sd$x_names
  colnames(X_test) <- sd$x_names
  L <- sd$L

  # Extract draws
  draws <- get_draws(fit)
  S <- dim(draws$f)[1]
  mu_pp_test <- list()
  mu_pp_train <- list()
  kl_divs <- rep(0, S)
  sigma_proj <- rep(0, S)
  counter <- 0
  weights <- list()
  terms_train <- list()
  terms_test <- list()
  mgcv_info <- list()
  sp <- list()

  for (s in seq_len(S)) {
    counter <- counter + 1
    perc <- round(100 * s / S, 2)
    if (counter >= S / 10) {
      counter <- 0
      cat(" * ", perc, "%", "\n", sep = "") # print progress
    }

    # Project parameters using training data
    pp <- project_draw_kopt(inds, X_train, sd$B, L, draws, s)

    # Predictions of projected model at test points
    yp <- rep(0, ncol(sd$X_test))
    mu_pp_test[[s]] <- yp
    mu_pp_train[[s]] <- pp$mu_proj
    sigma_proj[s] <- pp$sigma_proj
    kl_divs[s] <- pp$kl_div
    weights[[s]] <- 0 # pp$gam_fit$coefficients
  }
  ident <- function(x) {
    x
  }

  # Return
  list(
    sigma_proj = unlist(sigma_proj),
    mu_proj_test = sapply(mu_pp_test, ident),
    mu_proj_train = sapply(mu_pp_train, ident),
    kl_divs = kl_divs,
    weights = sapply(weights, ident)
  )
}

# Project a single draw
project_draw_kopt <- function(inds, X, B, L, draws, s) {
  # beta_ref <- as.numeric(draws$glm_b[s, ])
  mu_ref <- as.numeric(draws$f[s, ])
  sigma_ref <- draws$sigma[s]
  gam_fit <- project_gam(inds, mu_ref, X, B, L)
  mu_proj <- rep(0, length(mu_ref)) # gam_fit$fitted.values
  # sigma_proj_1 <- sqrt(p$sig2) # alt
  sigma_proj <- project_sigma(sigma_ref, mu_ref, mu_proj)
  # print(c(sigma_proj, sigma_proj_1))
  kl_div <- kl_divergence(mu_ref, mu_proj, sigma_ref, sigma_proj)

  list(
    kl_div = kl_div,
    mu_proj = mu_proj,
    sigma_proj = sigma_proj
  )
}
