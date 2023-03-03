# Compute half range of each column
half_ranges <- function(X1, X2) {
  D <- ncol(X1)
  X_hr <- rep(0, D)
  for (d in seq_len(D)) {
    tmp1 <- abs(X1[, d])
    tmp2 <- abs(X2[, d])
    X_hr[d] <- max(c(tmp1, tmp2))
  }
  return(X_hr)
}

# Create data
create_stan_data <- function(X_train, y_train, X_test, y_test, scale_bf, B) {
  N_train <- N_test <- nrow(X_test)
  J <- ncol(X_train)
  hranges <- half_ranges(X_train, X_test)
  L <- array(scale_bf * hranges)
  message(" * printing L\n")
  print(L)
  list(
    N_train = nrow(X_train),
    N_test = nrow(X_test),
    J = J,
    B = B,
    L = L,
    X_train = t(X_train),
    X_test = t(X_test),
    y_train = y_train,
    y_test = y_test
  )
}


# Create input data list
create_stan_input <- function(dat, split, scale_bf, B, x_cols = NULL) {
  if (is.null(x_cols)) {
    x_cols <- seq_len(ncol(dat$x))
  }
  X_train <- as.matrix(dat$x[split$train, x_cols])
  y_train <- dat$y[split$train]
  X_test <- as.matrix(dat$x[split$test, x_cols])
  y_test <- dat$y[split$test]
  sd <- create_stan_data(X_train, y_train, X_test, y_test, scale_bf, B)
  sd$split <- split
  sd$x_names <- colnames(dat$x)
  sd
}

# Run the GQ model
generate_phi <- function(sd, model_file) {
  g <- stan_model(model_file)
  message("Generating PHI (sampling with chains=1, iter=1, algorithm=Fixed_param).")
  fit <- sampling(g,
    iter = 1, chains = 1, data = sd, algorithm = "Fixed_param",
    refresh = 0
  )
  list(
    train = extract(fit, "PHI_train_gq")$PHI_train_gq[1, , ],
    test = extract(fit, "PHI_test_gq")$PHI_test_gq[1, , ]
  )
}

# Get alpha params
get_alpha_draws <- function(fit) {
  get_scalar_draws(fit, "alpha")
}

# Get draws of one component
get_f_comp_draws <- function(fit, idx, test = FALSE) {
  if (test) {
    nam <- "f_comps_test"
  } else {
    nam <- "f_comps_train"
  }
  out <- rstan::extract(fit, nam)[[nam]]
  if (is.null(idx)) {
    return(out)
  }
  out[, idx, ]
}

# Get draws of the sum f
get_f_draws <- function(fit, test = FALSE) {
  if (test) {
    nam <- "f_sum_test"
  } else {
    nam <- "f_sum_train"
  }
  rstan::extract(fit, nam)[[nam]]
}

# Get component estimate
get_f_comp_estimate <- function(fit, idx, est = mean, test = FALSE) {
  fc <- get_f_comp_draws(fit, idx, test = test)
  apply(fc, 2, est)
}

# Get component quantiles
get_f_comp_quantiles <- function(fit, idx, test = FALSE) {
  ff <- function(x) {
    quantile(x, probs = c(0.05, 0.5, 0.95))
  }
  fc <- get_f_comp_draws(fit, idx, test = test)
  apply(fc, 2, ff)
}

# Get draws of scalar parameter
get_scalar_draws <- function(fit, nam) {
  rstan::extract(fit, nam)[[nam]]
}

# Get draws of scalar parameter
get_mlpd_draws <- function(fit, test = FALSE) {
  if (test) {
    nam <- "loglik_test"
  } else {
    nam <- "loglik_train"
  }
  a <- rstan::extract(fit, nam)[[nam]]
  apply(a, 1, mean) # mean over data points
}

# Get draws of several parameters and quantities
get_draws <- function(fit) {
  list(
    f = get_f_draws(fit, test = FALSE), #  dim = c(num_draws, num_train)
    sigma = get_scalar_draws(fit, "sigma"),
    xi = extract(fit, "xi")$xi,
    alpha = extract(fit, "alpha")$alpha,
    ell = extract(fit, "ell")$ell,
    glm_b = extract(fit, "glm_b")$glm_b
  )
}
