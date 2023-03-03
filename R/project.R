# Column indices matrix corresponding to one additive term
inds_submatrix_single <- function(j, B) {
  is <- 1 + (j - 1) * B
  ie <- j * B
  is:ie
}

# Get column indices of a matrix corresponding to model specified by inds
inds_submatrix <- function(inds, B) {
  cinds <- c()
  for (i in inds) {
    cinds <- c(cinds, inds_submatrix_single(i, B))
  }
  cinds
}

# Design matrix of a submodel
design_submatrix <- function(PHI, B, inds) {
  cinds_phi <- inds_submatrix(inds, B)
  PHI[, cinds_phi]
}

# Build formula for GAM projection
project_gam.build_formula <- function(inds, B, L) {
  form <- "y ~ 1"
  bs_str <- getOption("bs.type", default = "hs")

  j <- 0
  for (idx in inds) {
    j <- j + 1
    f_j <- paste0(
      "s(x", idx, ", bs='", bs_str, "', k=", B,
      ", xt=list(dom_size=", L[idx], ")", ")"
    )
    form <- paste0(form, " + ", f_j)
  }
  as.formula(form)
}

# Project to a GAM using mgcv
project_gam.mgcv <- function(form, dat) {
  min_sp <- getOption("min.sp", default = NULL)
  if (!is.null(min_sp)) {
    D <- length(attr(terms(form), "term.labels"))
    min_sp <- rep(min_sp, D)
  }
  mgcv::gam(form, data = dat, min.sp = min_sp)
}

# Project to a submodel that is a GAM
project_gam <- function(inds, mu_ref, X, B, L) {
  y <- as.numeric(mu_ref)
  dat <- data.frame(cbind(y, X))
  form <- project_gam.build_formula(inds, B, L)
  project_gam.mgcv(form, dat)
}

# Project a single draw
project_draw <- function(inds, X, B, L, draws, s) {
  beta_ref <- as.numeric(draws$glm_b[s, ])
  mu_ref <- as.numeric(draws$f[s, ])
  sigma_ref <- draws$sigma[s]
  gam_fit <- project_gam(inds, mu_ref, X, B, L)
  mu_proj <- gam_fit$fitted.values
  # sigma_proj_1 <- sqrt(p$sig2) # alt
  sigma_proj <- project_sigma(sigma_ref, mu_ref, mu_proj)
  # print(c(sigma_proj, sigma_proj_1))
  kl_div <- kl_divergence(mu_ref, mu_proj, sigma_ref, sigma_proj)

  list(
    gam_fit = gam_fit,
    kl_div = kl_div,
    mu_proj = mu_proj,
    sigma_proj = sigma_proj
  )
}

# Get optimization info
get_mgcv_info <- function(x) {
  conv <- unlist(x$mgcv.conv)
  names <- c("sig2", "gcv.ubre.dev", "aic", "deviance", "boundary")
  x <- unlist(x[names])
  c(x, conv)
}

# Project all draws
project_draws <- function(inds, fit, sd) {
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
      # cat(" * ", perc, "%", "\n", sep = "") # print progress
    }

    # Project parameters using training data
    pp <- project_draw(inds, X_train, sd$B, L, draws, s)

    # Predictions of projected model at test points
    yp <- predict(pp$gam_fit, newdata = data.frame(X_test))
    mu_pp_test[[s]] <- yp
    mu_pp_train[[s]] <- pp$mu_proj
    sigma_proj[s] <- pp$sigma_proj
    kl_divs[s] <- pp$kl_div
    weights[[s]] <- pp$gam_fit$coefficients
    terms_train[[s]] <- predict(pp$gam_fit, newdata = data.frame(X_train), type = "terms")
    terms_test[[s]] <- predict(pp$gam_fit, newdata = data.frame(X_test), type = "terms")
    mgcv_info[[s]] <- get_mgcv_info(pp$gam_fit)
    sp[[s]] <- pp$gam_fit$sp
  }
  ident <- function(x) {
    x
  }
  mgcv_info <- sapply(mgcv_info, ident)
  mgcv_info <- tibble::as_tibble(t(mgcv_info))
  sp <- t(sapply(sp, ident))
  nams <- colnames(sp)
  sp <- data.frame(sp)
  colnames(sp) <- nams

  # Return
  list(
    sigma_proj = unlist(sigma_proj),
    mu_proj_test = sapply(mu_pp_test, ident),
    mu_proj_train = sapply(mu_pp_train, ident),
    kl_divs = kl_divs,
    weights = sapply(weights, ident),
    terms_train = terms_train,
    terms_test = terms_test,
    mgcv_info = mgcv_info,
    smoothing_params = sp
  )
}

# log likelihood at each point and draw
log_lik_gaussian <- function(y, mu, sigma) {
  S <- length(sigma)
  log_liks <- array(0, dim(mu))
  for (s in seq_len(S)) {
    log_dens <- stats::dnorm(y,
      mean = mu[, s],
      sd = sigma[s],
      log = TRUE
    )
    log_liks[, s] <- log_dens
  }
  log_liks
}

# Mean log predictive density
mlpd_gaussian <- function(y, mu, sigma) {
  log_liks <- log_lik_gaussian(y, mu, sigma)
  mean(colMeans(log_liks))
}

# Project to submodel
project_model <- function(inds, fit, sd, option = "gam") {
  str <- paste(inds, collapse = " ")
  cat("Projecting to submodel [", str, "] \n", sep = "")

  # Project draws
  pd <- project_draws(inds, fit, sd)

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

# Get draws of a projected component
get_projected_component <- function(subp, j, test = FALSE) {
  b <- subp$pd$terms_train
  if (test) {
    b <- subp$pd$terms_test
  }
  getter <- function(x) {
    x[, j]
  }
  sapply(b, getter)
}

# Visualize the projection
plot_projection <- function(inds, subp, sd, draws = TRUE) {
  didx <- NULL
  pcomp_plots <- list()
  NC <- length(inds)
  for (j in seq_len(NC)) {
    vidx <- inds[j]
    x1 <- sd$X_train[vidx, ]
    x2 <- sd$X_test[vidx, ]
    x_j <- c(x1, x2)

    # Additional df for plotting data
    f1 <- rep("train", length(x1))
    f2 <- rep("test", length(x2))
    type <- as.factor(c(f1, f2))
    y <- c(sd$y_train, sd$y_test)
    df_data <- data.frame(x_j, y, type)
    pc_j_train <- get_projected_component(subp, j)
    pc_j_test <- get_projected_component(subp, j, test = TRUE)
    pc <- rbind(pc_j_train, pc_j_test)
    plt_j <- plot_projected_component(x_j, pc, type, draws)
    plt_j <- plt_j + xlab(paste0("x", vidx)) +
      geom_point(
        data = df_data, mapping = aes(x = x_j, y = y, color = type),
        inherit.aes = FALSE
      )
    if (is.null(didx)) {
      plt_j <- plt_j + ylab("Posterior")
    } else {
      yl <- paste0("Projected f[", vidx, "](x", vidx, ")")
      plt_j <- plt_j + ylab(yl)
    }
    pcomp_plots[[j]] <- plt_j
  }
  pcomp_plots
}

# Plot projected component
plot_projected_component <- function(x, pcomp_draws, type, draws) {
  if (!draws) {
    out <- plot_projected_component_dist(x, pcomp_draws, type)
  } else {
    out <- plot_projected_component_draws(x, pcomp_draws, type)
  }
  out
}

# Plot projected component distribution
plot_projected_component_dist <- function(x, pcomp_draws, type) {
  mean <- apply(pcomp_draws, 1, mean)
  std <- apply(pcomp_draws, 1, stats::sd)
  df <- data.frame(x, mean, std, type)
  aesth <- aes(x = x, y = mean, ymin = mean - 2 * std, ymax = mean + 2 * std)
  ggplot(df, aesth) +
    geom_ribbon(fill = "steelblue3", alpha = 0.5) +
    geom_line(color = "steelblue4")
}

# Plot projected component (one draw)
plot_projected_component_draws <- function(x, y, type) {
  yy <- as.numeric(y)
  S <- ncol(y)
  xx <- rep(x, times = S)
  draw_idx <- as.factor(rep(1:S, each = length(x)))
  df <- data.frame(xx, yy, draw_idx)
  aesth <- aes(x = xx, y = yy, group = draw_idx)
  col <- "gray20"
  ggplot(df, aesth) +
    geom_line(color = col, alpha = 0.3)
}
