# Fit submodel
fit_sub_model <- function(model, ds, ms, ...) {
  model_str <- paste(model, collapse = ", ")
  message("Fitting submodel with variables {", model_str, "}.")
  sd <- create_stan_input(ds$dat, ds$split, ms$scale_bf, ms$B, x_cols = model)
  rstan::sampling(ms$stan_model, data = sd, refresh = 0, ...)
}

# Fit full model
fit_full_model <- function(ds, ms, nctot = NULL, ...) {
  D <- dim(ds$dat$x)[2]
  if (is.null(nctot)) {
    nctot <- D
    dat <- ds$dat
    dr <- NULL
  } else {
    message("Running supervised PCA.")
    dd <- create_data_frame(ds, test = FALSE)
    dr <- dimreduce::spca(x = dd[, 1:D], y = dd$y, nctot = D)
    new_x <- predict(dr, ds$dat$x)
    new_x <- new_x[, 1:nctot]
    dat <- list(y = ds$dat$y, x = new_x)
  }

  message("Creating input.")
  sd_spca <- create_stan_input(dat, ds$split, ms$scale_bf, ms$B)
  JB <- sd_spca$J * sd_spca$B
  N <- sd_spca$N_train
  if (JB > N) {
    stop("More total basis functions than data points!")
  }
  sd <- create_stan_input(dat, ds$split, ms$scale_bf, ms$B)

  # Sample full model
  message("Fitting full model.")
  fit <- rstan::sampling(ms$stan_model, data = sd_spca, refresh = 0, ...)
  message("Sampling done.")
  plt_comps <- plot_comps(fit, sd_spca, ncol = 4, nrow = 2)


  # Return
  list(
    sd = sd,
    sd_spca = sd_spca,
    fit = fit,
    plt_comps = plt_comps,
    ds = ds,
    ms = ms,
    dr = dr,
    J = nctot
  )
}


# Selection using projection predictive forward search
selection_pp <- function(fmf, path = NULL) {
  t_start <- Sys.time()

  # Forward search
  message("Performing forward search.")
  se <- pp_forward_search(fmf$fit, fmf$sd, option = "gam", path = path)

  # Return
  message("Search done.")
  t_end <- Sys.time()
  list(
    search = se,
    plt_kl = plot_pp_forward_search(se),
    plt_mlpd_test = plot_pp_forward_search_mlpd(se),
    plt_mlpd_train = plot_pp_forward_search_mlpd(se, F),
    time = t_end - t_start
  )
}

# Number of selected when using fv
compute_num_sel <- function(fv, thresh) {

}

# Selection using component variances
selection_fv <- function(ds, ms, full_fit, thresh = 0.95, median = FALSE, ...) {
  message("Ordering the variables.")
  # aa_est <- get_alpha_estimates(full_fit, median = median)
  J <- ncol(ds$dat$x)
  fv <- comp_vars(full_fit, J)
  order <- sort(fv, index.return = TRUE, decreasing = TRUE)$ix
  p_exp <- fv / sum(fv)
  r <- sort(p_exp, decreasing = TRUE, index.return = TRUE)
  num_sel <- length(which(cumsum(r$x) < thresh)) + 1

  # Return
  list(
    path = order,
    fv = fv,
    p_exp = p_exp,
    num_sel = num_sel,
    selected_model = order[1:num_sel]
  )
}
