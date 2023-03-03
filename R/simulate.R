library(MASS)

# Covariates
simulate_x <- function(n, rho, D = 4) {
  Sigma <- rho * matrix(1, D, D) + (1 - rho) * diag(nrow = D, ncol = D)
  mu <- rep(0, nrow(Sigma))
  x1 <- MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma)
  x2 <- MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma)
  cbind(x1, x2)
}

# Function components
simulate_f <- function(x) {
  f <- 0 * x
  for (j in 1:ncol(f)) {
    a <- 0.2 + 1.2 * runif(1)
    b <- pi * runif(1)
    f_j <- sin(a * x[, j] + b)
    f_j <- f_j - mean(f_j)
    f_j <- f_j / stats::sd(f_j)
    f[, j] <- f_j
  }
  f
}

# Multiply columns by relevance
decide_relevant <- function(f, rel) {
  N <- dim(f)[1]
  D <- dim(f)[2]
  f * matrix(rel, nrow = N, ncol = D, byrow = TRUE)
}

# Entire simulation
simulate <- function(n, rho, rel_true, sigma) {
  D <- length(rel_true) / 2
  checkmate::assert_integerish(D)
  x <- simulate_x(n, rho, D)
  # D <- length(rel_true)
  # x <- matrix(runif(n = n * D), n, D)
  colnames(x) <- paste0("x", seq_len(ncol(x)))
  f <- simulate_f(x)
  f <- decide_relevant(f, rel_true)
  f_tot <- rowSums(f)
  y <- f_tot + sigma * rnorm(n = length(f_tot))
  list(
    x = x,
    f = f,
    f_tot = f_tot,
    y = y
  )
}

# Create data split
create_split <- function(n, p_test) {
  checkmate::assert_integerish(n, lower = 1)
  checkmate::assert_number(p_test, lower = 0, upper = 1)
  p_train <- 1 - p_test
  n_train <- round(n * p_train)
  idx_train <- sort(sample(n, n_train))
  list(
    train = idx_train,
    test = sort(setdiff(seq_len(n), idx_train))
  )
}

# Take a subset of data
subset_data <- function(sim, inds) {
  list(
    x = sim$x[inds, ],
    f = sim$f[inds, ],
    f_tot = sim$f_tot[inds],
    y = sim$y[inds]
  )
}

# Create data frame for other methods
create_data_frame <- function(ds, test = FALSE) {
  df <- cbind(ds$dat$x, ds$dat$y)
  df <- as.data.frame(df)
  idx <- ds$split$train
  if (test) {
    idx <- ds$split$test
  }
  df <- df[idx, ]
  D <- dim(df)[2] - 1
  colnames(df) <- c(paste("x", seq_len(D), sep = ""), "y")
  df
}
