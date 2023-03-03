library(rstan)
library(posterior)
library(ggplot2)
library(mgcv)
library(dimreduce)
source("R/rstan.R")
source("R/search.R")
source("R/project.R")
source("R/project_old.R") # for project_sigma()
source("R/simulate.R")
source("R/plotting.R")
source("R/selection.R")
source("R/hs_smooth.R")
source("R/relevances.R")
options(projpred.extra_verbose = TRUE)

CHAINS <- 1
ITER <- 600
CN <- list(adapt_delta = 0.9)

# Simulation setup
rho <- 0.7
sigma <- 3
rel_true <- c(1, 1, 0, 0, 1, 0, 0, 0)
n <- 600
D <- length(rel_true)
set.seed(7899)

# Run sim and create data setup
dat <- simulate(n, rho, rel_true, sigma)
n <- length(dat$y)
splt <- create_split(n, 0.5)
ds <- list(dat = dat, split = splt)

# Model setup
# Create the Stan model
mod <- rstan::stan_model("stan/model_sample.stan")
ms <- list(
  stan_model = mod,
  B = 24,
  scale_bf = 1.5
)

dat_df <- create_data_frame(ds, test = FALSE)


fmf1 <- fit_full_model(ds, ms,
  chains = CHAINS, iter = ITER, control = CN,
  thin = 10
)

# Use SPCA reference model
# fmf2 <- fit_full_model(ds, ms,
#  chains = CHAINS, iter = ITER, control = CN,
#  thin = 10, nctot = 5
# )

# Rank
cv1 <- comp_vars(fmf1$fit, fmf1$J)
path <- sort(cv1, index.return = T, decreasing = T)$ix

# Run selection with projection predictive method
options(pp.threshold = Inf) # no stopping condition

# Use HS GP basis
options(bs.type = "hs")
options(hs.penalty.power = 4)
sel_pp_hs4 <- selection_pp(fmf1, path)

# Use thin plate spline basis
options(bs.type = "tp")
sel_pp_tp <- selection_pp(fmf1, path)

sels <- list(sel_pp_hs4, sel_pp_tp)
names(sels) <- c("hs4", "tp")
plot_mlpd <- function(sels, test = TRUE) {
  get_mlpd_test <- function(x) {
    x$search$history$mlpd_test
  }
  get_mlpd_train <- function(x) {
    x$search$history$mlpd_train
  }
  if (test) {
    mlpd <- data.frame(sapply(sels, get_mlpd_test))
    ylabb <- "Test MLPD"
  } else {
    mlpd <- data.frame(sapply(sels, get_mlpd_train))
    ylabb <- "Train MLPD"
  }

  mlpd$num_vars <- 0:8
  df <- reshape2::melt(mlpd, id.vars = "num_vars")
  ggplot(df, aes(x = num_vars, y = value, group = variable, color = variable)) +
    geom_line() +
    labs(color = "base") +
    ylab(ylabb) +
    geom_point() +
    theme(legend.position = "top")
}
paths <- sapply(sels, function(x) {
  x$search$path
})
plt1 <- plot_mlpd(sels, FALSE)
plt2 <- plot_mlpd(sels, TRUE)
plt_mlpd <- ggpubr::ggarrange(plt1, plt2)

# Run projection again for problematic model
model_name <- "correct-152"
problem_name <- "hsbad-at-correct"
a <- fmf1
subm <- c(1, 5, 2)
options(bs.type = "hs")
p1 <- project_model(subm, a$fit, a$sd)
options(bs.type = "tp")
p2 <- project_model(subm, a$fit, a$sd)

get_mt <- function(x) {
  x$metrics$mlpd_test
}
projs <- list(hs = p1, tp = p2)
v1 <- sapply(projs, get_mt) # check that these match earlier values
w <- p1$pd$weights
create_plt_z <- function(w, idx) {
  df_w <- data.frame(w = w[, idx], x = rownames(w))
  term <- sapply(strsplit(df_w$x, "[.]"), function(x) {
    x[1]
  })
  df_w$term <- as.factor(term)
  ggplot(df_w, aes(x = x, y = w, color = term)) +
    geom_point()
}
plt_z_a <- create_plt_z(w, 27) + ggtitle("Weights, draw_idx = 27") +
  theme(legend.position = "top")
plt_z_b <- create_plt_z(w, 1) + ggtitle("Weights, draw_idx = 1") +
  theme(legend.position = "top")
plt_z <- ggpubr::ggarrange(plt_z_a, plt_z_b)

# Plots
x_test <- a$ds$dat$x[a$ds$split$test, ]
x_train <- a$ds$dat$x[a$ds$split$train, ]
train_ranges <- apply(x_train, 2, range)
get_mpt <- function(x) {
  p <- x$pd$mu_proj_test
  lpd <- rowMeans(x$all_lpd_test)
  cbind(x_test, lpd, p)
}
mu_proj_test <- lapply(projs, get_mpt)

# Helper function
to_df <- function(mat) {
  a <- data.frame(mat)
  df <- reshape(a,
    direction = "long", varying = 10:ncol(a),
    sep = ""
  )
  colnames(df)[9:12] <- c("mlpd", "draw_idx", "mu", "data_idx")
  df
}
df <- lapply(mu_proj_test, to_df)
nrows <- sapply(df, nrow)
df <- rbind(df[[1]], df[[2]])
b1 <- rep(names(nrows)[1], nrows[1])
b2 <- rep(names(nrows)[2], nrows[2])
df$base <- as.factor(c(b1, b2))
rownames(df) <- NULL

# Plot to identify problematic draw/data point
plt <- ggplot(df, aes(
  x = draw_idx, y = mu, group = base, color = mlpd,
  pch = base
)) +
  geom_point(alpha = 0.75)

idx_max <- which(abs(df$mu) == max(abs(df$mu)))
idx_prob <- df[idx_max, ]$draw_idx
df_10 <- df[df$draw_idx == idx_prob, ]
#

create_4_plots <- function(df, cc, gg) {
  aes1 <- aes_string(
    x = "data_idx", y = "mu", group = gg, color = cc, pch = gg
  )
  aes2 <- aes_string(
    x = "x1", y = "mu", group = gg, color = cc, pch = gg
  )
  aes3 <- aes_string(
    x = "x2", y = "mu", group = gg, color = cc, pch = gg
  )
  aes4 <- aes_string(
    x = "x5", y = "mu", group = gg, color = cc, pch = gg
  )
  plt_10 <- ggplot(df, aes1) +
    geom_point(alpha = 0.75) +
    xlab("test_point_idx")

  plt_10_x1 <- ggplot(df, aes2) +
    geom_vline(xintercept = train_ranges[, 1], lty = 3) +
    geom_point(alpha = 0.75) +
    theme(legend.position = "top")
  plt_10_x2 <- ggplot(df, aes3) +
    geom_vline(xintercept = train_ranges[, 2], lty = 3) +
    geom_point(alpha = 0.75) +
    theme(legend.position = "top")
  plt_10_x5 <- ggplot(df, aes4) +
    geom_vline(xintercept = train_ranges[, 3], lty = 3) +
    geom_point(alpha = 0.75) +
    theme(legend.position = "top")

  pp <- ggpubr::ggarrange(plt_10_x1, plt_10_x2, plt_10_x5, nrow = 1, ncol = 3)
  ggpubr::ggarrange(plt_10, pp,
    nrow = 2, ncol = 1
  )
}

plt_a <- create_4_plots(df_10, "base", "base")
df_10_hs <- df_10[which(df_10$base == "hs"), ]
plt_b <- create_4_plots(df_10_hs, "mlpd", NULL)
plt_c <- ggplot(df_10_hs, aes(x = x1, y = x2, color = mlpd)) +
  geom_point() +
  geom_vline(xintercept = train_ranges[, 1], lty = 3) +
  geom_hline(yintercept = train_ranges[, 2], lty = 3) +
  ggtitle("MLPD (using HS basis) at test points",
    subtitle = "Dotted lines = training data range"
  )
fn1 <- paste0("problem-debug/base_", model_name, "_", problem_name, ".pdf")
fn2 <- paste0("problem-debug/mlpd_", model_name, "_", problem_name, ".pdf")
fn3 <- paste0("problem-debug/mlpd_2d", model_name, "_", problem_name, ".pdf")
fn4 <- paste0("problem-debug/weights", model_name, "_", problem_name, ".pdf")
fn5 <- paste0("problem-debug/path", model_name, "_", problem_name, ".pdf")
ggsave(plt_a, file = fn1, width = 7.6, height = 5.1)
ggsave(plt_b, file = fn2, width = 7.6, height = 5.1)
ggsave(plt_c, file = fn3, width = 7.6, height = 5.1)
ggsave(plt_z, file = fn4, width = 7.6, height = 5.1)
ggsave(plt_mlpd, file = fn5, width = 7.6, height = 5.1)
