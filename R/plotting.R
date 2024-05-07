# Plot
plot_comp <- function(fit, sd, idx, test = FALSE) {
  fq <- get_f_comp_quantiles(fit, idx, test = test)
  if (test) {
    x <- sd$X_test[idx, ]
    ttl <- "Test data"
  } else {
    x <- sd$X_train[idx, ]
    ttl <- "Train data"
  }
  df <- t(fq)
  colnames(df) <- c("lower", "median", "upper")
  df <- data.frame(x, df)
  ggplot(df, aes(x = x, y = median, ymin = lower, ymax = upper)) +
    geom_line() +
    geom_ribbon(col = "firebrick", fill = "firebrick", alpha = 0.5) +
    ylab(paste0("f", idx)) +
    xlab(paste0("x", idx)) +
    ggtitle(ttl)
}

# Plot
plot_comps <- function(fit, sd, ylims = c(-12, 12), test = FALSE, ...) {
  plts <- list()
  num_comps <- sd$J
  for (j in 1:num_comps) {
    plts[[j]] <- plot_comp(fit, sd, j, test = test) + ylim(ylims)
  }
  ggpubr::ggarrange(plotlist = plts, ...)
}

# Plot search path KL
plot_se_path <- function(se, log = FALSE) {
  if (!log) {
    plot(se$kl, type = "o", pch = 16, xlab = "#Terms", ylab = "KL")
  } else {
    # Plot on log scale
    plot(log(se$kl), type = "o", pch = 16, xlab = "#Terms", ylab = "log KL")
  }
  grid()
}

# Plot perc explained for each repetition
plot_pexp_paths <- function(res_pp, thresh = 0.95, alpha = 0.5) {
  L <- length(res_pp)
  pe <- c()
  for (j in seq_len(L)) {
    pe_j <- res_pp[[j]]$history$p_explained
    pe <- c(pe, pe_j)
  }
  J <- length(pe_j) - 1
  nv <- rep(0:J, times = L)
  idx <- rep(1:L, each = J + 1)
  df <- data.frame(as.factor(idx), nv, pe)
  colnames(df) <- c("rep_idx", "num_var", "p_exp")
  ggplot(df, aes(x = num_var, y = p_exp, group = rep_idx)) +
    geom_line(alpha = alpha) +
    ylab("Explained variation (%)") +
    geom_hline(yintercept = thresh, color = "firebrick3", lty = 2) +
    xlab("Number of variables/terms") +
    ggtitle("Forward search") +
    geom_point(alpha = alpha)
}
