library(R6)

# Base class
ForwardSearch <- R6Class(
  classname = "ForwardSearch",
  public = list(
    num_comps = NULL,
    verbosity = 0,
    score_name = "unknown",
    option = NULL,
    mlpd = list(),
    path = NULL, # possible predefined path

    # Init
    initialize = function(J, option, path = NULL, verbosity = 1) {
      self$num_comps <- J
      self$option <- option
      self$verbosity <- verbosity
      self$path <- path
    },

    # Score a model
    score = function(model, ...) {
      cat("* Inheriting class should override the score method!\n")
      return(runif(1))
    },

    # Compute score for each candidate
    step = function(model, candidates, ...) {
      scores <- rep(0, length(candidates))
      jj <- 0
      for (ca in candidates) {
        jj <- jj + 1
        c_model <- c(model, ca)
        scores[jj] <- self$score(c_model, ...)
      }
      scores
    },

    # Get candidate models at current step
    get_candidates = function(cur_model, J) {
      j <- length(cur_model) + 1
      if (!is.null(self$path)) {
        message("Using pre-defined search path")
        cands <- self$path[j] # use pre-defined path
      } else {
        message("Possibly many alternative models at this step")
        cands <- setdiff(seq_len(J), cur_model)
      }
    },

    # Perform search
    run = function(...) {
      J <- self$num_comps
      model <- NULL
      score <- NULL
      while (length(model) < J) {
        cands <- self$get_candidates(model, J)
        cand_scores <- self$step(model, cands, ...)
        i_best <- which(cand_scores == max(cand_scores))
        score <- c(score, cand_scores[i_best])
        model <- c(model, cands[i_best])
        if (self$verbosity > 0) {
          print(model)
        }
      }
      list(path = model, score = score)
    }
  )
)

# Forward search using projection
ProjectionForwardSearch <- R6Class(
  classname = "ProjectionForwardSearch",
  inherit = ForwardSearch,
  public = list(
    score_name = "negative_KL",

    # Score a model
    score = function(model, fit, sd) {
      res <- project_model(model, fit, sd, self$option)$metrics
      res$score <- -res$kl
      res
    },

    # Compute score for each candidate
    step = function(model, candidates, ...) {
      df <- NULL
      names <- NULL
      for (ca in candidates) {
        c_model <- c(model, ca)
        res <- self$score(c_model, ...)
        df <- rbind(df, res)
        names <- c(names, paste(c_model, collapse = ""))
      }
      rownames(df) <- names
      data.frame(df)
    },

    # Perform search
    run = function(...) {
      cat("Running a ProjectionForwardSearch with option =", self$option, "\n")
      J <- self$num_comps
      model <- NULL
      score <- NULL
      mlpd <- NULL
      history <- self$score(model, ...)
      history$p_explained <- 0.0
      kl0 <- history$kl
      rownames(history) <- "empty"
      thresh <- getOption("pp.threshold", default = 0.95)

      # Loop
      j <- 0
      while (length(model) < J) {
        j <- j + 1
        msg <- paste0("Step ", j, "/", J, ".")
        message(msg)
        cands <- self$get_candidates(model, J)
        step_res <- self$step(model, cands, ...)
        cand_scores <- step_res$score
        i_best <- which(cand_scores == max(cand_scores))
        new_row <- step_res[i_best, ]
        p_exp <- 1.0 - new_row$kl / kl0
        message(paste0(" * p_exp = ", p_exp))
        new_row$p_explained <- p_exp
        history <- rbind(history, new_row)
        model <- c(model, cands[i_best])
        if (self$verbosity > 0) {
          print(model)
        }
        if (p_exp >= thresh) {
          break
        }
      }

      # Return
      list(path = model, history = history)
    }
  )
)


# Projection predictive forward search
pp_forward_search <- function(fit, sd, option, path) {
  a <- ProjectionForwardSearch$new(sd$J, option, path)
  a$run(fit, sd)
}

# Plot result of above
plot_pp_forward_search <- function(res, explained = TRUE, thresh = 0.95) {
  J <- length(res$path)
  num_vars <- c(0, 1:J)
  kl_div <- res$history$kl
  p_exp <- res$history$p_explained
  df <- data.frame(num_vars, kl_div, p_exp)
  if (explained) {
    out <- ggplot(df, aes(x = num_vars, y = p_exp)) +
      ylab("Explained variation (%)") +
      geom_hline(yintercept = thresh, color = "firebrick3", lty = 2)
  } else {
    out <- ggplot(df, aes(x = num_vars, y = kl_div)) +
      ylab("KL divergence")
  }
  out <- out + geom_line() + geom_point() +
    xlab("Number of variables/terms") + ggtitle("Forward search")
  return(out)
}

# Plot result of above
plot_pp_forward_search_mlpd <- function(res, test = TRUE) {
  J <- length(res$path)
  num_vars <- c(0, 1:J)
  if (test) {
    ynam <- "Test MLPD"
    mlpd <- res$history$mlpd_test
  } else {
    ynam <- "Train MLPD"
    mlpd <- res$history$mlpd_train
  }

  df <- data.frame(num_vars, mlpd)
  out <- ggplot(df, aes(x = num_vars, y = mlpd)) +
    ylab(ynam)
  out <- out + geom_line() + geom_point() +
    xlab("Number of variables/terms") + ggtitle("Forward search")
  return(out)
}
