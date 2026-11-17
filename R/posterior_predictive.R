# Posterior Predictive Check Module for AnophelesEntoTrials
#
# Functions to compare the observed data with data replicated from the
# fitted model (y_rep), to check that the model reproduces the data.

#' Posterior predictive check
#'
#' Compares the observed outcome proportions of each arm with those of data
#' replicated from the fitted model. The replicated data are simulated in Stan
#' from each posterior draw, with the same number of mosquitoes per
#' observation as in the data. Available for the basic and semifield models.
#'
#' @param stan_fit A stanfit object from \code{\link{fit_ento_model}}, fitted
#'   with \code{save_y_rep = TRUE}.
#' @param data The data used to fit the model, same rows in the same order:
#'   formatted EHT data for the basic model, or the list with y0 and y1 for
#'   the semifield model (arms labelled "Control" and "Intervention"; the
#'   trap outcome T is left out when there is no trap).
#' @param n_draws Integer. Number of posterior draws to use (random subset).
#' @param probs Numeric vector. Quantiles of the replicated proportions.
#'   Default is c(0.025, 0.975) for a 95% interval.
#'
#' @return A list with:
#'   \itemize{
#'     \item summary: A tibble with one row per arm and outcome: treatment,
#'       outcome (UA, UD, FA, FD for the basic model; H1, H2, H3, H4, T, L for
#'       the semifield model), observed proportion, median, lower and
#'       upper quantiles of the replicated proportions, and p_value, the share
#'       of replicates greater than or equal to the observed proportion
#'       (values close to 0 or 1 indicate that the model does not reproduce
#'       the data)
#'     \item replicates: A tibble with the replicated proportion of each
#'       draw, arm and outcome
#'   }
#'   Both include vector_control_product if the data have this column.
#'
#' @examples
#' \dontrun{
#' fit <- fit_ento_model(eht_data, model_type = "basic", save_y_rep = TRUE)
#' ppc <- posterior_predictive_check(fit, eht_data)
#' ppc$summary
#' plot_posterior_predictive(ppc)
#' }
#'
#' @seealso \code{\link{plot_posterior_predictive}}
#'
#' @export
posterior_predictive_check <- function(stan_fit,
                                       data,
                                       n_draws = 200,
                                       probs = c(0.025, 0.975)) {

  # pars_oi: quantities saved in the fit
  saved <- stan_fit@sim$pars_oi
  semifield <- all(c("y0_rep", "y1_rep") %in% saved)
  if (!semifield && !"y_rep" %in% saved) {
    stop("The fit has no replicated data (y_rep). ",
         "Fit the model with fit_ento_model(..., save_y_rep = TRUE).",
         call. = FALSE)
  }

  mismatch <- "`data` does not match the data used to fit the model (different number of observations)."

  if (semifield) {
    # Both arms stacked: control nights, then intervention nights
    stan_data <- prepare_stan_data_semifield(data)
    y0_rep <- rstan::extract(stan_fit, pars = "y0_rep")$y0_rep  # draws x n x 6
    y1_rep <- rstan::extract(stan_fit, pars = "y1_rep")$y1_rep
    if (!identical(dim(y0_rep)[2:3], dim(stan_data$y0))) {
      stop(mismatch, call. = FALSE)
    }

    n <- stan_data$n
    y <- rbind(stan_data$y0, stan_data$y1)
    y_rep <- array(0L, dim = c(dim(y0_rep)[1], 2 * n, stan_data$K))
    y_rep[, seq_len(n), ] <- y0_rep
    y_rep[, n + seq_len(n), ] <- y1_rep
    treat <- rep(c(0, 1), each = n)
    outcomes <- c("H1", "H2", "H3", "H4", "T", "L")
    labels <- data.frame(treatment = c(0, 1),
                         vector_control_product = c("Control", "Intervention"))

    # Without a trap, T is always 0
    if (all(stan_data$beta0 == 0) && all(stan_data$beta1 == 0)) {
      keep <- outcomes != "T"
      y <- y[, keep, drop = FALSE]
      y_rep <- y_rep[, , keep, drop = FALSE]
      outcomes <- outcomes[keep]
    }
  } else {
    stan_data <- prepare_stan_data(data)
    y <- stan_data$y
    y_rep <- rstan::extract(stan_fit, pars = "y_rep")$y_rep  # draws x N x 4
    if (!identical(dim(y_rep)[2:3], dim(y))) {
      stop(mismatch, call. = FALSE)
    }
    treat <- stan_data$treat
    outcomes <- c("UA", "UD", "FA", "FD")

    # Product names per treatment, if available
    labels <- NULL
    if ("vector_control_product" %in% names(data)) {
      labels <- unique(data.frame(treatment = treat,
                                  vector_control_product = data$vector_control_product))
    }
  }

  # Random subset of draws
  n_total <- dim(y_rep)[1]
  draws <- if (n_draws < n_total) sample(n_total, n_draws) else seq_len(n_total)

  ppc_summarise(y, y_rep[draws, , , drop = FALSE], treat, probs = probs,
                labels = labels, outcomes = outcomes)
}

#' Summarise observed and replicated outcome proportions
#'
#' Internal helper of \code{\link{posterior_predictive_check}}.
#'
#' @param y Matrix of observed counts (observations x outcomes).
#' @param y_rep Array of replicated counts (draws x observations x outcomes).
#' @param treat Treatment code of each observation.
#' @param probs Numeric vector. Quantiles of the replicated proportions.
#' @param labels Optional data frame with treatment and vector_control_product.
#' @param outcomes Character vector. Names of the outcomes (columns of y).
#'
#' @return A list with the summary and replicates tibbles.
#'
#' @keywords internal
ppc_summarise <- function(y, y_rep, treat, probs = c(0.025, 0.975), labels = NULL,
                          outcomes = c("UA", "UD", "FA", "FD")) {

  # Proportion of each outcome per arm
  arm_props <- function(counts) {
    sums <- rowsum(matrix(counts, ncol = length(outcomes),
                          dimnames = list(NULL, outcomes)), treat)
    sums / rowSums(sums)
  }

  observed <- arm_props(y)                                    # arms x outcomes
  replicated <- vapply(seq_len(dim(y_rep)[1]),
                       function(d) arm_props(y_rep[d, , ]),
                       observed)                              # arms x outcomes x draws

  arms <- as.integer(rownames(observed))
  n_arms <- length(arms)

  summary <- data.frame(
    treatment = rep(arms, times = length(outcomes)),
    outcome = factor(rep(outcomes, each = n_arms), levels = outcomes),
    observed = as.vector(observed),
    median = as.vector(apply(replicated, c(1, 2), stats::median)),
    lower = as.vector(apply(replicated, c(1, 2), stats::quantile, probs = probs[1])),
    upper = as.vector(apply(replicated, c(1, 2), stats::quantile, probs = probs[2])),
    p_value = as.vector(apply(sweep(replicated, c(1, 2), observed, ">="), c(1, 2), mean))
  )

  replicates <- data.frame(
    draw = rep(seq_len(dim(replicated)[3]), each = n_arms * length(outcomes)),
    treatment = rep(arms, times = length(outcomes) * dim(replicated)[3]),
    outcome = factor(rep(rep(outcomes, each = n_arms), times = dim(replicated)[3]),
                     levels = outcomes),
    proportion = as.vector(replicated)
  )

  if (!is.null(labels)) {
    summary$vector_control_product <-
      labels$vector_control_product[match(summary$treatment, labels$treatment)]
    replicates$vector_control_product <-
      labels$vector_control_product[match(replicates$treatment, labels$treatment)]
  }

  list(summary = tibble::as_tibble(summary),
       replicates = tibble::as_tibble(replicates))
}

#' Plot a posterior predictive check
#'
#' Plots the replicated outcome proportions (histograms) with the observed
#' proportion (red line), for each arm and outcome. If the model reproduces
#' the data, the red lines fall within the histograms. Requires the ggplot2
#' package.
#'
#' @param ppc Output of \code{\link{posterior_predictive_check}}.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' ppc <- posterior_predictive_check(fit, eht_data)
#' plot_posterior_predictive(ppc)
#' }
#'
#' @seealso \code{\link{posterior_predictive_check}}
#'
#' @export
plot_posterior_predictive <- function(ppc) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is needed for this plot. ",
         "Install it with install.packages(\"ggplot2\").", call. = FALSE)
  }

  # Arms labelled by product name if available
  arm_label <- function(d) {
    if ("vector_control_product" %in% names(d)) {
      d$vector_control_product
    } else {
      paste("Treatment", d$treatment)
    }
  }
  replicates <- ppc$replicates
  replicates$arm <- arm_label(replicates)
  observed <- ppc$summary
  observed$arm <- arm_label(observed)

  ggplot2::ggplot(replicates, ggplot2::aes(x = .data$proportion)) +
    ggplot2::geom_histogram(bins = 30, fill = "grey70") +
    ggplot2::geom_vline(data = observed,
                        ggplot2::aes(xintercept = .data$observed),
                        colour = "red", inherit.aes = FALSE) +
    ggplot2::facet_grid(arm ~ outcome, scales = "free_x") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Posterior Predictive Check",
                  x = "Proportion of mosquitoes", y = "Replicates")
}
