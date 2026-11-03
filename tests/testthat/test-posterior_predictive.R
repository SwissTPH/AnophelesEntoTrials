# Tests for posterior predictive check functions

# Observed counts: 2 observations per arm (control 0, treatment 1)
ppc_test_y <- matrix(c(6, 1, 2, 1,
                       4, 2, 3, 1,
                       2, 5, 1, 2,
                       3, 4, 2, 1), ncol = 4, byrow = TRUE)
ppc_test_treat <- c(0, 0, 1, 1)

test_that("ppc_summarise returns observed and replicated proportions per arm", {
  # Two draws: one identical to the data, one with all mosquitoes unfed alive
  y_rep <- array(0, dim = c(2, 4, 4))
  y_rep[1, , ] <- ppc_test_y
  y_rep[2, , 1] <- rowSums(ppc_test_y)

  result <- ppc_summarise(ppc_test_y, y_rep, ppc_test_treat)

  expect_named(result, c("summary", "replicates"))
  expect_equal(nrow(result$summary), 8)          # 2 arms x 4 outcomes
  expect_equal(nrow(result$replicates), 16)      # 2 draws x 2 arms x 4 outcomes
  expect_equal(levels(result$summary$outcome), c("UA", "UD", "FA", "FD"))

  control <- result$summary[result$summary$treatment == 0, ]
  expect_equal(control$observed, c(10, 3, 5, 2) / 20)
  expect_equal(sum(control$observed), 1)

  # UA: replicates 0.5 (identical) and 1 (all unfed alive), both >= observed
  ua_control <- control[control$outcome == "UA", ]
  expect_equal(ua_control$median, 0.75)
  expect_equal(ua_control$p_value, 1)
  # UD: replicates 0.15 and 0, only one >= observed
  expect_equal(control$p_value[control$outcome == "UD"], 0.5)
})

test_that("ppc_summarise adds product names when given", {
  y_rep <- array(ppc_test_y, dim = c(1, 4, 4))
  y_rep[1, , ] <- ppc_test_y
  labels <- data.frame(treatment = c(0, 1),
                       vector_control_product = c("Control", "NewITN"))

  result <- ppc_summarise(ppc_test_y, y_rep, ppc_test_treat, labels = labels)

  expect_equal(unique(result$summary$vector_control_product), c("Control", "NewITN"))
  expect_equal(unique(result$replicates$vector_control_product), c("Control", "NewITN"))
})

test_that("plot_posterior_predictive returns a ggplot", {
  skip_if_not_installed("ggplot2")
  y_rep <- array(0, dim = c(1, 4, 4))
  y_rep[1, , ] <- ppc_test_y

  p <- plot_posterior_predictive(ppc_summarise(ppc_test_y, y_rep, ppc_test_treat))

  expect_s3_class(p, "ggplot")
})

test_that("ppc_summarise works with other outcomes (semifield)", {
  # 1 night per arm, outcomes H1, H2, H3, H4, L
  y <- matrix(c(4, 3, 2, 6, 5,
                1, 1, 1, 2, 15), ncol = 5, byrow = TRUE)
  y_rep <- array(0, dim = c(1, 2, 5))
  y_rep[1, , ] <- y
  outcomes <- c("H1", "H2", "H3", "H4", "L")

  result <- ppc_summarise(y, y_rep, treat = c(0, 1), outcomes = outcomes)

  expect_equal(nrow(result$summary), 10)         # 2 arms x 5 outcomes
  expect_equal(levels(result$summary$outcome), outcomes)
  expect_equal(result$summary$observed[result$summary$treatment == 1],
               c(1, 1, 1, 2, 15) / 20)
})

test_that("save_y_rep is not available for the hierarchical decay model", {
  decay_data <- data.frame(
    UA = rep(10, 4), UD = rep(2, 4), FA = rep(5, 4), FD = rep(1, 4),
    treatment = c(0, 0, 1, 1), day = c(1, 2, 1, 2), datej = c(1, 2, 1, 2),
    fed = rep(6, 4), total = rep(18, 4)
  )

  expect_error(
    fit_ento_model(decay_data, model_type = "hierarchical_decay", save_y_rep = TRUE,
                   verbose = FALSE),
    "only available for the basic and semifield models"
  )
})
