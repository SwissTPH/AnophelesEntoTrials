# Tests for model fitting functions

# Helper function to create test data
create_test_model_data <- function(n = 20) {
  set.seed(123)
  data.frame(
    UA = sample(10:50, n, replace = TRUE),
    UD = sample(0:10, n, replace = TRUE),
    FA = sample(5:30, n, replace = TRUE),
    FD = sample(0:5, n, replace = TRUE),
    treatment = rep(c(0, 1), each = n/2),
    vector_control_product = rep(c("Control", "Treatment_A"), each = n/2),
    stringsAsFactors = FALSE
  )
}

# ============================================================================
# Tests for prepare_stan_data()
# ============================================================================

test_that("prepare_stan_data creates valid Stan data list", {
  data <- create_test_model_data()
  result <- prepare_stan_data(data)

  expect_type(result, "list")
  expect_true("N" %in% names(result))
  expect_true("tr" %in% names(result))
  expect_true("treat" %in% names(result))
  expect_true("y" %in% names(result))
  expect_true("use_likelihood" %in% names(result))
})

test_that("prepare_stan_data sets correct N value", {
  data <- create_test_model_data(n = 20)
  result <- prepare_stan_data(data)

  expect_equal(result$N, 20)
})

test_that("prepare_stan_data creates correct y matrix dimensions", {
  data <- create_test_model_data(n = 20)
  result <- prepare_stan_data(data)

  expect_true(is.matrix(result$y))
  expect_equal(nrow(result$y), 20)
  expect_equal(ncol(result$y), 4)  # UA, UD, FA, FD
})

test_that("prepare_stan_data computes correct tr value", {
  data <- create_test_model_data()  # 2 treatments (control + 1)
  result <- prepare_stan_data(data)

  expect_equal(result$tr, 1)  # tr = n_treatments - 1 (excluding control)
})

test_that("prepare_stan_data handles multiple treatments", {
  data <- data.frame(
    UA = rep(10, 30),
    UD = rep(2, 30),
    FA = rep(5, 30),
    FD = rep(1, 30),
    treatment = rep(c(0, 1, 2), each = 10)
  )

  result <- prepare_stan_data(data)

  expect_equal(result$tr, 2)  # 3 treatments - 1 = 2
})

test_that("prepare_stan_data handles lowercase column names", {
  data <- data.frame(
    ua = rep(10, 10),
    ud = rep(2, 10),
    fa = rep(5, 10),
    fd = rep(1, 10),
    treatment = rep(c(0, 1), each = 5)
  )

  result <- prepare_stan_data(data)

  expect_equal(result$N, 10)
  expect_equal(ncol(result$y), 4)
})

test_that("prepare_stan_data errors on missing columns", {
  data <- data.frame(
    UA = rep(10, 10),
    UD = rep(2, 10),
    treatment = rep(c(0, 1), each = 5)
  )

  expect_error(prepare_stan_data(data), "Missing required endpoint columns")
})

test_that("prepare_stan_data errors on missing treatment", {
  data <- data.frame(
    UA = rep(10, 10),
    UD = rep(2, 10),
    FA = rep(5, 10),
    FD = rep(1, 10)
  )

  expect_error(prepare_stan_data(data), "Missing required column: treatment")
})

test_that("prepare_stan_data sets use_likelihood correctly", {
  data <- create_test_model_data()

  result_likelihood <- prepare_stan_data(data, use_likelihood = 1)
  expect_equal(result_likelihood$use_likelihood, 1L)

  result_prior <- prepare_stan_data(data, use_likelihood = 0)
  expect_equal(result_prior$use_likelihood, 0L)
})

test_that("prepare_stan_data preserves data order in y matrix", {
  data <- data.frame(
    UA = c(10, 20, 30),
    UD = c(1, 2, 3),
    FA = c(5, 10, 15),
    FD = c(0, 1, 2),
    treatment = c(0, 1, 1)
  )

  result <- prepare_stan_data(data)

  expect_equal(result$y[1, 1], 10)  # First row, UA
  expect_equal(result$y[2, 1], 20)  # Second row, UA
  expect_equal(result$y[3, 4], 2)   # Third row, FD
})

# ============================================================================
# Tests for generate_inits()
# ============================================================================

test_that("generate_inits returns a function", {
  init_fn <- generate_inits(2)
  expect_type(init_fn, "closure")
})

test_that("generate_inits function returns correct structure", {
  init_fn <- generate_inits(3)
  result <- init_fn()

  expect_type(result, "list")
  expect_true("InitialPostprandialkillingEfficacy" %in% names(result))
  expect_true("KillingDuringHostSeeking" %in% names(result))
  expect_true("InitialRepellencyRate" %in% names(result))
})

test_that("generate_inits creates correct length vectors", {
  nb_treat <- 3
  init_fn <- generate_inits(nb_treat)
  result <- init_fn()

  expect_length(result$InitialPostprandialkillingEfficacy, nb_treat)
  expect_length(result$KillingDuringHostSeeking, nb_treat)
  expect_length(result$InitialRepellencyRate, nb_treat)
})

test_that("generate_inits returns values in valid range", {
  init_fn <- generate_inits(5)
  result <- init_fn()

  # All efficacy/rate parameters should be between 0 and 1
  expect_true(all(result$InitialPostprandialkillingEfficacy >= 0))
  expect_true(all(result$InitialPostprandialkillingEfficacy <= 1))
  expect_true(all(result$KillingDuringHostSeeking >= 0))
  expect_true(all(result$KillingDuringHostSeeking <= 1))
  expect_true(all(result$InitialRepellencyRate >= 0))
  expect_true(all(result$InitialRepellencyRate <= 1))
})

# ============================================================================
# Tests for get_default_parameters()
# ============================================================================

test_that("get_default_parameters returns expected parameters", {
  result <- get_default_parameters()

  expect_type(result, "character")
  expect_true("InitialPostprandialkillingEfficacy" %in% result)
  expect_true("KillingDuringHostSeeking" %in% result)
  expect_true("InitialRepellencyRate" %in% result)
  expect_true("InitialPreprandialkillingEfficacy" %in% result)
  expect_true("alpha_0" %in% result)
  expect_true("mu_0" %in% result)
  expect_true("pc_0" %in% result)
})

# ============================================================================
# Tests for fit_ento_model() - basic checks without Stan
# ============================================================================

test_that("fit_ento_model errors on non-existent custom model", {
  data <- create_test_model_data()

  expect_error(
    fit_ento_model(data, custom_model = "/nonexistent/model.stan"),
    "not found"
  )
})

test_that("fit_ento_model errors on non-existent stan_model path", {
  data <- create_test_model_data()

  expect_error(
    fit_ento_model(data, stan_model = "/nonexistent/model.stan"),
    "not found"
  )
})

# ============================================================================
# Integration tests (skipped if Stan not available)
# ============================================================================

# These tests require Stan to be installed and compiled
# They are skipped in most CI environments

test_that("full workflow with simulated data works", {
  skip_on_cran()
  skip_if_not_installed("rstan")

  # Check if Stan model exists
  stan_file <- system.file("stan", "mosq_foraging_mcmc_full_fit_priorpred.stan",
                           package = "AnophelesEntoTrials")
  skip_if(stan_file == "", "Default Stan model not found")

  # Create simulated data similar to reference code
  set.seed(42)

  # True parameters
  mu_0 <- 0.14
  alpha_0 <- 0.15
  pc_0 <- 0.92
  pi <- 0.69     # repellency
  kappa <- 3.8   # pre-prandial killing
  xi <- 0.52     # post-prandial killing

  # Probability function
  return_proba <- function(alpha, mu, pc) {
    c(
      exp(-alpha - mu),                           # UA
      (1 - exp(-alpha - mu)) * mu / (alpha + mu), # UD
      (1 - exp(-alpha - mu)) * alpha / (alpha + mu) * pc,      # FA
      (1 - exp(-alpha - mu)) * alpha / (alpha + mu) * (1 - pc) # FD
    )
  }

  proba_control <- return_proba(alpha_0, mu_0, pc_0)
  proba_i <- return_proba(alpha_0 * (1 - pi), mu_0 + alpha_0 * kappa, pc_0 * (1 - xi))

  n_samples <- 50  # Small for testing

  sample_control <- t(rmultinom(n_samples, 10, proba_control))
  sample_i <- t(rmultinom(n_samples, 10, proba_i))

  simulated_data <- rbind(
    data.frame(sample_control, treatment = 0, vector_control_product = "control"),
    data.frame(sample_i, treatment = 1, vector_control_product = "newITN")
  )
  names(simulated_data)[1:4] <- c("UA", "UD", "FA", "FD")
  simulated_data$fed <- simulated_data$FA + simulated_data$FD
  simulated_data$total <- rowSums(simulated_data[, 1:4])

  # Prepare Stan data
  stan_data <- prepare_stan_data(simulated_data)

  expect_equal(stan_data$N, 100)
  expect_equal(stan_data$tr, 1)
  expect_equal(ncol(stan_data$y), 4)
})
