# Tests for data simulation functions

# ============================================================================
# Tests for eht_outcome_probs()
# ============================================================================

test_that("eht_outcome_probs returns named probabilities summing to 1", {
  probs <- eht_outcome_probs(alpha = 0.15, mu = 0.14, pC = 0.92)

  expect_named(probs, c("UA", "UD", "FA", "FD"))
  expect_true(all(probs >= 0 & probs <= 1))
  expect_equal(sum(probs), 1)
})

test_that("eht_outcome_probs matches the Markov chain formulas", {
  alpha <- 0.2
  mu <- 0.1
  pC <- 0.8
  probs <- eht_outcome_probs(alpha, mu, pC)

  p_ua <- exp(-(alpha + mu))
  expect_equal(unname(probs["UA"]), p_ua)
  expect_equal(unname(probs["UD"]), (1 - p_ua) * mu / (alpha + mu))
  expect_equal(unname(probs["FA"] + probs["FD"]), (1 - p_ua) * alpha / (alpha + mu))
  expect_equal(unname(probs["FA"] / (probs["FA"] + probs["FD"])), pC)
})

# ============================================================================
# Tests for simulate_eht_data()
# ============================================================================

test_that("simulate_eht_data returns the expected structure", {
  sim <- simulate_eht_data(n_obs = 20, n_mosquitoes = 10, seed = 1)

  expect_s3_class(sim, "tbl_df")
  expect_named(sim, c("UA", "UD", "FA", "FD", "treatment", "vector_control_product"))
  expect_equal(nrow(sim), 40)  # control + 1 treatment arm
  expect_identical(sort(unique(sim$treatment)), c(0L, 1L))
  expect_equal(unique(sim$vector_control_product), c("Control", "NewITN"))
  expect_true(all(rowSums(sim[, c("UA", "UD", "FA", "FD")]) == 10))
})

test_that("simulate_eht_data handles several treatment arms", {
  sim <- simulate_eht_data(
    n_obs = 5,
    treatments = data.frame(
      vector_control_product = c("ITN_A", "ITN_B"),
      HostAvailabilityReduction = c(0.4, 0.69),
      KillingDuringHostSeeking = c(1.5, 3.8),
      PostprandialKilling = c(0.3, 0.52)
    ),
    seed = 1
  )

  expect_equal(nrow(sim), 15)
  expect_identical(unique(sim$treatment), c(0L, 1L, 2L))
  expect_equal(unique(sim$vector_control_product), c("Control", "ITN_A", "ITN_B"))
})

test_that("simulate_eht_data is reproducible with a seed", {
  expect_identical(simulate_eht_data(n_obs = 10, seed = 42),
                   simulate_eht_data(n_obs = 10, seed = 42))
})

test_that("simulated proportions are close to the model probabilities", {
  alpha_0 <- 0.15
  mu_0 <- 0.14
  pC_0 <- 0.92
  sim <- simulate_eht_data(n_obs = 5000, n_mosquitoes = 10,
                           alpha_0 = alpha_0, mu_0 = mu_0, pC_0 = pC_0,
                           seed = 123)
  endpoints <- c("UA", "UD", "FA", "FD")

  control <- sim[sim$treatment == 0, endpoints]
  expect_equal(colSums(control) / sum(control),
               eht_outcome_probs(alpha_0, mu_0, pC_0),
               tolerance = 0.02)

  treated <- sim[sim$treatment == 1, endpoints]
  expect_equal(colSums(treated) / sum(treated),
               eht_outcome_probs((1 - 0.69) * alpha_0,
                                 mu_0 + 3.8 * alpha_0,
                                 pC_0 * (1 - 0.52)),
               tolerance = 0.02)
})

test_that("simulate_eht_data output can be formatted and prepared for Stan", {
  sim <- simulate_eht_data(n_obs = 10, seed = 1)
  formatted <- format_eht_data(sim)
  stan_data <- prepare_stan_data(formatted)

  expect_equal(stan_data$N, 20)
  expect_equal(stan_data$tr, 1)
})

test_that("simulate_eht_data validates its inputs", {
  expect_error(simulate_eht_data(n_obs = 0), "n_obs")
  expect_error(simulate_eht_data(n_mosquitoes = 2.5), "n_mosquitoes")
  expect_error(simulate_eht_data(alpha_0 = -1), "alpha_0")
  expect_error(simulate_eht_data(pC_0 = 1.5), "pC_0")
  expect_error(
    simulate_eht_data(treatments = data.frame(vector_control_product = "X")),
    "missing columns"
  )
  expect_error(
    simulate_eht_data(treatments = data.frame(
      vector_control_product = "X",
      HostAvailabilityReduction = 1.2,
      KillingDuringHostSeeking = 1,
      PostprandialKilling = 0.5
    )),
    "HostAvailabilityReduction"
  )
})
