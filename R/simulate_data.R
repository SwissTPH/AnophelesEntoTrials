# Data Simulation Module for AnophelesEntoTrials
#
# Functions to simulate entomological hut trial data with known parameter
# values, e.g. to test the models or check that a fit recovers the truth.

#' Outcome probabilities of the basic EHT model
#'
#' Probabilities of the four hut trial outcomes for one arm, using the same
#' Markov chain formulas as the basic Stan model
#' (\code{mosq_foraging_mcmc_full_fit_priorpred.stan}).
#'
#' @param alpha Numeric. Attack rate of the arm.
#' @param mu Numeric. Mortality rate of the arm.
#' @param pB Numeric. Probability of surviving after blood-feeding.
#'
#' @return Named numeric vector of probabilities for UA, UD, FA and FD.
#'
#' @keywords internal
eht_outcome_probs <- function(alpha, mu, pB) {
  p_ua <- exp(-(alpha + mu))
  p_fed <- (1 - p_ua) * alpha / (alpha + mu)
  p_ud <- (1 - p_ua) * mu / (alpha + mu)

  c(UA = p_ua,
    UD = p_ud,
    FA = p_fed * pB,
    FD = p_fed * (1 - pB))
}

#' Simulate entomological hut trial data from the basic model
#'
#' Simulates hut trial outcomes (UA, UD, FA, FD) for a control arm and one or
#' more treatment arms from the basic model, with known parameter values.
#' Parameters use the same names as the Stan output, so the values used here
#' can be compared directly with \code{\link{extract_stan_summary}}.
#' Only the basic model (\code{model_type = "basic"}) is supported: the data
#' have no time or day structure and cannot be used for the hierarchical
#' decay or semifield models.
#'
#' @param n_obs Integer. Number of observations (hut-nights) per arm.
#' @param n_mosquitoes Integer. Number of mosquitoes per observation.
#' @param alpha_0 Numeric. Baseline attack rate (> 0).
#' @param mu_0 Numeric. Baseline mortality rate (> 0).
#' @param pB_0 Numeric. Baseline probability of surviving after
#'   blood-feeding (between 0 and 1).
#' @param treatments Data frame with one row per treatment arm and columns:
#'   \itemize{
#'     \item vector_control_product: Name of the product
#'     \item InitialRepellencyRate: Repellency (between 0 and 1)
#'     \item KillingDuringHostSeeking: Additional mortality during host
#'       seeking, as a multiple of alpha_0 (>= 0)
#'     \item InitialPostprandialkillingEfficacy: Post-prandial killing
#'       (between 0 and 1)
#'   }
#' @param control_name Character. Product name of the control arm.
#' @param seed Integer or NULL. If given, passed to \code{set.seed()} before
#'   simulating (this sets the global random seed).
#'
#' @return A tibble with one row per observation and columns UA, UD, FA, FD,
#'   treatment (0 = control, 1, 2, ... in the order of \code{treatments}) and
#'   vector_control_product. Pass it to \code{\link{format_eht_data}} to add
#'   the derived columns.
#'
#' @details
#' For a treatment arm, the rates are modified as in the basic Stan model:
#' \itemize{
#'   \item attack rate: \code{(1 - InitialRepellencyRate) * alpha_0}
#'   \item mortality rate: \code{mu_0 + KillingDuringHostSeeking * alpha_0}
#'   \item survival after feeding:
#'     \code{pB_0 * (1 - InitialPostprandialkillingEfficacy)}
#' }
#' The control arm uses \code{alpha_0}, \code{mu_0} and \code{pB_0}. Each
#' observation is a multinomial draw of \code{n_mosquitoes} over the four
#' outcome probabilities.
#'
#' @examples
#' # One treatment arm
#' sim <- simulate_eht_data(seed = 42)
#' head(sim)
#'
#' # Two treatment arms
#' sim <- simulate_eht_data(
#'   treatments = data.frame(
#'     vector_control_product = c("Standard_ITN", "New_ITN"),
#'     InitialRepellencyRate = c(0.4, 0.69),
#'     KillingDuringHostSeeking = c(1.5, 3.8),
#'     InitialPostprandialkillingEfficacy = c(0.3, 0.52)
#'   ),
#'   seed = 42
#' )
#' eht_data <- format_eht_data(sim)
#'
#' @export
simulate_eht_data <- function(n_obs = 100,
                              n_mosquitoes = 10,
                              alpha_0 = 0.15,
                              mu_0 = 0.14,
                              pB_0 = 0.92,
                              treatments = data.frame(
                                vector_control_product = "NewITN",
                                InitialRepellencyRate = 0.69,
                                KillingDuringHostSeeking = 3.8,
                                InitialPostprandialkillingEfficacy = 0.52
                              ),
                              control_name = "Control",
                              seed = NULL) {

  # Validate inputs
  is_count <- function(x) {
    is.numeric(x) && length(x) == 1 && !is.na(x) && x >= 1 && x == round(x)
  }
  if (!is_count(n_obs)) {
    stop("`n_obs` must be a single positive integer.", call. = FALSE)
  }
  if (!is_count(n_mosquitoes)) {
    stop("`n_mosquitoes` must be a single positive integer.", call. = FALSE)
  }
  if (!is.numeric(alpha_0) || length(alpha_0) != 1 || !isTRUE(alpha_0 > 0)) {
    stop("`alpha_0` must be a single number > 0.", call. = FALSE)
  }
  if (!is.numeric(mu_0) || length(mu_0) != 1 || !isTRUE(mu_0 > 0)) {
    stop("`mu_0` must be a single number > 0.", call. = FALSE)
  }
  if (!is.numeric(pB_0) || length(pB_0) != 1 ||
      !isTRUE(pB_0 >= 0 && pB_0 <= 1)) {
    stop("`pB_0` must be a single number between 0 and 1.", call. = FALSE)
  }

  treatments <- as.data.frame(treatments)
  effect_cols <- c("InitialRepellencyRate",
                   "KillingDuringHostSeeking",
                   "InitialPostprandialkillingEfficacy")
  missing_cols <- setdiff(c("vector_control_product", effect_cols), names(treatments))
  if (length(missing_cols) > 0) {
    stop("`treatments` is missing columns: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  if (nrow(treatments) == 0) {
    stop("`treatments` must have at least one row.", call. = FALSE)
  }
  for (col in effect_cols) {
    if (!is.numeric(treatments[[col]]) || anyNA(treatments[[col]])) {
      stop("`treatments$", col, "` must be numeric without missing values.",
           call. = FALSE)
    }
  }
  if (any(treatments$InitialRepellencyRate < 0 | treatments$InitialRepellencyRate > 1)) {
    stop("`InitialRepellencyRate` must be between 0 and 1.", call. = FALSE)
  }
  if (any(treatments$KillingDuringHostSeeking < 0)) {
    stop("`KillingDuringHostSeeking` must be >= 0.", call. = FALSE)
  }
  if (any(treatments$InitialPostprandialkillingEfficacy < 0 |
          treatments$InitialPostprandialkillingEfficacy > 1)) {
    stop("`InitialPostprandialkillingEfficacy` must be between 0 and 1.",
         call. = FALSE)
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # One row per arm, control first with no intervention effect
  arms <- data.frame(
    vector_control_product = c(control_name,
                               as.character(treatments$vector_control_product)),
    InitialRepellencyRate = c(0, treatments$InitialRepellencyRate),
    KillingDuringHostSeeking = c(0, treatments$KillingDuringHostSeeking),
    InitialPostprandialkillingEfficacy = c(0, treatments$InitialPostprandialkillingEfficacy),
    stringsAsFactors = FALSE
  )

  sims <- lapply(seq_len(nrow(arms)), function(j) {
    probs <- eht_outcome_probs(
      alpha = (1 - arms$InitialRepellencyRate[j]) * alpha_0,
      mu = mu_0 + arms$KillingDuringHostSeeking[j] * alpha_0,
      pB = pB_0 * (1 - arms$InitialPostprandialkillingEfficacy[j])
    )
    counts <- t(stats::rmultinom(n_obs, n_mosquitoes, probs))
    colnames(counts) <- names(probs)

    data.frame(
      counts,
      treatment = j - 1L,
      vector_control_product = arms$vector_control_product[j],
      stringsAsFactors = FALSE
    )
  })

  tibble::as_tibble(do.call(rbind, sims))
}
