#' Entomological Trials Model Fitting Module
#'
#' Functions for Bayesian fitting of entomological trials(Ento-Trials) data using Stan.
#' Supports custom and pre-loaded Stan models and custom parameter specification
#'
#' @name fit_ento_model
#' @docType _PACKAGE
NULL

#' List available predefined Stan models
#'
#' Returns information about predefined Stan models included in the package.
#'
#' @param details Logical. If TRUE, returns detailed information about each model.
#'   If FALSE (default), returns just model names.
#'
#' @return If details = FALSE, a character vector of model names.
#'   If details = TRUE, a data frame with model information.
#'
#' @examples
#' # List model names
#' list_available_models()
#'
#' # Get detailed information
#' list_available_models(details = TRUE)
#'
#' @export
list_available_models <- function(details = FALSE) {

  models <- list(
    basic = list(
      name = "basic",
      file = "mosq_foraging_mcmc_full_fit_priorpred.stan",
      description = "Standard multinomial EHT model",
      required_data = c("UA", "UD", "FA", "FD", "treatment"),
      parameters = c("InitialPostprandialkillingEfficacy",
                     "InitialPreprandialkillingEfficacy",
                     "InitialRepellencyRate",
                     "KillingDuringHostSeeking",
                     "alpha_0", "mu_0", "pB_0"),
      use_case = "Standard experimental hut trial data with 4 outcome categories"
    ),
    hierarchical_decay = list(
      name = "hierarchical_decay",
      file = "hierarchial_time_decay.stan",
      description = "Hierarchical model with Weibull time decay",
      required_data = c("UA", "UD", "FA", "FD", "treatment", "day", "datej", "fed", "total"),
      parameters = c("InitialPostprandialkillingEfficacy",
                     "InitialPreprandialkillingEfficacy",
                     "InitialRepellencyRate",
                     "KillingDuringHostSeeking",
                     "beta", "kappa", "L",
                     "alpha_0", "mu_0"),
      use_case = "Longitudinal EHT data with temporal decay of intervention effects"
    ),
    semifield = list(
      name = "semifield",
      file = "13071_2020_4560_MOESM4_ESM.stan",
      description = "Continuous-time Markov chain model for semifield experiments",
      required_data = c("y0", "y1", "n", "K"),
      parameters = c("rho", "alpha_H0", "alpha_T0", "mu0",
                     "alpha_H1", "alpha_T1", "mu1"),
      use_case = "Semifield experiments with time-to-event outcomes (H1-H4, T, L)"
    )
  )

  if (!details) {
    return(names(models))
  }

  # Return as data frame
  model_df <- do.call(rbind, lapply(names(models), function(nm) {
    m <- models[[nm]]
    data.frame(
      model = nm,
      file = m$file,
      description = m$description,
      required_data = paste(m$required_data, collapse = ", "),
      use_case = m$use_case,
      stringsAsFactors = FALSE
    )
  }))

  return(tibble::as_tibble(model_df))
}

#' Get model configuration
#'
#' Returns configuration for a predefined model type.
#'
#' @param model_type Character. One of "basic", "hierarchical_decay", "semifield".
#'
#' @return List with model file path, default parameters, and data preparation function.
#'
#' @keywords internal
get_model_config <- function(model_type) {

  configs <- list(
    basic = list(
      file = "mosq_foraging_mcmc_full_fit_priorpred.stan",
      pars = c("InitialPostprandialkillingEfficacy",
               "KillingDuringHostSeeking",
               "InitialRepellencyRate",
               "InitialPreprandialkillingEfficacy",
               "alpha_0", "mu_0", "pB_0"),
      prepare_data = "prepare_stan_data",
      init_generator = "generate_inits_basic"
    ),
    hierarchical_decay = list(
      file = "hierarchial_time_decay.stan",
      pars = c("InitialPostprandialkillingEfficacy",
               "KillingDuringHostSeeking",
               "InitialRepellencyRate",
               "InitialPreprandialkillingEfficacy",
               "beta", "kappa", "L",
               "alpha_0", "mu_0"),
      prepare_data = "prepare_stan_data_decay",
      init_generator = "generate_inits_decay"
    ),
    semifield = list(
      file = "13071_2020_4560_MOESM4_ESM.stan",
      pars = c("rho", "alpha_H0", "alpha_T0", "mu0",
               "alpha_H1", "alpha_T1", "mu1"),
      prepare_data = "prepare_stan_data_semifield",
      init_generator = "generate_inits_semifield"
    )
  )

  if (!model_type %in% names(configs)) {
    stop("Unknown model type: ", model_type,
         ". Available: ", paste(names(configs), collapse = ", "),
         call. = FALSE)
  }

  return(configs[[model_type]])
}

#' Fit Bayesian Ento-Trials model
#'
#' High-level wrapper function to fit entomological hut trial data using Stan.
#' Estimates intervention effects on mosquito behavior (Post and Pre-prandial
#' killing efficacy, repellency, mortality during host seeking).
#' Supports predefined models and custom Stan models.
#'
#' @param data Data frame. Formatted data with endpoints and treatment column.
#'   Required columns depend on model_type (see \code{\link{list_available_models}}).
#' @param model_type Character. Predefined model type: "basic" (default),
#'   "hierarchical_decay", or "semifield". Use \code{list_available_models()}
#'   to see available options. Ignored if custom_model is provided.
#' @param custom_model Character or NULL. Path to custom Stan model file.
#'   If provided, overrides model_type.
#' @param stan_data List or NULL. Pre-prepared Stan data list. If NULL,
#'   data is prepared automatically based on model_type.
#' @param iter Integer. Total iterations per chain (default: 6000).
#' @param warmup Integer. Warmup/burn-in iterations (default: 3000).
#' @param chains Integer. Number of MCMC chains (default: 3).
#' @param control List. Stan control parameters. See \code{rstan::stan}.
#' @param init Character, function, or list. Initialization method.
#'   "auto" generates model-appropriate inits. Can also pass a function
#'   or list for custom initialization.
#' @param pars Character vector. Parameters to save. If NULL, uses model defaults.
#' @param prior_settings List. Model-specific prior settings (for supported models).
#' @param output_path Character. Path to save fitted model (.rds).
#' @param seed Integer. Random seed for reproducibility.
#' @param verbose Logical. Print fitting progress? (default: TRUE)
#' @param ... Additional arguments passed to \code{rstan::stan}.
#'
#' @return A \code{stanfit} object with MCMC samples. Parameters depend on model:
#'
#' **Basic model:**
#' \itemize{
#'   \item InitialPostprandialkillingEfficacy: Post-feeding killing
#'   \item InitialPreprandialkillingEfficacy: Pre-feeding killing
#'   \item InitialRepellencyRate: Repellency/deterrence
#'   \item KillingDuringHostSeeking: Mortality during approach
#'   \item alpha_0, mu_0, pB_0: Baseline parameters (attack rate, mortality, blood-feeding survival)
#' }
#'
#' **Hierarchical decay model (additional):**
#' \itemize{
#'   \item beta, kappa: Weibull decay parameters
#'   \item L: Half-life of intervention effect
#' }
#'
#' @details
#' ## Predefined Models
#'
#' Use \code{model_type} to select from predefined models:
#' \itemize{
#'   \item \code{"basic"}: Standard multinomial model for EHT data
#'   \item \code{"hierarchical_decay"}: Includes temporal Weibull decay
#'   \item \code{"semifield"}: For semifield experiments (advanced)
#' }
#'
#' ## Custom Stan Models
#'
#' For custom models, provide path via \code{custom_model} parameter.
#' You must also provide \code{stan_data} as a properly formatted list,
#' or ensure your data has columns matching your model's requirements.
#'
#' ## Custom Parameters
#'
#' Use \code{pars} to specify which parameters to extract. For custom models,
#' this should match your model's parameter block.
#'
#' @examples
#' \dontrun{
#' # Basic usage with default model
#' fit <- fit_ento_model(data = my_data)
#'
#' # Use hierarchical decay model for longitudinal data
#' fit_decay <- fit_ento_model(
#'   data = longitudinal_data,
#'   model_type = "hierarchical_decay"
#' )
#'
#' # With custom Stan model and custom parameters
#' fit_custom <- fit_ento_model(
#'   data = my_data,
#'   custom_model = "path/to/my_model.stan",
#'   stan_data = my_custom_stan_list,
#'   pars = c("my_param1", "my_param2"),
#'   init = my_init_function
#' )
#'
#' # Adjust MCMC settings
#' fit <- fit_ento_model(
#'   data = data,
#'   control = list(adapt_delta = 0.95, max_treedepth = 12)
#' )
#' }
#'
#' @seealso
#' \code{\link{list_available_models}}, \code{\link{prepare_stan_data}},
#' \code{\link{extract_stan_summary}}, \code{\link{extract_stan_posteriors}}
#'
#' @export
fit_ento_model <- function(data,
                          model_type = "basic",
                          custom_model = NULL,
                          stan_data = NULL,
                          iter = 6000,
                          warmup = 3000,
                          chains = 3,
                          control = list(adapt_delta = 0.8),
                          init = "auto",
                          pars = NULL,
                          prior_settings = NULL,
                          output_path = NULL,
                          seed = NULL,
                          verbose = TRUE,
                          ...) {

  # Determine model configuration
  if (!is.null(custom_model)) {
    # Custom model path provided
    if (!file.exists(custom_model)) {
      stop("Custom Stan model file not found: ", custom_model, call. = FALSE)
    }
    stan_file <- custom_model
    model_config <- NULL

    if (verbose) {
      cli::cli_alert_info("Using custom Stan model: {.file {basename(custom_model)}}")
    }
  } else {
    # Use predefined model
    model_config <- get_model_config(model_type)
    stan_file <- system.file("stan", model_config$file,
                             package = "AnophelesEntoTrials")

    if (stan_file == "") {
      stop("Model file not found: ", model_config$file,
           ". Package may not be installed correctly.", call. = FALSE)
    }

    if (verbose) {
      cli::cli_alert_info("Using predefined model: {model_type}")
    }
  }

 # Prepare Stan data if not provided
  if (is.null(stan_data)) {
    if (!is.null(model_config)) {
      # Use model-specific data preparation
      stan_data <- switch(model_type,
        "basic" = prepare_stan_data(data),
        "hierarchical_decay" = prepare_stan_data_decay(data, prior_settings),
        "semifield" = prepare_stan_data_semifield(data, prior_settings),
        prepare_stan_data(data)  # fallback
      )
    } else {
      # Custom model without stan_data - try basic preparation
      stan_data <- prepare_stan_data(data)
      if (verbose) {
        cli::cli_alert_warning(
          "Using basic data preparation for custom model. ",
          "Consider providing stan_data directly."
        )
      }
    }
  }

  # Get parameters to save
  if (is.null(pars)) {
    if (!is.null(model_config)) {
      pars <- model_config$pars
    } else {
      pars <- get_default_parameters()
    }
  }

  # Generate initial values
  if (identical(init, "auto")) {
    if (!is.null(model_config)) {
      init_fun <- switch(model_type,
        "basic" = generate_inits_basic(stan_data$tr + 1),
        "hierarchical_decay" = generate_inits_decay(
          nb_treat = stan_data$tr + 1,
          nb_days = stan_data$d
        ),
        "semifield" = generate_inits_semifield(stan_data$n),
        generate_inits_basic(stan_data$tr + 1)  # fallback
      )
    } else {
      # Custom model - use basic inits or let Stan decide
      if ("tr" %in% names(stan_data)) {
        init_fun <- generate_inits_basic(stan_data$tr + 1)
      } else {
        init_fun <- "random"
      }
    }
  } else if (is.function(init)) {
    init_fun <- init
  } else {
    init_fun <- init
  }

  # Fit model
  if (verbose) {
    cli::cli_alert_info("Fitting Stan model...")
    cli::cli_alert_info("Iterations: {iter} | Warmup: {warmup} | Chains: {chains}")
  }

  stan_fit <- tryCatch(
    {
      rstan::stan(
        file = stan_file,
        data = stan_data,
        pars = pars,
        include = TRUE,
        iter = iter,
        init = init_fun,
        warmup = warmup,
        chains = chains,
        control = control,
        seed = seed,
        verbose = FALSE,
        ...
      )
    },
    error = function(e) {
      stop("Stan model fitting failed: ", e$message, call. = FALSE)
    }
  )

  # Check for warnings
  check_stan_warnings(stan_fit, verbose = verbose)

  # Save if requested
  if (!is.null(output_path)) {
    output_dir <- dirname(output_path)
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    saveRDS(stan_fit, output_path)
    if (verbose) {
      cli::cli_alert_success("Model saved to: {.file {output_path}}")
    }
  }

  if (verbose) {
    cli::cli_alert_success("Model fitting complete")
  }

  return(stan_fit)
}

#' Prepare data for Stan
#'
#' Converts formatted EHT data to Stan-compatible list.
#'
#' @param data Formatted EHT data frame. Can have uppercase (UA, UD, FA, FD)
#'   or lowercase (ua, ud, fa, fd) endpoint columns.
#' @param use_likelihood Integer. 1 for normal fitting (default),
#'   0 for prior predictive checks.
#'
#' @return List with elements:
#'   \itemize{
#'     \item N: Number of observations
#'     \item tr: Number of treatments (excluding control)
#'     \item treat: Treatment indicators (integer vector)
#'     \item y: N×4 matrix of endpoints (UA, UD, FA, FD)
#'     \item use_likelihood: Flag for likelihood computation
#'   }
#'
#' @examples
#' \dontrun{
#' stan_data <- prepare_stan_data(formatted_data)
#' }
#'
#' @export
prepare_stan_data <- function(data, use_likelihood = 1) {

  # Get column names (handle both upper and lowercase)
  col_names <- names(data)
  col_names_lower <- tolower(col_names)

  # Find endpoint columns
  ua_col <- col_names[col_names_lower == "ua"]
  ud_col <- col_names[col_names_lower == "ud"]
  fa_col <- col_names[col_names_lower == "fa"]
  fd_col <- col_names[col_names_lower == "fd"]

  # Validate required columns exist
  missing <- c()
  if (length(ua_col) == 0) missing <- c(missing, "UA")
  if (length(ud_col) == 0) missing <- c(missing, "UD")
  if (length(fa_col) == 0) missing <- c(missing, "FA")
  if (length(fd_col) == 0) missing <- c(missing, "FD")

  if (length(missing) > 0) {
    stop("Missing required endpoint columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }

  # Get treatment column
  treat_col <- col_names[col_names_lower == "treatment"]
  if (length(treat_col) == 0) {
    stop("Missing required column: treatment", call. = FALSE)
  }

  # Number of treatments (excluding control)
  nb_treat <- length(unique(data[[treat_col]]))

  # Create outcome matrix
  y_matrix <- matrix(
    c(data[[ua_col]], data[[ud_col]], data[[fa_col]], data[[fd_col]]),
    ncol = 4
  )

  # Ensure treatment is integer
  treatment_vec <- as.integer(data[[treat_col]])

  stan_data <- list(
    N = nrow(data),
    tr = nb_treat - 1,  # Excluding control
    treat = treatment_vec,
    y = y_matrix,
    use_likelihood = as.integer(use_likelihood)
  )

  return(stan_data)
}

#' Prepare data for hierarchical decay Stan model
#'
#' Converts longitudinal EHT data to Stan-compatible list for the
#' hierarchical time decay model.
#'
#' @param data Formatted EHT data frame with required columns:
#'   UA, UD, FA, FD, treatment, day, datej, fed, total.
#' @param prior_settings List. Optional prior hyperparameters:
#'   \itemize{
#'     \item priorsigma_mean_logrates: Prior sigma for log rates (default: 6)
#'     \item hierarchy: Variance for hierarchical effects (default: 1)
#'   }
#'
#' @return List with elements for hierarchical decay Stan model.
#'
#' @export
prepare_stan_data_decay <- function(data, prior_settings = NULL) {

  # Get column names (case-insensitive)
  col_names <- names(data)
  col_names_lower <- tolower(col_names)

  # Find required columns
  get_col <- function(target) {
    col_names[col_names_lower == tolower(target)]
  }

  # Validate required columns
  required <- c("UA", "UD", "FA", "FD", "treatment", "day", "datej", "fed", "total")
  missing <- c()
  for (req in required) {
    if (length(get_col(req)) == 0) missing <- c(missing, req)
  }

  if (length(missing) > 0) {
    stop("Missing required columns for hierarchical_decay model: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }

  # Extract data
  ua_col <- get_col("UA")
  ud_col <- get_col("UD")
  fa_col <- get_col("FA")
  fd_col <- get_col("FD")

  # Number of treatments and days
  nb_treat <- length(unique(data[[get_col("treatment")]])) - 1
  nb_days <- length(unique(data[[get_col("day")]]))

  # Create y matrix (UA, F=FA+FD, UD) - 3 columns for this model
  y_matrix <- matrix(
    c(data[[ua_col]],
      data[[get_col("FA")]] + data[[get_col("FD")]],
      data[[ud_col]]),
    ncol = 3
  )

  # Default prior settings
  defaults <- list(
    priorsigma_mean_logrates = 6,
    hierarchy = 1
  )

  if (!is.null(prior_settings)) {
    defaults[names(prior_settings)] <- prior_settings
  }

  stan_data <- list(
    N = nrow(data),
    tr = nb_treat,
    d = nb_days,
    total = as.integer(data[[get_col("total")]]),
    datej = as.integer(data[[get_col("datej")]]),
    day = as.integer(data[[get_col("day")]]),
    treat = as.integer(data[[get_col("treatment")]]),
    fed = as.integer(data[[get_col("fed")]]),
    FA = as.integer(data[[get_col("FA")]]),
    y = y_matrix,
    priorsigma_mean_logrates = defaults$priorsigma_mean_logrates,
    hierarchy = defaults$hierarchy
  )

  return(stan_data)
}

#' Prepare data for semifield Stan model
#'
#' Converts semifield experiment data to Stan-compatible list.
#'
#' @param data Data frame or list with semifield experiment results.
#'   Should contain y0 (control) and y1 (intervention) outcome matrices,
#'   or columns H1, H2, H3, H4, T, L for each arm.
#' @param prior_settings List. Optional prior hyperparameters:
#'   \itemize{
#'     \item priorsigma_mean_logrates: Prior sigma (default: 2)
#'     \item priorsigma_pikappa: Prior sigma for effect params (default: 1)
#'     \item hierarchy: Variance for hierarchical effects (default: 1)
#'   }
#'
#' @return List with elements for semifield Stan model.
#'
#' @export
prepare_stan_data_semifield <- function(data, prior_settings = NULL) {

  # Check if data is already in matrix format
  if (is.list(data) && all(c("y0", "y1") %in% names(data))) {
    y0 <- as.matrix(data$y0)
    y1 <- as.matrix(data$y1)
  } else {
    stop("Semifield model requires data with y0 and y1 matrices. ",
         "See ?prepare_stan_data_semifield for details.", call. = FALSE)
  }

  n <- nrow(y0)
  K <- ncol(y0)

  # Default prior settings
  defaults <- list(
    priorsigma_mean_logrates = 2,
    priorsigma_pikappa = 1,
    hierarchy = 1
  )

  if (!is.null(prior_settings)) {
    defaults[names(prior_settings)] <- prior_settings
  }

  stan_data <- list(
    n = n,
    K = K,
    y0 = y0,
    y1 = y1,
    priorsigma_mean_logrates = defaults$priorsigma_mean_logrates,
    priorsigma_pikappa = defaults$priorsigma_pikappa,
    hierarchy = defaults$hierarchy
  )

  return(stan_data)
}

#' Generate initial values for basic Stan model
#'
#' Creates random initial values for MCMC chains for the basic model.
#' This function returns an init function that Stan will call once per chain.
#'
#' @param nb_treat Number of treatments (including control). The basic model
#'   estimates treatment-level parameters as vectors of length `nb_treat`.
#'
#' @return A function that generates a named list of initial values when called.
#'   The returned function takes no arguments and produces:
#'   \itemize{
#'     \item \code{InitialPostprandialkillingEfficacy}: Vector of length nb_treat,
#'       initialized to random values in (0.1, 0.9)
#'     \item \code{KillingDuringHostSeeking}: Vector of length nb_treat,
#'       initialized to random values in (0.1, 0.9)
#'     \item \code{InitialRepellencyRate}: Vector of length nb_treat,
#'       initialized to random values in (0.1, 0.9)
#'     \item \code{a}: Scalar, log-scale baseline attack rate, initialized near 0
#'     \item \code{m}: Scalar, log-scale baseline mortality rate, initialized near 0
#'     \item \code{b}: Scalar, logit-scale baseline blood-feeding survival, initialized near 0
#'   }
#'
#' @details
#' ## How Stan Initialization Works
#'
#' When you call `fit_ento_model()` with `init = "auto"`, this function is used
#' to create an init function. Stan calls this function once per chain to get
#' starting values for the MCMC sampler.
#'
#' **Why random initialization matters:**
#' - Different chains start from different points in parameter space

#' - This helps diagnose convergence: if chains converge to similar values
#'   from different starting points, we have more confidence in the results
#' - Avoids getting stuck in local modes
#'
#' **What happens if parameters are not initialized:**
#' - Stan uses its default initialization: uniform random on (-2, 2) for
#'   unconstrained parameters, then transforms to constrained space
#' - This usually works but can be inefficient for complex models
#'
#' ## Customizing Initialization
#'
#' You can provide your own init function to `fit_ento_model()`:
#' \preformatted{
#' my_inits <- function() {
#'   list(
#'     InitialPostprandialkillingEfficacy = rep(0.5, nb_treat),
#'     KillingDuringHostSeeking = rep(0.3, nb_treat),
#'     InitialRepellencyRate = rep(0.4, nb_treat),
#'     a = 0, m = 0, b = 0
#'   )
#' }
#' fit <- fit_ento_model(data, init = my_inits)
#' }
#'
#' @examples
#' \dontrun{
#' # Create init function for 2 treatments (control + 1 intervention)
#' init_fn <- generate_inits_basic(nb_treat = 2)
#'
#' # See what one set of inits looks like
#' init_fn()
#'
#' # Use in model fitting
#' fit <- fit_ento_model(data, init = init_fn)
#' }
#'
#' @seealso
#' \code{\link{generate_inits_decay}} for hierarchical decay model,
#' \code{\link{generate_inits_semifield}} for semifield model,
#' \code{\link{fit_ento_model}} for the main fitting function
#'
#' @export
generate_inits_basic <- function(nb_treat) {
  function() {
    list(
      # Treatment-level efficacy parameters (vectors)
      InitialPostprandialkillingEfficacy = runif(nb_treat, 0.1, 0.9),
      KillingDuringHostSeeking = runif(nb_treat, 0.1, 0.9),
      InitialRepellencyRate = runif(nb_treat, 0.1, 0.9),
      # Scalar baseline parameters (log/logit scale)
      a = rnorm(1, 0, 1),  # log(alpha_0) - baseline attack rate
      m = rnorm(1, 0, 1),  # log(mu_0) - baseline mortality
      b = rnorm(1, 0, 1)   # logit(pB_0) - baseline blood-feeding survival
    )
  }
}

#' Generate initial values for hierarchical decay model
#'
#' Creates random initial values for MCMC chains for the hierarchical decay model.
#' This model includes Weibull decay parameters for temporal effects.
#'
#' @param nb_treat Number of treatments (including control). Parameters are
#'   vectors of length `nb_treat`.
#' @param nb_days Number of days in the study. Used to size hierarchical
#'   random effects. If NULL, defaults to 10.
#'
#' @return A function that generates a named list of initial values:
#'   \itemize{
#'     \item \code{InitialPostprandialkillingEfficacy}: Vector (nb_treat)
#'     \item \code{KillingDuringHostSeeking}: Vector (nb_treat)
#'     \item \code{InitialRepellencyRate}: Vector (nb_treat)
#'     \item \code{beta}: Weibull scale parameter for decay (nb_treat)
#'     \item \code{kappa}: Weibull shape parameter for decay (nb_treat)
#'     \item \code{a}, \code{m}: Log-scale baseline rates
#'     \item \code{r_prob_surviving_feeding_d}: Logit-scale survival
#'     \item \code{sigma_a}, \code{sigma_m}, \code{sigma_prob_surviving_feeding_d}:
#'       Hierarchical standard deviations
#'     \item \code{phi}, \code{psi}, \code{omega_prob_surviving_feeding_d}:
#'       Day-level random effects (nb_days)
#'   }
#'
#' @details
#' The hierarchical decay model adds:
#' - **Weibull decay**: Effect decays as `exp(-(beta * t)^kappa)` where t is
#'   days since treatment. Half-life L = (log(2)^(1/kappa)) / beta.
#' - **Day-level hierarchy**: Random effects phi, psi, omega capture day-to-day
#'   variation in baseline rates.
#'
#' @seealso \code{\link{generate_inits_basic}}, \code{\link{fit_ento_model}}
#'
#' @export
generate_inits_decay <- function(nb_treat, nb_days = NULL) {
  if (is.null(nb_days)) nb_days <- 10  # default
  function() {
    list(
      # Treatment-level efficacy parameters
      InitialPostprandialkillingEfficacy = runif(nb_treat, 0.1, 0.9),
      KillingDuringHostSeeking = runif(nb_treat, 0.1, 0.9),
      InitialRepellencyRate = runif(nb_treat, 0.1, 0.9),
      # Weibull decay parameters
      beta = runif(nb_treat, 0.001, 0.1),
      kappa = runif(nb_treat, 0.5, 2),
      # Baseline rate parameters (log/logit scale)
      a = rnorm(1, 0, 1),
      m = rnorm(1, 0, 1),
      r_prob_surviving_feeding_d = rnorm(1, 0, 1),
      # Hierarchical standard deviations
      sigma_a = abs(rnorm(1, 0, 0.5)),
      sigma_m = abs(rnorm(1, 0, 0.5)),
      sigma_prob_surviving_feeding_d = abs(rnorm(1, 0, 0.5)),
      # Day-level random effects (standard normal)
      phi = rnorm(nb_days, 0, 1),
      psi = rnorm(nb_days, 0, 1),
      omega_prob_surviving_feeding_d = rnorm(nb_days, 0, 1)
    )
  }
}

#' Generate initial values for semifield model
#'
#' Creates random initial values for MCMC chains for the semifield model.
#' The semifield model uses continuous-time Markov chains for time-to-event data.
#'
#' @param n Number of experiments (replicates). Used to size hierarchical
#'   random effects.
#'
#' @return A function that generates a named list of initial values:
#'   \itemize{
#'     \item \code{rho}: Effect multiplier for intervention (log-normal prior)
#'     \item \code{a}, \code{b}, \code{m}: Log-scale baseline rates
#'     \item \code{sigma_a}, \code{sigma_b}, \code{sigma_m}: Hierarchical SDs
#'     \item \code{phi}, \code{eta}, \code{psi}: Experiment-level random effects (n)
#'   }
#'
#' @details
#' The semifield model estimates:
#' - `alpha_H`: Rate of host-seeking success
#' - `alpha_T`: Rate of trap entry
#' - `mu`: Mortality rate
#'
#' The intervention effect is modeled as `alpha_T1 = rho * alpha_H0`.
#'
#' @seealso \code{\link{generate_inits_basic}}, \code{\link{fit_ento_model}}
#'
#' @export
generate_inits_semifield <- function(n) {
  function() {
    list(
      # Effect parameter
      rho = runif(1, 0.5, 2),
      # Log-scale baseline rates
      a = rnorm(1, 0, 0.5),
      b = rnorm(1, 0, 0.5),
      m = rnorm(1, 0, 0.5),
      # Hierarchical standard deviations
      sigma_a = abs(rnorm(1, 0, 0.3)),
      sigma_b = abs(rnorm(1, 0, 0.3)),
      sigma_m = abs(rnorm(1, 0, 0.3)),
      # Experiment-level random effects
      phi = rnorm(n, 0, 1),
      eta = rnorm(n, 0, 1),
      psi = rnorm(n, 0, 1)
    )
  }
}

# Backward compatibility alias
generate_inits <- generate_inits_basic

#' Get default parameters to extract
#'
#' Returns standard parameters for EHT model.
#'
#' @keywords internal
get_default_parameters <- function() {
  c(
    "InitialPostprandialkillingEfficacy",
    "KillingDuringHostSeeking",
    "InitialRepellencyRate",
    "InitialPreprandialkillingEfficacy",
    "alpha_0",
    "mu_0",
    "pB_0"
  )
}

#' Check Stan warnings
#'
#' Checks for convergence issues and prints warnings.
#'
#' @param stan_fit Fitted Stan model
#' @param verbose Print warnings?
#'
#' @keywords internal
check_stan_warnings <- function(stan_fit, verbose = TRUE) {
  # Check for divergences
  sampler_params <- rstan::get_sampler_params(stan_fit, inc_warmup = FALSE)
  n_divergent <- sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))

  # MAYbe too specific
  if (n_divergent > 0 && verbose) {
    cli::cli_alert_warning("{n_divergent} divergent transitions detected")
    cli::cli_alert_info("Consider increasing adapt_delta")
  }

  # Check R-hat
  summary_fit <- rstan::summary(stan_fit)$summary
  max_rhat <- max(summary_fit[, "Rhat"], na.rm = TRUE)

  if (max_rhat > 1.1 && verbose) {
    cli::cli_alert_warning("Maximum Rhat = {round(max_rhat, 3)} > 1.1")
    cli::cli_alert_info("Model may not have converged")
  }

  invisible(TRUE)
}

#' Extract summary statistics from Stan fit
#'
#' Extracts parameter summaries (mean, credible intervals) from a fitted Stan model.
#' Joins with treatment labels if data with vector_control_product is provided.
#'
#' @param stan_fit A stanfit object from \code{fit_ento_model}.
#' @param data Optional. Original data frame with vector_control_product and treatment
#'   columns for labeling.
#' @param pars Character vector. Parameters to extract. If NULL, uses default
#'   intervention parameters.
#' @param probs Numeric vector. Quantiles for credible intervals.
#'   Default is c(0.025, 0.975) for 95% CI.
#' @param digits Integer. Number of decimal places to round to.
#'
#' @return A tibble with columns:
#'   \itemize{
#'     \item param: Parameter name
#'     \item treatment: Treatment index (1, 2, ...)
#'     \item mean: Posterior mean
#'     \item lower: Lower credible interval bound
#'     \item upper: Upper credible interval bound
#'     \item vector_control_product: Treatment label (if data provided)
#'   }
#'
#' @examples
#' \dontrun{
#' # Extract summaries
#' summaries <- extract_stan_summary(fit)
#'
#' # With treatment labels
#' summaries <- extract_stan_summary(fit, data = original_data)
#' }
#'
#' @export
extract_stan_summary <- function(stan_fit,
                                  data = NULL,
                                  pars = NULL,
                                  probs = c(0.025, 0.975),
                                  digits = 3) {

  # Default parameters
  if (is.null(pars)) {
    pars <- c("InitialPostprandialkillingEfficacy",
              "InitialPreprandialkillingEfficacy",
              "InitialRepellencyRate",
              "KillingDuringHostSeeking")
  }

  # Check which parameters exist in the model
  model_pars <- stan_fit@model_pars
  pars_exist <- pars[pars %in% model_pars]

  if (length(pars_exist) == 0) {
    stop("None of the requested parameters found in model. Available: ",
         paste(model_pars, collapse = ", "), call. = FALSE)
  }

  # Get summary
  fit_summary <- rstan::summary(stan_fit, pars = pars_exist, probs = probs)$summary

  # Convert to data frame
  summary_df <- as.data.frame(fit_summary)
  summary_df$index <- rownames(fit_summary)

  # Parse parameter names and indices
  # Pattern: paramname[index] or paramname
  summary_df <- summary_df %>%
    dplyr::mutate(
      param = gsub("\\[.*\\]", "", index),
      treatment = as.integer(gsub(".*\\[([0-9]+)\\].*", "\\1", index))
    ) %>%
    dplyr::mutate(
      treatment = ifelse(grepl("\\[", index), treatment, NA_integer_)
    )

  # Filter to treatment effects only (index > 1, since 1 is control)
  summary_df <- summary_df %>%
    dplyr::filter(is.na(treatment) | treatment > 1) %>%
    dplyr::mutate(treatment = treatment - 1)  # Convert to 1-indexed treatments

  # Select and rename columns
  lower_col <- paste0(probs[1] * 100, "%")
  upper_col <- paste0(probs[2] * 100, "%")

  result <- summary_df %>%
    dplyr::select(
      param,
      treatment,
      mean = mean,
      lower = dplyr::all_of(lower_col),
      upper = dplyr::all_of(upper_col)
    ) %>%
    dplyr::mutate(
      mean = round(mean, digits),
      lower = round(lower, digits),
      upper = round(upper, digits)
    )

  # Join with vector control product names if data provided

  if (!is.null(data) && "vector_control_product" %in% names(data)) {
    treatment_labels <- unique(data[, c("vector_control_product", "treatment")])
    treatment_labels <- treatment_labels[treatment_labels$treatment > 0, ]

    result <- result %>%
      dplyr::left_join(treatment_labels, by = "treatment")
  }

  return(tibble::as_tibble(result))
}

#' Extract posterior samples from Stan fit
#'
#' Extracts posterior samples or MAP (Maximum A Posteriori) estimates
#' from a fitted Stan model.
#'
#' @param stan_fit A stanfit object from \code{fit_ento_model}.
#' @param data Optional. Original data frame with vector_control_product and treatment
#'   columns for labeling.
#' @param pars Character vector. Parameters to extract. If NULL, uses default
#'   intervention parameters.
#' @param n_samples Integer. Number of posterior samples to return.
#'   If NULL, returns all samples.
#' @param map_only Logical. If TRUE, returns only the MAP estimate (sample with
#'   highest log posterior). Default is FALSE.
#'
#' @return A tibble with columns:
#'   \itemize{
#'     \item sample: Sample index
#'     \item param: Parameter name
#'     \item treatment: Treatment index
#'     \item value: Parameter value
#'     \item vector_control_product: Treatment label (if data provided)
#'   }
#'
#' @examples
#' \dontrun{
#' # Extract all posterior samples
#' samples <- extract_stan_posteriors(fit)
#'
#' # Extract MAP estimate only
#' map_estimate <- extract_stan_posteriors(fit, map_only = TRUE)
#'
#' # Extract 100 random samples
#' samples <- extract_stan_posteriors(fit, n_samples = 100)
#' }
#'
#' @export
extract_stan_posteriors <- function(stan_fit,
                                     data = NULL,
                                     pars = NULL,
                                     n_samples = NULL,
                                     map_only = FALSE) {

  # Default parameters (include lp__ for MAP)
  if (is.null(pars)) {
    pars <- c("InitialPostprandialkillingEfficacy",
              "InitialPreprandialkillingEfficacy",
              "InitialRepellencyRate",
              "KillingDuringHostSeeking",
              "alpha_0", "mu_0")
  }

  # Always include lp__ for MAP
  pars_with_lp <- unique(c(pars, "lp__"))

  # Check which parameters exist
  model_pars <- stan_fit@model_pars
  pars_exist <- pars_with_lp[pars_with_lp %in% c(model_pars, "lp__")]

  # Extract posterior samples
  fit_extract <- rstan::extract(stan_fit, pars = pars_exist, inc_warmup = FALSE)

  # Get total number of samples
  total_samples <- length(fit_extract$lp__)

  if (map_only) {
    # Find MAP estimate
    id_MAP <- which.max(fit_extract$lp__)
    sample_ids <- id_MAP
  } else if (!is.null(n_samples) && n_samples < total_samples) {
    # Random subset of samples
    sample_ids <- sample(1:total_samples, n_samples)
  } else {
    sample_ids <- 1:total_samples
  }

  # Build result data frame
  results <- list()

  for (par_name in pars) {
    if (!(par_name %in% names(fit_extract))) next

    par_samples <- fit_extract[[par_name]]

    if (is.matrix(par_samples)) {
      # Parameter with indices (e.g., efficacy[treatment])
      n_indices <- ncol(par_samples)

      for (idx in 1:n_indices) {
        for (s in sample_ids) {
          results[[length(results) + 1]] <- data.frame(
            sample = s,
            param = par_name,
            treatment = idx - 1,  # 0-indexed (0 = control)
            value = par_samples[s, idx],
            stringsAsFactors = FALSE
          )
        }
      }
    } else {
      # Scalar parameter
      for (s in sample_ids) {
        results[[length(results) + 1]] <- data.frame(
          sample = s,
          param = par_name,
          treatment = 0,
          value = par_samples[s],
          stringsAsFactors = FALSE
        )
      }
    }
  }

  result <- dplyr::bind_rows(results)

  # Filter to treatment effects only (exclude control, treatment = 0)
  result <- result %>%
    dplyr::filter(treatment > 0)

  # Join with vector control product names if data provided
  if (!is.null(data) && "vector_control_product" %in% names(data)) {
    treatment_labels <- unique(data[, c("vector_control_product", "treatment")])
    treatment_labels <- treatment_labels[treatment_labels$treatment > 0, ]

    result <- result %>%
      dplyr::left_join(treatment_labels, by = "treatment")
  }

  return(tibble::as_tibble(result))
}

#' Extract posteriorMax estimates
#'
#' Convenience function to extract only the posteriorMax estimate from a fitted model.
#' This is the posterior sample with the highest log-posterior density.
#'
#' @param stan_fit A stanfit object from \code{fit_ento_model}.
#' @param data Optional. Original data frame for treatment labels.
#' @param pars Character vector. Parameters to extract.
#'
#' @return A tibble with posteriorMax estimates for each parameter and treatment.
#'
#' @examples
#' \dontrun{
#' posteriorMax <- extract_posteriorMax(fit, data = original_data)
#' }
#'
#' @export
extract_posteriorMax <- function(stan_fit, data = NULL, pars = NULL) {
  extract_stan_posteriors(stan_fit, data = data, pars = pars, map_only = TRUE)
}
