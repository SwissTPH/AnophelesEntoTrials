# Fit Bayesian Ento-Trials model

High-level wrapper function to fit entomological hut trial data using
Stan. Estimates intervention effects on mosquito behavior (Post and
Pre-prandial killing efficacy, reduction in host availability, mortality
during host seeking). Supports predefined models and custom Stan models.

## Usage

``` r
fit_ento_model(
  data,
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
  save_y_rep = FALSE,
  seed = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- data:

  Data frame. Formatted data with endpoints and treatment column.
  Required columns depend on model_type (see
  [`list_available_models`](https://swisstph.github.io/AnophelesEntoTrials/reference/list_available_models.md)).

- model_type:

  Character. Predefined model type: "basic" (default) or "semifield".
  Use
  [`list_available_models()`](https://swisstph.github.io/AnophelesEntoTrials/reference/list_available_models.md)
  to see available options. Ignored if custom_model is provided.

- custom_model:

  Character or NULL. Path to custom Stan model file. If provided,
  overrides model_type.

- stan_data:

  List or NULL. Pre-prepared Stan data list. If NULL, data is prepared
  automatically based on model_type.

- iter:

  Integer. Total iterations per chain (default: 6000).

- warmup:

  Integer. Warmup/burn-in iterations (default: 3000).

- chains:

  Integer. Number of MCMC chains (default: 3).

- control:

  List. Stan control parameters. See
  [`rstan::stan`](https://mc-stan.org/rstan/reference/stan.html).

- init:

  Character, function, or list. Initialization method. "auto" generates
  model-appropriate inits. Can also pass a function or list for custom
  initialization.

- pars:

  Character vector. Parameters to save. If NULL, uses model defaults, or
  saves all quantities for a custom model.

- prior_settings:

  List. Model-specific prior settings (for supported models).

- output_path:

  Character. Path to save fitted model (.rds).

- save_y_rep:

  Logical. If TRUE, also save the replicated data `y_rep` used by
  [`posterior_predictive_check`](https://swisstph.github.io/AnophelesEntoTrials/reference/posterior_predictive_check.md)
  (default: FALSE, as it can make the fit large).

- seed:

  Integer. Random seed for reproducibility.

- verbose:

  Logical. Print fitting progress? (default: TRUE)

- ...:

  Additional arguments passed to
  [`rstan::stan`](https://mc-stan.org/rstan/reference/stan.html).

## Value

A `stanfit` object with MCMC samples. Parameters depend on model:

**Basic model:**

- PostprandialKilling: Post-feeding killing

- PreprandialKilling: Pre-feeding killing

- HostAvailabilityReduction: Reduction in host availability

- KillingDuringHostSeeking: Mortality during approach, relative to the
  attack rate alpha_0

- KillingDuringHostSeekingAbsolute: The same as an absolute additional
  mortality rate (alpha_0 \* KillingDuringHostSeeking)

- alpha_0, mu_0, pC_0: Baseline parameters (attack rate, mortality,
  blood-feeding survival)

## Details

### Predefined Models

Use `model_type` to select from predefined models:

- `"basic"`: Standard multinomial model for EHT data

- `"semifield"`: For matched semifield experiments with a spatial
  repellent (advanced)

### Custom Stan Models

For custom models, provide path via `custom_model` parameter. You must
also provide `stan_data` as a properly formatted list, or ensure your
data has columns matching your model's requirements.

### Custom Parameters

Use `pars` to specify which parameters to extract. For custom models,
this should match your model's parameter block.

## See also

[`list_available_models`](https://swisstph.github.io/AnophelesEntoTrials/reference/list_available_models.md),
[`prepare_stan_data`](https://swisstph.github.io/AnophelesEntoTrials/reference/prepare_stan_data.md),
[`extract_stan_summary`](https://swisstph.github.io/AnophelesEntoTrials/reference/extract_stan_summary.md),
[`extract_stan_posteriors`](https://swisstph.github.io/AnophelesEntoTrials/reference/extract_stan_posteriors.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage with default model
fit <- fit_ento_model(data = my_data)

# With custom Stan model and custom parameters
fit_custom <- fit_ento_model(
  data = my_data,
  custom_model = "path/to/my_model.stan",
  stan_data = my_custom_stan_list,
  pars = c("my_param1", "my_param2"),
  init = my_init_function
)

# Adjust MCMC settings
fit <- fit_ento_model(
  data = data,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)
} # }
```
