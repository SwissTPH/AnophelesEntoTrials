# Extract posterior samples from Stan fit

Extracts posterior samples or the posteriorMax estimate from a fitted
Stan model. The posteriorMax estimate is the posterior draw with the
highest log posterior density (`lp__`) among the samples.

## Usage

``` r
extract_stan_posteriors(
  stan_fit,
  data = NULL,
  pars = NULL,
  n_samples = NULL,
  posteriorMax_only = FALSE,
  map_only = NULL
)
```

## Arguments

- stan_fit:

  A stanfit object from `fit_ento_model`.

- data:

  Optional. Original data frame with vector_control_product and
  treatment columns for labeling.

- pars:

  Character vector. Parameters to extract. If NULL, uses the
  intervention parameters and the baseline parameters alpha_0, mu_0 and
  pC_0 (see
  [`get_default_parameters()`](https://swisstph.github.io/AnophelesEntoTrials/reference/get_default_parameters.md)).
  Parameters not in the model are skipped.

- n_samples:

  Integer. Number of posterior samples to return. If NULL, returns all
  samples.

- posteriorMax_only:

  Logical. If TRUE, returns only the posteriorMax estimate (sample with
  highest log posterior). Default is FALSE.

- map_only:

  Deprecated, use `posteriorMax_only`.

## Value

A tibble with columns:

- sample: Sample index

- param: Parameter name

- treatment: Treatment index (NA for scalar parameters such as alpha_0,
  mu_0 and pC_0)

- value: Parameter value

- vector_control_product: Treatment label (if data provided)

## Examples

``` r
if (FALSE) { # \dontrun{
# Extract all posterior samples
samples <- extract_stan_posteriors(fit)

# Extract posteriorMax estimate only
posteriorMax <- extract_stan_posteriors(fit, posteriorMax_only = TRUE)

# Extract 100 random samples
samples <- extract_stan_posteriors(fit, n_samples = 100)
} # }
```
