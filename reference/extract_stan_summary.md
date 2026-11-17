# Extract summary statistics from Stan fit

Extracts parameter summaries (mean, credible intervals) from a fitted
Stan model. Joins with treatment labels if data with
vector_control_product is provided.

## Usage

``` r
extract_stan_summary(
  stan_fit,
  data = NULL,
  pars = NULL,
  probs = c(0.025, 0.975),
  digits = 3
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

- probs:

  Numeric vector. Quantiles for credible intervals. Default is c(0.025,
  0.975) for 95% CI.

- digits:

  Integer. Number of decimal places to round to.

## Value

A tibble with columns:

- param: Parameter name

- treatment: Treatment index (1, 2, ...; NA for scalar parameters such
  as alpha_0, mu_0 and pC_0)

- mean: Posterior mean

- lower: Lower credible interval bound

- upper: Upper credible interval bound

- vector_control_product: Treatment label (if data provided)

## Examples

``` r
if (FALSE) { # \dontrun{
# Extract summaries
summaries <- extract_stan_summary(fit)

# With treatment labels
summaries <- extract_stan_summary(fit, data = original_data)
} # }
```
