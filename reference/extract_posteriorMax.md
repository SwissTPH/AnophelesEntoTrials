# Extract posteriorMax estimates

Convenience function to extract only the posteriorMax estimate from a
fitted model. This is the posterior sample with the highest
log-posterior density.

## Usage

``` r
extract_posteriorMax(stan_fit, data = NULL, pars = NULL)
```

## Arguments

- stan_fit:

  A stanfit object from `fit_ento_model`.

- data:

  Optional. Original data frame for treatment labels.

- pars:

  Character vector. Parameters to extract. If NULL, uses the same
  defaults as
  [`extract_stan_posteriors`](https://swisstph.github.io/AnophelesEntoTrials/reference/extract_stan_posteriors.md).

## Value

A tibble with posteriorMax estimates for each parameter and treatment.

## Examples

``` r
if (FALSE) { # \dontrun{
posteriorMax <- extract_posteriorMax(fit, data = original_data)
} # }
```
