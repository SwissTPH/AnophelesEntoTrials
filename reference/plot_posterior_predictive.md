# Plot a posterior predictive check

Plots the replicated outcome proportions (histograms) with the observed
proportion (red line), for each arm and outcome. If the model reproduces
the data, the red lines fall within the histograms. Requires the ggplot2
package.

## Usage

``` r
plot_posterior_predictive(ppc)
```

## Arguments

- ppc:

  Output of
  [`posterior_predictive_check`](https://swisstph.github.io/AnophelesEntoTrials/reference/posterior_predictive_check.md).

## Value

A ggplot object.

## See also

[`posterior_predictive_check`](https://swisstph.github.io/AnophelesEntoTrials/reference/posterior_predictive_check.md)

## Examples

``` r
if (FALSE) { # \dontrun{
ppc <- posterior_predictive_check(fit, eht_data)
plot_posterior_predictive(ppc)
} # }
```
