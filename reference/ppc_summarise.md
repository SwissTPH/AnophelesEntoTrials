# Summarise observed and replicated outcome proportions

Internal helper of
[`posterior_predictive_check`](https://swisstph.github.io/AnophelesEntoTrials/reference/posterior_predictive_check.md).

## Usage

``` r
ppc_summarise(
  y,
  y_rep,
  treat,
  probs = c(0.025, 0.975),
  labels = NULL,
  outcomes = c("UA", "UD", "FA", "FD")
)
```

## Arguments

- y:

  Matrix of observed counts (observations x outcomes).

- y_rep:

  Array of replicated counts (draws x observations x outcomes).

- treat:

  Treatment code of each observation.

- probs:

  Numeric vector. Quantiles of the replicated proportions.

- labels:

  Optional data frame with treatment and vector_control_product.

- outcomes:

  Character vector. Names of the outcomes (columns of y).

## Value

A list with the summary and replicates tibbles.
