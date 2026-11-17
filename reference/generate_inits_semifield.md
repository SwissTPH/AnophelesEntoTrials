# Generate initial values for semifield model

Creates random initial values for MCMC chains for the semifield model.
The semifield model uses continuous-time Markov chains for time-to-event
data.

## Usage

``` r
generate_inits_semifield(n)
```

## Arguments

- n:

  Number of experiments (replicates). Used to size hierarchical random
  effects.

## Value

A function that generates a named list of initial values:

- `HostAvailabilityReduction`: Reduction in host availability (pi in the
  article)

- `KillingDuringHostSeeking`: Killing/disarming during host seeking
  (kappa in the article)

- `a`, `m`: Log-scale baseline rates

- `sigma_a`, `sigma_m`: Hierarchical SDs

- `phi`, `psi`: Experiment-level random effects (n)

## Details

The semifield model estimates, per night:

- `alpha_H`: Rate of host-seeking success

- `mu`: Mortality rate

The intervention effects are modeled as
`alpha_H1 = (1 - HostAvailabilityReduction) * alpha_H0` and
`mu1 = mu0 + KillingDuringHostSeeking * alpha_H0`.

## See also

[`generate_inits_basic`](https://swisstph.github.io/AnophelesEntoTrials/reference/generate_inits_basic.md),
[`fit_ento_model`](https://swisstph.github.io/AnophelesEntoTrials/reference/fit_ento_model.md)
