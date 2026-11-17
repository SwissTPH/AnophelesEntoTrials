# Generate initial values for basic Stan model

Creates random initial values for MCMC chains for the basic model. This
function returns an init function that Stan will call once per chain.

## Usage

``` r
generate_inits_basic(nb_treat)
```

## Arguments

- nb_treat:

  Number of treatments (including control). The basic model estimates
  treatment-level parameters as vectors of length `nb_treat`.

## Value

A function that generates a named list of initial values when called.
The returned function takes no arguments and produces:

- `PostprandialKilling`: Vector of length nb_treat, initialized to
  random values in (0.1, 0.9)

- `KillingDuringHostSeeking`: Vector of length nb_treat, initialized to
  random values in (0.1, 0.9)

- `HostAvailabilityReduction`: Vector of length nb_treat, initialized to
  random values in (0.1, 0.9)

- `a`: Scalar, log-scale baseline attack rate, initialized near 0

- `m`: Scalar, log-scale baseline mortality rate, initialized near 0

- `b`: Scalar, logit-scale baseline blood-feeding survival, initialized
  near 0

## Details

### How Stan Initialization Works

When you call
[`fit_ento_model()`](https://swisstph.github.io/AnophelesEntoTrials/reference/fit_ento_model.md)
with `init = "auto"`, this function is used to create an init function.
Stan calls this function once per chain to get starting values for the
MCMC sampler.

**Why random initialization matters:**

- Different chains start from different points in parameter space

- This helps diagnose convergence: if chains converge to similar values
  from different starting points, we have more confidence in the results

- Avoids getting stuck in local modes

**What happens if parameters are not initialized:**

- Stan uses its default initialization: uniform random on (-2, 2) for
  unconstrained parameters, then transforms to constrained space

- This usually works but can be inefficient for complex models

### Customizing Initialization

You can provide your own init function to
[`fit_ento_model()`](https://swisstph.github.io/AnophelesEntoTrials/reference/fit_ento_model.md):


    my_inits <- function() {
      list(
        PostprandialKilling = rep(0.5, nb_treat),
        KillingDuringHostSeeking = rep(0.3, nb_treat),
        HostAvailabilityReduction = rep(0.4, nb_treat),
        a = 0, m = 0, b = 0
      )
    }
    fit <- fit_ento_model(data, init = my_inits)

## See also

[`generate_inits_semifield`](https://swisstph.github.io/AnophelesEntoTrials/reference/generate_inits_semifield.md)
for semifield model,
[`fit_ento_model`](https://swisstph.github.io/AnophelesEntoTrials/reference/fit_ento_model.md)
for the main fitting function

## Examples

``` r
if (FALSE) { # \dontrun{
# Create init function for 2 treatments (control + 1 intervention)
init_fn <- generate_inits_basic(nb_treat = 2)

# See what one set of inits looks like
init_fn()

# Use in model fitting
fit <- fit_ento_model(data, init = init_fn)
} # }
```
