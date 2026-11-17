# Prepare data for Stan

Converts formatted EHT data to Stan-compatible list.

## Usage

``` r
prepare_stan_data(data, use_likelihood = 1)
```

## Arguments

- data:

  Formatted EHT data frame. Can have uppercase (UA, UD, FA, FD) or
  lowercase (ua, ud, fa, fd) endpoint columns. Treatment must be coded 0
  (control), 1, ..., tr (see
  [`format_eht_data`](https://swisstph.github.io/AnophelesEntoTrials/reference/format_eht_data.md)).

- use_likelihood:

  Integer. 1 for normal fitting (default), 0 for prior predictive
  checks.

## Value

List with elements:

- N: Number of observations

- tr: Number of treatments (excluding control)

- treat: Treatment indicators (integer vector)

- y: N×4 matrix of endpoints (UA, UD, FA, FD)

- use_likelihood: Flag for likelihood computation

## Examples

``` r
if (FALSE) { # \dontrun{
stan_data <- prepare_stan_data(formatted_data)
} # }
```
