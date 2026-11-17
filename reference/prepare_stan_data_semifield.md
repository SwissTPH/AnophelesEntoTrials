# Prepare data for semifield Stan model

Converts semifield experiment data to Stan-compatible list for the
matched semifield model for a spatial repellent (Denz et al. 2021,
[doi:10.1186/s13071-020-04560-x](https://doi.org/10.1186/s13071-020-04560-x)
).

## Usage

``` r
prepare_stan_data_semifield(data, prior_settings = NULL)
```

## Arguments

- data:

  List with the outcome matrices of matched experiments (one row per
  night, same night order in both arms):

  - y0: Control arm

  - y1: Intervention arm

  Columns are H1, H2, H3, H4 (human landing catches in hours 1, 2, 3 and
  4-12), T (trap) and L (lost). Without a trap, the T column can be left
  out (5 columns) and is then set to 0. Optionally, beta0 and beta1:
  trap availability per night in each arm (default 0, i.e. no trap, for
  a spatial repellent only). If the compartments contained a trap (e.g.
  unbaited), pass its availability, as in Denz et al. (2021), rather
  than leaving it at 0.

- prior_settings:

  List. Optional prior hyperparameters (defaults as in Denz et al.
  2021):

  - priorsigma_mean_logrates: Prior sigma for log rates (default: 6)

  - priorsigma_pikappa: Prior sigma for effect params (default: 5)

  - hierarchy: Variance for hierarchical effects (default: 1)

## Value

List with elements for semifield Stan model.
