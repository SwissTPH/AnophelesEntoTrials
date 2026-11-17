# Simulate entomological hut trial data from the basic model

Simulates hut trial outcomes (UA, UD, FA, FD) for a control arm and one
or more treatment arms from the basic model, with known parameter
values. Parameters use the same names as the Stan output, so the values
used here can be compared directly with
[`extract_stan_summary`](https://swisstph.github.io/AnophelesEntoTrials/reference/extract_stan_summary.md).
Only the basic model (`model_type = "basic"`) is supported: the data
have no time or day structure and cannot be used for the semifield
model.

## Usage

``` r
simulate_eht_data(
  n_obs = 100,
  n_mosquitoes = 10,
  alpha_0 = 0.15,
  mu_0 = 0.14,
  pC_0 = 0.92,
  treatments = data.frame(vector_control_product = "NewITN", HostAvailabilityReduction =
    0.69, KillingDuringHostSeeking = 3.8, PostprandialKilling = 0.52),
  control_name = "Control",
  seed = NULL
)
```

## Arguments

- n_obs:

  Integer. Number of observations (hut-nights) per arm.

- n_mosquitoes:

  Integer. Number of mosquitoes per observation.

- alpha_0:

  Numeric. Baseline attack rate (\> 0).

- mu_0:

  Numeric. Baseline mortality rate (\> 0).

- pC_0:

  Numeric. Baseline probability of surviving after blood-feeding
  (between 0 and 1).

- treatments:

  Data frame with one row per treatment arm and columns:

  - vector_control_product: Name of the product

  - HostAvailabilityReduction: Reduction in host availability (between 0
    and 1)

  - KillingDuringHostSeeking: Additional mortality during host seeking,
    as a multiple of alpha_0 (\>= 0)

  - PostprandialKilling: Post-prandial killing (between 0 and 1)

- control_name:

  Character. Product name of the control arm.

- seed:

  Integer or NULL. If given, passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before simulating
  (this sets the global random seed).

## Value

A tibble with one row per observation and columns UA, UD, FA, FD,
treatment (0 = control, 1, 2, ... in the order of `treatments`) and
vector_control_product. Pass it to
[`format_eht_data`](https://swisstph.github.io/AnophelesEntoTrials/reference/format_eht_data.md)
to add the derived columns.

## Details

For a treatment arm, the rates are modified as in the basic Stan model:

- attack rate: `(1 - HostAvailabilityReduction) * alpha_0`

- mortality rate: `mu_0 + KillingDuringHostSeeking * alpha_0`

- survival after feeding: `pC_0 * (1 - PostprandialKilling)`

The control arm uses `alpha_0`, `mu_0` and `pC_0`. Each observation is a
multinomial draw of `n_mosquitoes` over the four outcome probabilities.

## Examples

``` r
# One treatment arm
sim <- simulate_eht_data(seed = 42)
head(sim)
#> # A tibble: 6 × 6
#>      UA    UD    FA    FD treatment vector_control_product
#>   <int> <int> <int> <int>     <int> <chr>                 
#> 1     6     3     1     0         0 Control               
#> 2     6     2     2     0         0 Control               
#> 3     7     0     3     0         0 Control               
#> 4     7     1     2     0         0 Control               
#> 5     5     2     3     0         0 Control               
#> 6     5     5     0     0         0 Control               

# Two treatment arms
sim <- simulate_eht_data(
  treatments = data.frame(
    vector_control_product = c("Standard_ITN", "New_ITN"),
    HostAvailabilityReduction = c(0.4, 0.69),
    KillingDuringHostSeeking = c(1.5, 3.8),
    PostprandialKilling = c(0.3, 0.52)
  ),
  seed = 42
)
eht_data <- format_eht_data(sim)
```
