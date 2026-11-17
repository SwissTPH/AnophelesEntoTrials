# Check treatment coding expected by the Stan models

The Stan models index treatment-level parameters with `treat + 1` and
treat `0` as the control arm, so treatment must be coded as consecutive
integers 0, 1, ..., tr.

## Usage

``` r
check_treatment_coding(treatment)
```

## Arguments

- treatment:

  Vector of treatment codes.

## Value

The treatment codes as an integer vector.
