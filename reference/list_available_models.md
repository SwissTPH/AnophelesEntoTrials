# List available predefined Stan models

Returns information about predefined Stan models included in the
package.

## Usage

``` r
list_available_models(details = FALSE)
```

## Arguments

- details:

  Logical. If TRUE, returns detailed information about each model. If
  FALSE (default), returns just model names.

## Value

If details = FALSE, a character vector of model names. If details =
TRUE, a data frame with model information.

## Examples

``` r
# List model names
list_available_models()
#> [1] "basic"     "semifield"

# Get detailed information
list_available_models(details = TRUE)
#> # A tibble: 2 × 5
#>   model     file                     description          required_data use_case
#>   <chr>     <chr>                    <chr>                <chr>         <chr>   
#> 1 basic     eht_basic.stan           Standard multinomia… UA, UD, FA, … Standar…
#> 2 semifield semifield_repellent.stan Matched semifield m… y0, y1        Semifie…
```
