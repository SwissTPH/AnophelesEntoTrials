# Format Variable Names for Entomological Data

Applies consistent naming rules for imported datasets:

- Endpoints (`UA`, `UD`, `FA`, `FD`) are returned in **upper case**

- All other variables are converted to **sentence case** (first letter
  capitalised, rest lower case)

## Usage

``` r
format_varnames(names)
```

## Arguments

- names:

  A character vector of variable names to format.

## Value

A character vector with formatted variable names.

## Details

This function is used internally when importing raw EHT data to ensure
consistent column naming across heterogeneous input files.

## Examples

``` r
format_varnames(c("ua", "treatment", "Wash_STATUS"))
#>            ua     treatment   Wash_STATUS 
#>          "UA"   "Treatment" "Wash_status" 
# Returns: "UA", "Treatment", "Wash_status"
```
