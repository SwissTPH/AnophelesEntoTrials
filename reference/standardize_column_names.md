# Standardize column names to expected EHT format

Maps common column name variations to standard EHT format. Handles
variations in naming conventions across different datasets.

## Usage

``` r
standardize_column_names(data, mapping = NULL)
```

## Arguments

- data:

  A data frame or tibble.

- mapping:

  List. Custom column name mapping to override or extend defaults. Each
  element should be named with the target column name and contain a
  character vector of possible source names.

## Value

Data frame with standardized column names.

## Details

Default mappings include:

- Unfed alive: UA, unfedalive, UF.alive, u_alive, unfed_alive, etc.

- Unfed dead: UD, unfeddead, UF.dead, u_dead, unfed_dead, etc.

- Fed alive: FA, fedalive, F.alive, f_alive, fed_alive, etc.

- Fed dead: FD, feddead, F.dead, f_dead, fed_dead, etc.

- treatment: treatment, treat, nettype, treatment_type

- vector_control_product: vector_control_product, insecticide, net_type,
  ingredient, lntype

- wash_status: wash_status, washed, condition, wash

## Examples

``` r
if (FALSE) { # \dontrun{
# Standardize column names
data_std <- standardize_column_names(raw_data)

# With custom mapping
data_std <- standardize_column_names(raw_data,
  mapping = list(my_col = c("old_name1", "old_name2")))
} # }
```
