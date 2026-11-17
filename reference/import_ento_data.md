# Import entomological hut trial data

Universal import function that reads entomologica data from CSV or Excel
files and returns a standardized tibble with appropriate data types.

## Usage

``` r
import_ento_data(
  file_path,
  required_vars = c("UA", "UD", "FA", "FD", "treatment"),
  treatment_var = "treatment",
  product_var = NULL,
  quiet = FALSE
)
```

## Arguments

- file_path:

  Character. Path to the data file.

- required_vars:

  Character vector. Required columns to check. `"treatment"` is replaced
  by `treatment_var`, which is always required.

- treatment_var:

  Character. Name of the treatment column in the file
  (case-insensitive). It is renamed to `treatment`, the name used by
  [`format_eht_data`](https://swisstph.github.io/AnophelesEntoTrials/reference/format_eht_data.md)
  and the model functions.

- product_var:

  Character or NULL. Name of the column with the vector control product
  names (case-insensitive). If given, it is required and renamed to
  `vector_control_product`, which the extract functions use to label
  results. If NULL (default), no product column is required.

- quiet:

  Logical. If FALSE, print import status.

## Value

A tibble with standardized column names and data types.

- Endpoints (ua, ud, fa, fd) are numeric.

- Treatment variable: renamed to `treatment`, values preserved as-is

- Product variable (if `product_var` is given): renamed to
  `vector_control_product`, values preserved as-is

- Other variables: preserved as-is

## Details

The import function:

- Validates file integrity

- Reads data with appropriate reader function

- Converts to tibble format

- Standardizes column names to lowercase

- Renames the treatment column to `treatment`

- Renames the product column to `vector_control_product` (if
  `product_var` is given)

- Coerces endpoints to numeric

- Returns as `tibble`

## Examples

``` r
if (FALSE) { # \dontrun{
# Import CSV file
data <- import_ento_data("ITN_data.csv")

# Import Excel file
data <- import_ento_data("trial_data.xlsx")

# Treatment stored in a column called "Intervention"
data <- import_ento_data("trial_data.csv", treatment_var = "Intervention")

# Treatment codes in "arm", product names in "net_name"
data <- import_ento_data("trial_data.csv",
  treatment_var = "arm", product_var = "net_name")
} # }
```
