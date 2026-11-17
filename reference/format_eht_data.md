# Format EHT data to standard structure

Transforms raw EHT data into standardized format required for analysis.
Handles column renaming, derived variable calculation, and type
coercion. Stops if a treatment code has more than one
vector_control_product, since the extract functions would then duplicate
rows when labelling results. The same product may appear under several
codes (e.g. washed and unwashed).

## Usage

``` r
format_eht_data(
  data,
  calculate_derived = TRUE,
  coerce_types = TRUE,
  control_value = NULL,
  mapping = NULL
)
```

## Arguments

- data:

  Data frame or tibble. Raw EHT data.

- calculate_derived:

  Logical. If TRUE, calculate fed, dead, total columns.

- coerce_types:

  Logical. If TRUE, coerce columns to appropriate types.

- control_value:

  Value that identifies the control group in the treatment column. If
  NULL, attempts to detect automatically (looks for 0 or "control").
  Numeric treatment columns are only re-encoded when `control_value` is
  given (e.g. `control_value = 1` for data coded 1, 2, 3).

- mapping:

  List. Custom column name mapping passed to
  [`standardize_column_names`](https://swisstph.github.io/AnophelesEntoTrials/reference/standardize_column_names.md),
  e.g. `list(UA = "alive", FA = "alivefed")`. Each element is named with
  the standard column name and contains the source name(s) in the data.

## Value

Formatted tibble with standardized structure:

- Endpoints (UA, UD, FA, FD) as numeric

- treatment as integer (0 = control, 1+ = treatments)

- vector_control_product as character (if present; otherwise created
  from text treatment labels before they are encoded)

- Derived variables: fed, dead, total (if calculate_derived = TRUE)

## Examples

``` r
if (FALSE) { # \dontrun{
# Format raw data
data_formatted <- format_eht_data(raw_data)

# Specify control group
data_formatted <- format_eht_data(raw_data, control_value = "untreated")

# Map non-standard column names
data_formatted <- format_eht_data(raw_data,
  mapping = list(UA = "alive", UD = "dead", FA = "alivefed", FD = "deadfed"))
} # }
```
