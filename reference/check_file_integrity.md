# Check file integrity and required columns

Validates file existence and format, checks for required columns in the
data, and optionally previews data structure.

## Usage

``` r
check_file_integrity(
  file_path,
  required_vars = c("UA", "UD", "FA", "FD", "treatment"),
  verbose = TRUE
)
```

## Arguments

- file_path:

  Character. Path to the data file (CSV, Excel).

- required_vars:

  Character vector. Names of required columns. Default includes standard
  EHT endpoints (UA, UD, FA, FD) and key variables.

- verbose:

  Logical. If TRUE, print validation status messages.

## Value

Returns the detected file format string ("csv" or "excel") invisibly if
validation passes. Throws an error if validation fails.

## Details

The function performs the following checks:

- File exists and is readable

- File format is supported (.csv, .xlsx, .xls)

- All required columns are present

## Examples

``` r
if (FALSE) { # \dontrun{
# Check a CSV file
check_file_integrity("data.csv")

# Check with custom required variables
check_file_integrity("data.csv",
  required_vars = c("UA", "UD", "FA", "FD", "treatment", "my_new_variable"))
} # }
```
