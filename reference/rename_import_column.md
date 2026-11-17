# Rename an imported column to its standard name

Used by
[`import_ento_data`](https://swisstph.github.io/AnophelesEntoTrials/reference/import_ento_data.md)
to rename the user's treatment and product columns. Column names are
expected to be lowercase already.

## Usage

``` r
rename_import_column(data, from, to)
```

## Arguments

- data:

  Data frame with lowercase column names.

- from:

  Character. Name of the column in the file (case-insensitive).

- to:

  Character. Standard column name.

## Value

`data` with the column renamed.
