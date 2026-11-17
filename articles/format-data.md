# Formatting entomological data with format_eht_data()

``` r

library(AnophelesEntoTrials)
```

## Formatting of entomological data

Entomological hut trial (EHT) data often comes from various sources with
inconsistent naming. Different teams may use:

- Different column names for the same variables (e.g., `UnfedAlive`,
  `UF.alive`, `u_alive`, `UA`)
- Different cases (e.g., `Treatment`, `treatment`, `TREATMENT`)
- Different treatment encodings (e.g., character labels vs. numeric
  codes)

The
[`format_eht_data()`](https://swisstph.github.io/AnophelesEntoTrials/reference/format_eht_data.md)
function standardizes these data inconsistencies.

## What does format_eht_data() function do?

1.  **Standardizes column names**: Maps common variations to standard
    names (UA, UD, FA, FD, treatment, etc.)
2.  **Coerces data types**: Converts endpoints to numeric
3.  **Encodes treatments**: Converts treatment labels to numeric
    encoding (0 = control, 1+ = treatments), keeping the labels in
    `vector_control_product` to label the results
4.  **Calculates variables**: Adds `fed`, `dead`, and `total` columns

## Example

``` r

# A messy dataset with inconsistent naming conventions
messy_data <- data.frame(
  # Different naming for endpoints
  unfedalive = c(25, 30, 22, 28, 35, 15, 18, 12, 20, 14),
  UF.dead    = c(3, 2, 4, 3, 2, 8, 10, 12, 9, 11),
  fed_alive  = c(18, 15, 20, 17, 12, 10, 8, 6, 9, 7),
  F.dead     = c(2, 1, 3, 2, 1, 7, 9, 10, 8, 8),

  # Treatment as character with various labels
  nettype = c("Control", "Control", "Control", "Control", "Control",
              "IG2", "IG2", "IG2", "IG2", "IG2"),

  # Additional metadata
  hut_id = c(1, 2, 3, 4, 5, 1, 2, 3, 4, 5),
  night = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2)
)

print(messy_data)
#>    unfedalive UF.dead fed_alive F.dead nettype hut_id night
#> 1          25       3        18      2 Control      1     1
#> 2          30       2        15      1 Control      2     1
#> 3          22       4        20      3 Control      3     1
#> 4          28       3        17      2 Control      4     1
#> 5          35       2        12      1 Control      5     1
#> 6          15       8        10      7     IG2      1     2
#> 7          18      10         8      9     IG2      2     2
#> 8          12      12         6     10     IG2      3     2
#> 9          20       9         9      8     IG2      4     2
#> 10         14      11         7      8     IG2      5     2
str(messy_data)
#> 'data.frame':    10 obs. of  7 variables:
#>  $ unfedalive: num  25 30 22 28 35 15 18 12 20 14
#>  $ UF.dead   : num  3 2 4 3 2 8 10 12 9 11
#>  $ fed_alive : num  18 15 20 17 12 10 8 6 9 7
#>  $ F.dead    : num  2 1 3 2 1 7 9 10 8 8
#>  $ nettype   : chr  "Control" "Control" "Control" "Control" ...
#>  $ hut_id    : num  1 2 3 4 5 1 2 3 4 5
#>  $ night     : num  1 1 1 1 1 2 2 2 2 2
```

Now let’s format this data:

``` r

# Apply format_eht_data()
clean_data <- format_eht_data(messy_data, control_value = "Control")

print(clean_data)
#> # A tibble: 10 × 11
#>       UA    UD    FA    FD treatment hut_id night vector_control_product   fed
#>    <dbl> <dbl> <dbl> <dbl>     <int>  <dbl> <dbl> <chr>                  <dbl>
#>  1    25     3    18     2         0      1     1 Control                   20
#>  2    30     2    15     1         0      2     1 Control                   16
#>  3    22     4    20     3         0      3     1 Control                   23
#>  4    28     3    17     2         0      4     1 Control                   19
#>  5    35     2    12     1         0      5     1 Control                   13
#>  6    15     8    10     7         1      1     2 IG2                       17
#>  7    18    10     8     9         1      2     2 IG2                       17
#>  8    12    12     6    10         1      3     2 IG2                       16
#>  9    20     9     9     8         1      4     2 IG2                       17
#> 10    14    11     7     8         1      5     2 IG2                       15
#> # ℹ 2 more variables: dead <dbl>, total <dbl>
str(clean_data)
#> tibble [10 × 11] (S3: tbl_df/tbl/data.frame)
#>  $ UA                    : num [1:10] 25 30 22 28 35 15 18 12 20 14
#>  $ UD                    : num [1:10] 3 2 4 3 2 8 10 12 9 11
#>  $ FA                    : num [1:10] 18 15 20 17 12 10 8 6 9 7
#>  $ FD                    : num [1:10] 2 1 3 2 1 7 9 10 8 8
#>  $ treatment             : int [1:10] 0 0 0 0 0 1 1 1 1 1
#>  $ hut_id                : num [1:10] 1 2 3 4 5 1 2 3 4 5
#>  $ night                 : num [1:10] 1 1 1 1 1 2 2 2 2 2
#>  $ vector_control_product: chr [1:10] "Control" "Control" "Control" "Control" ...
#>  $ fed                   : num [1:10] 20 16 23 19 13 17 17 16 17 15
#>  $ dead                  : num [1:10] 5 3 7 5 3 15 19 22 17 19
#>  $ total                 : num [1:10] 48 48 49 50 50 40 45 40 46 40
```

## Handling Different Naming Conventions

The
[`standardize_column_names()`](https://swisstph.github.io/AnophelesEntoTrials/reference/standardize_column_names.md)
function (called internally) recognizes many common variations:

| Standard Name | Recognized Variations |
|----|----|
| UA | `unfedalive`, `UF.alive`, `u_alive`, `unfed_alive`, `gaufalive24` |
| UD | `unfeddead`, `UF.dead`, `u_dead`, `unfed_dead`, `gaufdead24` |
| FA | `fedalive`, `F.alive`, `f_alive`, `fed_alive`, `gafedalive24` |
| FD | `feddead`, `F.dead`, `f_dead`, `fed_dead`, `gafeddead24` |
| treatment | `treat`, `nettype`, `treatment_type` |
| vector_control_product | `insecticide_name`, `insecticide`, `net_type`, `treatment_name`, `ingredient`, `lntype` |
| wash_status | `washed`, `condition`, `wash` |

Other names can be mapped with the `mapping` argument (standard name =
name in the data):

``` r

custom_data <- data.frame(
  alive = c(30, 15), dead = c(2, 10), alivefed = c(15, 8), deadfed = c(1, 9),
  arm = c("Control", "IG2")
)

format_eht_data(custom_data, control_value = "Control",
                mapping = list(UA = "alive", UD = "dead", FA = "alivefed",
                               FD = "deadfed", treatment = "arm"))
#> # A tibble: 2 × 9
#>      UA    UD    FA    FD treatment vector_control_product   fed  dead total
#>   <dbl> <dbl> <dbl> <dbl>     <int> <chr>                  <dbl> <dbl> <dbl>
#> 1    30     2    15     1         0 Control                   16     3    48
#> 2    15    10     8     9         1 IG2                       17    19    42
```

## Specifying the control group

By default,
[`format_eht_data()`](https://swisstph.github.io/AnophelesEntoTrials/reference/format_eht_data.md)
attempts to auto-detect the control group by looking for: 1. A treatment
value of `0` 2. A treatment label containing “control”
(case-insensitive)

If neither is found, it uses the first unique value and issues a
warning. You can explicitly specify the control (also for numeric codes,
e.g. `control_value = 1` for data coded 1, 2, 3):

``` r

# Data with non-obvious control labeling
data_odd_labels <- data.frame(
  UA = c(30, 28, 15, 12),
  UD = c(2, 3, 10, 11),
  FA = c(15, 17, 8, 6),
  FD = c(1, 2, 9, 10),
  treatment = c("Untreated_Net", "Untreated_Net", "PBO_Net", "PBO_Net")
)

# Explicitly specify control
formatted <- format_eht_data(data_odd_labels, control_value = "Untreated_Net")
print(formatted[, c("treatment", "vector_control_product", "UA", "UD", "FA", "FD")])
#> # A tibble: 4 × 6
#>   treatment vector_control_product    UA    UD    FA    FD
#>       <int> <chr>                  <dbl> <dbl> <dbl> <dbl>
#> 1         0 Untreated_Net             30     2    15     1
#> 2         0 Untreated_Net             28     3    17     2
#> 3         1 PBO_Net                   15    10     8     9
#> 4         1 PBO_Net                   12    11     6    10
```

Each treatment must have a single product name in
`vector_control_product`; otherwise
[`format_eht_data()`](https://swisstph.github.io/AnophelesEntoTrials/reference/format_eht_data.md)
stops, as the results could not be labelled.

## Washed and unwashed nets

To fit washed and unwashed nets as separate arms, combine the net and
its wash status into the treatment label before formatting:

``` r

wash_data <- data.frame(
  UA = c(30, 15, 20), UD = c(2, 10, 6), FA = c(15, 8, 11), FD = c(1, 9, 5),
  net = c("Control", "IG2", "IG2"),
  wash = c("unwashed", "unwashed", "washed")
)
wash_data$treatment <- ifelse(wash_data$net == "Control", "Control",
                              paste0(wash_data$net, " (", wash_data$wash, ")"))

format_eht_data(wash_data, control_value = "Control")[, c("treatment", "vector_control_product")]
#> # A tibble: 3 × 2
#>   treatment vector_control_product
#>       <int> <chr>                 
#> 1         0 Control               
#> 2         1 IG2 (unwashed)        
#> 3         2 IG2 (washed)
```

All arms then share the control arm, so they should come from the same
trial.

## Skipping Derived Variables

If you don’t need the derived variables (`fed`, `dead`, `total`), you
can skip their calculation:

``` r

minimal_format <- format_eht_data(messy_data,
                                   calculate_derived = FALSE,
                                   control_value = "Control")
names(minimal_format)
#> [1] "UA"                     "UD"                     "FA"                    
#> [4] "FD"                     "treatment"              "hut_id"                
#> [7] "night"                  "vector_control_product"
```

## Integration with the analysis pipeline

``` r

# 1. Import raw data
raw_data <- import_ento_data("my_trial_data.csv")
# or, if the treatment and product names are in other columns:
# raw_data <- import_ento_data("my_trial_data.csv",
#                              treatment_var = "arm", product_var = "net_name")

# 2. Format the data
formatted_data <- format_eht_data(raw_data, control_value = "Control")

# 3. Fit the model
fit <- fit_ento_model(
  data = formatted_data,
  model_type = "basic",
  iter = 6000,
  warmup = 3000,
  chains = 3
)
```
