# Tests for data import functions

# Helper function to create test data
create_test_data <- function(n = 10, col_case = "upper") {
  set.seed(123)

  if (col_case == "upper") {
    data.frame(
      UA = sample(10:50, n, replace = TRUE),
      UD = sample(0:10, n, replace = TRUE),
      FA = sample(5:30, n, replace = TRUE),
      FD = sample(0:5, n, replace = TRUE),
      treatment = rep(c(0, 1), each = n/2),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      ua = sample(10:50, n, replace = TRUE),
      ud = sample(0:10, n, replace = TRUE),
      fa = sample(5:30, n, replace = TRUE),
      fd = sample(0:5, n, replace = TRUE),
      treatment = rep(c(0, 1), each = n/2),
      stringsAsFactors = FALSE
    )
  }
}

# ============================================================================
# Tests for standardize_column_names()
# ============================================================================

test_that("standardize_column_names handles standard column names", {
  data <- create_test_data(col_case = "upper")
  result <- standardize_column_names(data)

  expect_true("UA" %in% names(result))
  expect_true("UD" %in% names(result))
  expect_true("FA" %in% names(result))
  expect_true("FD" %in% names(result))
  expect_true("treatment" %in% names(result))
})

test_that("standardize_column_names handles lowercase column names", {
  data <- create_test_data(col_case = "lower")
  result <- standardize_column_names(data)

  expect_true("UA" %in% names(result))
  expect_true("UD" %in% names(result))
  expect_true("FA" %in% names(result))
  expect_true("FD" %in% names(result))
})

test_that("standardize_column_names handles alternative column name formats", {
  data <- data.frame(
    unfedalive = c(10, 20),
    unfeddead = c(1, 2),
    fedalive = c(5, 10),
    feddead = c(0, 1),
    treat = c(0, 1)
  )

  result <- standardize_column_names(data)

  expect_true("UA" %in% names(result))
  expect_true("UD" %in% names(result))
  expect_true("FA" %in% names(result))
  expect_true("FD" %in% names(result))
  expect_true("treatment" %in% names(result))
})

test_that("standardize_column_names handles treatment_name column", {
  data <- data.frame(
    UA = c(10, 20),
    UD = c(1, 2),
    FA = c(5, 10),
    FD = c(0, 1),
    treatment = c(0, 1),
    treatment_name = c("Control", "Treatment A")
  )

  result <- standardize_column_names(data)
  expect_true("vector_control_product" %in% names(result))
})

test_that("standardize_column_names handles custom mapping", {
  data <- data.frame(
    my_custom_ua = c(10, 20),
    UD = c(1, 2),
    FA = c(5, 10),
    FD = c(0, 1),
    treatment = c(0, 1)
  )

  result <- standardize_column_names(data, mapping = list(UA = "my_custom_ua"))
  expect_true("UA" %in% names(result))
})

test_that("standardize_column_names preserves unrecognized columns", {
  data <- data.frame(
    UA = c(10, 20),
    UD = c(1, 2),
    FA = c(5, 10),
    FD = c(0, 1),
    treatment = c(0, 1),
    my_extra_column = c("a", "b")
  )

  result <- standardize_column_names(data)
  expect_true("my_extra_column" %in% names(result))
})

# ============================================================================
# Tests for format_eht_data()
# ============================================================================

test_that("format_eht_data returns tibble with required columns", {
  data <- create_test_data()
  result <- format_eht_data(data)

  expect_s3_class(result, "tbl_df")
  expect_true("UA" %in% names(result))
  expect_true("UD" %in% names(result))
  expect_true("FA" %in% names(result))
  expect_true("FD" %in% names(result))
  expect_true("treatment" %in% names(result))
})

test_that("format_eht_data calculates derived columns correctly", {
  data <- data.frame(
    UA = c(10, 20),
    UD = c(2, 4),
    FA = c(5, 10),
    FD = c(1, 2),
    treatment = c(0, 1)
  )

  result <- format_eht_data(data)
#
  expect_equal(result$fed, c(6, 12))    # FA + FD
expected_dead <- c(3, 6)   # UD + FD
  expect_equal(result$dead, expected_dead)
  expect_equal(result$total, c(18, 36)) # UA + UD + FA + FD
})

test_that("format_eht_data skips derived columns when requested", {
  data <- create_test_data()
  result <- format_eht_data(data, calculate_derived = FALSE)

  expect_false("fed" %in% names(result))
  expect_false("dead" %in% names(result))
  expect_false("total" %in% names(result))
})

test_that("format_eht_data coerces endpoints to numeric", {
  data <- data.frame(
    UA = c("10", "20"),
    UD = c("2", "4"),
    FA = c("5", "10"),
    FD = c("1", "2"),
    treatment = c(0, 1),
    stringsAsFactors = FALSE
  )

  result <- format_eht_data(data)

  expect_type(result$UA, "double")
  expect_type(result$UD, "double")
  expect_type(result$FA, "double")
  expect_type(result$FD, "double")
})

test_that("format_eht_data handles character treatment values", {
  data <- data.frame(
    UA = c(10, 20, 30),
    UD = c(2, 4, 6),
    FA = c(5, 10, 15),
    FD = c(1, 2, 3),
    treatment = c("control", "ITN_A", "ITN_B"),
    stringsAsFactors = FALSE
  )

  result <- format_eht_data(data)

  expect_type(result$treatment, "integer")
  expect_equal(result$treatment[result$treatment == 0], 0)  # control = 0
})

test_that("format_eht_data uses specified control_value", {
  data <- data.frame(
    UA = c(10, 20, 30),
    UD = c(2, 4, 6),
    FA = c(5, 10, 15),
    FD = c(1, 2, 3),
    treatment = c("untreated", "ITN_A", "ITN_B"),
    stringsAsFactors = FALSE
  )

  result <- format_eht_data(data, control_value = "untreated")

  # Control should be 0
  expect_equal(result$treatment[1], 0)
})

test_that("format_eht_data re-encodes numeric treatment when control_value is given", {
  data <- data.frame(
    UA = c(10, 20, 30),
    UD = c(2, 4, 6),
    FA = c(5, 10, 15),
    FD = c(1, 2, 3),
    treatment = c(1, 2, 3)
  )

  result <- format_eht_data(data, control_value = 1)

  expect_identical(result$treatment, c(0L, 1L, 2L))
})

test_that("format_eht_data leaves numeric treatment unchanged without control_value", {
  data <- data.frame(
    UA = c(10, 20),
    UD = c(2, 4),
    FA = c(5, 10),
    FD = c(1, 2),
    treatment = c(0, 1)
  )

  result <- format_eht_data(data)

  expect_equal(result$treatment, c(0, 1))
})

test_that("format_eht_data errors when control_value is not in the data", {
  data <- data.frame(
    UA = c(10, 20),
    UD = c(2, 4),
    FA = c(5, 10),
    FD = c(1, 2),
    treatment = c("A", "B")
  )

  expect_error(
    format_eht_data(data, control_value = "control"),
    "control_value 'control' not found"
  )
})

test_that("format_eht_data keeps text treatment labels as product names", {
  data <- data.frame(
    UA = c(10, 20, 30),
    UD = c(2, 4, 6),
    FA = c(5, 10, 15),
    FD = c(1, 2, 3),
    treatment = c("control", "ITN_A", "ITN_B")
  )

  result <- format_eht_data(data)

  expect_identical(result$treatment, c(0L, 1L, 2L))
  expect_equal(result$vector_control_product, c("control", "ITN_A", "ITN_B"))
})

test_that("format_eht_data does not overwrite an existing product column", {
  data <- data.frame(
    UA = c(10, 20),
    UD = c(2, 4),
    FA = c(5, 10),
    FD = c(1, 2),
    treatment = c("control", "arm_1"),
    vector_control_product = c("Untreated", "PBO_Net")
  )

  result <- format_eht_data(data)

  expect_equal(result$vector_control_product, c("Untreated", "PBO_Net"))
})

test_that("format_eht_data applies a custom column mapping", {
  data <- data.frame(
    alive = c(10, 20),
    dead = c(2, 4),
    alivefed = c(5, 10),
    deadfed = c(1, 2),
    treatment = c(0, 1)
  )

  result <- format_eht_data(data, mapping = list(
    UA = "alive", UD = "dead", FA = "alivefed", FD = "deadfed"
  ))

  expect_true(all(c("UA", "UD", "FA", "FD") %in% names(result)))
  expect_equal(result$total, c(18, 36))
})

test_that("format_eht_data errors on missing required columns", {
  data <- data.frame(
    UA = c(10, 20),
    UD = c(2, 4),
    treatment = c(0, 1)
  )

  expect_error(format_eht_data(data), "Missing required endpoint columns")
})


# ============================================================================
# Tests for check_file_integrity() - file-based tests
# ============================================================================

test_that("check_file_integrity errors on non-existent file", {
  expect_error(
    check_file_integrity("/nonexistent/path/file.csv"),
    "File does not exist"
  )
})

test_that("check_file_integrity errors on unsupported file format", {
  # Create a temp file with unsupported extension
  temp_file <- tempfile(fileext = ".txt")
  writeLines("test", temp_file)
  on.exit(unlink(temp_file))

  expect_error(
    check_file_integrity(temp_file),
    "Unsupported file format"
  )
})

test_that("check_file_integrity validates CSV files", {
  # Create a temp CSV file with required columns
  temp_file <- tempfile(fileext = ".csv")
  data <- create_test_data()
  write.csv(data, temp_file, row.names = FALSE)
  on.exit(unlink(temp_file))

  result <- check_file_integrity(temp_file, verbose = FALSE)
  expect_equal(result, "csv")
})

test_that("check_file_integrity detects missing columns", {
  temp_file <- tempfile(fileext = ".csv")
  data <- data.frame(UA = 1, UD = 2)  # Missing FA, FD, treatment
  write.csv(data, temp_file, row.names = FALSE)
  on.exit(unlink(temp_file))

  expect_error(
    check_file_integrity(temp_file, verbose = FALSE),
    "Missing required columns"
  )
})

# ============================================================================
# Tests for import_ento_data() - file-based tests
# ============================================================================

test_that("import_ento_data imports CSV correctly", {
  temp_file <- tempfile(fileext = ".csv")
  data <- create_test_data()
  write.csv(data, temp_file, row.names = FALSE)
  on.exit(unlink(temp_file))

  result <- import_ento_data(temp_file, quiet = TRUE)

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 10)
})

test_that("import_ento_data converts endpoints to numeric", {
  temp_file <- tempfile(fileext = ".csv")
  data <- create_test_data()
  write.csv(data, temp_file, row.names = FALSE)
  on.exit(unlink(temp_file))

  result <- import_ento_data(temp_file, quiet = TRUE)

  expect_type(result$ua, "double")
  expect_type(result$ud, "double")
  expect_type(result$fa, "double")
  expect_type(result$fd, "double")
})

test_that("import_ento_data handles case-insensitive column names", {
  temp_file <- tempfile(fileext = ".csv")
  data <- data.frame(
    ua = c(10, 20),
    UD = c(2, 4),
    Fa = c(5, 10),
    fD = c(1, 2),
    TREATMENT = c(0, 1)
  )
  write.csv(data, temp_file, row.names = FALSE)
  on.exit(unlink(temp_file))

  result <- import_ento_data(temp_file, quiet = TRUE)

  expect_true("treatment" %in% names(result))
})

test_that("import_ento_data renames treatment_var to treatment", {
  temp_file <- tempfile(fileext = ".csv")
  data <- data.frame(
    UA = c(10, 20),
    UD = c(2, 4),
    FA = c(5, 10),
    FD = c(1, 2),
    Intervention = c("control", "ITN")
  )
  write.csv(data, temp_file, row.names = FALSE)
  on.exit(unlink(temp_file))

  result <- import_ento_data(temp_file, treatment_var = "Intervention", quiet = TRUE)

  expect_true("treatment" %in% names(result))
  expect_false("intervention" %in% names(result))
  expect_equal(result$treatment, c("control", "ITN"))
})

test_that("import_ento_data errors if treatment_var would duplicate treatment", {
  temp_file <- tempfile(fileext = ".csv")
  data <- data.frame(
    UA = c(10, 20),
    UD = c(2, 4),
    FA = c(5, 10),
    FD = c(1, 2),
    treatment = c(0, 1),
    Intervention = c("control", "ITN")
  )
  write.csv(data, temp_file, row.names = FALSE)
  on.exit(unlink(temp_file))

  expect_error(
    import_ento_data(temp_file, treatment_var = "Intervention", quiet = TRUE),
    "already has a 'treatment' column"
  )
})

test_that("import_ento_data errors if treatment_var is missing", {
  temp_file <- tempfile(fileext = ".csv")
  data <- create_test_data()
  write.csv(data, temp_file, row.names = FALSE)
  on.exit(unlink(temp_file))

  expect_error(
    import_ento_data(temp_file, treatment_var = "Intervention", quiet = TRUE),
    "Missing required columns: intervention"
  )
})

test_that("import_ento_data renames product_var to vector_control_product", {
  temp_file <- tempfile(fileext = ".csv")
  data <- data.frame(
    UA = c(10, 20),
    UD = c(2, 4),
    FA = c(5, 10),
    FD = c(1, 2),
    arm = c(0, 1),
    Net_Name = c("Untreated", "PBO_Net")
  )
  write.csv(data, temp_file, row.names = FALSE)
  on.exit(unlink(temp_file))

  result <- import_ento_data(temp_file, treatment_var = "arm",
                             product_var = "Net_Name", quiet = TRUE)

  expect_equal(result$treatment, c(0, 1))
  expect_equal(result$vector_control_product, c("Untreated", "PBO_Net"))
})

test_that("import_ento_data errors if product_var is missing", {
  temp_file <- tempfile(fileext = ".csv")
  data <- create_test_data()
  write.csv(data, temp_file, row.names = FALSE)
  on.exit(unlink(temp_file))

  expect_error(
    import_ento_data(temp_file, product_var = "net_name", quiet = TRUE),
    "Missing required columns: net_name"
  )
})

test_that("import_ento_data allows custom required_vars", {
  temp_file <- tempfile(fileext = ".csv")
  data <- data.frame(
    UA = c(10, 20),
    UD = c(2, 4),
    FA = c(5, 10),
    FD = c(1, 2),
    treatment = c(0, 1),
    custom_col = c("a", "b")
  )
  write.csv(data, temp_file, row.names = FALSE)
  on.exit(unlink(temp_file))

  result <- import_ento_data(
    temp_file,
    required_vars = c("UA", "UD", "FA", "FD", "treatment", "custom_col"),
    quiet = TRUE
  )

  expect_true("custom_col" %in% names(result))
})
