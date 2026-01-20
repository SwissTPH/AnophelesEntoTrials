#' Data Import Module for AnophelesInterventions
#'
#' Functions to import entomological hut trial data from multiple file formats
#' (CSV, Excel) with automatic format detection and preprocessing.
#'
#' @name data_import
#' @docType _PACKAGE
NULL

#' Check file integrity and required columns
#'
#' Validates file existence and format, checks for required columns in the data,
#' and optionally previews data structure.
#'
#' @param file_path Character. Path to the data file (CSV, Excel).
#' @param required_vars Character vector. Names of required columns.
#'   Default includes standard EHT endpoints (UA, UD, FA, FD) and key variables.
#' @param verbose Logical. If TRUE, print validation status messages.
#'
#' @return Returns the detected file format string ("csv" or "excel") invisibly
#'   if validation passes. Throws an error if validation fails.
#'
#' @details
#' The function performs the following checks:
#' \itemize{
#'   \item File exists and is readable
#'   \item File format is supported (.csv, .xlsx, .xls)
#'   \item All required columns are present
#' }
#'
#' @examples
#' \dontrun{
#' # Check a CSV file
#' check_file_integrity("data.csv")
#'
#' # Check with custom required variables
#' check_file_integrity("data.csv",
#'   required_vars = c("UA", "UD", "FA", "FD", "treatment", "my_new_variable"))
#' }
#'
#' @export
check_file_integrity <- function(file_path,
                                 required_vars = c("UA", "UD", "FA", "FD", "treatment"),
                                 verbose = TRUE) {

  # Check file exists
  if (!file.exists(file_path)) {
    stop("File does not exist: ", file_path, call. = FALSE)
  }

  # Detect file format
  ext <- tolower(tools::file_ext(file_path))
  file_format <- switch(ext,
                        "csv" = "csv",
                        "xlsx" = "excel",
                        "xls" = "excel",
                        stop("Unsupported file format: .", ext,
                             " (Only .csv, .xlsx, .xls supported)", call. = FALSE)
  )

  # Read preview of data to check columns
  preview <- tryCatch(
    {
      if (file_format == "csv") {
        readr::read_csv(file_path, n_max = 0, show_col_types = FALSE)
      } else if (file_format == "excel") {
        readxl::read_excel(file_path, n_max = 0)
      }
    },
    error = function(e) {
      stop("Error reading file preview: ", e$message, call. = FALSE)
    }
  )

  #check for required columns
  preview_names_lower <- tolower(names(preview))
  required_vars_lower <- tolower(required_vars)
  missing_cols <- setdiff(required_vars_lower, preview_names_lower)

  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "),
         call. = FALSE)
  }

  if (verbose) {
    cli::cli_alert_success("File validation passed: {.file {file_path}}")
    cli::cli_alert_info("File format: {file_format}")
    cli::cli_alert_info("Found columns: {paste(names(preview), collapse = ', ')}")
  }

  invisible(file_format)
}

#' Import entomological hut trial data
#'
#' Universal import function that reads entomologica data
#' from CSV or Excel files
#' and returns a standardized tibble with appropriate data types.
#'
#' @param file_path Character. Path to the data file.
#' @param required_vars Character vector. Required columns to check.
#' @param treatment_var Character. Name of the treatment variable.
#' @param quiet Logical. If FALSE, print import status.
#'
#' @return A tibble with standardized column names and data types.
#'   - Endpoints (UA, UD, FA, FD) are numeric.
#'   - Treatment variable: factor (with levels sorted
#'   alphabetically, control group as reference)
#'   - Other variables: preserved as-is
#'
#' @details
#' The import function:
#' \itemize{
#'   \item Validates file integrity
#'   \item Reads data with appropriate reader function
#'   \item Converts to tibble format
#'   \item Standardizes column names to lowercase
#'   \item Coerces endpoints to numeric
#'   \item Returns as \code{tibble}
#' }
#'
#' @examples
#' \dontrun{
#' # Import CSV file
#' data <- import_ento_data("ITN_data.csv")
#'
#' # Import Excel file
#' data <- import_ento_data("trial_data.xlsx")
#' }
#'
#' @export
import_ento_data <- function(file_path,
                            required_vars = c("UA", "UD", "FA", "FD", "treatment"),
                            treatment_var = "treatment",
                            quiet = FALSE) {

  # Validate file
  file_format <- check_file_integrity(file_path, required_vars, verbose = !quiet)

  # Import data
  data <- tryCatch(
    {
      if (file_format == "csv") {
        readr::read_csv(file_path, show_col_types = FALSE)
      } else if (file_format == "excel") {
        readxl::read_excel(file_path)
      }
    },
    error = function(e) {
      stop("Error importing file: ", e$message, call. = FALSE)
    }
  )

  # Clean names
  names(data) <- tolower(names(data))

  # Ensure Treatment variable
  if (!treatment_var %in% names(data)) {
    # Try to find treatment variable with different case
    possible_names <- names(data)[tolower(names(data)) == tolower(treatment_var)]
    if (length(possible_names) > 0) {
      names(data)[names(data) == possible_names[1]] <- treatment_var
    } else {
      stop("Treatment variable '", treatment_var, "' not found in data. Available columns: ",
           paste(names(data), collapse = ", "), call. = FALSE)
    }
  }

  # Convert endpoints to numeric
  endpoint_cols <- c("ua", "ud", "fa", "fd")
  for (col in endpoint_cols[endpoint_cols %in% names(data)]) {
    data[[col]] <- as.numeric(data[[col]])
  }

  # Convert to tibble
  data <- tibble::as_tibble(data)

  if (!quiet) {
    cli::cli_alert_success("Data imported successfully")
    cli::cli_alert_info("Dimensions: {nrow(data)} rows × {ncol(data)} columns")
  }

  return(data)
}

#' Format Variable Names for Entomological Data
#'
#' Applies consistent naming rules for imported datasets:
#' \itemize{
#'   \item Endpoints (`UA`, `UD`, `FA`, `FD`) are returned in **upper case**
#'   \item All other variables are converted to **sentence case**
#'         (first letter capitalised, rest lower case)
#' }
#'
#' This function is used internally when importing raw EHT data to ensure
#' consistent column naming across heterogeneous input files.
#'
#' @param names A character vector of variable names to format.
#'
#' @return A character vector with formatted variable names.
#'
#' @examples
#' format_varnames(c("ua", "treatment", "Wash_STATUS"))
#' # Returns: "UA", "Treatment", "Wash_status"
#'
#' @keywords internal
#' @export
format_varnames <- function(names) {

  # Input validation
  if (!is.character(names)) {
    stop("`names` must be a character vector.", call. = FALSE)
  }

  # Variables that must remain uppercase
  upper_vars <- c("UA", "UD", "FA", "FD")

  # Transformation
  vapply(
    names,
    FUN = function(nm) {

      if (is.na(nm) || nm == "") {
        return(nm)  # keep NA or empty strings unchanged
      }

      nm_clean <- toupper(nm)

      # Keep UA/UD/FA/FD uppercase
      if (nm_clean %in% upper_vars) {
        return(nm_clean)
      }

      # Sentence case for all other variables
      nm_lower <- tolower(nm)
      paste0(
        toupper(substr(nm_lower, 1, 1)),
        substr(nm_lower, 2, nchar(nm_lower))
      )
    },
    FUN.VALUE = character(1)
  )
}

#' Standardize column names to expected EHT format
#'
#' Maps common column name variations to standard EHT format.
#' Handles variations in naming conventions across different datasets.
#'
#' @param data A data frame or tibble.
#' @param mapping List. Custom column name mapping to override or extend defaults.
#'   Each element should be named with the target column name and contain a
#'   character vector of possible source names.
#'
#' @return Data frame with standardized column names.
#'
#' @details
#' Default mappings include:
#' \itemize{
#'   \item Unfed alive: UA, unfedalive, UF.alive, u_alive, unfed_alive, etc.
#'   \item Unfed dead: UD, unfeddead, UF.dead, u_dead, unfed_dead, etc.
#'   \item Fed alive: FA, fedalive, F.alive, f_alive, fed_alive, etc.
#'   \item Fed dead: FD, feddead, F.dead, f_dead, fed_dead, etc.
#'   \item treatment: treatment, treat, nettype, treatment_type
#'   \item insecticide_name: insecticide_name, insecticide, net_type, ingredient, lntype
#'   \item wash_status: wash_status, washed, condition, wash
#' }
#'
#' @examples
#' \dontrun{
#' # Standardize column names
#' data_std <- standardize_column_names(raw_data)
#'
#' # With custom mapping
#' data_std <- standardize_column_names(raw_data,
#'   mapping = list(my_col = c("old_name1", "old_name2")))
#' }
#'
#' @export
standardize_column_names <- function(data, mapping = NULL) {

  # Default mapping for common variations
  default_mapping <- list(
    # Unfed alive variations
    UA = c("UA", "unfedalive", "UF.alive", "u_alive",
           "unfed_alive", "unfedalive24", "unfedalive72",
           "gaufalive24", "gaufalive72"),
    # Unfed dead variations
    UD = c("UD", "unfeddead", "UF.dead", "u_dead",
           "unfed_dead", "unfeddead24", "unfeddead72",
           "gaufdead24", "gaufdead72"),
    # Fed alive variations
    FA = c("FA", "fedalive", "F.alive", "f_alive",
           "fed_alive", "fedalive24", "fedalive72",
           "gafedalive24", "gafedalive72"),
    # Fed dead variations
    FD = c("FD", "feddead", "F.dead", "f_dead",
           "fed_dead", "feddead24", "feddead72",
           "gafeddead24", "gafeddead72"),
    # Treatment variations
    treatment = c("treatment", "treat", "nettype", "treatment_type"),
    # Insecticide name variations
    insecticide_name = c("insecticide_name", "insecticide", "net_type",
                         "treatment_name", "ingredient", "lntype"),
    # Wash status variations
    wash_status = c("wash_status", "washed", "condition", "wash")
  )

  # Merge with user mapping (user mapping takes precedence)
  if (!is.null(mapping)) {
    for (nm in names(mapping)) {
      default_mapping[[nm]] <- mapping[[nm]]
    }
  }

  # Create inverse mapping: old_name -> new_name
  col_map <- data.frame(
    old = unlist(default_mapping),
    new = rep(names(default_mapping), lengths(default_mapping)),
    stringsAsFactors = FALSE
  )

  # Match current column names (case-insensitive)
  current_names <- names(data)
  current_names_lower <- tolower(current_names)

  for (i in seq_along(current_names)) {
    match_idx <- which(tolower(col_map$old) == current_names_lower[i])
    if (length(match_idx) > 0) {
      names(data)[i] <- col_map$new[match_idx[1]]
    }
  }

  return(data)
}

#' Format EHT data to standard structure
#'
#' Transforms raw EHT data into standardized format required for analysis.
#' Handles column renaming, derived variable calculation, and type coercion.
#'
#' @param data Data frame or tibble. Raw EHT data.
#' @param calculate_derived Logical. If TRUE, calculate fed, dead, total columns.
#' @param coerce_types Logical. If TRUE, coerce columns to appropriate types.
#' @param control_value Value that identifies the control group in the treatment
#'   column. If NULL, attempts to detect automatically (looks for 0 or "control").
#'
#' @return Formatted tibble with standardized structure:
#'   \itemize{
#'     \item Endpoints (UA, UD, FA, FD) as numeric
#'     \item treatment as integer (0 = control, 1+ = treatments)
#'     \item insecticide_name as character (if present)
#'     \item Derived variables: fed, dead, total (if calculate_derived = TRUE)
#'   }
#'
#' @examples
#' \dontrun{
#' # Format raw data
#' data_formatted <- format_eht_data(raw_data)
#'
#' # Specify control group
#' data_formatted <- format_eht_data(raw_data, control_value = "untreated")
#' }
#'
#' @export
format_eht_data <- function(data,
                            calculate_derived = TRUE,
                            coerce_types = TRUE,
                            control_value = NULL) {

  # Standardize column names first
  data <- standardize_column_names(data)

  # Define required endpoints
  required_endpoints <- c("UA", "UD", "FA", "FD")

  # Check if required endpoints exist
  missing_endpoints <- setdiff(required_endpoints, names(data))

  if (length(missing_endpoints) > 0) {
    stop("Missing required endpoint columns: ",
         paste(missing_endpoints, collapse = ", "),
         call. = FALSE)
  }

  # Coerce endpoints to numeric
  if (coerce_types) {
    for (col in required_endpoints) {
      data[[col]] <- suppressWarnings(as.numeric(data[[col]]))
    }
  }

  # Ensure treatment column exists and is properly encoded
  if ("treatment" %in% names(data)) {
    # Convert treatment to numeric encoding
    if (!is.numeric(data$treatment)) {
      unique_treats <- unique(data$treatment)

      # Detect control value
      if (is.null(control_value)) {
        if (0 %in% unique_treats) {
          control_value <- 0
        } else if ("control" %in% tolower(unique_treats)) {
          control_value <- unique_treats[tolower(unique_treats) == "control"][1]
        } else {
          # Use first unique value as control
          control_value <- unique_treats[1]
          warning("Control group not detected. Using '", control_value, "' as control.",
                  call. = FALSE)
        }
      }

      # Create numeric encoding with control = 0
      treat_levels <- c(control_value, setdiff(unique_treats, control_value))
      data$treatment <- match(data$treatment, treat_levels) - 1
    }
  }

  # Calculate derived variables
  if (calculate_derived) {
    data$fed <- data$FA + data$FD
    data$dead <- data$UD + data$FD
    data$total <- data$UA + data$UD + data$FA + data$FD
  }

  # Convert to tibble
  data <- tibble::as_tibble(data)

  return(data)
}
