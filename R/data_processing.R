#' Process Officer-Involved Shootings Data
#'
#' @title Process OIS Data for Specification Curve Analysis
#' @description Processes officer-involved shootings data from Excel format,
#' aggregates to yearly level, and calculates derived metrics.
#'
#' @param file_path Character. Path to the Excel file containing OIS data.
#'
#' @return Data table with processed OIS data including:
#'   \itemize{
#'     \item \code{stateid} - State identifier
#'     \item \code{year} - Year
#'     \item \code{ois} - Number of officer-involved shootings
#'     \item \code{unarmed} - Number of unarmed shootings
#'     \item \code{mh} - Mental health related shootings
#'     \item \code{officers} - Number of officers
#'     \item \code{pop} - Population
#'     \item \code{oisp} - OIS per million population
#'     \item \code{oiso} - OIS per thousand officers
#'     \item \code{ofp} - Officers per population
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Process OIS data
#' ois_data <- process_ois_data("path/to/ois_data.xlsx")
#' }
process_ois_data <- function(file_path) {
  # Load dataset
  dataset <- read_xlsx(file_path)
  dataset <- data.table(dataset)
  
  # Data preprocessing
  dataset[is.na(ois), ois := 0]
  dataset[is.na(unarmed), unarmed := 0]
  dataset[is.na(mh), mh := 0]
  
  # Create date variables and convert to year level
  dataset[, mdate := as.Date(modate, format = "%Y-%m-%d")]
  dataset[, year := year(mdate)]
  dataset[, month := month(mdate)]
  
  # Aggregate to year level
  dataset <- dataset[, .(
    ois = sum(ois),
    unarmed = sum(unarmed),
    mh = sum(mh),
    officers = mean(officers),
    pop = mean(pop)
  ), by = c("stateid", "year")]
  
  # Calculate derived metrics
  dataset[, oisp := (ois / pop) * 1000000]  # OIS per million population
  dataset[, oiso := (ois / officers) * 1000] # OIS per thousand officers
  dataset[, ofp := officers / pop]          # Officers per population
  
  return(dataset)
}

#' Process Homicide Data
#'
#' @title Process Homicide Data for Specification Curve Analysis
#' @description Processes homicide data, calculates rates, and cleans unit identifiers.
#'
#' @param file_path Character. Path to the CSV file containing homicide data.
#'
#' @return Data table with processed homicide data including:
#'   \itemize{
#'     \item \code{ori9} - Cleaned unit identifier
#'     \item \code{num_homicide} - Number of homicides
#'     \item \code{hr_rate} - Homicide rate (homicides per population)
#'     \item \code{population} - Population
#'     \item Additional variables from original dataset
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Process homicide data
#' homicide_data <- process_homicide_data("path/to/homicide_data.csv")
#' }
process_homicide_data <- function(file_path) {
  # Load dataset
  dataset <- fread(file_path)
  
  # Preprocessing
  dataset[, hr_rate := num_homicide/population]
  dataset[, ori9 := gsub(",| |-|\\.", "_", ori9)]
  dataset[, ori9 := gsub("_+", "_", ori9)]
  dataset[, num_homicide := as.numeric(num_homicide)]
  dataset[, hr_rate := as.numeric(hr_rate)]
  
  return(dataset)
}

#' Create Analysis Configuration
#'
#' @title Create Configuration for Multiple Dataset Analysis
#' @description Creates standardized configuration objects for running
#' specification curve analysis across multiple datasets.
#'
#' @param analysis_type Character. Type of analysis ("homicide" or "ois").
#' @param custom_params List. Custom parameters to override defaults.
#'
#' @return List containing configuration parameters for the specified analysis type.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create homicide analysis configuration
#' homicide_config <- create_analysis_config("homicide")
#'
#' # Create OIS analysis configuration with custom parameters
#' ois_config <- create_analysis_config("ois", 
#'   custom_params = list(treated_period = 2021))
#' }
create_analysis_config <- function(analysis_type = c("homicide", "ois"), custom_params = list()) {
  analysis_type <- match.arg(analysis_type)
  
  # Default configurations
  base_configs <- list(
    homicide = list(
      outcomes = c("num_homicide", "hr_rate"),
      col_name_unit_name = "ori9",
      name_treated_unit = "PAPEP0000",
      covagg = list(
        basic_homicide = list(
          label = "Basic Homicide",
          var = c("num_homicide", "hr_rate"),
          each_period = TRUE
        ),
        extended_homicide = list(
          label = "Extended Homicide", 
          var = c("num_homicide", "cleared_cases", "hr_rate"),
          each_period = TRUE
        ),
        avg_homicide = list(
          label = "Average Homicide",
          var = "hr_rate",
          average = "full_pre"
        )
      ),
      treated_period = 2015,
      min_period = 2010,
      end_period = 2019,
      col_name_period = "year",
      feature_weights = c("uniform"),
      donor_sample = c("all", "most_similar"),
      outcome_models = c("none", "augsynth", "lasso", "ridge", "ols"),
      constraints = list(
        list(name = "simplex"),
        list(name = "lasso")
      )
    ),
    
    ois = list(
      outcomes = c("oiso", "oisp"),
      col_name_unit_name = "stateid",
      name_treated_unit = "CA",
      covagg = list(
        basic_ois = list(
          label = "Basic OIS",
          var = c("oisp", "oiso"),
          each_period = TRUE
        ),
        extended_ois = list(
          label = "Extended OIS",
          var = c("oisp", "oiso", "ofp"),
          each_period = TRUE
        ),
        avg_ois = list(
          label = "Average OIS",
          var = "oiso", 
          average = "full_pre"
        )
      ),
      treated_period = 2020,
      min_period = 2015,
      end_period = 2022,
      col_name_period = "year",
      feature_weights = c("uniform", "optimize"),
      donor_sample = c("all", "most_similar"),
      outcome_models = c("none", "augsynth", "lasso", "ridge", "ols"),
      constraints = list(
        list(name = "simplex"),
        list(name = "lasso")
      )
    )
  )
  
  # Get base configuration
  config <- base_configs[[analysis_type]]
  
  # Override with custom parameters
  if (length(custom_params) > 0) {
    for (param_name in names(custom_params)) {
      config[[param_name]] <- custom_params[[param_name]]
    }
  }
  
  return(config)
}