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
          operations = list(
            list(var = "num_homicide", partition_periods = list(type = "by_period")),
            list(var = "hr_rate", partition_periods = list(type = "by_period"))
          )
        ),
        extended_homicide = list(
          label = "Extended Homicide", 
          operations = list(
            list(var = "num_homicide", partition_periods = list(type = "by_period")),
            list(var = "cleared_cases", partition_periods = list(type = "by_period")),
            list(var = "hr_rate", partition_periods = list(type = "by_period"))
          )
        ),
        avg_homicide = list(
          label = "Average Homicide",
          operations = list(
            list(var = "hr_rate", compute = "mean")
          )
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
          operations = list(
            list(var = "oisp", partition_periods = list(type = "by_period")),
            list(var = "oiso", partition_periods = list(type = "by_period"))
          )
        ),
        extended_ois = list(
          label = "Extended OIS",
          operations = list(
            list(var = "oisp", partition_periods = list(type = "by_period")),
            list(var = "oiso", partition_periods = list(type = "by_period")),
            list(var = "ofp", partition_periods = list(type = "by_period"))
          )
        ),
        avg_ois = list(
          label = "Average OIS",
          operations = list(
            list(var = "oiso", compute = "mean")
          )
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
