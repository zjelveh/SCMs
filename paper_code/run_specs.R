# Load required packages
library(data.table)
library(readxl)
library(ggplot2)
library(devtools)

# Data preparation functions
library(data.table)
library(readxl)

# Create a properly structured config
config <- list(
  output_dir = "~/Dropbox/research/synth_control_paper/shap_exploration/data/",
  analyses = list(
    homicide = list(
      outcomes = c("num_homicide", 'hr_rate'),
      col_name_unit_name = "ori9",
      name_treated_unit = "PAPEP0000",
      covagg = list(
        list(every_period = c("num_homicide", "hr_rate")),
        list(every_period = c("num_homicide", "cleared_cases", "hr_rate"))
      ),
      treated_period = 2015,
      min_period = 2010,
      end_period = 2019,
      col_name_period = "year",
      feature_weights = c("uniform", "optimize"),
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
        list(every_period = c("oisp", "oiso")),
        list(every_period = c("oisp", "oiso", "ofp"))
      ),
      treated_period = 2020,
      min_period = 2015,
      end_period = 2022,
      col_name_period = "year",
      num_pre_period_years = c(5),
      feature_weights = c("uniform", "optimize"),
      donor_sample = c("all", "most_similar"),
      outcome_models = c("none", "augsynth", "lasso", "ridge", "ols"),
      constraints = list(
        list(name = "simplex"), 
        list(name = "lasso")
      )
    )
  ),
  
  # Add visualization settings
  visualization = list(
    homicide = list(
      rmse_threshold = Inf,
      normalize_outcomes = TRUE
    ),
    ois = list(
      rmse_threshold = Inf,
      normalize_outcomes = TRUE
    )
  )
) 


# Data preparation functions
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



# Analysis functions
run_spec_curve <- function(dataset, outcomes = NULL, params, output_dir = NULL) {
  # Check if outcomes is provided as an argument or in params
  if (is.null(outcomes)) {
    # Use outcomes from params if not explicitly provided
    all_params <- params
    actual_outcomes <- params$outcomes
  } else {
    # If outcomes is explicitly provided, override the one in params
    params_copy <- params
    params_copy$outcomes <- NULL  # Remove outcomes from params to avoid duplication
    
    all_params <- c(
      list(
        dataset = dataset,
        outcomes = outcomes
      ),
      params_copy
    )
    actual_outcomes <- outcomes
  }
  
  # Add dataset if not already in the parameters
  if (!"dataset" %in% names(all_params)) {
    all_params$dataset <- dataset
  }
  
  # Call spec_curve with all parameters
  results <- do.call(spec_curve, all_params)
  
  # Save results if output_dir is provided
  if (!is.null(output_dir)) {
    save_path <- file.path(output_dir, paste0(params$name_treated_unit, "_", 
                                              paste(actual_outcomes, collapse="_"), "_sc.rdata"))
    save(results, file = save_path)
  }
  
  return(results)
}



create_spec_curve_plots <- function(results, viz_config) {
  plots <- list()
  
  for (outcome in names(results)) {
    print(outcome)
    
    # Find the treated unit name - try different methods
    treated_unit <- NULL
    
    # Try to find the first valid specification to get the treated unit name
    for (const in names(results[[outcome]])) {
      for (fw in names(results[[outcome]][[const]])) {
        for (feat in names(results[[outcome]][[const]][[fw]])) {
          for (ds in names(results[[outcome]][[const]][[fw]][[feat]])) {
            for (ny in names(results[[outcome]][[const]][[fw]][[feat]][[ds]])) {
              model_spec <- results[[outcome]][[const]][[fw]][[feat]][[ds]][[ny]]
              if (!is.null(model_spec) && !is.null(model_spec$estimate)) {
                treated_unit <- model_spec$estimate$name_treated_unit
                break
              }
            }
            if (!is.null(treated_unit)) break
          }
          if (!is.null(treated_unit)) break
        }
        if (!is.null(treated_unit)) break
      }
      if (!is.null(treated_unit)) break
    }
    
    if (is.null(treated_unit)) {
      warning(paste("Could not find treated unit name for outcome:", outcome))
      next
    }
    
    plot <- plot_spec_curve(
      results,
      outcomes = outcome,
      name_treated_unit = treated_unit,
      rmse_threshold = viz_config$rmse_threshold,
      normalize_outcomes = viz_config$normalize_outcomes
    )
    
    plots[[outcome]] <- plot
  }
  
  return(plots)
}


extract_results <- function(results) {
  sc_results_list <- list()
  
  # Loop through all model specifications and compile results
  for(outcome in names(results)) {
    for(const in names(results[[outcome]])) {
      for(fw in names(results[[outcome]][[const]])) {
        for(feat in names(results[[outcome]][[const]][[fw]])) {
          for(ds in names(results[[outcome]][[const]][[fw]][[feat]])) {
            for(ny in names(results[[outcome]][[const]][[fw]][[feat]][[ds]])) {
              # Get the current result
              current_result <- results[[outcome]][[const]][[fw]][[feat]][[ds]][[ny]]
              
              # Skip if NULL
              if (is.null(current_result)) next
              
              # Process inference results
              tryCatch({
                # Check if infer component exists and has the right structure
                if (!is.null(current_result$infer) && 
                    is.list(current_result$infer) && 
                    "inference.results" %in% names(current_result$infer) &&
                    "rmse" %in% names(current_result$infer)) {
                  
                  # Get treatment effects data
                  tau_data <- current_result$infer$inference.results
                  
                  # Get RMSE data
                  rmse_data <- current_result$infer$rmse
                  
                  # Check if both datasets exist
                  if (!is.null(tau_data) && !is.null(rmse_data)) {
                    # Make sure they're data.tables
                    if (!is.data.table(tau_data)) tau_data <- as.data.table(tau_data)
                    if (!is.data.table(rmse_data)) rmse_data <- as.data.table(rmse_data)
                    
                    # Merge the datasets on the common keys
                    merged_data <- merge(
                      tau_data,
                      rmse_data,
                      by = c("unit_name", "unit_type", "outcome_model"),
                      all.x = TRUE
                    )
                    
                    # Add metadata
                    merged_data$outcome <- outcome
                    merged_data$const <- const
                    merged_data$fw <- fw
                    merged_data$data_sample <- ds
                    merged_data$feat <- feat
                    merged_data$num_pre_period_years <- ny
                    
                    # Add to results list
                    sc_results_list[[length(sc_results_list) + 1]] <- merged_data
                  }
                }
              }, error = function(e) {
                warning("Error processing result: ", e$message, 
                        " in spec: ", outcome, "/", const, "/", fw, "/", feat, "/", ds, "/", ny)
              })
            }
          }
        }
      }
    }
  }
  
  # Combine all results into a single dataframe
  if (length(sc_results_list) > 0) {
    sc_results_df <- rbindlist(sc_results_list, fill = TRUE)
    return(sc_results_df)
  } else {
    warning("No valid results found to extract")
    return(NULL)
  }
}



save_results <- function(results_df, base_path) {
  # Clean up results if needed
  results_df[, unit_type := NULL]
  
  # Save to CSV
  fwrite(results_df, paste0(base_path, ".csv"))
}

config$analyses$homicide$cores=32
config$analyses$ois$cores=32

# hogan
hogan_path = 'Dropbox/research/synth_control_paper/repo/variation_scm/data/hogan_dataset_processed.csv'
dataset = process_homicide_data(hogan_path)

# Load custom SCM package

load_all('~/Dropbox/research/synth_control_paper/repo/SCMs/SCMs/')

# Then run with the properly structured config
results = run_spec_curve(
  dataset = dataset,
  params = config$analyses$homicide,
  output_dir = config$output_dir
)

save(results, file='~/Dropbox/research/synth_control_paper/shap_exploration/data/hogan_sc.rdata')

extracted_results <- extract_results(results)

save_results(extracted_results, '~/Dropbox/research/synth_control_paper/shap_exploration/data/spec_curve_hogan_results')



###kaplan
homi_path = 'Dropbox/research/synth_control_paper/repo/variation_scm/data/kaplan_dataset_processed.csv'
dataset = process_homicide_data(homi_path)

# Then run with the properly structured config
results = run_spec_curve(
  dataset = dataset,
  params = config$analyses$homicide,
  output_dir = config$output_dir
)

save(results, file='~/Dropbox/research/synth_control_paper/shap_exploration/data/homi_sc.rdata')

extracted_results <- extract_results(results)

save_results(extracted_results, '~/Dropbox/research/synth_control_paper/shap_exploration/data/spec_curve_homi_results')


######OIS
ois_path = '~/Dropbox/research/synth_control_paper/data/OIS_cali_replication_data/fillinpanel3.xlsx'
dataset = process_ois_data(ois_path)

# Then run with the properly structured config
results = run_spec_curve(
  dataset = dataset,
  params = config$analyses$ois,
  output_dir = config$output_dir
)

save(results, file='~/Dropbox/research/synth_control_paper/shap_exploration/data/ois_sc.rdata')

extracted_results <- extract_results(results)

save_results(extracted_results, '~/Dropbox/research/synth_control_paper/shap_exploration/data/spec_curve_ois_results')




viz_config = config$visualization$homicide

save_plots <- function(plots, base_path) {
  for (plot_name in names(plots)) {
    plot <- plots[[plot_name]]
    file_path <- paste0(base_path, "_", plot_name, ".pdf")
    ggsave(file_path, plot, width = 10, height = 8)
  }
}




# Create visualization of specification curve
p <- plot_spec_curve(sc, outcomes=c('oiso', 'oisp'), name_treated_unit = 'CA')
p
print(p)


p <- plot_spec_curve(sc, outcomes=c('num_homicide'), name_treated_unit = 'PAPEP0000')
p

p <- plot_spec_curve(results,normalize_outcomes = TRUE, rmse_threshold = 20,
                     outcomes=c( 'num_homicide'),
                     name_treated_unit = 'PAPEP0000')
p


hogan_results = fread('~/Dropbox/research/synth_control_paper/shap_exploration/data/spec_curve_hogan_results.csv')
kaplan_results = fread('~/Dropbox/research/synth_control_paper/shap_exploration/data/spec_curve_HOMI_results.csv')
ois_results = fread('~/Dropbox/research/synth_control_paper/shap_exploration/data/spec_curve_OIS_results.csv')

names(hogan_results)
names(kaplan_results)
names(ois_results)


