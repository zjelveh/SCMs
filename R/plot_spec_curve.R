#' @title Plot Specification Curve for Synthetic Control Method Results
#' @description This function creates a specification curve plot based on the results from multiple Synthetic Control Method (SCM) estimations across different specifications.
#'
#' @param spec_curve_results List. Results from the spec_curve function.
#' @param outcomes Character vector. Names of the outcome variables to plot.
#' @param name_treated_unit Character. Name of the treated unit.
#' @param normalize_outcomes Logical. Whether to normalize outcomes. Default is FALSE.
#' @param rmse_threshold Numeric. Threshold for root mean square error to filter results. Default is Inf.
#' @param file_path_save Character or NA. File path to save the plot. If NA, plot is not saved. Default is NA.
#' @param width Numeric. Width of the saved plot in inches. Default is 6.
#' @param height Numeric. Height of the saved plot in inches. Default is 10.
#'
#' @return A ggplot object representing the specification curve.
#'
#'
#' @export
#'
#' @examples
#' # Example usage (replace with actual example when available)
#' # plot <- plot_spec_curve(spec_curve_results = my_results, outcomes = c("gdp", "unemployment"),
#' #                         name_treated_unit = "California", normalize_outcomes = TRUE)
plot_spec_curve <- function(
    spec_curve_results,
    outcomes,
    name_treated_unit,
    normalize_outcomes=FALSE,
    rmse_threshold=Inf,
    file_path_save=NA,
    width=6,
    height=10
){
  # Copy input to avoid modifying original data
  sc = copy(spec_curve_results)
  
  # Initialize list to store results
  sc_results_list = list()
  
  # Iterate through all specifications and collect results
  for(outcome in names(sc)){
    if(outcome %in% outcomes){
      for(const in names(sc[[outcome]])){
        for(fw in names(sc[[outcome]][[const]])){
          for(feat in names(sc[[outcome]][[const]][[fw]])){
            for(ds in names(sc[[outcome]][[const]][[fw]][[feat]])){
              for(ny in names(sc[[outcome]][[const]][[fw]][[feat]][[ds]])){
                # Get the current result
                model_specification = sc[[outcome]][[const]][[fw]][[feat]][[ds]][[ny]]
                
                # Skip if NULL
                if (is.null(model_specification)) next
                
                tryCatch({
                  # Check for the new structure
                  if (!is.null(model_specification$infer) && 
                      is.list(model_specification$infer) && 
                      "inference.results" %in% names(model_specification$infer) &&
                      "rmse" %in% names(model_specification$infer)) {
                    
                    # Get treatment effects data
                    tau_data <- model_specification$infer$inference.results
                    
                    # Get RMSE data
                    rmse_data <- model_specification$infer$rmse
                    
                    # Make sure they're data.tables
                    if (!is.data.table(tau_data)) tau_data <- as.data.table(tau_data)
                    if (!is.data.table(rmse_data)) rmse_data <- as.data.table(rmse_data)
                    
                    # Merge the datasets on the common keys
                    combined <- merge(
                      tau_data,
                      rmse_data,
                      by = c("unit_name", "unit_type", "outcome_model"),
                      all.x = TRUE
                    )
                    
                    estee = model_specification$estimate
                    combined$outcome = outcome
                    combined$const = const
                    combined$fw = fw
                    combined$data_sample = ds
                    combined$feat = feat
                    combined$num_pre_period_years = ny
                    
                    # Use pre_rmse from the merged data instead of calculating it
                    # But keep the old way as a fallback
                    if ("pre_rmse" %in% names(combined)) {
                      combined$rmse <- combined$pre_rmse
                    } else {
                      combined$rmse = sqrt(mean((estee$data$Y.pre-estee$est.results$Y.pre.fit)^2))
                    }
                    
                    # Normalize outcomes if requested
                    if(normalize_outcomes){
                      treat_period = estee$treated_period
                      min_period = estee$min_period
                      col_name_period = estee$col_name_period
                      sd_outcome = sd(estee$data$Y.donors)
                      combined[, tau:=tau/sd_outcome] 
                    }
                    
                    sc_results_list[[length(sc_results_list) + 1]] = combined
                  } else {
                    # Old structure fallback
                    # Try to handle old structure if the new one isn't found
                    combined = rbindlist(model_specification$infer)
                    
                    estee = model_specification$estimate
                    combined$outcome = outcome
                    combined$const = const
                    combined$fw = fw
                    combined$data_sample = ds
                    combined$feat = feat
                    combined$num_pre_period_years = ny
                    combined$rmse = sqrt(mean((estee$data$Y.pre-estee$est.results$Y.pre.fit)^2))
                    
                    # Normalize outcomes if requested
                    if(normalize_outcomes){
                      treat_period = estee$treated_period
                      min_period = estee$min_period
                      col_name_period = estee$col_name_period
                      sd_outcome = sd(estee$data$Y.donors)
                      combined[, tau:=tau/sd_outcome] 
                    }
                    
                    sc_results_list[[length(sc_results_list) + 1]] = combined
                  }
                }, error = function(e) {
                  warning(sprintf("Error processing specification %s/%s/%s/%s/%s/%s: %s", 
                                  outcome, const, fw, feat, ds, ny, e$message))
                })
              } 
            } 
          }
        }
      }
    }  
  }
  
  # Combine results and filter based on RMSE threshold
  if (length(sc_results_list) == 0) {
    stop("No valid results to plot. Check your results structure or adjust your thresholds.")
  }
  
  sc_results_df = rbindlist(sc_results_list, fill = TRUE)
  sc_results_df = sc_results_df[rmse<rmse_threshold]
  
  if(nrow(sc_results_df)==0){
    stop("Not enough specifications to plot after RMSE filtering")
  }
  
  # Check if post_period column exists, if not create it
  if (!"post_period" %in% names(sc_results_df)) {
    warning("post_period column missing - this may indicate an issue with donor_sample='all'. Creating post_period column based on time periods.")
    
    # Create post_period column based on the time column and treated period
    # We need to infer the treated period from the data
    if ("time" %in% names(sc_results_df)) {
      # Find the treated period from the first specification's data
      treated_periods <- unique(sc_results_df[unit_type == "treated"]$time)
      if (length(treated_periods) > 0) {
        treated_start <- min(treated_periods)
        sc_results_df[, post_period := time >= treated_start]
      } else {
        # Fallback: assume post-period is latter half of time periods
        all_times <- sort(unique(sc_results_df$time))
        mid_point <- all_times[ceiling(length(all_times)/2)]
        sc_results_df[, post_period := time >= mid_point]
      }
    } else {
      stop("Cannot create post_period column - no time variable found")
    }
  }
  
  # Calculate average effects AND keep RMSE values
  average_effect_df = sc_results_df[
    post_period==TRUE,
    .(tau=mean(tau), rmse=mean(rmse)), # Keep RMSE by taking mean per specification
    by=c('unit_name', 'outcome', 'outcome_model', 'data_sample',
         'const', 'unit_type', 'fw', 'feat', 'num_pre_period_years')]
  
  # Sort and number specifications
  df_treated_sorted = average_effect_df[unit_name==name_treated_unit]
  df_treated_sorted = df_treated_sorted[order(tau)]
  df_treated_sorted[, new_spec_number:=1:.N]
  
  # Keep these columns from the sorted data
  keep_cols = c('outcome', 'outcome_model', 'const', 'fw', 'feat', 'data_sample', 
                'num_pre_period_years', 'new_spec_number')
  df_treated_sorted = df_treated_sorted[, ..keep_cols]
  
  # Merge and prepare final dataset for plotting
  match_cols = c('outcome', 'outcome_model', 'const', 'fw', 'feat', 'data_sample', 'num_pre_period_years')
  average_effect_df_all = merge(average_effect_df, df_treated_sorted, by=match_cols)
  average_effect_df_all = average_effect_df_all[order(new_spec_number,-unit_type)]
  
  # Rename and recode variables for clarity
  average_effect_df_all[const=='simplex', const:="SCM Weights - Original"]
  average_effect_df_all[const=='lasso', const:="SCM Weights - Penalty Lasso"]
  average_effect_df_all[const=='ridge', const:="SCM Weights - Penalty Ridge"]
  average_effect_df_all[const=='ols', const:="OLS (unconstrated) Weights"]
  

  # Use setnames individually to avoid length mismatch errors
  setnames(average_effect_df_all, 'const', 'Weight_Method')
  setnames(average_effect_df_all, 'tau', 'Estimate')
  setnames(average_effect_df_all, 'outcome', 'Outcome')
  setnames(average_effect_df_all, 'outcome_model', 'Outcome\nModel')
  setnames(average_effect_df_all, 'Weight_Method', 'Weight\nMethod')
  setnames(average_effect_df_all, 'fw', 'V Weights')
  setnames(average_effect_df_all, 'feat', 'Matching\nVars')
  setnames(average_effect_df_all, 'data_sample', 'Donor\nPool')
  setnames(average_effect_df_all, 'num_pre_period_years', 'Pre Length')
  setnames(average_effect_df_all, 'unit_name', 'Unit Name')
  setnames(average_effect_df_all, 'unit_type', 'Unit Type')
  setnames(average_effect_df_all, 'rmse', 'RMSE')
  setnames(average_effect_df_all, 'new_spec_number', 'Specification')
  
  average_effect_df_all[, Type:='Average']
  
  estimates = copy(average_effect_df_all)
  
  spec_cols = c('Outcome', 'Outcome\nModel',
                'Weight\nMethod', 'V Weights', 'Matching\nVars', 
                'Donor\nPool', 'Pre Length')
  
  df = copy(estimates[!is.na(Estimate)])
  var = df$Estimate
  
  group = NULL
  choices = c(spec_cols)
  desc = FALSE
  null = 0
  
  value <- key <- NULL
  var <- enquo(var)
  group <- enquo(group)
  
  
  # Clean up variable names
  df$`Matching\nVars` = gsub('every_period__', '', df$`Matching\nVars`)
  df$`Matching\nVars` = gsub('c\\(|\\)', '', df$`Matching\nVars`)
  
  df2 = df[!duplicated(Specification)]  
  var = df2$Estimate
  
  
  # Create specification plot
  p2 = df2 %>%
    format_results(var = var, group = group, null = null, desc = desc) %>%
    tidyr::gather(key, value, choices) %>%
    dplyr::mutate(key = factor(.data$key, levels = choices)) %>%
    ggplot(aes(x = .data$Specification,
               y = .data$value,
               color = .data$color)) +
    geom_point(aes(x = .data$Specification,
                   y = .data$value),
               shape = 124,
               size = 3.35) +
    scale_color_identity() +
    theme_minimal() +
    facet_grid(.data$key~1, scales = "free_y", space = "free_y") +
    theme(
      axis.line = element_line("black", size = .5),
      legend.position = "none",
      panel.spacing = unit(.75, "lines"),
      axis.text = element_text(colour = "black"),
      strip.text.x = element_blank(),
      strip.text.y = element_text(angle = 90, hjust = 0.5)) +
    labs(x = "", y = "")
  
  var = df$Estimate
  
  # Create main effect plot
  p1inf <- df %>%
    format_results(var = var, group = group, null = null, desc = desc) %>%
    ggplot(aes(x = .data$RMSE,
               y = .data$Estimate,
               fill = `Unit Type`,
               color = `Unit Type`)
           ) + 
    geom_point(aes(alpha = `Unit Type`)) +
    theme_minimal() +
    scale_alpha_manual(values = c("treated" = 1, "control" = .25)) +
    scale_size_continuous(name = "Model Fit", 
                          guide = guide_legend(title = "Pre-period Fit",
                                               reverse = TRUE,
                                               override.aes = list(alpha = 1))) +
    geom_hline(yintercept=0, alpha=.5, linetype='dashed') + 
    theme(strip.text = element_blank(),
          axis.line = element_line("black", size = .5),
          legend.position = "right",
          panel.spacing = unit(.75, "lines"),
          axis.text = element_text(colour = "black")) +
    labs(x = "RMSE", y='Average Annual Treatment Effect')
  
  # Create main effect plot with point size varying by RMSE
  p1 <- df %>%
    format_results(var = var, group = group, null = null, desc = desc) %>%
    ggplot(aes(x = .data$Specification,
               y = .data$Estimate,
               fill = `Unit Type`,
               color = `Unit Type`,
               size = (.data$RMSE))) +  # Inverse RMSE for size (larger = better fit)
    geom_point(aes(alpha = `Unit Type`)) +
    theme_minimal() +
    scale_alpha_manual(values = c("treated" = 1, "control" = .25)) +
    scale_size_continuous(name = "Model Fit", 
                          guide = guide_legend(title = "Pre-period Fit",
                                               reverse = TRUE,
                                               override.aes = list(alpha = 1))) +
    geom_hline(yintercept=0, alpha=.5, linetype='dashed') + 
    theme(strip.text = element_blank(),
          axis.line = element_line("black", size = .5),
          legend.position = "right",
          panel.spacing = unit(.75, "lines"),
          axis.text = element_text(colour = "black")) +
    labs(x = "Specification #", y='Average Annual Treatment Effect')
  
  
  # Combine plots
  p = cowplot::plot_grid(p1,
                         p2,
                         labels = c('A','B'),
                         align = "v",
                         axis = "rbl",
                         rel_heights = c(.5,.5),
                         ncol = 1
  )
  
  # Save plot if file path is provided
  if(!is.na(file_path_save)){
    ggsave(file_path_save, plot = p, width=width, height = height)
  }
  
  return(p)
}