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
    rmse_threshold= Inf,
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
                model_specification = sc[[outcome]][[const]][[fw]][[feat]][[ds]][[ny]]
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
            } 
          }
        }
      }
    }  
  }
  
  # Combine results and filter based on RMSE threshold
  sc_results_df = rbindlist(sc_results_list)
  sc_results_df = sc_results_df[rmse<rmse_threshold]

  if(nrow(sc_results_df)==0){
    stop("Not enough specifications to plot")
  }

  # Calculate average effects
  average_effect_df = sc_results_df[
    post_period==TRUE,
    .(tau=mean(tau)),
    by=c('unit_name', 'outcome', 'outcome_model', 'data_sample',
         'const', 'unit_type', 'fw', 'feat', 'num_pre_period_years')]
  
  # Sort and number specifications
  df_treated_sorted = average_effect_df[unit_name==name_treated_unit]
  df_treated_sorted = df_treated_sorted[order(tau)]
  df_treated_sorted[, new_spec_number:=1:.N]
  df_treated_sorted[, tau:=NULL]
  df_treated_sorted[, unit_name:=NULL]
  df_treated_sorted[, unit_type:=NULL]
  
  # Merge and prepare final dataset for plotting
  match_cols = c('outcome', 'outcome_model', 'const', 'fw', 'feat', 'data_sample', 'num_pre_period_years')
  average_effect_df_all = merge(average_effect_df, df_treated_sorted, by=match_cols)
  average_effect_df_all = average_effect_df_all[order(new_spec_number,-unit_type)]
  
  # Rename and recode variables for clarity
  average_effect_df_all[const=='simplex', const:="SCM Weights - Original"]
  average_effect_df_all[const=='lasso', const:="SCM Weights - Penalty Lasso"]
  average_effect_df_all[const=='ridge', const:="SCM Weights - Penalty Ridge"]
  average_effect_df_all[const=='ols', const:="OLS (unconstrated) Weights"]
  
  setnames(average_effect_df_all, 'const', 'Weight_Method')
  
  setnames(average_effect_df_all, c('Outcome', 'Outcome\nModel',
                             'Weight\nMethod', 'V Weights', 'Matching\nVars', 
                             'Donor\nPool', 'Pre Length', 'Unit Name', 'Unit Type', 'Estimate', 
                             'Specification'
  ))

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
      strip.text.x = element_blank()) +
    labs(x = "", y = "")
  
  var = df$Estimate
  
  # Create main effect plot
  p1 <- df %>%
    format_results(var = var, group = group, null = null, desc = desc) %>%
    ggplot(aes(x = .data$Specification,
               y = .data$Estimate,
               fill=`Unit Type`,
               color = `Unit Type`)) +
    geom_point(aes(alpha = `Unit Type`),
               size = 1) +
    theme_minimal() +
    scale_alpha_manual(values = c("treated" = 1, "control" = .25)) +
    geom_hline(yintercept=0, alpha=.5, linetype='dashed') + 
    theme(strip.text = element_blank(),
          axis.line = element_line("black", size = .5),
          legend.position = "none",
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