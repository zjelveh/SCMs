#' Plot Shapley Value Distributions
#'
#' @title Visualize SHAP Value Distributions with Treated Unit Highlighted
#' @description Creates histograms of SHAP values across all units with the
#' treated unit's values highlighted as vertical lines.
#'
#' @param xgb_results List. Results from \code{run_xgboost_shap_analysis}.
#' @param metric Character. Which SHAP metric to plot ("mean_shap", "sd_shap", or "mean_abs_shap").
#' @param ncol Integer. Number of columns for facet layout.
#'
#' @return ggplot2 object showing SHAP value distributions.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Plot mean absolute SHAP values
#' p <- plot_shapley_distributions(xgb_results, metric = "mean_abs_shap")
#' print(p)
#' }
plot_shapley_distributions <- function(xgb_results, metric = "mean_abs_shap", ncol = 2) {
  # Extract results data table
  results_dt <- xgb_results$results
  
  # Get the treated unit name from config
  treated_unit_name <- xgb_results$config$treated_unit_name
  
  if(nrow(results_dt) == 0 || is.null(results_dt)) {
    stop("No results available to plot")
  }
  
  # Get values for the treated unit
  treated_unit_values <- results_dt[unit == treated_unit_name, ]
  
  if(nrow(treated_unit_values) == 0) {
    warning(paste("Treated unit", treated_unit_name, "not found in results"))
    return(NULL)
  }
  
  # Create a data frame for the vertical lines
  if(metric == "mean_shap") {
    vertical_lines <- data.table(
      feature = treated_unit_values$feature,
      value = treated_unit_values$mean_shap
    )
    x_label <- "Mean SHAP Value"
  } else if(metric == "sd_shap") {
    vertical_lines <- data.table(
      feature = treated_unit_values$feature,
      value = treated_unit_values$sd_shap
    )
    x_label <- "Standard Deviation of SHAP Values"
  } else {
    vertical_lines <- data.table(
      feature = treated_unit_values$feature,
      value = treated_unit_values$mean_abs_shap
    )
    x_label <- "Mean Absolute SHAP Value"
  }
  
  # Create the plot
  p <- ggplot(results_dt, aes_string(x = metric)) +
    geom_histogram(bins = 30, fill = "skyblue", color = "darkblue", alpha = 0.7) + 
    # Add vertical lines for treated unit values
    geom_vline(data = vertical_lines, 
               aes(xintercept = value), 
               color = "red", 
               linetype = "dashed", 
               size = 1) +
    # Add labels for the vertical lines
    geom_text(data = vertical_lines,
              aes(x = value, y = 0, label = treated_unit_name),
              color = "red", 
              angle = 90, 
              vjust = -0.5,
              hjust = -0.1,
              size = 3) +
    facet_wrap(~feature, scales = "free", ncol = ncol) +
    labs(title = paste("Distribution of", x_label, "with Treated Unit Values Highlighted"),
         subtitle = paste("Treated Unit:", treated_unit_name),
         x = x_label,
         y = "Count") +
    theme_minimal() +
    theme(
      strip.background = element_rect(fill = "lightgray", color = NA),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  return(p)
}

#' Calculate Specification Interactions  
#'
#' @title Calculate Variance Explained by Specification Feature Interactions
#' @description Computes how much variance in treatment effects is explained
#' by individual specification features and their pairwise interactions.
#'
#' @param xgb_results List. Results from \code{run_xgboost_shap_analysis}.
#' @param dataset_name Character. Name for the dataset (for labeling).
#'
#' @return List containing:
#'   \itemize{
#'     \item \code{individual_summary} - Variance explained by individual features
#'     \item \code{interaction_summary} - Variance explained by feature pairs
#'     \item \code{detailed_effects} - Detailed effects for each feature combination
#'     \item \code{dataset} - Dataset name
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Calculate specification interactions
#' interactions <- calculate_specification_interactions(xgb_results, "homicide_study")
#' }
calculate_specification_interactions <- function(xgb_results, dataset_name) {
  # Check if results exist
  if(is.null(xgb_results) || is.null(xgb_results$predictions) || 
     nrow(xgb_results$predictions) == 0) {
    warning("Missing prediction data for dataset: ", dataset_name)
    return(NULL)
  }

  # Get treated unit name
  treated_unit <- xgb_results$config$treated_unit_name

  # Extract predictions data
  preds_data <- xgb_results$predictions

  # Get specification features
  spec_features <- c("outcome_model", "const", "fw", "feat", "data_sample")
  available_features <- intersect(spec_features, names(preds_data))

  # Process interactions for the treated unit
  treated_data <- preds_data[unit == treated_unit]

  if(nrow(treated_data) == 0) {
    warning("No data for treated unit in dataset: ", dataset_name)
    return(NULL)
  }

  # Compute interaction effect sizes
  # For each pair of features
  interaction_dt <- data.table()

  for(i in 1:(length(available_features)-1)) {
    feat_i <- available_features[i]
    
    for(j in (i+1):length(available_features)) {
      feat_j <- available_features[j]
      
      # Get unique combinations of these features
      combos <- unique(treated_data[, c(feat_i, feat_j), with=FALSE])
      
      if(nrow(combos) <= 1) {
        # Skip if there's only one unique combination
        next
      }
      
      # Calculate effect for each combination
      combo_effects <- data.table()
      
      for(row in 1:nrow(combos)) {
        val_i <- combos[[feat_i]][row]
        val_j <- combos[[feat_j]][row]
        
        # Get effects for this combination - using subsetting for safety
        subset_data <- treated_data[treated_data[[feat_i]] == val_i & 
                                      treated_data[[feat_j]] == val_j, ]
        
        if(nrow(subset_data) > 0) {
          effects <- data.table(
            mean_effect = mean(subset_data$actual, na.rm=TRUE),
            sd_effect = sd(subset_data$actual, na.rm=TRUE),
            count = nrow(subset_data)
          )
          
          combo_effects <- rbind(combo_effects, data.table(
            feature1 = feat_i,
            feature2 = feat_j,
            value1 = val_i,
            value2 = val_j,
            mean_effect = effects$mean_effect,
            sd_effect = effects$sd_effect,
            count = effects$count
          ))
        }
      }
      
      # Calculate overall variance explained by this interaction
      if(nrow(combo_effects) > 1) {
        # Calculate total variance
        total_variance <- var(treated_data$actual, na.rm=TRUE)
        
        # Calculate weighted between-group variance for this interaction
        mean_overall <- mean(treated_data$actual, na.rm=TRUE)
        weighted_between_var <- sum(combo_effects$count * 
                                      (combo_effects$mean_effect - mean_overall)^2) / 
          sum(combo_effects$count)
        
        # Add interaction summary
        interaction_dt <- rbind(interaction_dt, data.table(
          feature1 = feat_i,
          feature2 = feat_j,
          variance_explained = weighted_between_var,
          variance_explained_pct = weighted_between_var / total_variance * 100,
          total_variance = total_variance,
          n_combinations = nrow(combo_effects)
        ))
      }
    }
  }

  # Also compute individual feature variance explained
  individual_dt <- data.table()
  detailed_effects <- list()

  for(feature in available_features) {
    # Skip if feature doesn't exist in data
    if(!(feature %in% names(treated_data))) {
      warning(paste("Feature column", feature, "not found in treated_data"))
      next
    }
    
    # Use data.table's by capability to safely compute by group
    value_effects <- treated_data[, .(
      mean_effect = mean(actual, na.rm = TRUE),
      sd_effect = sd(actual, na.rm = TRUE),
      count = .N
    ), by = feature]
    
    # Only proceed if we have more than one unique value
    if(nrow(value_effects) > 1) {
      # Rename the grouping column to 'value' for consistency
      setnames(value_effects, feature, "value")
      
      # Add feature column
      value_effects[, feature := feature]
      
      # Calculate total variance
      total_variance <- var(treated_data$actual, na.rm=TRUE)
      
      # Calculate weighted between-group variance
      mean_overall <- mean(treated_data$actual, na.rm=TRUE)
      weighted_between_var <- sum(value_effects$count * 
                                    (value_effects$mean_effect - mean_overall)^2) / 
        sum(value_effects$count)
      
      # Add summary
      individual_dt <- rbind(individual_dt, data.table(
        feature = feature,
        variance_explained = weighted_between_var,
        variance_explained_pct = weighted_between_var / total_variance * 100,
        total_variance = total_variance,
        n_values = nrow(value_effects)
      ))
      
      # Store detailed value data
      detailed_effects[[feature]] <- value_effects
    }
  }

  return(list(
    individual_summary = individual_dt,
    interaction_summary = interaction_dt,
    detailed_effects = detailed_effects,
    dataset = dataset_name
  ))
}

#' Plot Specification Interaction Heatmap
#'
#' @title Visualize Feature Interactions as Heatmap
#' @description Creates a heatmap showing variance explained by specification
#' feature interactions, with individual feature effects on the diagonal.
#'
#' @param interaction_data List. Results from \code{calculate_specification_interactions}.
#'
#' @return ggplot2 object showing interaction heatmap.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create interaction heatmap
#' p <- plot_interaction_heatmap(interaction_results)
#' print(p)
#' }
plot_interaction_heatmap <- function(interaction_data) {
  if(is.null(interaction_data) || is.null(interaction_data$interaction_summary) || 
     nrow(interaction_data$interaction_summary) == 0) {
    return(ggplot() + labs(title = paste("No interaction data for", 
                                         interaction_data$dataset)))
  }

  # Create symmetrical interaction matrix
  int_matrix <- interaction_data$interaction_summary

  # Create complete matrix with both feature1-feature2 and feature2-feature1
  reversed <- int_matrix[, .(
    feature1 = feature2,
    feature2 = feature1,
    variance_explained = variance_explained,
    variance_explained_pct = variance_explained_pct
  )]

  full_matrix <- rbind(int_matrix, reversed, fill=TRUE)

  # Add individual feature effects for diagonal
  if(!is.null(interaction_data$individual_summary) && 
     nrow(interaction_data$individual_summary) > 0) {
    diag_entries <- interaction_data$individual_summary[, .(
      feature1 = feature,
      feature2 = feature,
      variance_explained = variance_explained,
      variance_explained_pct = variance_explained_pct
    )]
    
    full_matrix <- rbind(full_matrix, diag_entries, fill=TRUE)
  }

  # Create heatmap
  # Get unique features in order of importance
  features <- unique(c(full_matrix$feature1, full_matrix$feature2))

  # Reorder to put most important features first
  if(!is.null(interaction_data$individual_summary) && 
     nrow(interaction_data$individual_summary) > 0) {
    ordered_features <- interaction_data$individual_summary[
      order(-variance_explained_pct), feature]
    features <- union(ordered_features, features)  # preserves order of ordered_features
  }

  # Create the heatmap
  p <- ggplot(full_matrix, aes(x = factor(feature1, levels = features), 
                               y = factor(feature2, levels = rev(features)), 
                               fill = variance_explained_pct)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = round(variance_explained_pct, 1)), 
              color = "black", size = 3) +
    scale_fill_gradient2(low = "white", mid = "lightblue", high = "darkblue",
                         midpoint = max(full_matrix$variance_explained_pct, na.rm=TRUE) / 2,
                         name = "% Variance\nExplained") +
    labs(title = paste("Specification Feature Interactions:", interaction_data$dataset),
         subtitle = "Diagonal = Individual Effects, Off-diagonal = Interaction Effects",
         x = "Feature 1", y = "Feature 2") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))

  return(p)
}

#' Save Visualization Plots
#'
#' @title Save Multiple Plots to Files
#' @description Convenience function to save multiple ggplot objects to PDF files.
#'
#' @param plots Named list of ggplot objects.
#' @param base_path Character. Base path for output files (without extension).
#' @param width Numeric. Plot width in inches.
#' @param height Numeric. Plot height in inches.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Save plots
#' plots <- list(
#'   shap_dist = plot_shapley_distributions(results),
#'   interactions = plot_interaction_heatmap(interaction_data)
#' )
#' save_visualization_plots(plots, "output/analysis_plots")
#' }
save_visualization_plots <- function(plots, base_path, width = 10, height = 8) {
  for (plot_name in names(plots)) {
    plot <- plots[[plot_name]]
    file_path <- paste0(base_path, "_", plot_name, ".pdf")
    ggsave(file_path, plot, width = width, height = height)
  }
}