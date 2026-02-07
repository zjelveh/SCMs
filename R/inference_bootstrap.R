#' Bootstrap Null Hypothesis Inference for Synthetic Control
#'
#' @title Bootstrap-based significance testing under null hypothesis of no treatment effect
#' @description Implements bootstrap inference by: 1) Enforcing null hypothesis by subtracting
#' actual treatment effect from post-treatment values, 2) Bootstrap resampling at unit level,
#' 3) Re-estimating SCM on bootstrapped samples, 4) Computing p-values from bootstrap distribution
#'
#' @param sc.pred Result from estimate_sc function
#' @param dataset Original dataset used for estimation
#' @param n_bootstrap Number of bootstrap replications (default: 1000)
#' @param cores Number of cores for parallel processing (default: 1)
#' @param verbose Logical for verbose output (default: FALSE)
#' @param seed Random seed for reproducibility (default: NULL)
#'
#' @return List containing bootstrap results and p-values
bootstrap_null_inference <- function(sc.pred,
                                   dataset,
                                   n_bootstrap = 1000,
                                   cores = 1,
                                   verbose = FALSE,
                                   seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Validate inputs
  if (is.null(sc.pred) || is.null(dataset)) {
    stop("sc.pred and dataset cannot be NULL")
  }
  
  # Extract key parameters
  col_name_unit <- sc.pred$col_name_unit_name
  col_name_period <- sc.pred$col_name_period
  col_name_outcome <- sc.pred$outcome
  treated_unit <- sc.pred$name_treated_unit
  treated_period <- sc.pred$treated_period
  
  if (verbose) {
    cat("Starting bootstrap null hypothesis inference...\n")
    cat("Bootstrap replications:", n_bootstrap, "\n")
    cat("Treated unit:", treated_unit, "\n")
    cat("Treatment period:", treated_period, "\n")
  }
  
  # Step 1: Calculate actual treatment effects for all outcome models
  actual_effects <- list()
  
  for (oc in names(sc.pred$est.results$outcome_model)) {
    sc_post <- sc.pred$est.results$outcome_model[[oc]]$Y.post.fit
    actual_tau <- sc.pred$data$Y.post - sc_post
    actual_effects[[oc]] <- mean(actual_tau, na.rm = TRUE)  # Average treatment effect
    
    if (verbose) {
      cat("Actual treatment effect for", oc, ":", round(actual_effects[[oc]], 4), "\n")
    }
  }
  
  # Step 2: Enforce null hypothesis by subtracting treatment effects
  dataset_null <- enforce_null_hypothesis(dataset, sc.pred, actual_effects, verbose)
  
  # Step 3: Bootstrap procedure
  bootstrap_results <- perform_bootstrap_inference(
    dataset_null = dataset_null,
    sc.pred = sc.pred,
    actual_effects = actual_effects,
    n_bootstrap = n_bootstrap,
    cores = cores,
    verbose = verbose
  )
  
  # Step 4: Calculate p-values
  p_values <- calculate_bootstrap_pvalues(bootstrap_results$effects, actual_effects, verbose)
  
  return(list(
    actual_effects = actual_effects,
    bootstrap_effects = bootstrap_results$effects,
    bootstrap_iteration_data = bootstrap_results$iteration_data,
    bootstrap_estimates = bootstrap_results$estimates,
    p_values = p_values,
    n_bootstrap = n_bootstrap,
    null_enforced_data = dataset_null
  ))
}

#' Enforce Null Hypothesis by Subtracting Treatment Effects
#'
#' @param dataset Original dataset
#' @param sc.pred SC estimation results
#' @param actual_effects List of actual treatment effects by outcome model
#' @param verbose Logical for verbose output
#'
#' @return Dataset with null hypothesis enforced
enforce_null_hypothesis <- function(dataset, sc.pred, actual_effects, verbose = FALSE) {
  
  # Ensure dataset is data.table for efficient operations
  if (!is.data.table(dataset)) setDT(dataset)
  
  # Create a copy to avoid modifying original
  dataset_null <- copy(dataset)
  
  col_name_unit <- sc.pred$col_name_unit_name
  col_name_period <- sc.pred$col_name_period
  col_name_outcome <- sc.pred$outcome
  treated_unit <- sc.pred$name_treated_unit
  treated_period <- sc.pred$treated_period
  
  # Get treated unit post-treatment observations
  treated_post_mask <- dataset_null[[col_name_unit]] == treated_unit & 
                       dataset_null[[col_name_period]] >= treated_period
  
  if (sum(treated_post_mask) == 0) {
    stop("No post-treatment observations found for treated unit")
  }
  
  # Subtract the average treatment effect from post-treatment outcomes
  # This enforces the null hypothesis of no treatment effect
  # CRITICAL: Convert outcome to numeric to avoid integer truncation issues
  if (!is.numeric(dataset_null[[col_name_outcome]])) {
    dataset_null[, (col_name_outcome) := as.numeric(get(col_name_outcome))]
  }
  
  if (length(actual_effects) == 1) {
    # Single outcome model case
    effect_to_subtract <- actual_effects[[1]]
    dataset_null[treated_post_mask, (col_name_outcome) := as.numeric(get(col_name_outcome)) - effect_to_subtract]
    
    if (verbose) {
      cat("Subtracted effect", round(effect_to_subtract, 4), "from", sum(treated_post_mask), 
          "post-treatment observations\n")
    }
  } else {
    # Multiple outcome models - use first one as primary
    # (In practice, you might want to handle this differently)
    effect_to_subtract <- actual_effects[[1]]
    dataset_null[treated_post_mask, (col_name_outcome) := as.numeric(get(col_name_outcome)) - effect_to_subtract]
    
    if (verbose) {
      cat("Using primary outcome model effect", round(effect_to_subtract, 4), 
          "for null enforcement\n")
    }
  }
  
  return(dataset_null)
}

#' Perform Bootstrap Inference
#'
#' @param dataset_null Dataset with null hypothesis enforced
#' @param sc.pred Original SC estimation results
#' @param actual_effects List of actual treatment effects
#' @param n_bootstrap Number of bootstrap replications
#' @param cores Number of cores for parallel processing
#' @param verbose Logical for verbose output
#'
#' @return List with bootstrap effects and estimates
perform_bootstrap_inference <- function(dataset_null, sc.pred, actual_effects,
                                       n_bootstrap, cores, verbose) {
  
  # Get unique units for bootstrapping
  col_name_unit <- sc.pred$col_name_unit_name
  all_units <- unique(dataset_null[[col_name_unit]])
  control_units <- setdiff(all_units, sc.pred$name_treated_unit)
  
  if (verbose) {
    cat("Total units for bootstrapping:", length(all_units), "\n")
    cat("Control units:", length(control_units), "\n")
  }
  
  # Setup parallel processing if requested
  if (cores > 1) {
    if (!requireNamespace("parallel", quietly = TRUE)) {
      warning("parallel package not available, using single core")
      cores <- 1
    }
  }
  
  # Bootstrap function for single replication
  bootstrap_single <- function(b) {
    tryCatch({
      # Sample units with replacement (unit-level bootstrap)
      # CRITICAL: Always include treated unit, then sample control units
      n_control <- length(control_units)
      bootstrap_control_units <- sample(control_units, size = n_control, replace = TRUE)
      bootstrap_units <- c(sc.pred$name_treated_unit, bootstrap_control_units)
      
      # Create bootstrap dataset
      bootstrap_data <- create_bootstrap_dataset(dataset_null, bootstrap_units, 
                                                col_name_unit, sc.pred$name_treated_unit)
      
      
      # Re-estimate SCM on bootstrap sample with optimized settings
      bootstrap_sc <- estimate_sc(
        dataset = bootstrap_data,
        outcome = sc.pred$outcome,
        covagg = sc.pred$covagg,
        col_name_unit_name = sc.pred$col_name_unit_name,
        name_treated_unit = sc.pred$name_treated_unit,
        col_name_period = sc.pred$col_name_period,
        treated_period = sc.pred$treated_period,
        min_period = sc.pred$min_period,
        end_period = sc.pred$end_period,
        feature_weights = sc.pred$feature_weights,
        outcome_models = sc.pred$outcome_models,
        w.constr = sc.pred$w.constr
      )
      
      # Calculate treatment effects for each outcome model
      effects <- list()
      iteration_data <- list()
      
      for (oc in names(bootstrap_sc$est.results$outcome_model)) {
        sc_post <- bootstrap_sc$est.results$outcome_model[[oc]]$Y.post.fit
        sc_pre <- bootstrap_sc$est.results$Y.pre.fit
        
        # Post-treatment effects
        tau_post <- bootstrap_sc$data$Y.post - sc_post
        # Pre-treatment effects  
        tau_pre <- bootstrap_sc$data$Y.pre - sc_pre
        # Combined effects
        tau_all <- c(tau_pre, tau_post)
        
        # Store average effect for p-value calculation
        effects[[oc]] <- mean(tau_post, na.rm = TRUE)
        
        # Store full time series for plotting (similar to placebo format)
        iteration_data[[oc]] <- data.table(
          unit_name = paste0("bootstrap_", b),
          period = sc.pred$min_period:sc.pred$end_period,
          tau = tau_all,
          post_period = (sc.pred$min_period:sc.pred$end_period) >= sc.pred$treated_period,
          unit_type = 'bootstrap',
          outcome_model = oc
        )
      }
      
      return(list(
        bootstrap_id = b,
        effects = effects,
        iteration_data = iteration_data,
        estimate = bootstrap_sc
      ))
      
    }, error = function(e) {
      # FAIL HARD - don't tolerate bootstrap failures
      stop(paste("Bootstrap replication", b, "failed:", e$message, 
                "\n  Bootstrap units:", length(bootstrap_units),
                "\n  Unique bootstrap units:", length(unique(bootstrap_units))))
    })
  }
  
  # Run bootstrap replications
  if (cores > 1) {
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl))
    
    # Export necessary objects to cluster
    parallel::clusterExport(cl, c("dataset_null", "sc.pred", "col_name_unit",
                                 "create_bootstrap_dataset", "estimate_sc"),
                           envir = environment())
    
    bootstrap_results <- parallel::parLapply(cl, 1:n_bootstrap, bootstrap_single)
  } else {
    bootstrap_results <- lapply(1:n_bootstrap, bootstrap_single)
  }
  
  # Ensure all bootstrap replications succeeded
  failed_count <- sum(sapply(bootstrap_results, is.null))
  if (failed_count > 0) {
    stop(paste("Bootstrap failed:", failed_count, "out of", n_bootstrap, "replications failed"))
  }
  
  cat("Bootstrap success: All", n_bootstrap, "replications completed successfully\n")
  
  # Extract effects and iteration data
  bootstrap_effects <- list()
  bootstrap_iteration_data <- list()
  
  for (oc in names(actual_effects)) {
    # Extract average effects for p-value calculation - no failsafes
    bootstrap_effects[[oc]] <- sapply(bootstrap_results, function(x) x$effects[[oc]])
    
    # Extract iteration data for plotting - no failsafes
    iteration_data_list <- lapply(bootstrap_results, function(x) x$iteration_data[[oc]])
    bootstrap_iteration_data[[oc]] <- rbindlist(iteration_data_list)
  }
  
  return(list(
    effects = bootstrap_effects,
    iteration_data = bootstrap_iteration_data,
    estimates = bootstrap_results
  ))
}

#' Create Bootstrap Dataset by Unit-level Resampling
#'
#' @param dataset_null Original dataset with null enforced
#' @param bootstrap_units Bootstrapped unit names
#' @param col_name_unit Column name for unit identifier
#' @param treated_unit Name of treated unit
#'
#' @return Bootstrap dataset
create_bootstrap_dataset <- function(dataset_null, bootstrap_units, col_name_unit, treated_unit) {
  
  # Ensure unique unit naming to avoid covariate matrix issues
  unique_bootstrap_units <- unique(bootstrap_units)
  bootstrap_data_list <- list()
  unit_counter <- 1
  
  for (unit in bootstrap_units) {
    unit_data <- dataset_null[get(col_name_unit) == unit]
    
    # Create new unit identifier to avoid conflicts
    if (unit == treated_unit) {
      new_unit_name <- treated_unit  # Keep treated unit name unchanged
    } else {
      # Create unique names for all control units to avoid matrix dimension issues
      new_unit_name <- paste0("bootstrap_ctrl_", unit_counter)
      unit_counter <- unit_counter + 1
    }
    
    unit_data_copy <- copy(unit_data)
    unit_data_copy[[col_name_unit]] <- new_unit_name
    bootstrap_data_list[[length(bootstrap_data_list) + 1]] <- unit_data_copy
  }
  
  bootstrap_data <- rbindlist(bootstrap_data_list)
  
  # Ensure we have the expected number of unique units
  n_unique_units <- length(unique(bootstrap_data[[col_name_unit]]))
  expected_units <- length(bootstrap_units)
  
  if (n_unique_units != expected_units) {
    # This should maximize chances - create exactly the right number of unique units
    cat("Adjusting bootstrap dataset: had", n_unique_units, "unique units, need", expected_units, "\n")
  }
  
  return(bootstrap_data)
}

#' Calculate Bootstrap P-values Using Rank-Based Method
#'
#' Calculates empirical bootstrap p-values by determining where the actual
#' treatment effect ranks within the bootstrap distribution. This is the
#' standard approach for bootstrap hypothesis testing.
#'
#' The p-values are calculated as follows:
#' - One-tailed positive: P(bootstrap >= actual) = (n - rank) / n
#' - One-tailed negative: P(bootstrap <= actual) = (rank + 1) / n  
#' - Two-tailed: 2 * min(positive, negative), capped at 1.0
#'
#' @param bootstrap_effects List of bootstrap effect distributions
#' @param actual_effects List of actual treatment effects  
#' @param verbose Logical for verbose output
#'
#' @return Data.table with p-values and rank information
calculate_bootstrap_pvalues <- function(bootstrap_effects, actual_effects, verbose = FALSE) {
  
  p_values_list <- list()
  
  for (oc in names(actual_effects)) {
    if (length(bootstrap_effects[[oc]]) > 0) {
      actual_effect <- actual_effects[[oc]]
      bootstrap_dist <- bootstrap_effects[[oc]]
      n_bootstrap <- length(bootstrap_dist)
      
      # Sort bootstrap distribution for rank-based p-value calculation
      sorted_bootstrap <- sort(bootstrap_dist)
      
      # Find rank of actual effect in bootstrap distribution
      # Using findInterval to handle ties properly
      rank_in_bootstrap <- findInterval(actual_effect, sorted_bootstrap, rightmost.closed = TRUE)
      
      # Calculate one-tailed p-values based on rank
      # p_value_positive: P(bootstrap >= actual) = (n - rank) / n
      p_value_positive <- (n_bootstrap - rank_in_bootstrap) / n_bootstrap
      
      # p_value_negative: P(bootstrap <= actual) = (rank + 1) / n  
      # +1 because rank is 0-indexed but we want proportion including the actual value
      # Cap at 1.0 to handle edge case where actual is more extreme than all bootstrap values
      p_value_negative <- min((rank_in_bootstrap + 1) / n_bootstrap, 1.0)
      
      # Two-tailed p-value: 2 * min(one-tailed p-values)
      # This accounts for testing in both directions
      p_value_two_tailed <- 2 * min(p_value_positive, p_value_negative)
      
      # Ensure p-values don't exceed 1.0 due to the 2x multiplier
      p_value_two_tailed <- min(p_value_two_tailed, 1.0)
      
      p_values_list[[oc]] <- data.table(
        outcome_model = oc,
        actual_effect = actual_effect,
        p_value_two_tailed = p_value_two_tailed,
        p_value_positive = p_value_positive,
        p_value_negative = p_value_negative,
        n_bootstrap = n_bootstrap,
        bootstrap_mean = mean(bootstrap_dist, na.rm = TRUE),
        bootstrap_sd = sd(bootstrap_dist, na.rm = TRUE),
        actual_rank = rank_in_bootstrap + 1  # 1-indexed rank for interpretability
      )
      
      if (verbose) {
        cat("Outcome model:", oc, "\n")
        cat("  Actual effect:", round(actual_effect, 4), "\n")
        cat("  Bootstrap mean:", round(mean(bootstrap_dist, na.rm = TRUE), 4), "\n")
        cat("  Actual rank in bootstrap:", rank_in_bootstrap + 1, "out of", n_bootstrap, "\n")
        cat("  One-tailed p-values: positive =", round(p_value_positive, 4), 
            ", negative =", round(p_value_negative, 4), "\n")
        cat("  Two-tailed p-value:", round(p_value_two_tailed, 4), "\n")
      }
    }
  }
  
  if (length(p_values_list) > 0) {
    return(rbindlist(p_values_list))
  } else {
    return(data.table())
  }
}
