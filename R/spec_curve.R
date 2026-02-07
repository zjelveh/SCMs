#' Generate Specification Curve (Core Computational Engine)
#'
#' @title **CORE ENGINE** - Specification Curve Computation
#' @description **INTERNAL COMPUTATIONAL CORE** - This function performs the core computation
#' for specification curve analysis by systematically varying modeling choices and estimating
#' treatment effects across all combinations. Called internally by `run_spec_curve_analysis()`.
#' For user-facing specification curve analysis, use `run_spec_curve_analysis()` instead.
#'
#' @param dataset Data frame containing panel data in long format with units, time periods, and outcomes.
#' @param outcomes Character vector of outcome variable names to analyze. Multiple outcomes will be processed separately.
#' @param col_name_unit_name Character. Name of column containing unit identifiers (e.g., "state", "country").
#' @param name_treated_unit Character. Identifier of the treated unit as it appears in the data.
#' @param covagg List of covariate specifications for pre-treatment matching.
#'   Each top-level entry must be a named specification with optional \code{label}
#'   and required \code{operations}, where \code{operations} is a non-empty list
#'   of operation signatures:
#'   \code{list(var, select_periods, partition_periods, compute)}.
#'   Use \code{"outcome_var"} in \code{var} to reference the current outcome.
#'   Default is empty list.
#' @param treated_period Numeric. First time period when treatment is active for the treated unit.
#' @param min_period Numeric. Earliest time period available in the dataset.
#' @param end_period Numeric. Latest time period available in the dataset.
#' @param col_name_period Character. Name of column containing time period identifiers.
#' @param feature_weights Character vector of feature weighting methods:
#'   \itemize{
#'     \item \code{"uniform"} - Equal weights for all pre-treatment periods
#'     \item \code{"optimize"} - Data-driven optimization of feature weights  
#'   }
#'   Default is \code{c('uniform', 'optimize')}.
#' @param outcome_models Character vector of outcome modeling approaches:
#'   \itemize{
#'     \item \code{"none"} - Standard synthetic control (no outcome model)
#'     \item \code{"augsynth"} - Augmented synthetic control method
#'     \item \code{"ridge"} - Ridge regression outcome model
#'     \item \code{"lasso"} - Lasso regression outcome model  
#'     \item \code{"ols"} - OLS regression outcome model
#'   }
#'   Default is \code{c('none', 'augsynth', 'ridge', 'lasso', 'ols')}.
#' @param donor_sample Character vector specifying donor pool selection:
#'   \itemize{
#'     \item \code{"all"} - Use all available control units
#'     \item \code{"most_similar"} - Use only most similar control units
#'   }
#'   Default is \code{c('all', 'most_similar')}.
#' @param sim_function Function to determine similarity for donor selection when \code{donor_sample} includes "most_similar".
#'   Default uses built-in similarity function.
#' @param constraints List of constraint specifications for weight optimization. Each element should be a list
#'   with constraint parameters (e.g., \code{list(name = "simplex")}, \code{list(name = "lasso", Q = 0.1)}).
#' @param cores Integer. Number of CPU cores for parallel processing. Default is 1 (sequential).
#'   Higher values significantly speed up computation but require more memory.
#' @param verbose Logical. Whether to display progress messages and warnings. Default is TRUE.
#' @param inference_type Character. Type of inference to perform:
#'   \itemize{
#'     \item \code{"placebo"} - Abadie-style placebo inference only (default)
#'     \item \code{"bootstrap"} - Bootstrap null hypothesis inference only
#'     \item \code{"all"} - Both placebo and bootstrap inference
#'   }
#' @param inference_config List. Configuration settings for inference methods:
#'   \itemize{
#'     \item \code{bootstrap_n_replications} - Number of bootstrap replications (default: 1000)
#'     \item \code{bootstrap_cores} - Number of cores for bootstrap processing (default: 1)
#'   }
#'
#' @details
#' The function creates a Cartesian product of all specification choices and estimates synthetic controls
#' for each combination. This implementation includes:
#' \itemize{
#'   \item \strong{Parallel Processing}: Distributes computations across multiple cores for efficiency
#'   \item \strong{Error Handling}: Gracefully handles failed estimations and continues processing
#'   \item \strong{Progress Tracking}: Provides real-time updates on estimation progress
#'   \item \strong{Memory Management}: Efficiently handles large specification spaces
#' }
#'
#' The resulting specification curve can reveal:
#' \itemize{
#'   \item Sensitivity of treatment effect estimates to modeling choices
#'   \item Distribution of effect sizes across specifications  
#'   \item Identification of modeling choices that drive results
#'   \item Assessment of overall robustness of findings
#' }
#'
#' @return A list structure containing:
#' \itemize{
#'   \item \code{results} - Data.table in long format with all specification results (treatment effects, unit data, RMSE)
#'   \item \code{abadie_inference} - List containing Abadie placebo inference results:
#'     \itemize{
#'       \item \code{p_values} - Specification-level p-values and significance tests
#'       \item \code{post_pre_ratios} - Unit-level post/pre RMSPE ratios for ranking
#'       \item \code{filtered_results} - Filtered significance results (if available)
#'     }
#'   \item \code{bootstrap_inference} - List containing bootstrap inference results:
#'     \itemize{
#'       \item \code{p_values} - Specification-level bootstrap p-values and statistics
#'     }
#' }
#'
#' The \code{results} data.table contains treatment effects (tau), unit information, and specification
#' metadata for all units (treated, control placebo, bootstrap iterations) across all specifications.
#' Inference results are separated to avoid duplication and provide cleaner access to statistical tests.
#'
#' @references
#' \itemize{
#'   \item Simonsohn, U., Simmons, J. P., & Nelson, L. D. (2020). Specification curve analysis.
#'         \emph{Nature Human Behaviour}, 4(11), 1208-1214.
#'   \item Huntington-Klein, N. et al. (2021). The influence of hidden researcher decisions in applied microeconomics.
#'         \emph{Economic Inquiry}, 59(3), 944-960.
#' }
#'
#' @examples
#' \dontrun{
#' # Basic specification curve analysis
#' results <- spec_curve(
#'   dataset = state_panel,
#'   outcomes = "gdp_per_capita",
#'   col_name_unit_name = "state",
#'   name_treated_unit = "California",
#'   covagg = list(
#'     "Outcome Only" = list(
#'       label = "Outcome Only",
#'       operations = list(
#'         list(
#'           var = "outcome_var",
#'           partition_periods = list(type = "by_period")
#'         )
#'       )
#'     ),
#'     "Outcome + Means" = list(
#'       label = "Outcome + Means",
#'       operations = list(
#'         list(
#'           var = "outcome_var",
#'           partition_periods = list(type = "by_period")
#'         ),
#'         list(var = "population", compute = "mean"),
#'         list(var = "income", compute = "mean")
#'       )
#'     )
#'   ),
#'   treated_period = 2000,
#'   min_period = 1990,
#'   end_period = 2010,
#'   col_name_period = "year",
#'   constraints = list(
#'     list(name = "simplex"),
#'     list(name = "lasso", Q = 0.1),
#'     list(name = "ridge", Q = 0.1)
#'   ),
#'   cores = 4
#' )
#'
#' # Plot results
#' plot_spec_curve(results)
#' }
#' 

#' @export
spec_curve <- function(
    dataset,
    outcomes,
    col_name_unit_name,
    name_treated_unit,
    covagg = list(),
    treated_period,
    min_period,
    end_period,
    col_name_period,
    feature_weights = c('uniform', 'optimize'),
    outcome_models = c('none', 'augsynth', 'ridge', 'lasso', 'ols'),
    donor_sample = c('all', 'most_similar'),
    sim_function = most_similar,
    constraints = list(
      list(name='ols'),
      list(name='simplex'),
      list(name='lasso'),
      list(name='ridge'),
      list(name='L1-L2')
    ),
    constants = c(FALSE, TRUE),
    cores = 1,
    verbose = TRUE,
    inference_type = "placebo",
    inference_config = list(),
    expected_direction = "negative"
) {
  # Input validation
  if (!any(c("data.frame", "data.table") %in% class(dataset))) {
    stop("dataset must be a data.frame or data.table")
  }
  
  validate_character(outcomes, "outcomes", length = length(outcomes))
  validate_character(col_name_unit_name, "col_name_unit_name")
  validate_character(col_name_period, "col_name_period")
  validate_character(name_treated_unit, "name_treated_unit")
  validate_numeric(treated_period, "treated_period")
  validate_numeric(min_period, "min_period")
  validate_numeric(end_period, "end_period")
  validate_period_relationships(min_period, treated_period, end_period)
  
  validate_dataframe(
    dataset,
    arg_name = "dataset",
    allow_empty = FALSE,
    required_cols = c(outcomes, col_name_unit_name, col_name_period)
  )
  
  validate_values_in_column(dataset, col_name_unit_name, name_treated_unit, "name_treated_unit")

  validate_character(feature_weights, "feature_weights", length = length(feature_weights),
                     choices = c("uniform", "optimize"))
  validate_character(outcome_models, "outcome_models", length = length(outcome_models),
                     choices = c("none", "augsynth", "ridge", "lasso", "ols"))

  # Validate inference_type parameter
  valid_inference_types <- c("placebo", "bootstrap", "all")
  if (!inference_type %in% valid_inference_types) {
    stop("inference_type must be one of: ", paste(valid_inference_types, collapse = ", "))
  }
  
  # Validate donor_sample options early to avoid undefined behavior downstream
  validate_character(donor_sample, "donor_sample", length = length(donor_sample),
                     choices = c("all", "most_similar"))
  
  # Validate cores parameter
  if (!is.numeric(cores) || length(cores) != 1 || cores < 1) {
    stop("cores must be >= 1")
  }


  # Validate that we have at least some specifications
  if (length(covagg) == 0) {
    stop("No covariate specifications provided. Must provide covagg specifications.")
  }

  # Validate covagg spec structure (new API):
  # covagg = list(
  #   "Spec Name" = list(
  #     label = "Spec Name",
  #     operations = list(
  #       list(var = "outcome_var", ...),
  #       list(var = "population", ...)
  #     )
  #   )
  # )
  covagg_names <- names(covagg)
  if (is.null(covagg_names) || any(is.na(covagg_names) | covagg_names == "")) {
    stop("covagg entries must be named specifications")
  }

  for (i in seq_along(covagg)) {
    spec <- covagg[[i]]
    if (!is.list(spec)) {
      stop("covagg specification '", covagg_names[[i]], "' must be a list")
    }
    if (is.null(spec$operations) || !is.list(spec$operations) || length(spec$operations) == 0) {
      stop("covagg specification '", covagg_names[[i]], "' must define a non-empty operations list")
    }
    # Validate operations now against the available pre-period range
    compile_covagg_signature(
      covagg = spec$operations,
      outcome_var = outcomes[[1]],
      period_pre = min_period:(treated_period - 1),
      require_named = FALSE,
      context = paste0("spec_curve covagg specification '", covagg_names[[i]], "'")
    )
  }

  # Set default inference config values
  default_config <- list(
    bootstrap_n_replications = 1000,
    bootstrap_cores = 1
  )

  # Merge user config with defaults
  inference_config <- utils::modifyList(default_config, inference_config)

  # Always return long format (no more nested format)
  collect_long_format <- TRUE

  # Ensure dataset is a data.table for faster operations
  if (!data.table::is.data.table(dataset)) {
    dataset <- data.table::as.data.table(dataset)
  }


  # Pre-compute donor samples to avoid redundancy
  donor_samples <- list()
  if ("most_similar" %in% donor_sample) {
    for (outc in outcomes) {
      if (verbose) message(sprintf("Pre-computing donor sample for outcome: %s", outc))
      donor_samples[[outc]] <- sim_function(dataset,
                                            outc,
                                            col_name_unit_name,
                                            name_treated_unit,
                                            treated_period,
                                            min_period,
                                            col_name_period)
    }
  }

  # Generate parameter grid for all combinations
  create_param_grid <- function() {
    param_list <- list(
      outcome = outcomes,
      constraint = lapply(constraints, function(x) x$name),
      feature_weight = feature_weights,
      covariate_agg = seq_along(covagg),
      donor = donor_sample,
      constant = constants
    )

    # Create all combinations
    expand.grid(
      outcome = param_list$outcome,
      constraint_idx = seq_along(param_list$constraint),
      feature_weight = param_list$feature_weight,
      covariate_agg = param_list$covariate_agg,
      donor = param_list$donor,
      constant = param_list$constant,
      stringsAsFactors = FALSE
    )
  }

  # Function to process a single specification
  process_spec <- function(
    outc, const_idx, fw, ca_idx, ds, constant, spec_number,
    constraints, covagg, donor_samples, dataset,
    min_period, treated_period, end_period,
    col_name_unit_name, name_treated_unit, col_name_period,
    outcome_models, verbose
  ) {
    const <- constraints[[const_idx]]
    ca <- covagg[[ca_idx]]

    # Create feature names - Use label field or covagg name
    covagg_name <- names(covagg)[ca_idx]

    # Check if specification has a label field (preferred approach)
    if (is.list(ca) && !is.null(ca$label) && is.character(ca$label) && length(ca$label) == 1) {
      feature_names <- ca$label
    } else {
      # Use covagg name directly, but skip dummy names
      if (!is.null(covagg_name) && nchar(covagg_name) > 0 && !grepl("^spec_\\d+$", covagg_name)) {
        feature_names <- covagg_name
      } else {
        # FAIL HARD - no fallback logic
        stop(paste("Covariate specification at index", ca_idx, "must have a descriptive name or label field."))
      }
    }

    # Use the provided min_period directly
    local_min_period <- min_period

    # Select donor sample (avoid unnecessary copy for 'all')
    if (ds == 'all') {
      dataset2 <- dataset  # No copy needed if not modifying
    } else if (ds == 'most_similar') {
      dataset2 <- donor_samples[[outc]]
    }

    # Create specification identifier for logging
    spec_name <- paste(outc, const$name, fw, feature_names, ds, sep = "-")

    if (verbose) message(sprintf("Processing: %s", spec_name))

    if (!is.list(ca) || is.null(ca$operations) || !is.list(ca$operations)) {
      stop("Invalid covagg specification: each entry must be a list with operations")
    }
    ca_formatted <- ca$operations

    sc.pred <- tryCatch({
      estimate_sc(
        dataset2,
        outc,
        ca_formatted,
        col_name_unit_name,
        name_treated_unit,
        col_name_period,
        treated_period,
        local_min_period,
        end_period,
        outcome_models = outcome_models,
        feature_weights = fw,
        w.constr = const,
        constant = constant
      )
    }, error = function(e) {
      message(sprintf("ERROR in estimate_sc for spec %s: %s", spec_name, e$message))
      return(NULL)
    })

    if (!is.null(sc.pred)) {

      # Initialize inference structure
      sc.infer <- list()

      # Run placebo inference if requested
      if (inference_type %in% c("placebo", "all")) {
        placebo_results <- tryCatch({
          inference_sc(
            sc.pred,
            dataset2,
            cores = 1,  # No nested parallelism
            verbose = FALSE,
            expected_direction = expected_direction
          )
        }, error = function(e) {
          warning(sprintf("Error in placebo inference for spec '%s': %s", spec_name, e$message))
          return(NULL)
        })


        if (!is.null(placebo_results)) {
          # Standardize structure to match bootstrap mode expectations
          sc.infer <- list(
            inference.results = placebo_results$inference.results,
            rmse = placebo_results$rmse,
            abadie_significance = placebo_results$abadie_significance
          )

        }
      }

      # Run bootstrap inference if requested
      if (inference_type %in% c("bootstrap", "all")) {
        bootstrap_results <- tryCatch({
          bootstrap_null_inference(
            sc.pred = sc.pred,
            dataset = dataset2,
            n_bootstrap = inference_config$bootstrap_n_replications,
            cores = inference_config$bootstrap_cores,
            verbose = verbose
          )
        }, error = function(e) {
          if (verbose) {
            warning(sprintf("Error in bootstrap inference for spec '%s': %s", spec_name, e$message))
          }
          return(NULL)
        })

        # Add bootstrap results to inference structure
        if (!is.null(bootstrap_results)) {
          if (inference_type == "bootstrap") {
            # Bootstrap-only mode: create complete inference structure with bootstrap data

            # Create treated unit data (same as placebo mode)
            treated_taus_list <- list()
            treated_rmse_list <- list()

            for(oc in names(sc.pred$est.results$outcome_model)){
              sc_post = sc.pred$est.results$outcome_model[[oc]]$Y.post.fit
              sc_pre  = sc.pred$est.results$Y.pre.fit

              treated_tau = sc.pred$data$Y.post - sc_post
              treated_pre_tau = sc.pred$data$Y.pre - sc_pre
              treated_all_tau = c(treated_pre_tau, treated_tau)
              treated_pre_rmse = sqrt(mean(treated_pre_tau^2))

              treated_rmse_list[[oc]] = data.table(
                unit_name = sc.pred$name_treated_unit,
                pre_rmse = treated_pre_rmse,
                unit_type = 'treated',
                outcome_model = oc
              )

              treated_taus_list[[oc]] = data.table(
                unit_name = rep(sc.pred$name_treated_unit, length(treated_all_tau)),
                period = sc.pred$min_period:sc.pred$end_period,
                tau = treated_all_tau,
                post_period = (sc.pred$min_period:sc.pred$end_period) >= sc.pred$treated_period,
                unit_type = 'treated',
                outcome_model = oc
              )
            }

            # Combine treated data with bootstrap iteration data
            all_taus_list <- treated_taus_list
            all_rmse_list <- treated_rmse_list

            # Add bootstrap iteration data as "control" units
            for(oc in names(bootstrap_results$bootstrap_iteration_data)) {
              bootstrap_data <- bootstrap_results$bootstrap_iteration_data[[oc]]
              if (nrow(bootstrap_data) > 0) {
                # Add to tau data
                all_taus_list[[paste0("bootstrap_", oc)]] <- bootstrap_data

                # Create RMSE data for bootstrap iterations (calculate from pre-period)
                bootstrap_rmse <- bootstrap_data[post_period == FALSE, .(
                  pre_rmse = sqrt(mean(tau^2))
                ), by = unit_name]
                bootstrap_rmse[, `:=`(unit_type = 'bootstrap', outcome_model = oc)]
                all_rmse_list[[paste0("bootstrap_rmse_", oc)]] <- bootstrap_rmse
              }
            }

            # Create inference structure compatible with existing extraction code
            sc.infer <- list(
              inference.results = rbindlist(all_taus_list, fill = TRUE),
              rmse = rbindlist(all_rmse_list, fill = TRUE),
              bootstrap_significance = bootstrap_results
            )

          } else {
            # "all" mode: add bootstrap results to existing placebo inference
            sc.infer$bootstrap_significance <- bootstrap_results

              # 2. Add the raw bootstrap iteration data to the main tables
            for(oc in names(bootstrap_results$bootstrap_iteration_data)) {
                bootstrap_data <- bootstrap_results$bootstrap_iteration_data[[oc]]
                if (nrow(bootstrap_data) > 0) {
                    # Add bootstrap tau data to the main inference table
                    sc.infer$inference.results <- rbindlist(
                      list(sc.infer$inference.results, bootstrap_data),
                      fill = TRUE
                    )

                    # Create and add bootstrap RMSE data
                    bootstrap_rmse <- bootstrap_data[post_period == FALSE, .(
                        pre_rmse = sqrt(mean(tau^2))
                    ), by = unit_name]
                    bootstrap_rmse[, `:=`(unit_type = 'bootstrap', outcome_model = oc)]

                    sc.infer$rmse <- rbindlist(
                      list(sc.infer$rmse, bootstrap_rmse),
                      fill = TRUE
                    )
                }
            }
          }
        }
      }
      result <- list(
        estimate = sc.pred,
        infer = sc.infer
      )

    } else {
      result <- NULL
    }

    return(list(
      outc = outc,
      spec_number = spec_number,
      const_name = const$name,
      fw = fw,
      feature_names = feature_names,
      ds = ds,
      constant = constant,
      result = result
    ))
  }

  # Create parameter grid
  param_grid <- create_param_grid()
  total_specs <- nrow(param_grid)

  # Add specification numbers early - each row is one unique specification
  param_grid$spec_number <- 1:total_specs


  if (verbose) {
    message(sprintf("Running %d specifications with %d cores", total_specs, cores))
  }

  # Run specifications (parallel or sequential)
  if (cores > 1 && total_specs > 1) {
    # Limit cores to reasonable maximum and available cores
    max_cores <- min(cores, parallel::detectCores(), total_specs)

    if (verbose) {
      cat("Setting up parallel processing with", max_cores, "cores\n")
    }

    # Setup parallel backend with error handling
    tryCatch({
      # Ensure required packages are loaded
      if (!requireNamespace("doParallel", quietly = TRUE)) {
        stop("doParallel package not available")
      }
      if (!requireNamespace("foreach", quietly = TRUE)) {
        stop("foreach package not available")
      }

      cl <- parallel::makeCluster(max_cores)
      doParallel::registerDoParallel(cl)
      on.exit(parallel::stopCluster(cl))

      if (verbose) {
        cat("Parallel backend registered successfully\n")
      }
    }, error = function(e) {
      warning("Failed to setup parallel processing: ", e$message, ". Falling back to sequential processing.")
      cores <<- 1  # Force sequential processing
    })

  }

  # Execute specifications
  if (cores > 1 && total_specs > 1) {
    # Packages for parallel execution imported via NAMESPACE

    # Use foreach for parallel execution
    results_list <- foreach::foreach(
      i = 1:nrow(param_grid),
      .packages = c('data.table', 'SCMs', 'optimx', 'kernlab'),
      .noexport = c('data.table', '.SD', '.N', '.BY')
    ) %dopar% {
      # SCMs package functions available in parallel worker

      p <- param_grid[i, ]
      process_spec(
        outc = p$outcome,
        const_idx = p$constraint_idx,
        fw = p$feature_weight,
        ca_idx = p$covariate_agg,
        ds = p$donor,
        constant = p$constant,
        spec_number = p$spec_number,
        constraints = constraints,
        covagg = covagg,
        donor_samples = donor_samples,
        dataset = dataset,
        min_period = min_period,
        treated_period = treated_period,
        end_period = end_period,
        col_name_unit_name = col_name_unit_name,
        name_treated_unit = name_treated_unit,
        col_name_period = col_name_period,
        outcome_models = outcome_models,
        verbose = FALSE  # Disable verbose in parallel workers
      )
    }
  } else {
    # Sequential processing
    results_list <- vector("list", total_specs)
    for (i in 1:nrow(param_grid)) {
      p <- param_grid[i, ]

      results_list[[i]] <- process_spec(
        outc = p$outcome,
        const_idx = p$constraint_idx,
        fw = p$feature_weight,
        ca_idx = p$covariate_agg,
        ds = p$donor,
        constant = p$constant,
        spec_number = p$spec_number,
        constraints = constraints,
        covagg = covagg,
        donor_samples = donor_samples,
        dataset = dataset,
        min_period = min_period,
        treated_period = treated_period,
        end_period = end_period,
        col_name_unit_name = col_name_unit_name,
        name_treated_unit = name_treated_unit,
        col_name_period = col_name_period,
        outcome_models = outcome_models,
        verbose = verbose
      )

      # Force garbage collection periodically
      if (i %% 20 == 0) gc()
    }
  }

  # Process results into separate data structures
  long_results_list <- list()
  abadie_results_list <- list()
  bootstrap_results_list <- list()

  for (res in results_list) {
    if (!is.null(res) && !is.null(res$result)) {

      # Extract data from the result structure
      if (!is.null(res$result$infer) &&
          is.list(res$result$infer) &&
          "inference.results" %in% names(res$result$infer) &&
          "rmse" %in% names(res$result$infer)) {

        # Get treatment effects data
        tau_data <- res$result$infer$inference.results
        rmse_data <- res$result$infer$rmse

        # Make sure they're data.tables (using setDT for efficiency)
        if (!data.table::is.data.table(tau_data)) data.table::setDT(tau_data)
        if (!data.table::is.data.table(rmse_data)) data.table::setDT(rmse_data)

        # Merge the main datasets on the common keys
        combined <- merge(
          tau_data,
          rmse_data,
          by = c("unit_name", "unit_type", "outcome_model"),
          all.x = TRUE
        )

        # Add specification metadata to main results
        combined$outcome <- res$outc
        combined$const <- res$const_name
        combined$fw <- res$fw
        combined$data_sample <- res$ds
        combined$feat <- res$feature_names
        combined$spec_number <- res$spec_number
        combined$constant <- res$constant
        combined$rmse <- combined$pre_rmse

        ## FIX: Create a new, fully unique specification ID
        combined[, full_spec_id := paste(spec_number, outcome_model, sep = "_")]

        long_results_list[[length(long_results_list) + 1]] <- combined

        # Process Abadie inference results separately
        if ("abadie_significance" %in% names(res$result$infer)) {
          abadie_data <- res$result$infer$abadie_significance

          # Create Abadie results for this specification
          abadie_spec_results <- list()

          # Process each test statistic type
          test_stats <- c("rmse_ratio", "treatment_effect", "normalized_te")

          for (test_stat in test_stats) {
            if (!is.null(abadie_data[[test_stat]])) {
              test_stat_data <- abadie_data[[test_stat]]

              # Add p-values for this test statistic
              if (!is.null(test_stat_data$p_values) && nrow(test_stat_data$p_values) > 0) {
                p_values_data <- copy(test_stat_data$p_values)
                p_values_data[, `:=`(
                  full_spec_id = paste(res$spec_number, outcome_model, sep = "_"),
                  inference_type = paste0("abadie_pvalue_", test_stat)
                )]
                abadie_spec_results[[paste0("p_values_", test_stat)]] <- p_values_data
              }

              # Add test statistics for this test statistic type
              if (!is.null(test_stat_data$test_statistics) && nrow(test_stat_data$test_statistics) > 0) {
                test_stats_data <- copy(test_stat_data$test_statistics)
                test_stats_data[, `:=`(
                  full_spec_id = paste(res$spec_number, outcome_model, sep = "_"),
                  inference_type = paste0("abadie_teststat_", test_stat)
                )]
                abadie_spec_results[[paste0("test_statistics_", test_stat)]] <- test_stats_data
              }

              # Add filtered results (only for rmse_ratio)
              if (test_stat == "rmse_ratio" &&
                  !is.null(test_stat_data$filtered_results) &&
                  nrow(test_stat_data$filtered_results) > 0) {
                filtered_data <- copy(test_stat_data$filtered_results)
                filtered_data[, `:=`(
                  full_spec_id = paste(res$spec_number, outcome_model, sep = "_"),
                  inference_type = "abadie_filtered_rmse_ratio"
                )]
                abadie_spec_results[["filtered_results"]] <- filtered_data
              }
            }
          }


          abadie_results_list[[length(abadie_results_list) + 1]] <- abadie_spec_results
        }

        # Process bootstrap inference results separately
        if ("bootstrap_significance" %in% names(res$result$infer)) {
          bootstrap_data <- res$result$infer$bootstrap_significance

          # Create bootstrap results for this specification
          bootstrap_spec_results <- list()

          # Add bootstrap p-values (specification-level)
          if (!is.null(bootstrap_data$p_values) && nrow(bootstrap_data$p_values) > 0) {
            bootstrap_p_data <- copy(bootstrap_data$p_values)
            # Only add essential merge key and inference type - no redundant spec metadata
            bootstrap_p_data[, `:=`(
              full_spec_id = paste(res$spec_number, outcome_model, sep = "_"),
              inference_type = "bootstrap_pvalue"
            )]
            bootstrap_spec_results[["p_values"]] <- bootstrap_p_data
          }

          bootstrap_results_list[[length(bootstrap_results_list) + 1]] <- bootstrap_spec_results
        }
      } else {
      }
    }
  }

  # Always return structured results with separate inference tables
  if (length(long_results_list) > 0) {
    # Combine all long format results (spec_number already included from process_spec)
    long_df <- data.table::rbindlist(long_results_list, fill = TRUE)

    # Create return structure
    final_results <- list(
      results = long_df,
      expected_direction = expected_direction
    )

    # Process and add Abadie inference results if available
    if (length(abadie_results_list) > 0) {
      abadie_combined <- list()

      # Initialize lists for each test statistic type
      test_stats <- c("rmse_ratio", "treatment_effect", "normalized_te")

      for (test_stat in test_stats) {
        p_values_key <- paste0("p_values_", test_stat)
        test_stats_key <- paste0("test_statistics_", test_stat)

        p_values_list <- list()
        test_stats_list <- list()

        for (spec_results in abadie_results_list) {
          if (p_values_key %in% names(spec_results)) {
            p_values_list[[length(p_values_list) + 1]] <- spec_results[[p_values_key]]
          }
          if (test_stats_key %in% names(spec_results)) {
            test_stats_list[[length(test_stats_list) + 1]] <- spec_results[[test_stats_key]]
          }
        }

        if (length(p_values_list) > 0) {
          abadie_combined[[p_values_key]] <- data.table::rbindlist(p_values_list, fill = TRUE)
        }
        if (length(test_stats_list) > 0) {
          abadie_combined[[test_stats_key]] <- data.table::rbindlist(test_stats_list, fill = TRUE)
        }
      }

      # Handle filtered results (only for rmse_ratio)
      abadie_filtered <- list()
      for (spec_results in abadie_results_list) {
        if ("filtered_results" %in% names(spec_results)) {
          abadie_filtered[[length(abadie_filtered) + 1]] <- spec_results$filtered_results
        }
      }

      if (length(abadie_filtered) > 0) {
        abadie_combined$filtered_results <- data.table::rbindlist(abadie_filtered, fill = TRUE)
      }


      final_results$abadie_inference <- abadie_combined
    }

    # Process and add bootstrap inference results if available
    if (length(bootstrap_results_list) > 0) {
      bootstrap_combined <- list()

      # Combine p-values across specifications
      bootstrap_p_values <- list()

      for (spec_results in bootstrap_results_list) {
        if ("p_values" %in% names(spec_results)) {
          bootstrap_p_values[[length(bootstrap_p_values) + 1]] <- spec_results$p_values
        }
      }

      if (length(bootstrap_p_values) > 0) {
        bootstrap_combined$p_values <- data.table::rbindlist(bootstrap_p_values, fill = TRUE)
      }

      final_results$bootstrap_inference <- bootstrap_combined
    }

  } else {
    stop("No valid results found")
  }

  # Assign S3 class for proper method dispatch
  class(final_results) <- c("spec_curve", "list")
  
  if (verbose) message("Specification curve generation complete")
  return(final_results)
}
