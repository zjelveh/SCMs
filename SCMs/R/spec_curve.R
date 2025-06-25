#' Generate Specification Curve for Synthetic Control Method
#'
#' This function performs comprehensive specification curve analysis for Synthetic Control Methods
#' by systematically varying multiple modeling choices and estimating treatment effects across
#' all combinations. This approach helps assess the robustness of results to modeling assumptions
#' and identifies the sensitivity of estimates to researcher degrees of freedom.
#'
#' @param dataset Data frame containing panel data in long format with units, time periods, and outcomes.
#' @param outcomes Character vector of outcome variable names to analyze. Multiple outcomes will be processed separately.
#' @param col_name_unit_name Character. Name of column containing unit identifiers (e.g., "state", "country").
#' @param name_treated_unit Character. Identifier of the treated unit as it appears in the data.
#' @param covagg List of covariate specifications for pre-treatment matching. Each element should be a vector of variable names.
#' @param treated_period Numeric. First time period when treatment is active for the treated unit.
#' @param min_period Numeric. Earliest time period available in the dataset.
#' @param end_period Numeric. Latest time period available in the dataset.
#' @param col_name_period Character. Name of column containing time period identifiers.
#' @param feature_weights Character vector of feature weighting methods:
#'   \itemize{
#'     \item \code{"uniform"} - Equal weights for all pre-treatment periods
#'     \item \code{"optimize"} - Data-driven optimization of feature weights  
#'   }
#'   Default is \code{c('uniform', 'optimize')}.
#' @param num_pre_period_years Numeric or NA. Number of pre-treatment periods to include in estimation.
#'   If NA, uses all available pre-treatment periods. Default is NA.
#' @param outcome_models Character vector of outcome modeling approaches:
#'   \itemize{
#'     \item \code{"none"} - Standard synthetic control (no outcome model)
#'     \item \code{"augsynth"} - Augmented synthetic control method
#'     \item \code{"ridge"} - Ridge regression outcome model
#'     \item \code{"lasso"} - Lasso regression outcome model  
#'     \item \code{"ols"} - OLS regression outcome model
#'   }
#'   Default is \code{c('none', 'augsynth', 'ridge', 'lasso', 'ols')}.
#' @param donor_sample Character vector specifying donor pool selection:
#'   \itemize{
#'     \item \code{"all"} - Use all available control units
#'     \item \code{"most_similar"} - Use only most similar control units
#'   }
#'   Default is \code{c('all', 'most_similar')}.
#' @param sim_function Function to determine similarity for donor selection when \code{donor_sample} includes "most_similar".
#'   Default uses built-in similarity function.
#' @param constraints List of constraint specifications for weight optimization. Each element should be a list
#'   with constraint parameters (e.g., \code{list(name = "simplex")}, \code{list(name = "lasso", Q = 0.1)}).
#' @param cores Integer. Number of CPU cores for parallel processing. Default is 1 (sequential).
#'   Higher values significantly speed up computation but require more memory.
#' @param use_cache Logical. Whether to cache intermediate results to avoid recomputation. Default is FALSE.
#'   Useful for large specification spaces or when iterating on analysis.
#' @param cache_dir Character. Directory path for storing cache files when \code{use_cache = TRUE}.
#'   Default is temporary directory from \code{tempdir()}.
#' @param verbose Logical. Whether to display progress messages and warnings. Default is TRUE.
#'
#' @details
#' The function creates a Cartesian product of all specification choices and estimates synthetic controls
#' for each combination. This implementation includes:
#' \itemize{
#'   \item \strong{Parallel Processing}: Distributes computations across multiple cores for efficiency
#'   \item \strong{Error Handling}: Gracefully handles failed estimations and continues processing
#'   \item \strong{Caching System}: Optionally stores results to disk using digest-based keys
#'   \item \strong{Progress Tracking}: Provides real-time updates on estimation progress
#'   \item \strong{Memory Management}: Efficiently handles large specification spaces
#' }
#'
#' The resulting specification curve can reveal:
#' \itemize{
#'   \item Sensitivity of treatment effect estimates to modeling choices
#'   \item Distribution of effect sizes across specifications  
#'   \item Identification of modeling choices that drive results
#'   \item Assessment of overall robustness of findings
#' }
#'
#' @return A nested list structure containing:
#' \itemize{
#'   \item \code{results} - List of estimation results for each specification
#'   \item \code{specifications} - Data frame describing each specification combination
#'   \item \code{summary_stats} - Summary statistics across all specifications
#'   \item \code{failed_runs} - Information about any failed estimations
#'   \item \code{metadata} - Information about the analysis (timing, parameters, etc.)
#' }
#' 
#' Each element in \code{results} contains the full output from \code{scest()} and \code{inference_sc()}
#' for that specification, allowing detailed post-hoc analysis.
#'
#' @references 
#' \itemize{
#'   \item Simonsohn, U., Simmons, J. P., & Nelson, L. D. (2020). Specification curve analysis. 
#'         \emph{Nature Human Behaviour}, 4(11), 1208-1214.
#'   \item Huntington-Klein, N. et al. (2021). The influence of hidden researcher decisions in applied microeconomics. 
#'         \emph{Economic Inquiry}, 59(3), 944-960.
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
#'   covagg = list(c("population", "income")),
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
#' # Multi-outcome analysis with caching
#' results_multi <- spec_curve(
#'   dataset = policy_data,
#'   outcomes = c("unemployment", "gdp", "investment"),
#'   col_name_unit_name = "country",
#'   name_treated_unit = "Germany", 
#'   covagg = list(
#'     c("population", "gdp_lag"),
#'     c("population", "gdp_lag", "trade_openness")
#'   ),
#'   treated_period = 2005,
#'   min_period = 1995,
#'   end_period = 2015,
#'   col_name_period = "year",  
#'   use_cache = TRUE,
#'   cache_dir = "./cache",
#'   cores = 8
#' )
#' 
#' # Plot results
#' plot_spec_curve(results)
#' }
#' @export
spec_curve <- function(
    dataset,
    outcomes,
    col_name_unit_name,
    name_treated_unit,
    covagg,
    treated_period,
    min_period,
    end_period, 
    col_name_period,
    feature_weights = c('uniform', 'optimize'),
    num_pre_period_years = NA,
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
    cores = 1,
    use_cache = FALSE,
    cache_dir = tempdir(),
    verbose = TRUE
) {
  # Input validation
  if (!any(c("data.frame", "data.table") %in% class(dataset))) {
    stop("dataset must be a data.frame or data.table")
  }
  
  # Ensure dataset is a data.table for faster operations
  if (!data.table::is.data.table(dataset)) {
    dataset <- data.table::as.data.table(dataset)
  }
  
  # Setup cache if requested
  if (use_cache) {
    if (!dir.exists(cache_dir)) {
      dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    # Function to generate cache key
    get_cache_key <- function(...) {
      args <- list(...)
      # Use hash of serialized arguments as key
      digest::digest(args, algo = "md5")
    }
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
      pre_period = if(all(is.na(num_pre_period_years))) NA else num_pre_period_years
    )
    
    # Create all combinations
    expand.grid(
      outcome = param_list$outcome,
      constraint_idx = seq_along(param_list$constraint),
      feature_weight = param_list$feature_weight,
      covariate_agg = param_list$covariate_agg,
      donor = param_list$donor,
      pre_period = if(all(is.na(param_list$pre_period))) 1 else param_list$pre_period,
      stringsAsFactors = FALSE
    )
  }
  
  # Function to process a single specification
  process_spec <- function(
    outc, const_idx, fw, ca_idx, ds, ny, 
    constraints, covagg, donor_samples, dataset, 
    min_period, treated_period, end_period,
    col_name_unit_name, name_treated_unit, col_name_period,
    outcome_models, use_cache, cache_dir, verbose
  ) {
    const <- constraints[[const_idx]]
    ca <- covagg[[ca_idx]]
    
    # Create feature names - handle both old and new covagg formats
    # Check if user provided a custom label first
    if (!is.null(ca$label) && is.character(ca$label) && length(ca$label) == 1) {
      feature_names <- ca$label
    } else if (is.list(ca) && "var" %in% names(ca)) {
      # NEW FORMAT: Create descriptive name from specification
      spec_parts <- c()
      spec_parts <- c(spec_parts, ca$var)
      
      if (!is.null(ca$each_period) && ca$each_period) {
        spec_parts <- c(spec_parts, "each_period")
      } else if (!is.null(ca$average)) {
        spec_parts <- c(spec_parts, paste0("avg_", ca$average))
      } else if (!is.null(ca$every_n)) {
        spec_parts <- c(spec_parts, paste0("every_", ca$every_n))
      } else if (!is.null(ca$periods)) {
        spec_parts <- c(spec_parts, paste0("periods_", paste(ca$periods, collapse="-")))
      } else if (!is.null(ca$first)) {
        spec_parts <- c(spec_parts, paste0("first_", ca$first))
      } else if (!is.null(ca$last)) {
        spec_parts <- c(spec_parts, paste0("last_", ca$last))
      } else if (!is.null(ca$rolling)) {
        spec_parts <- c(spec_parts, paste0("roll_", ca$rolling))
      } else if (!is.null(ca$growth)) {
        spec_parts <- c(spec_parts, paste0("growth_", ca$growth))
      } else if (!is.null(ca$volatility)) {
        spec_parts <- c(spec_parts, paste0("vol_", ca$volatility))
      }
      
      feature_names <- paste(spec_parts, collapse = "_")
    } else {
      # OLD FORMAT: Original logic
      # Check for custom label in old format too
      if (!is.null(ca$label) && is.character(ca$label) && length(ca$label) == 1) {
        feature_names <- ca$label
      } else {
        # Remove label from the list before processing if it exists
        ca_for_name <- ca
        if (!is.null(ca_for_name$label)) {
          ca_for_name$label <- NULL
        }
        feature_names <- paste0(
          paste(names(ca_for_name), collapse='_'), '__',
          paste0(ca_for_name, collapse='_')
        )
      }
    }
    
    # Adjust min_period if needed
    local_min_period <- min_period
    if (!is.na(ny)) {
      local_min_period <- treated_period - ny
    }
    
    # Select donor sample
    if (ds == 'all') {
      dataset2 <- data.table::copy(dataset)
    } else if (ds == 'most_similar') {
      dataset2 <- donor_samples[[outc]]
    }
    
    # Create specification identifier for logging/caching
    spec_name <- paste(outc, const$name, fw, feature_names, ds, 
                       if(is.na(ny)) paste0('n_pp_years_', treated_period-min_period) 
                       else paste0('n_pp_years_', ny), sep = "-")
    
    if (verbose) message(sprintf("Processing: %s", spec_name))
    
    # Check cache if enabled
    cache_hit <- FALSE
    if (use_cache) {
      cache_key <- get_cache_key(outc, const$name, fw, feature_names, ds, ny, 
                                 local_min_period, outcome_models)
      cache_file <- file.path(cache_dir, paste0(cache_key, ".rds"))
      
      if (file.exists(cache_file)) {
        if (verbose) message("  Using cached result")
        result <- readRDS(cache_file)
        cache_hit <- TRUE
      }
    }
    
    if (!use_cache || !cache_hit) {
      # Estimate Synthetic Control
      # Ensure ca is in proper format for estimate_sc
      if (is.list(ca) && "var" %in% names(ca)) {
        # NEW FORMAT: Wrap single spec in a list with a name
        ca_formatted <- list()
        ca_formatted[[names(covagg)[ca_idx]]] <- ca
      } else {
        # OLD FORMAT: Use as-is
        ca_formatted <- ca
      }
      
      sc.pred <- 
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
          w.constr = const                    
        )
      
      
      if (!is.null(sc.pred)) {
        
        # Perform inference
        sc.infer <- tryCatch({
          inference_sc(
            sc.pred,
            dataset2,
            inference_type = 'placebo',
            P = NULL,
            u.missp = TRUE,
            u.sigma = "HC1",
            u.order = 1,
            u.lags = 0,
            u.design = NULL,
            u.alpha = 0.05,
            e.method = "all",
            e.order = 1,
            e.lags = 0,
            e.design = NULL,
            e.alpha = 0.05,
            sims = 200,
            rho = NULL,
            rho.max = 0.2,
            lgapp = "generalized",
            cores = 1,  # No nested parallelism
            w.bounds = NULL,
            e.bounds = NULL,
            verbose = FALSE
          )
        }, error = function(e) {
          warning(sprintf("Error in inference_sc for spec '%s': %s", spec_name, e$message))
          return(NULL)
        })
        
        result <- list(
          estimate = sc.pred,
          infer = sc.infer
        )
        
        # Cache result if enabled
        if (use_cache) {
          saveRDS(result, cache_file)
        }
      } else {
        result <- NULL
      }
    }
    
    return(list(
      outc = outc,
      const_name = const$name,
      fw = fw,
      feature_names = feature_names,
      ds = ds,
      pre_period_label = if(is.na(ny)) paste0('n_pp_years_', treated_period-min_period) 
      else paste0('n_pp_years_', ny),
      result = result
    ))
  }
  
  # Create parameter grid
  param_grid <- create_param_grid()
  total_specs <- nrow(param_grid)


  if (verbose) {
    message(sprintf("Running %d specifications with %d cores", total_specs, cores))
  }
  
  # Run specifications (parallel or sequential)
  if (cores > 1 && total_specs > 1) {
    # Setup parallel backend
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
    
    
    # Use foreach for parallel execution
    results_list <- foreach::foreach(
      i = 1:nrow(param_grid),
      .packages = c('data.table', 'digest', 'devtools', 'SCMs'),  # Include SCMs package
      .export =  c('estimate_sc', 'inference_sc', 'create_scm_dataset', 
                   'sim_function', 'most_similar', 'process_covariates',
                   'process_each_period_dt', 'process_specific_periods_dt',
                   'process_every_n_periods_dt', 'process_averages_dt',
                   'process_first_n_dt', 'process_last_n_dt')  # Export new functions
    ) %dopar% {
      # Load SCMs package (use library instead of load_all for installed package)
      library(SCMs)
      
      p <- param_grid[i, ]
      process_spec(
        outc = p$outcome, 
        const_idx = p$constraint_idx, 
        fw = p$feature_weight, 
        ca_idx = p$covariate_agg,
        ds = p$donor, 
        ny = if(p$pre_period == 1 && all(is.na(num_pre_period_years))) NA else p$pre_period,
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
        use_cache = use_cache,
        cache_dir = cache_dir,
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
        ny = if(p$pre_period == 1 && all(is.na(num_pre_period_years))) NA else p$pre_period,
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
        use_cache = use_cache,
        cache_dir = cache_dir,
        verbose = verbose
      )
      
      # Force garbage collection periodically 
      if (i %% 20 == 0) gc()
    }
  }
  
  # Convert flat results list to nested structure
  nested_results <- list()
  for (res in results_list) {
    if (!is.null(res) && !is.null(res$result)) {
      # Initialize nested structure if needed
      if (!res$outc %in% names(nested_results)) {
        nested_results[[res$outc]] <- list()
      }
      if (!res$const_name %in% names(nested_results[[res$outc]])) {
        nested_results[[res$outc]][[res$const_name]] <- list()
      }
      if (!res$fw %in% names(nested_results[[res$outc]][[res$const_name]])) {
        nested_results[[res$outc]][[res$const_name]][[res$fw]] <- list()
      }
      if (!res$feature_names %in% names(nested_results[[res$outc]][[res$const_name]][[res$fw]])) {
        nested_results[[res$outc]][[res$const_name]][[res$fw]][[res$feature_names]] <- list()
      }
      if (!res$ds %in% names(nested_results[[res$outc]][[res$const_name]][[res$fw]][[res$feature_names]])) {
        nested_results[[res$outc]][[res$const_name]][[res$fw]][[res$feature_names]][[res$ds]] <- list()
      }
      
      # Assign result
      nested_results[[res$outc]][[res$const_name]][[res$fw]][[res$feature_names]][[res$ds]][[res$pre_period_label]] <- res$result
    }
  }
  
  if (verbose) message("Specification curve generation complete")
  return(nested_results)
}

