#' Generate Specification Curve for Synthetic Control Method
#'
#' This function generates a specification curve by running multiple variations of 
#' Synthetic Control Method (SCM) estimations across different specifications.
#'
#' @param dataset A data frame containing the panel data.
#' @param outcomes Character vector. Names of the outcome variables to analyze.
#' @param col_name_unit_name Character. Column name in dataset containing unit names.
#' @param name_treated_unit Character. Name of the treated unit.
#' @param covagg List of covariates used for matching.
#' @param treated_period Numeric. Time period when treatment starts for the treated unit.
#' @param min_period Numeric. Earliest time period in the dataset.
#' @param end_period Numeric. Latest time period in the dataset.
#' @param col_name_period Character. Column name in dataset containing time periods.
#' @param feature_weights Character vector. Methods for assigning weights to predictors. Default is c('uniform', 'optimize').
#' @param num_pre_period_years Numeric or NA. Number of pre-treatment periods to use. If NA, uses all available pre-treatment periods.
#' @param outcome_models Character vector. Outcome models to fit. Default is c('none', 'augsynth', 'ridge', 'lasso', 'ols').
#' @param donor_sample Character vector. Method for selecting donor units. Default is c('all', 'most_similar').
#' @param sim_function Function. Function to select most similar donor units. Default is most_similar.
#' @param constraints List of lists. Each inner list specifies a constraint type for the SCM estimation.
#' @param cores Integer. Number of cores to use for parallel processing. Default is 1 (sequential).
#' @param use_cache Logical. Whether to cache results to disk. Default is FALSE.
#' @param cache_dir Character. Directory to store cache files if use_cache is TRUE. Default is tempdir().
#' @param verbose Logical. Whether to display progress messages. Default is TRUE.
#'
#' @return A nested list containing results for each combination of specifications.
#'
#' @examples
#' \dontrun{
#' results <- spec_curve(dataset = my_data, outcomes = c("gdp"),
#'                       col_name_unit_name = "state", name_treated_unit = "California",
#'                       covagg = list(c("population")), treated_period = 2000,
#'                       min_period = 1990, end_period = 2010, col_name_period = "year",
#'                       cores = 2)
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
    
    # Create feature names
    feature_names <- paste0(
      names(ca), '__',
      paste0(ca, collapse='_')
    )
    
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
      sc.pred <- 
        estimate_sc(
          dataset2,
          outc,
          ca,
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
      .packages = c('data.table', 'digest', 'devtools'),  # Add your package here
      .export =  c('estimate_sc', 'inference_sc', 'create_scm_dataset', 
                   'sim_function', 'most_similar')  # Explicitly export functions
    ) %dopar% {
      # Load custom SCM package
      load_all('~/Dropbox/research/synth_control_paper/repo/SCMs/SCMs/')
      
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

