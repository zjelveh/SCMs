#' Performance Optimization Utilities
#' 
#' @title Performance Optimization Functions
#' @description Internal utilities for optimizing computational performance in synthetic control
#' estimation and specification curve analysis.
#' 
#' @name performance-utils
#' @keywords internal
NULL

#' Cached Matrix Computation Manager
#' 
#' @description Manages caching of expensive matrix computations to avoid redundant calculations
#' in specification curve analysis where many similar matrices are processed.
#' 
#' @keywords internal
MatrixCache <- R6::R6Class("MatrixCache",
  public = list(
    #' @field cache_store Internal hash table for storing computed matrices
    cache_store = NULL,
    
    #' @description Initialize the cache
    initialize = function() {
      self$cache_store <- new.env(hash = TRUE)
    },
    
    #' @description Generate cache key from matrix characteristics
    #' @param matrix Matrix to generate key for
    #' @param operation Character string describing the operation
    get_cache_key = function(matrix, operation) {
      # Create deterministic key based on matrix dimensions and checksum
      dims <- paste(dim(matrix), collapse = "x")
      checksum <- digest::digest(matrix, algo = "md5")
      return(paste(operation, dims, checksum, sep = "_"))
    },
    
    #' @description Store computation result in cache
    #' @param key Cache key
    #' @param result Computation result to store
    store_result = function(key, result) {
      self$cache_store[[key]] <- result
    },
    
    #' @description Retrieve cached result
    #' @param key Cache key
    get_result = function(key) {
      self$cache_store[[key]]
    },
    
    #' @description Check if result exists in cache
    #' @param key Cache key
    has_result = function(key) {
      exists(key, envir = self$cache_store)
    },
    
    #' @description Clear cache to free memory
    clear_cache = function() {
      rm(list = ls(self$cache_store), envir = self$cache_store)
    },
    
    #' @description Get cache statistics
    get_stats = function() {
      list(
        entries = length(ls(self$cache_store)),
        memory_usage = object.size(self$cache_store)
      )
    }
  )
)

#' Optimized Matrix Operations
#' 
#' @description High-performance implementations of common matrix operations
#' with caching and numerical stability improvements.
#' 
#' @param A Matrix A
#' @param B Matrix B
#' @param V Weighting matrix V
#' @param cache MatrixCache object (optional)
#' @keywords internal
optimized_matrix_operations <- function(A, B, V, cache = NULL) {
  
  # Check for valid inputs
  if (!is.matrix(A) || !is.matrix(B) || !is.matrix(V)) {
    stop("All inputs must be matrices")
  }
  
  if (nrow(A) != nrow(B) || nrow(A) != nrow(V) || ncol(V) != nrow(V)) {
    stop("Matrix dimensions incompatible")
  }
  
  # Use cache if provided
  use_cache <- !is.null(cache) && inherits(cache, "MatrixCache")
  
  # Optimized weighted squared error computation: ||A - B*w||_V^2
  compute_weighted_sse <- function(A, B, w, V) {
    residual <- A - B %*% w
    return(as.numeric(t(residual) %*% V %*% residual))
  }
  
  # Optimized Cholesky-based solve for positive definite systems
  compute_cholesky_solve <- function(M, rhs) {
    # Add regularization for numerical stability
    reg_param <- 1e-8
    M_reg <- M + reg_param * diag(nrow(M))
    
    tryCatch({
      chol_decomp <- chol(M_reg)
      return(backsolve(chol_decomp, forwardsolve(t(chol_decomp), rhs)))
    }, error = function(e) {
      # Fall back to standard solve if Cholesky fails
      return(solve(M_reg, rhs))
    })
  }
  
  # Pre-compute expensive operations if caching enabled
  if (use_cache) {
    # Cache V decomposition
    v_key <- cache$get_cache_key(V, "V_chol")
    if (cache$has_result(v_key)) {
      V_chol <- cache$get_result(v_key)
    } else {
      V_chol <- tryCatch({
        chol(V + 1e-8 * diag(nrow(V)))
      }, error = function(e) {
        NULL  # Use V directly if Cholesky fails
      })
      cache$store_result(v_key, V_chol)
    }
    
    # Cache B'VB computation for efficiency
    bvb_key <- cache$get_cache_key(B, paste("BVB", digest::digest(V), sep = "_"))
    if (cache$has_result(bvb_key)) {
      BVB <- cache$get_result(bvb_key)
    } else {
      BVB <- t(B) %*% V %*% B
      cache$store_result(bvb_key, BVB)
    }
  } else {
    # Compute without caching
    V_chol <- tryCatch({
      chol(V + 1e-8 * diag(nrow(V)))
    }, error = function(e) {
      NULL
    })
    BVB <- t(B) %*% V %*% B
  }
  
  return(list(
    compute_weighted_sse = compute_weighted_sse,
    compute_cholesky_solve = compute_cholesky_solve,
    V_chol = V_chol,
    BVB = BVB
  ))
}

#' Memory-Efficient Specification Curve Processing
#' 
#' @description Processes specification curves in batches to manage memory usage
#' for large specification spaces.
#' 
#' @param spec_combinations Data.table of specification combinations
#' @param batch_size Integer, number of specifications to process per batch
#' @param process_function Function to process each batch
#' @param ... Additional arguments passed to process_function
#' 
#' @return Combined results from all batches
#' @keywords internal
process_specifications_in_batches <- function(spec_combinations, batch_size = 100, 
                                             process_function, ...) {
  
  if (!is.data.table(spec_combinations)) {
    stop("spec_combinations must be a data.table")
  }
  
  if (!is.function(process_function)) {
    stop("process_function must be a function")
  }
  
  n_specs <- nrow(spec_combinations)
  n_batches <- ceiling(n_specs / batch_size)
  
  results <- vector("list", n_batches)
  
  for (i in seq_len(n_batches)) {
    start_idx <- (i - 1) * batch_size + 1
    end_idx <- min(i * batch_size, n_specs)
    
    batch_specs <- spec_combinations[start_idx:end_idx, ]
    
    # Process batch with error handling
    batch_result <- tryCatch({
      process_function(batch_specs, ...)
    }, error = function(e) {
      warning(paste("Batch", i, "failed:", e$message))
      return(NULL)
    })
    
    results[[i]] <- batch_result
    
    # Optional garbage collection after each batch
    if (i %% 10 == 0) {
      gc()
    }
  }
  
  # Combine non-NULL results
  valid_results <- results[!sapply(results, is.null)]
  
  if (length(valid_results) == 0) {
    stop("All specification batches failed")
  }
  
  return(rbindlist(valid_results, fill = TRUE))
}

#' Optimized Data Preparation for Large Panels
#' 
#' @description Efficient data processing for large panel datasets using data.table operations.
#' 
#' @param data Data.table with panel data
#' @param id_col Character, name of unit identifier column
#' @param time_col Character, name of time column  
#' @param outcome_col Character, name of outcome column
#' @param treated_unit Character, identifier of treated unit
#' @param pre_periods Numeric vector of pre-treatment periods
#' @param post_periods Numeric vector of post-treatment periods
#' 
#' @return List with processed panel data components
#' @keywords internal
optimize_panel_data_prep <- function(data, id_col, time_col, outcome_col, 
                                    treated_unit, pre_periods, post_periods) {
  
  # Convert to data.table if needed
  if (!is.data.table(data)) {
    data <- as.data.table(data)
  }
  
  # Pre-filter relevant periods for efficiency
  relevant_periods <- c(pre_periods, post_periods)
  data_subset <- data[get(time_col) %in% relevant_periods]
  
  # Create efficient matrices using data.table operations
  # Pre-treatment outcomes for treated unit
  Y_pre_treated <- data_subset[
    get(id_col) == treated_unit & get(time_col) %in% pre_periods,
    get(outcome_col)
  ]
  
  # Post-treatment outcomes for treated unit  
  Y_post_treated <- data_subset[
    get(id_col) == treated_unit & get(time_col) %in% post_periods,
    get(outcome_col)
  ]
  
  # Donor outcomes matrix (periods x donors)
  donor_units <- setdiff(unique(data_subset[[id_col]]), treated_unit)
  
  Y_donors_pre <- data_subset[
    get(id_col) %in% donor_units & get(time_col) %in% pre_periods
  ][
    order(get(time_col), get(id_col))
  ] |> 
  dcast(formula = as.formula(paste(time_col, "~", id_col)), 
        value.var = outcome_col)
  
  # Convert to matrix (remove time column)
  time_col_idx <- which(names(Y_donors_pre) == time_col)
  Y_donors_matrix <- as.matrix(Y_donors_pre[, -time_col_idx, with = FALSE])
  
  # Post-treatment donor matrix
  Y_donors_post <- data_subset[
    get(id_col) %in% donor_units & get(time_col) %in% post_periods
  ][
    order(get(time_col), get(id_col))
  ] |>
  dcast(formula = as.formula(paste(time_col, "~", id_col)),
        value.var = outcome_col)
  
  Y_donors_post_matrix <- as.matrix(Y_donors_post[, -time_col_idx, with = FALSE])
  
  return(list(
    Y_pre_treated = Y_pre_treated,
    Y_post_treated = Y_post_treated,
    Y_donors_pre = Y_donors_matrix,
    Y_donors_post = Y_donors_post_matrix,
    donor_units = donor_units,
    n_pre_periods = length(pre_periods),
    n_post_periods = length(post_periods)
  ))
}

#' Get Optimal Solver for Problem Characteristics
#' 
#' @description Selects the most efficient CVXR solver based on problem size and constraints.
#' 
#' @param n_vars Number of variables
#' @param constraint_type Type of constraint ("simplex", "lasso", "ridge", "ols")
#' @param available_solvers Character vector of available solvers
#' 
#' @return Character, recommended solver name
#' @keywords internal
get_optimal_solver <- function(n_vars, constraint_type, available_solvers = CVXR::installed_solvers()) {
  
  # Default preferences based on problem characteristics
  if (constraint_type == "lasso") {
    # OSQP is often better for L1 problems
    if ("OSQP" %in% available_solvers) {
      return("OSQP")
    }
  }
  
  if (constraint_type == "simplex" && n_vars <= 100) {
    # ECOS is good for small to medium simplex problems
    if ("ECOS" %in% available_solvers) {
      return("ECOS")
    }
  }
  
  if (constraint_type == "ridge" || constraint_type == "ols") {
    # SCS handles quadratic problems well
    if ("SCS" %in% available_solvers) {
      return("SCS") 
    }
  }
  
  # For large problems, prefer scalable solvers
  if (n_vars > 500) {
    if ("GUROBI" %in% available_solvers) {
      return("GUROBI")  # Commercial solver, best for large problems
    }
    if ("OSQP" %in% available_solvers) {
      return("OSQP")    # Good open-source option for large problems
    }
  }
  
  # Fallback to first available solver
  return(available_solvers[1])
}