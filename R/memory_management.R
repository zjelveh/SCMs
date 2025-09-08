#' Memory Management for Large Specification Curves
#' 
#' @title Memory-Efficient Specification Curve Storage and Processing
#' @description Implements efficient storage and retrieval systems for large specification curve
#' analyses that may exceed memory limits. Uses on-disk storage, compression, and lazy loading.
#'
#' @name memory-management
#' @keywords internal
NULL

#' Specification Results Storage Manager
#' 
#' @description R6 class for managing large specification curve results with disk-based storage,
#' compression, and memory-efficient retrieval patterns.
#' 
#' @export
SpecResultsManager <- R6::R6Class("SpecResultsManager",
  public = list(
    #' @field storage_dir Directory for storing specification results
    storage_dir = NULL,
    
    #' @field compression_level Compression level (0-9, higher = more compression)
    compression_level = 6,
    
    #' @field max_memory_mb Maximum memory usage in MB before offloading to disk
    max_memory_mb = 1000,
    
    #' @field current_memory_usage Current memory usage tracking
    current_memory_usage = 0,
    
    #' @description Initialize storage manager
    #' @param storage_dir Character, directory path for storage
    #' @param max_memory_mb Numeric, maximum memory in MB
    #' @param compression_level Integer 0-9, compression level
    initialize = function(storage_dir = tempdir(), max_memory_mb = 1000, compression_level = 6) {
      self$storage_dir <- storage_dir
      self$max_memory_mb <- max_memory_mb
      self$compression_level <- compression_level
      
      # Create storage directory if needed
      if (!dir.exists(self$storage_dir)) {
        dir.create(self$storage_dir, recursive = TRUE)
      }
      
      # Initialize metadata file
      metadata_file <- file.path(self$storage_dir, "spec_metadata.rds")
      if (!file.exists(metadata_file)) {
        metadata <- data.table(
          spec_id = character(0),
          file_path = character(0),
          size_bytes = numeric(0),
          created = character(0),
          n_specifications = numeric(0),
          outcome_vars = character(0)
        )
        saveRDS(metadata, metadata_file, compress = TRUE)
      }
    },
    
    #' @description Store specification results efficiently
    #' @param spec_results List or data.table with specification results
    #' @param spec_id Character, unique identifier for this specification set
    #' @param metadata List with additional metadata
    store_results = function(spec_results, spec_id, metadata = list()) {
      
      # Generate file path
      file_name <- paste0("spec_", spec_id, ".rds")
      file_path <- file.path(self$storage_dir, file_name)
      
      # Compress and save
      saveRDS(spec_results, file_path, compress = TRUE)
      
      # Update metadata
      self$update_metadata(spec_id, file_path, spec_results, metadata)
      
      # Check memory usage and trigger cleanup if needed
      self$check_memory_usage()
      
      return(invisible(file_path))
    },
    
    #' @description Retrieve specification results
    #' @param spec_id Character, specification set identifier
    #' @param lazy Logical, whether to use lazy loading for large results
    retrieve_results = function(spec_id, lazy = TRUE) {
      
      metadata <- self$get_metadata()
      file_info <- metadata[spec_id == spec_id]
      
      if (nrow(file_info) == 0) {
        stop(paste("Specification set", spec_id, "not found"))
      }
      
      file_path <- file_info$file_path[1]
      
      if (!file.exists(file_path)) {
        stop(paste("Results file not found:", file_path))
      }
      
      if (lazy) {
        # Return a lazy loader function instead of the data
        return(function() readRDS(file_path))
      } else {
        return(readRDS(file_path))
      }
    },
    
    #' @description Get filtered subset of results without loading all data
    #' @param spec_id Character, specification set identifier  
    #' @param filter_function Function to apply for filtering
    #' @param select_cols Character vector of columns to select
    get_filtered_results = function(spec_id, filter_function = NULL, select_cols = NULL) {
      
      # Load results
      results <- self$retrieve_results(spec_id, lazy = FALSE)
      
      if (is.data.table(results$results)) {
        dt <- results$results
      } else if (is.data.frame(results$results)) {
        dt <- as.data.table(results$results)
      } else {
        stop("Results format not supported for filtering")
      }
      
      # Apply filter
      if (!is.null(filter_function)) {
        dt <- dt[filter_function(dt)]
      }
      
      # Select columns
      if (!is.null(select_cols)) {
        missing_cols <- setdiff(select_cols, names(dt))
        if (length(missing_cols) > 0) {
          warning(paste("Missing columns:", paste(missing_cols, collapse = ", ")))
        }
        available_cols <- intersect(select_cols, names(dt))
        dt <- dt[, ..available_cols]
      }
      
      return(dt)
    },
    
    #' @description Update metadata for stored results
    #' @param spec_id Character, specification identifier
    #' @param file_path Character, path to results file
    #' @param spec_results The results object
    #' @param additional_metadata List with additional metadata
    update_metadata = function(spec_id, file_path, spec_results, additional_metadata = list()) {
      
      metadata_file <- file.path(self$storage_dir, "spec_metadata.rds")
      metadata <- readRDS(metadata_file)
      
      # Get file size
      file_size <- file.info(file_path)$size
      
      # Extract relevant info from results
      n_specs <- 0
      outcome_vars <- ""
      
      if (is.list(spec_results) && "results" %in% names(spec_results)) {
        if (is.data.table(spec_results$results) || is.data.frame(spec_results$results)) {
          n_specs <- nrow(spec_results$results)
          if ("outcome" %in% names(spec_results$results)) {
            outcome_vars <- paste(unique(spec_results$results$outcome), collapse = ",")
          }
        }
      }
      
      # Create new metadata row
      new_row <- data.table(
        spec_id = spec_id,
        file_path = file_path,
        size_bytes = file_size,
        created = as.character(Sys.time()),
        n_specifications = n_specs,
        outcome_vars = outcome_vars
      )
      
      # Remove existing entry for this spec_id if exists
      metadata <- metadata[spec_id != spec_id]
      
      # Add new row
      metadata <- rbind(metadata, new_row, fill = TRUE)
      
      # Save updated metadata
      saveRDS(metadata, metadata_file, compress = TRUE)
    },
    
    #' @description Get metadata for all stored results
    get_metadata = function() {
      metadata_file <- file.path(self$storage_dir, "spec_metadata.rds")
      return(readRDS(metadata_file))
    },
    
    #' @description Check current memory usage and trigger cleanup if needed
    check_memory_usage = function() {
      
      # Get current R memory usage
      mem_usage <- as.numeric(object.size(.GlobalEnv)) / 1024^2  # Convert to MB
      
      if (mem_usage > self$max_memory_mb) {
        self$cleanup_memory()
      }
    },
    
    #' @description Clean up memory by removing old cached objects
    cleanup_memory = function() {
      
      # Force garbage collection
      gc()
      
      # Clean up temporary objects in global environment if they look like spec curve objects
      global_objects <- ls(.GlobalEnv)
      spec_objects <- grep("spec_|result_|cache_", global_objects, value = TRUE)
      
      for (obj in spec_objects) {
        if (exists(obj, .GlobalEnv)) {
          obj_size <- object.size(get(obj, .GlobalEnv)) / 1024^2
          if (obj_size > 100) {  # Remove objects larger than 100MB
            rm(list = obj, envir = .GlobalEnv)
          }
        }
      }
      
      gc()
    },
    
    #' @description Get summary statistics about stored specifications
    get_storage_summary = function() {
      metadata <- self$get_metadata()
      
      if (nrow(metadata) == 0) {
        return(list(
          n_specification_sets = 0,
          total_specifications = 0,
          total_size_mb = 0,
          unique_outcomes = 0
        ))
      }
      
      # Calculate summary statistics
      total_size_bytes <- sum(metadata$size_bytes, na.rm = TRUE)
      total_specs <- sum(metadata$n_specifications, na.rm = TRUE)
      
      # Get unique outcomes
      all_outcomes <- unlist(strsplit(metadata$outcome_vars, ","))
      unique_outcomes <- length(unique(all_outcomes[all_outcomes != ""]))
      
      return(list(
        n_specification_sets = nrow(metadata),
        total_specifications = total_specs,
        total_size_mb = round(total_size_bytes / 1024^2, 2),
        unique_outcomes = unique_outcomes,
        storage_directory = self$storage_dir,
        oldest_entry = min(metadata$created, na.rm = TRUE),
        newest_entry = max(metadata$created, na.rm = TRUE)
      ))
    },
    
    #' @description Clean up storage directory
    #' @param older_than_days Remove results older than this many days
    cleanup_storage = function(older_than_days = 30) {
      
      metadata <- self$get_metadata()
      
      if (nrow(metadata) == 0) return(invisible(NULL))
      
      # Convert creation dates
      creation_dates <- as.POSIXct(metadata$created)
      cutoff_date <- Sys.time() - as.difftime(older_than_days, units = "days")
      
      # Identify old files
      old_files <- metadata[creation_dates < cutoff_date]
      
      if (nrow(old_files) > 0) {
        # Remove files
        for (i in seq_len(nrow(old_files))) {
          file_path <- old_files$file_path[i]
          if (file.exists(file_path)) {
            file.remove(file_path)
          }
        }
        
        # Update metadata
        updated_metadata <- metadata[creation_dates >= cutoff_date]
        metadata_file <- file.path(self$storage_dir, "spec_metadata.rds")
        saveRDS(updated_metadata, metadata_file, compress = TRUE)
        
        cat(paste("Removed", nrow(old_files), "old specification sets\n"))
      }
    },
    
    #' @description Merge multiple specification result sets efficiently
    #' @param spec_ids Character vector of specification IDs to merge
    #' @param output_spec_id Character, ID for merged results
    merge_specification_sets = function(spec_ids, output_spec_id) {
      
      if (length(spec_ids) < 2) {
        stop("Need at least 2 specification sets to merge")
      }
      
      # Load all results
      all_results <- list()
      for (i in seq_along(spec_ids)) {
        results <- self$retrieve_results(spec_ids[i], lazy = FALSE)
        all_results[[i]] <- results
      }
      
      # Merge results tables
      merged_results_data <- rbindlist(
        lapply(all_results, function(x) {
          if (is.data.table(x$results)) {
            return(x$results)
          } else if (is.data.frame(x$results)) {
            return(as.data.table(x$results))
          } else {
            return(NULL)
          }
        }),
        fill = TRUE
      )
      
      # Merge inference results if present
      merged_abadie <- NULL
      merged_bootstrap <- NULL
      
      if (all(sapply(all_results, function(x) "abadie_inference" %in% names(x)))) {
        # Merge Abadie inference
        abadie_results <- lapply(all_results, function(x) x$abadie_inference)
        # Merge p_values if present
        if (all(sapply(abadie_results, function(x) "p_values" %in% names(x)))) {
          merged_abadie <- list(
            p_values = rbindlist(lapply(abadie_results, function(x) x$p_values), fill = TRUE)
          )
        }
      }
      
      if (all(sapply(all_results, function(x) "bootstrap_inference" %in% names(x)))) {
        # Merge bootstrap inference
        bootstrap_results <- lapply(all_results, function(x) x$bootstrap_inference)
        if (all(sapply(bootstrap_results, function(x) "p_values" %in% names(x)))) {
          merged_bootstrap <- list(
            p_values = rbindlist(lapply(bootstrap_results, function(x) x$p_values), fill = TRUE)
          )
        }
      }
      
      # Create merged specification result
      merged_spec <- list(
        results = merged_results_data,
        expected_direction = all_results[[1]]$expected_direction
      )
      
      if (!is.null(merged_abadie)) {
        merged_spec$abadie_inference <- merged_abadie
      }
      
      if (!is.null(merged_bootstrap)) {
        merged_spec$bootstrap_inference <- merged_bootstrap
      }
      
      # Store merged results
      self$store_results(merged_spec, output_spec_id, 
                        metadata = list(merged_from = paste(spec_ids, collapse = ",")))
      
      return(output_spec_id)
    }
  )
)