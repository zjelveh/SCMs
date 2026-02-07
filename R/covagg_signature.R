#' Compile Covagg Signature to Internal Covariate Specs
#'
#' @description Internal utility that compiles a parametric covagg signature
#' into the internal per-variable/per-period specifications used by SCMs.
#'
#' @param covagg List of operation specifications.
#' @param outcome_var Character. Outcome variable name for resolving "outcome_var".
#' @param period_pre Numeric/integer vector of pre-treatment periods.
#' @param require_named Logical. Whether top-level entries must be named.
#' @param context Character. Context string used in error messages.
#'
#' @return Named list of compiled internal covariate specs with fields
#'   `var`, `periods`, and `agg_fun`.
#' @keywords internal
compile_covagg_signature <- function(covagg,
                                     outcome_var,
                                     period_pre,
                                     require_named = FALSE,
                                     context = "covagg") {
  if (!is.list(covagg)) {
    stop(context, " must be a list of operation specifications")
  }

  period_pre <- sort(unique(period_pre))
  if (length(period_pre) == 0) {
    stop("period_pre must contain at least one period")
  }

  covagg_names <- names(covagg)
  if (require_named && (is.null(covagg_names) || any(is.na(covagg_names) | covagg_names == ""))) {
    stop(context, " entries must be named")
  }

  compiled <- list()

  for (i in seq_along(covagg)) {
    op <- covagg[[i]]
    op_name <- if (!is.null(covagg_names) && length(covagg_names) >= i && !is.na(covagg_names[[i]]) && nzchar(covagg_names[[i]])) {
      covagg_names[[i]]
    } else {
      paste0("op_", i)
    }

    validate_covagg_operation(op, op_name)

    var_name <- if (identical(op$var, "outcome_var")) outcome_var else op$var

    selected_periods <- resolve_covagg_select_periods(op$select_periods, period_pre, op_name)
    groups <- resolve_covagg_partitions(op$partition_periods, selected_periods, op_name)
    compute <- resolve_covagg_compute(op$compute, op_name)

    for (g in seq_along(groups)) {
      group_periods <- sort(unique(groups[[g]]))
      if (length(group_periods) == 0) next

      spec_name <- make.unique(c(names(compiled), paste0(op_name, "__g", g)))[length(compiled) + 1]

      compiled[[spec_name]] <- list(
        var = var_name,
        periods = group_periods,
        agg_fun = compute$raw
      )
    }
  }

  compiled
}

#' @keywords internal
validate_covagg_operation <- function(op, op_name) {
  if (!is.list(op)) {
    stop("Operation '", op_name, "' must be a list")
  }

  allowed <- c("var", "select_periods", "partition_periods", "compute", "label")
  unknown <- setdiff(names(op), allowed)
  if (length(unknown) > 0) {
    stop("Operation '", op_name, "' has unsupported fields: ", paste(unknown, collapse = ", "))
  }

  if (is.null(op$var) || !is.character(op$var) || length(op$var) != 1) {
    stop("Operation '", op_name, "' must define var as a single character value")
  }
}

#' @keywords internal
resolve_covagg_select_periods <- function(select_spec, period_pre, op_name) {
  if (is.null(select_spec)) {
    return(period_pre)
  }

  if (is.numeric(select_spec) || is.integer(select_spec)) {
    return(validate_selected_periods(select_spec, period_pre, op_name))
  }

  if (is.function(select_spec)) {
    out <- select_spec(period_pre)
    return(resolve_selector_output(out, period_pre, op_name))
  }

  if (!is.list(select_spec)) {
    stop("Operation '", op_name, "': select_periods must be NULL, numeric vector, function, or list")
  }

  type <- select_spec$type %||% "all_pre"

  if (identical(type, "all_pre")) {
    return(period_pre)
  }

  if (identical(type, "explicit")) {
    if (is.null(select_spec$periods)) {
      stop("Operation '", op_name, "': select_periods$type='explicit' requires periods")
    }
    return(validate_selected_periods(select_spec$periods, period_pre, op_name))
  }

  if (identical(type, "every_n")) {
    n <- select_spec$n %||% NA
    offset <- select_spec$offset %||% 1

    if (!is.numeric(n) || length(n) != 1 || n <= 0 || n != as.integer(n)) {
      stop("Operation '", op_name, "': select_periods$type='every_n' requires integer n >= 1")
    }
    if (!is.numeric(offset) || length(offset) != 1 || offset <= 0 || offset != as.integer(offset)) {
      stop("Operation '", op_name, "': select_periods$type='every_n' requires integer offset >= 1")
    }

    idx <- seq(from = offset, to = length(period_pre), by = n)
    return(period_pre[idx])
  }

  if (identical(type, "predicate")) {
    fn <- select_spec$fn
    if (!is.function(fn)) {
      stop("Operation '", op_name, "': select_periods$type='predicate' requires function fn")
    }
    out <- fn(period_pre)
    return(resolve_selector_output(out, period_pre, op_name))
  }

  stop("Operation '", op_name, "': unsupported select_periods type '", type, "'")
}

#' @keywords internal
resolve_selector_output <- function(out, period_pre, op_name) {
  if (is.logical(out)) {
    if (length(out) != length(period_pre)) {
      stop("Operation '", op_name, "': predicate selector returned logical vector of wrong length")
    }
    selected <- period_pre[out]
  } else if (is.numeric(out) || is.integer(out)) {
    selected <- out
  } else {
    stop("Operation '", op_name, "': selector output must be logical or numeric")
  }

  validate_selected_periods(selected, period_pre, op_name)
}

#' @keywords internal
validate_selected_periods <- function(selected, period_pre, op_name) {
  if (!is.numeric(selected) && !is.integer(selected)) {
    stop("Operation '", op_name, "': selected periods must be numeric")
  }

  selected <- sort(unique(selected))
  if (length(selected) == 0) {
    stop("Operation '", op_name, "': period selection is empty")
  }

  missing <- setdiff(selected, period_pre)
  if (length(missing) > 0) {
    stop("Operation '", op_name, "': selected periods are outside pre-period window: ",
         paste(missing, collapse = ", "))
  }

  selected
}

#' @keywords internal
resolve_covagg_partitions <- function(partition_spec, selected_periods, op_name) {
  if (is.null(partition_spec)) {
    return(list(selected_periods))
  }

  if (is.function(partition_spec)) {
    groups <- partition_spec(selected_periods)
    return(validate_partition_groups(groups, selected_periods, op_name))
  }

  if (is.character(partition_spec) && length(partition_spec) == 1) {
    partition_spec <- list(type = partition_spec)
  }

  if (!is.list(partition_spec)) {
    stop("Operation '", op_name, "': partition_periods must be NULL, function, character, or list")
  }

  type <- partition_spec$type %||% "all"

  if (identical(type, "all")) {
    return(list(selected_periods))
  }

  if (identical(type, "by_period")) {
    return(lapply(selected_periods, function(p) p))
  }

  if (identical(type, "fixed_width")) {
    width <- partition_spec$width %||% NA
    step <- partition_spec$step %||% width

    if (!is.numeric(width) || length(width) != 1 || width <= 0 || width != as.integer(width)) {
      stop("Operation '", op_name, "': partition_periods$type='fixed_width' requires integer width >= 1")
    }
    if (!is.numeric(step) || length(step) != 1 || step <= 0 || step != as.integer(step)) {
      stop("Operation '", op_name, "': partition_periods$type='fixed_width' requires integer step >= 1")
    }

    starts <- seq(1, length(selected_periods), by = step)
    groups <- lapply(starts, function(i) selected_periods[i:min(i + width - 1, length(selected_periods))])
    return(validate_partition_groups(groups, selected_periods, op_name))
  }

  if (identical(type, "custom")) {
    groups <- partition_spec$groups
    return(validate_partition_groups(groups, selected_periods, op_name))
  }

  stop("Operation '", op_name, "': unsupported partition_periods type '", type, "'")
}

#' @keywords internal
validate_partition_groups <- function(groups, selected_periods, op_name) {
  if (!is.list(groups) || length(groups) == 0) {
    stop("Operation '", op_name, "': partition groups must be a non-empty list")
  }

  norm_groups <- lapply(groups, function(g) sort(unique(g)))

  for (i in seq_along(norm_groups)) {
    g <- norm_groups[[i]]
    if (!is.numeric(g) || length(g) == 0) {
      stop("Operation '", op_name, "': partition group ", i, " must be a non-empty numeric vector")
    }
    outside <- setdiff(g, selected_periods)
    if (length(outside) > 0) {
      stop("Operation '", op_name, "': partition group ", i,
           " contains periods outside selected set: ", paste(outside, collapse = ", "))
    }
  }

  norm_groups
}

#' @keywords internal
resolve_covagg_compute <- function(compute_spec, op_name) {
  raw <- compute_spec

  if (is.null(raw)) {
    raw <- "mean"
  }

  if (is.list(raw)) {
    raw <- raw$fn
  }

  if (is.null(raw)) {
    raw <- "mean"
  }

  fn <- if (is.character(raw)) match.fun(raw) else raw

  if (!is.function(fn)) {
    stop("Operation '", op_name, "': compute must be a function or function name")
  }

  # basic contract check
  test <- tryCatch(fn(c(1, 2, 3)), error = function(e) e)
  if (inherits(test, "error") || length(test) != 1 || !is.numeric(test)) {
    stop("Operation '", op_name, "': compute must return a single numeric value")
  }

  list(fn = fn, raw = raw)
}

#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
