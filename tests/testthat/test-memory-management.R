library(testthat)
library(data.table)

test_that("SpecResultsManager stores and retrieves results by spec_id", {
  storage_dir <- tempfile("scm_storage_")
  dir.create(storage_dir, recursive = TRUE)

  mgr <- SpecResultsManager$new(storage_dir = storage_dir, max_memory_mb = 1e6)

  res_a <- list(results = data.table(
    full_spec_id = c("s1", "s2"),
    outcome = "y",
    tau = c(1.0, 2.0),
    unit_type = "treated"
  ))
  res_b <- list(results = data.table(
    full_spec_id = c("s3", "s4", "s5"),
    outcome = "y",
    tau = c(-1.0, -2.0, -3.0),
    unit_type = "treated"
  ))

  mgr$store_results(res_a, "A")
  mgr$store_results(res_b, "B")

  metadata <- mgr$get_metadata()
  expect_setequal(metadata$spec_id, c("A", "B"))
  expect_equal(metadata[spec_id == "A", n_specifications], 2)
  expect_equal(metadata[spec_id == "B", n_specifications], 3)

  loaded_a <- mgr$retrieve_results("A", lazy = FALSE)
  loaded_b <- mgr$retrieve_results("B", lazy = FALSE)
  expect_equal(loaded_a$results$tau, c(1.0, 2.0))
  expect_equal(loaded_b$results$tau, c(-1.0, -2.0, -3.0))

  lazy_loader <- mgr$retrieve_results("A", lazy = TRUE)
  expect_true(is.function(lazy_loader))
  expect_equal(lazy_loader()$results$tau, c(1.0, 2.0))
})

test_that("SpecResultsManager updates an existing spec_id without duplication", {
  storage_dir <- tempfile("scm_storage_")
  dir.create(storage_dir, recursive = TRUE)

  mgr <- SpecResultsManager$new(storage_dir = storage_dir, max_memory_mb = 1e6)

  initial <- list(results = data.table(
    full_spec_id = c("s1", "s2"),
    outcome = "y",
    tau = c(1.0, 2.0),
    unit_type = "treated"
  ))
  updated <- list(results = data.table(
    full_spec_id = c("s1", "s2", "s3"),
    outcome = "y",
    tau = c(10.0, 20.0, 30.0),
    unit_type = "treated"
  ))

  mgr$store_results(initial, "A")
  mgr$store_results(updated, "A")

  metadata <- mgr$get_metadata()
  expect_equal(nrow(metadata[spec_id == "A"]), 1)
  expect_equal(metadata[spec_id == "A", n_specifications], 3)

  loaded <- mgr$retrieve_results("A", lazy = FALSE)
  expect_equal(loaded$results$tau, c(10.0, 20.0, 30.0))
})

test_that("SpecResultsManager filtering and missing-spec errors behave as expected", {
  storage_dir <- tempfile("scm_storage_")
  dir.create(storage_dir, recursive = TRUE)

  mgr <- SpecResultsManager$new(storage_dir = storage_dir, max_memory_mb = 1e6)

  res <- list(results = data.table(
    full_spec_id = c("s1", "s2", "s3"),
    outcome = "y",
    tau = c(1.0, -1.0, 3.0),
    unit_type = c("treated", "control", "treated")
  ))
  mgr$store_results(res, "A")

  filtered <- NULL
  expect_warning(
    filtered <- mgr$get_filtered_results(
      "A",
      filter_function = function(dt) dt$unit_type == "treated",
      select_cols = c("full_spec_id", "tau", "missing_col")
    ),
    "Missing columns"
  )

  expect_equal(nrow(filtered), 2)
  expect_equal(names(filtered), c("full_spec_id", "tau"))

  expect_error(
    mgr$retrieve_results("does_not_exist", lazy = FALSE),
    "not found"
  )
})
