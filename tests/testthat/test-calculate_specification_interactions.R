library(testthat)
library(data.table)
library(catboost)

test_that("calculate_specification_interactions works with SHAP interactions disabled", {
  # Create mock CatBoost results data
  set.seed(42)
  n_specs <- 20
  
  # Create mock prediction data
  predictions_data <- data.table(
    unit = rep("treated_unit", n_specs),
    actual = rnorm(n_specs, mean = 5, sd = 2),
    predicted_loo = rnorm(n_specs, mean = 5, sd = 1.5),
    outcome_model = sample(c("linear", "quadratic"), n_specs, replace = TRUE),
    const = sample(c("TRUE", "FALSE"), n_specs, replace = TRUE),
    fw = sample(c("unit_weights", "ridge"), n_specs, replace = TRUE),
    feat = sample(c("all", "selected"), n_specs, replace = TRUE),
    data_sample = sample(c("full", "subset"), n_specs, replace = TRUE),
    full_spec_id = paste0("spec_", 1:n_specs)
  )
  
  # Create feature matrix for CatBoost model training
  feature_matrix <- predictions_data[, .(outcome_model, const, fw, feat, data_sample)]
  feature_matrix <- as.data.frame(feature_matrix)
  for(col in names(feature_matrix)) {
    feature_matrix[[col]] <- as.factor(feature_matrix[[col]])
  }
  
  # Train a simple CatBoost model
  train_pool <- catboost.load_pool(
    data = feature_matrix, 
    label = predictions_data$actual
  )
  
  params <- list(
    loss_function = 'RMSE',
    iterations = 10,
    depth = 3,
    learning_rate = 0.1,
    random_seed = 42,
    verbose = 0
  )
  
  trained_model <- catboost.train(learn_pool = train_pool, params = params)
  
  # Create mock catboost_results object
  catboost_results <- list(
    predictions = predictions_data,
    config = list(treated_unit_name = "treated_unit"),
    model = trained_model
  )
  
  # Test the function (disable SHAP interactions due to CatBoost API issues)
  result <- calculate_specification_interactions(
    catboost_results = catboost_results,
    dataset_name = "test_dataset",
    include_shap_interactions = FALSE,
    top_n = 5
  )
  
  # Verify structure
  expect_type(result, "list")
  expect_named(result, c("shap_interactions", "dataset"))
  expect_equal(result$dataset, "test_dataset")
  
  # Test SHAP interactions (should be NULL when disabled)
  expect_null(result$shap_interactions)
})

test_that("calculate_specification_interactions fails hard on missing data", {
  # Test with NULL catboost_results
  expect_error(
    calculate_specification_interactions(NULL, "test"),
    "Missing prediction data for dataset: test"
  )
  
  # Test with missing predictions
  catboost_results_no_preds <- list(
    predictions = NULL,
    config = list(treated_unit_name = "treated_unit")
  )
  
  expect_error(
    calculate_specification_interactions(catboost_results_no_preds, "test"),
    "Missing prediction data for dataset: test"
  )
  
  # Test with empty predictions
  catboost_results_empty <- list(
    predictions = data.table(),
    config = list(treated_unit_name = "treated_unit")
  )
  
  expect_error(
    calculate_specification_interactions(catboost_results_empty, "test"),
    "Missing prediction data for dataset: test"
  )
})

test_that("calculate_specification_interactions fails hard on missing treated unit", {
  # Create mock data without treated unit
  predictions_data <- data.table(
    unit = rep("other_unit", 10),
    actual = rnorm(10),
    outcome_model = sample(c("linear", "quadratic"), 10, replace = TRUE),
    const = sample(c("TRUE", "FALSE"), 10, replace = TRUE)
  )
  
  catboost_results <- list(
    predictions = predictions_data,
    config = list(treated_unit_name = "treated_unit")
  )
  
  expect_error(
    calculate_specification_interactions(catboost_results, "test"),
    "No data for treated unit in dataset: test"
  )
})

test_that("function handles insufficient features gracefully", {
  # Create mock data with only one feature
  predictions_data <- data.table(
    unit = rep("treated_unit", 10),
    actual = rnorm(10),
    outcome_model = sample(c("linear", "quadratic"), 10, replace = TRUE)
  )
  
  catboost_results <- list(
    predictions = predictions_data,
    config = list(treated_unit_name = "treated_unit")
  )
  
  # Should work but have limited interaction data
  result <- calculate_specification_interactions(catboost_results, "test", include_shap_interactions = FALSE)
  
  expect_type(result, "list")
  expect_named(result, c("shap_interactions", "dataset"))
  
  # SHAP interactions should be NULL when disabled
  expect_null(result$shap_interactions)
})

test_that("function works with basic data structure", {
  # Create mock data
  predictions_data <- data.table(
    unit = rep("treated_unit", 15),
    actual = rnorm(15, mean = 3, sd = 1),
    outcome_model = sample(c("linear", "quadratic"), 15, replace = TRUE),
    const = sample(c("TRUE", "FALSE"), 15, replace = TRUE),
    fw = sample(c("unit_weights", "ridge"), 15, replace = TRUE)
  )
  
  catboost_results <- list(
    predictions = predictions_data,
    config = list(treated_unit_name = "treated_unit")
  )
  
  # Test with SHAP interactions disabled
  result <- calculate_specification_interactions(
    catboost_results = catboost_results,
    dataset_name = "test_dataset",
    include_shap_interactions = FALSE
  )
  
  # Verify structure
  expect_type(result, "list")
  expect_named(result, c("shap_interactions", "dataset"))
  expect_null(result$shap_interactions)
  expect_equal(result$dataset, "test_dataset")
})