make_rejection_state <- function(path, model) {
  state <- new.env(parent = globalenv())
  state$model_def <- NULL
  state$model_list <- list(m1 = model)
  state$prior_dist <- list(m1 = list(c("x", "unif", 0, 1)))
  state$ss_obs <- 0
  state$tmp_object_store_root <- file.path(dirname(path), "objects")
  state$store_summaries <- "none"
  state$store_outputs <- "none"
  state$dist_names <- "dist1"
  state$model_names <- "m1"
  state$param_names <- "x"
  state$column_names <- c(
    "attempt_id", "job_id", "accepted", "retained", "model", "x", "dist1"
  )
  state$thresholds <- 1
  save(list = ls(state), envir = state, file = path)
}

test_that("a batch publishes its control file only after its particles", {
  root <- file.path(tempdir(), paste0("batch-", sample.int(1e8, 1)))
  dir.create(root, recursive = TRUE)
  state_path <- file.path(root, "state.RData")
  make_rejection_state(state_path, function(x, obs) 0)
  spec <- BRREWABC:::newBatchSpec(2, 0, 1, 3, seed = 42)

  result <- BRREWABC:::runABCBatch(
    "rejection", spec, state_path, file.path(root, "batch-2")
  )
  particles <- readRDS(result$particles_path)

  expect_identical(result$status, "ok")
  expect_equal(result$n_attempted, 3)
  expect_equal(result$n_accepted, 3)
  expect_true(file.exists(file.path(root, "batch-2", "result.rds")))
  expect_equal(length(unique(particles$attempt_id)), 3)
  expect_match(particles$attempt_id, "^g0000-b00000002-a")
})

test_that("a failed batch writes an error control record", {
  root <- file.path(tempdir(), paste0("batch-error-", sample.int(1e8, 1)))
  dir.create(root, recursive = TRUE)
  state_path <- file.path(root, "state.RData")
  make_rejection_state(state_path, function(x, obs) stop("model failed"))
  spec <- BRREWABC:::newBatchSpec(1, 0, 1, 2, seed = 42)
  batch_dir <- file.path(root, "batch-1")

  result <- BRREWABC:::runABCBatch(
    "rejection", spec, state_path, batch_dir
  )

  expect_identical(result$status, "error")
  expect_match(result$message, "model failed")
  expect_true(file.exists(file.path(batch_dir, "result.rds")))
  expect_false(file.exists(file.path(batch_dir, "particles.rds")))
})

test_that("batch size must be a positive integer", {
  expect_equal(BRREWABC:::validateBatchSize(50), 50L)
  expect_error(BRREWABC:::validateBatchSize(0), "positive integer")
  expect_error(BRREWABC:::validateBatchSize(1.5), "positive integer")
  expect_equal(
    BRREWABC:::validatePositiveInteger(2, "workers"), 2L
  )
  expect_error(
    BRREWABC:::validatePositiveInteger(0, "workers"), "workers.*positive"
  )
})
