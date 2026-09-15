test_that("legacy and structured model results are accepted", {
  legacy <- BRREWABC:::normalizeModelResult(c(0.2), "dist1")
  expect_equal(legacy$distances, c(dist1 = 0.2))
  expect_length(legacy$summaries, 0)
  expect_length(legacy$outputs, 0)

  trajectory <- data.frame(timestep = 1:2, S = c(10, 9), I = c(0, 1))
  structured <- BRREWABC:::normalizeModelResult(
    list(
      distances = c(dist1 = 0.3),
      summaries = list(trajectory = trajectory, peak = 1),
      outputs = list(state = trajectory)
    ),
    "dist1"
  )
  expect_equal(structured$summaries$trajectory, trajectory)
  expect_equal(structured$outputs$state, trajectory)
})

test_that("vectors and tables are normalized for storage", {
  vector <- BRREWABC:::asStoredTable(c(a = 2, b = 3))
  expect_equal(vector$element_index, 1:2)
  expect_equal(vector$element_name, c("a", "b"))
  expect_equal(vector$value, c(2, 3))

  table <- data.frame(timestep = 1:2, pop_id = c("A", "B"), S = c(10, 20))
  expect_equal(BRREWABC:::asStoredTable(table), table)
})

test_that("Parquet data can be listed, filtered, and read", {
  staging <- file.path(tempdir(), paste0("stage-", sample.int(1e8, 1)))
  storage <- file.path(tempdir(), paste0("store-", sample.int(1e8, 1)))
  trajectory <- data.frame(timestep = 1:2, pop_id = c("A", "B"), S = c(10, 20))

  BRREWABC:::stageStoredCollection(
    list(trajectory = trajectory), "summaries", "all", TRUE,
    "a1", 1, 1, staging
  )
  BRREWABC:::stageStoredCollection(
    list(trajectory = trajectory), "summaries", "all", FALSE,
    "a2", 1, 1, staging
  )
  BRREWABC:::stageStoredCollection(
    list(state = trajectory), "outputs", "all", TRUE,
    "a1", 1, 1, staging
  )
  BRREWABC:::persistStoredGeneration(staging, storage, 1, "a1", "all", "all")
  result <- list(storage = list(path = storage))

  expect_equal(sort(list_abc_stored_data(result)$kind), c("outputs", "summaries"))
  expect_equal(unique(read_summary_statistics(result, status = "retained")$attempt_id), "a1")
  expect_equal(unique(read_summary_statistics(result, status = "rejected")$attempt_id), "a2")
  expect_equal(unique(read_summary_statistics(result, attempt_id = "a2")$attempt_id), "a2")
  expect_equal(unique(read_model_outputs(result)$stored_name), "state")
})

test_that("a changing table schema is rejected", {
  staging <- file.path(tempdir(), paste0("stage-schema-", sample.int(1e8, 1)))
  storage <- file.path(tempdir(), paste0("store-schema-", sample.int(1e8, 1)))
  BRREWABC:::stageStoredCollection(
    list(trajectory = data.frame(timestep = 1, S = 10)),
    "summaries", "all", TRUE, "a1", 1, 1, staging
  )
  BRREWABC:::stageStoredCollection(
    list(trajectory = data.frame(timestep = 1, I = 2)),
    "summaries", "all", TRUE, "a2", 1, 1, staging
  )
  expect_error(
    BRREWABC:::persistStoredGeneration(staging, storage, 1, c("a1", "a2"), "all", "none"),
    "schema"
  )
})
