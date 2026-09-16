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

test_that("Parquet fragments become visible only after a complete write", {
  directory <- file.path(tempdir(), paste0("atomic-", sample.int(1e8, 1)))
  dir.create(directory, recursive = TRUE)
  path <- file.path(directory, "fragment.parquet")

  BRREWABC:::writeParquetAtomically(data.frame(value = 1:3), path)

  expect_true(file.exists(path))
  expect_gt(file.info(path)$size, 0)
  expect_equal(as.data.frame(arrow::read_parquet(path))$value, 1:3)
  expect_length(list.files(directory, pattern = "\\.tmp-"), 0)
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
  BRREWABC:::stageStoredCollection(
    list(trajectory = trajectory), "summaries", "all", TRUE,
    "orphan", 1, 1, staging
  )
  BRREWABC:::persistStoredGeneration(
    staging, storage, 1, c("a1", "a2"), "a1", "all", "all"
  )
  result <- list(storage = list(path = storage))

  expect_equal(sort(list_abc_stored_data(result)$kind), c("outputs", "summaries"))
  expect_equal(unique(read_summary_statistics(result, status = "retained")$attempt_id), "a1")
  expect_equal(unique(read_summary_statistics(result, status = "rejected")$attempt_id), "a2")
  expect_equal(unique(read_summary_statistics(result, attempt_id = "a2")$attempt_id), "a2")
  expect_false("orphan" %in% read_summary_statistics(result)$attempt_id)
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
    BRREWABC:::persistStoredGeneration(
      staging, storage, 1, c("a1", "a2"), c("a1", "a2"), "all", "none"
    ),
    "schema"
  )
})

test_that("batch fragments are bounded, persisted, and consolidated", {
  staging <- file.path(tempdir(), paste0("stage-batch-", sample.int(1e8, 1)))
  storage <- file.path(tempdir(), paste0("store-batch-", sample.int(1e8, 1)))
  spec <- BRREWABC:::newBatchSpec(1, 2, 1, 3, seed = 42)
  buffer <- BRREWABC:::newStoredBuffer(
    spec, staging, chunk_rows = 2, chunk_mb = 1
  )
  for (attempt in 1:3) {
    BRREWABC:::appendStoredCollection(
      buffer,
      list(trajectory = data.frame(time = 1:2, value = attempt + 1:2)),
      "summaries", "all", attempt != 3,
      paste0("a", attempt), 2, 1, attempt
    )
  }
  fragments <- BRREWABC:::flushStoredBuffer(buffer)

  expect_equal(nrow(fragments), 3)
  expect_true(all(file.exists(fragments$file)))
  BRREWABC:::persistStoredGeneration(
    staging, storage, 2, paste0("a", 1:3), "a1",
    "all", "none", fragments
  )
  result <- list(storage = list(path = storage))
  manifest <- list_abc_stored_data(result)

  expect_equal(nrow(manifest), 3)
  expect_true(all(manifest$layout == "fragmented"))
  expect_false(any(file.exists(fragments$file)))
  expect_equal(
    unique(read_summary_statistics(result, status = "retained")$attempt_id),
    "a1"
  )
  expect_equal(
    unique(read_summary_statistics(result, attempt_id = "a3")$attempt_id),
    "a3"
  )
  expect_equal(
    nrow(read_summary_statistics(result, attempt_id = "missing")),
    0
  )

  plan <- consolidate_abc_storage(result, dry_run = TRUE)
  expect_equal(plan$source_files, 3)
  consolidated <- consolidate_abc_storage(result)
  expect_true(file.exists(consolidated$output_file))
  expect_equal(
    nrow(as.data.frame(arrow::read_parquet(consolidated$output_file))),
    sum(manifest$rows)
  )
})

test_that("replacement consolidation updates the active manifest", {
  storage <- file.path(tempdir(), paste0("replace-store-", sample.int(1e8, 1)))
  data_dir <- file.path(storage, "summaries", "value", "generation_0001")
  dir.create(data_dir, recursive = TRUE)
  paths <- file.path(data_dir, paste0("part-", 1:2, ".parquet"))
  for (index in 1:2) {
    arrow::write_parquet(data.frame(
      attempt_id = paste0("a", index), generation = 1L, job_id = 1L,
      accepted = TRUE, retained = TRUE, value = index
    ), paths[[index]])
  }
  schema <- BRREWABC:::schemaSignature(data.frame(value = 1L))
  manifest <- data.frame(
    format_version = 2L, kind = "summaries", name = "value",
    generation = 1L, layout = "fragmented", batch_id = 1:2,
    part = 1L, file = file.path(
      "summaries", "value", "generation_0001", basename(paths)
    ), rows = 1L, bytes = file.info(paths)$size, schema = schema,
    active = TRUE
  )
  arrow::write_parquet(manifest, file.path(storage, "manifest.parquet"))
  result <- list(storage = list(path = storage))

  report <- consolidate_abc_storage(
    result, mode = "replace", keep_fragments = FALSE
  )
  active <- list_abc_stored_data(result)

  expect_true(report$replaced)
  expect_equal(nrow(active), 1)
  expect_equal(active$layout, "consolidated")
  expect_false(any(file.exists(paths)))
  expect_equal(nrow(read_summary_statistics(result)), 2)
})
