# Store model summaries and outputs in Parquet

Model functions can return summary statistics and detailed outputs
alongside the distances used by the ABC algorithm. The historical
numeric distance vector remains supported.

``` r

epidemic_model <- function(x, ss_obs) {
  trajectory <- data.frame(
    timestep = rep(1:100, each = 2),
    pop_id = rep(c("A", "B"), 100),
    S = simulate_S(x),
    E = simulate_E(x),
    I = simulate_I(x),
    R = simulate_R(x)
  )

  list(
    distances = c(dist1 = compute_distance(trajectory, ss_obs)),
    summaries = list(
      epidemic_trajectory = trajectory,
      peak_infectious = max(trajectory$I)
    ),
    outputs = list(
      final_state = subset(trajectory, timestep == max(timestep))
    )
  )
}
```

Each element of `summaries` and `outputs` must have a unique name and
contain an atomic vector, a matrix, or a data frame. A tabular object’s
column names and types must remain fixed across simulations. The number
of rows may vary.

Use separate policies to select which objects are kept:

``` r

result <- abcsmc(
  model_list = list(m1 = epidemic_model),
  prior_dist = prior_dist,
  ss_obs = ss_obs,
  store_summaries = "all",
  store_outputs = "retained",
  storage_chunk_rows = 1000000,
  storage_chunk_mb = 128
)
```

Available policies are `"none"`, `"retained"`, `"accepted"`, and
`"all"`. Data are stored below `res/parquet`, separately for each object
and generation. Workers buffer records and publish bounded Parquet
fragments atomically. The row and approximate memory limits are
independent of `batch_size`, which only controls the number of
simulations assigned to a worker process.

The manifest lists the available datasets:

``` r

list_abc_stored_data(result)
```

Readers can restrict the data by object name, generation, attempt
identifier, or ABC status:

``` r

accepted_trajectory <- read_summary_statistics(
  result,
  name = "epidemic_trajectory",
  generation = 5,
  status = "accepted"
)

rejected_trajectory <- read_summary_statistics(
  result,
  name = "epidemic_trajectory",
  status = "rejected"
)

selected_output <- read_model_outputs(
  result,
  name = "final_state",
  attempt_id = c("g0005-b00000001-a0003"),
  status = "all"
)
```

The returned tables include `attempt_id`, `generation`, `job_id`,
`accepted`, and `retained`. The `stored_name` column identifies the
requested summary or output when several objects are read together.

## Consolidating fragments on demand

The fragmented layout is intended for robust computation and selective
reads. For transfer or archival,
[`consolidate_abc_storage()`](https://gaelbn.github.io/BRREWABC/reference/consolidate_abc_storage.md)
creates one Parquet file per named object and generation:

``` r

plan <- consolidate_abc_storage(
  result,
  kind = c("summaries", "outputs"),
  generation = 1:5,
  dry_run = TRUE
)

exported <- consolidate_abc_storage(
  result,
  generation = 1:5
)
```

The default `mode = "export"` writes below `res/parquet/consolidated`
and does not alter the active storage manifest. To make compacted files
active:

``` r

consolidate_abc_storage(
  result,
  kind = "outputs",
  mode = "replace",
  keep_fragments = TRUE
)
```

Replacement is transactional: each compacted file is finalized before
the manifest is updated. Keeping source fragments provides an additional
recovery path; set `keep_fragments = FALSE` only when disk space is more
important.
