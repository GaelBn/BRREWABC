# Consolidate stored ABC summaries and outputs

Consolidates the active Parquet fragments for each selected named object
and generation into one Parquet file. By default this creates a separate
export and leaves the operational fragmented storage unchanged.

## Usage

``` r
consolidate_abc_storage(
  x,
  kind = c("summaries", "outputs"),
  name = NULL,
  generation = NULL,
  mode = c("export", "replace"),
  output_dir = NULL,
  keep_fragments = TRUE,
  compression = "zstd",
  overwrite = FALSE,
  dry_run = FALSE
)
```

## Arguments

- x:

  an object returned by \[abcsmc()\] or \[abcrejection()\], or an
  experiment/storage directory.

- kind:

  one or both of \`"summaries"\` and \`"outputs"\`.

- name:

  optional stored object names.

- generation:

  optional generation numbers.

- mode:

  \`"export"\` to leave the active manifest unchanged, or \`"replace"\`
  to make consolidated files the active storage layout.

- output_dir:

  destination directory. Defaults to \`consolidated\` for an export and
  \`compacted\` for replacement, below the storage root.

- keep_fragments:

  whether source fragments should remain on disk after a successful
  replacement.

- compression:

  Parquet compression codec.

- overwrite:

  whether existing consolidated files may be replaced.

- dry_run:

  return the consolidation plan without writing files.

## Value

a data frame describing the consolidation performed or planned.
