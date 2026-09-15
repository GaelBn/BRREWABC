# Read model outputs stored in Parquet files

Read model outputs stored in Parquet files

## Usage

``` r
read_model_outputs(
  x,
  name = NULL,
  generation = NULL,
  attempt_id = NULL,
  status = c("all", "accepted", "retained", "rejected")
)
```

## Arguments

- x:

  an object returned by \[abcsmc()\] or \[abcrejection()\], or the path
  to an experiment directory or its \`res/parquet\` directory.

- name:

  optional names of summary statistics to read.

- generation:

  optional generation numbers.

- attempt_id:

  optional attempt identifiers.

- status:

  one of \`"all"\`, \`"accepted"\`, \`"retained"\`, or \`"rejected"\`.

## Value

a data frame containing the selected model outputs.
