# Read summary statistics stored in Parquet files

Read summary statistics stored in Parquet files

## Usage

``` r
read_summary_statistics(
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

a data frame containing the selected summary statistics.
