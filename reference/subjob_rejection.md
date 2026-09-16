# Run a batch task of ABC rejection

This compatibility wrapper runs a batch specification stored in
\`batch_manifest_path\` inside the saved ABC state.

## Usage

``` r
subjob_rejection(job_id, path_to_abc_state)
```

## Arguments

- job_id:

  array-task index in the batch manifest.

- path_to_abc_state:

  path to the saved ABC state.

## Value

the batch control object, invisibly.
