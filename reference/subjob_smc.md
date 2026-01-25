# Run a subtask of the ABC-SMC. Shouldn't have to be used by the user, this function is visible so that it can be used on cluster by the main script

Run a subtask of the ABC-SMC. Shouldn't have to be used by the user,
this function is visible so that it can be used on cluster by the main
script

## Usage

``` r
subjob_smc(job_id, path_to_abc_state)
```

## Arguments

- job_id:

  id of the current job

- path_to_abc_state:

  path to the .Rdata of the current experiment

## Value

nothing, write results in a tmp file
