# Run ABC-SMC inference in parallel

Run ABC-SMC inference in parallel

## Usage

``` r
abcsmc(
  model_list = list(),
  model_def = NULL,
  prior_dist = list(),
  ss_obs = NA,
  max_number_of_gen = 15,
  nb_acc_prtcl_per_gen = 1000,
  var_prtrbtn_krnl = NA,
  var_prtrbtn_krnl_min = 0.01,
  scaling_sd_prtrbtn_krnl = 1,
  model_jump_prob = 0.1,
  nb_threshold = 1,
  new_threshold_quantile = 0.8,
  distance_threshold_min = 0.01,
  max_attempts = 1e+05,
  acceptance_rate_min = 0.01,
  use_lhs_for_first_iter = TRUE,
  experiment_folderpath = "./",
  on_cluster = FALSE,
  cluster_type = NULL,
  slurm_script_template =
    "#!/bin/bash\n#SBATCH --ntasks=1\n#SBATCH --ntasks-per-node=1\n#SBATCH --hint=nomultithread\n#SBATCH --time=24:00:00\n# THE FOLLOWING SECTION SHOULD NOT BE MODIFIED\n#SBATCH --job-name=job-array_%%A_%%a   # nom du job\n#SBATCH --array=%s-%s%%%d\noutput_fpath=%s\nerror_fpath=%s\n#SBATCH --output=$output_fpath/output_%%A_%%a.out\n#SBATCH --error=$error_fpath/error_%%A_%%a.out\nmkdir -p $output_fpath\nmkdir -p $error_fpath\nRscript %s $SLURM_ARRAY_TASK_ID\n",
  sge_script_template =
    "#!/bin/bash\n#$ -S /bin/bash\n#$ -N subjob_abcsmc_prlll\n# #$ -q \"short.q|long.q\"\n# THE FOLLOWING SECTION SHOULD NOT BE MODIFIED\n#$ -cwd\n#$ -V\n#$ -t %s-%s\n#$ -tc %d\n#$ -o /dev/null\n#$ -e /dev/null\noutput_fpath=%s\nerror_fpath=%s\nmkdir -p $output_fpath\nmkdir -p $error_fpath\nRscript %s $SGE_TASK_ID >$output_fpath/subjob.${SGE_TASK_ID}.out 2>$error_fpath/subjob.${SGE_TASK_ID}.err\n",
  max_concurrent_jobs = 1,
  previous_gens = NA,
  previous_epsilons = NA,
  store_summaries = c("none", "retained", "accepted", "all"),
  store_outputs = c("none", "retained", "accepted", "all"),
  verbose = FALSE,
  progressbar = FALSE
)
```

## Arguments

- model_list:

  a list linking model name ( character string) to associated function

- model_def:

  a R file containing only the model(s) function(s)

- prior_dist:

  a list linking model name (character string) to a list describing the
  prior distribution of each parameter to be estimated

- ss_obs:

  the observed summary statistics

- max_number_of_gen:

  the maximum number of generations to be performed

- nb_acc_prtcl_per_gen:

  the number of particles (per model, the total number corresponding to
  this number multiplied by the number of models) to be accepted during
  a generation before moving on to the next one

- var_prtrbtn_krnl:

  the standard deviation to be used for the perturbation kernel, if NA
  the empirical standard deviation will be used instead

- var_prtrbtn_krnl_min:

  the minimum standard deviation to be used for the perturbation kernel
  (must be greater than zero)

- scaling_sd_prtrbtn_krnl:

  scaling parameter for increasing or decreasing the value of the
  standard deviation used for the perturbation kernel

- model_jump_prob:

  probability of changing model when creating a new particle

- nb_threshold:

  number of thresholds used, corresponding to the number of distances
  returned by the function(s) stipulated in the model list (model_list)

- new_threshold_quantile:

  the quantile to be used to decrease the threshold(s) during iterations

- distance_threshold_min:

  the minimum threshold below which the procedure stops (in order to
  prevent excessively long computations)

- max_attempts:

  the maximum number of particles to be tested during an iteration,
  beyond which the procedure stops (in order to prevent excessively long
  computations)

- acceptance_rate_min:

  the acceptance rate below which the procedure stops (in order to
  prevent excessively long computations)

- use_lhs_for_first_iter:

  whether or not to use a LHS sampling method for the first iteration

- experiment_folderpath:

  the folder in which to carry out the estimation procedure and save the
  results

- on_cluster:

  whether or not the procedure is run on a computation cluster

- cluster_type:

  cluster type used (sge and slurm currently supported)

- slurm_script_template:

  script used to launch jobs on a slurm cluster

- sge_script_template:

  script used to launch jobs on a sge cluster

- max_concurrent_jobs:

  maximum number of jobs/tasks run in parallel

- previous_gens:

  an object (dataframe) containing previous results (set of iterations),
  in order to start from the last iteration performed

- previous_epsilons:

  an object (dataframe) containing previous results (set of thresholds),
  in order to start from the last iteration performed

- store_summaries:

  which summary statistics to store in Parquet files: \`"none"\`,
  \`"retained"\`, \`"accepted"\`, or \`"all"\`

- store_outputs:

  which model outputs to store in Parquet files, using the same policies
  as \`store_summaries\`

- verbose:

  whether or not to display specific information

- progressbar:

  whether or not to display progressbar

## Value

a list containing accepted particles, thresholds, and a \`storage\`
descriptor for summary statistics and model outputs stored in Parquet.

## Examples

``` r
library(BRREWABC)

tmp_dir <- tempdir()

# model definition
compute_dist = function(x, ss_obs){
    ss_sim = c( x[["alpha"]] + x[["beta"]] + rnorm(1,0,0.1),
           x[["alpha"]] * x[["beta"]] + rnorm(1,0,0.1) ) # a simple toy model
    dist = sum((ss_sim-ss_obs)^2)
    return(c(dist))
}

MODEL_LIST <- list("m1" = compute_dist)
PRIOR_DIST <- list("m1" = list(c('alpha', 'unif', 0, 4),
                               c('beta', 'unif', 0, 1)))

# create a reference trajectory
sum_stat_obs = c(2.0,0.75)

# run abc smc procedure
res = abcsmc(model_list = MODEL_LIST, prior_dist = PRIOR_DIST,
ss_obs = sum_stat_obs, max_number_of_gen = 20, nb_acc_prtcl_per_gen = 2000,
new_threshold_quantile = 0.8, experiment_folderpath = tmp_dir,
max_concurrent_jobs = 2, verbose = FALSE)
#> Warning: EOF within quoted string
#> Warning: EOF within quoted string
#> Warning: EOF within quoted string
#> Warning: EOF within quoted string

# get results and plots
all_accepted_particles = res$particles
all_thresholds = res$thresholds
plot_abcsmc_res(data = all_accepted_particles, prior = PRIOR_DIST, colorpal = "YlOrBr", filename = file.path(tmp_dir, "abcsmc_results.png"))
#> [1] "Number of generations exceed the threshold (15) allowed by ggpairs, it may cause long processing times. You may (re)define the iter argument to choose which generations to plot."
#> [1] "Plot saved as 'png'."
plot_thresholds(data = all_thresholds, nb_threshold = 1, colorpal = "YlOrBr", filename = file.path(tmp_dir, "thresholds.png"))
#> [1] "Plot saved as 'png'."
plot_ess(data = all_accepted_particles, colorpal = "YlOrBr", filename = file.path(tmp_dir, "ess.png"))
#> [1] "Plot saved as 'png'."
#>    gen      ess
#> 1    1 2000.000
#> 2    2 1875.993
#> 3    3 1914.225
#> 4    4 1915.111
#> 5    5 1934.219
#> 6    6 1930.498
#> 7    7 1928.733
#> 8    8 1920.193
#> 9    9 1922.517
#> 10  10 1914.963
#> 11  11 1899.364
#> 12  12 1872.430
#> 13  13 1885.557
#> 14  14 1835.226
#> 15  15 1867.778
#> 16  16 1820.865
#> 17  17 1816.429
#> 18  18 1801.023
#> 19  19 1763.685
#> 20  20 1752.913
plot_densityridges(data = all_accepted_particles, prior = PRIOR_DIST, colorpal = "YlOrBr", filename = file.path(tmp_dir, "densityridges.png"))
#> [1] "Plot saved as 'png'."
```
