
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BRREWABC <img src="man/figures/icon.png" align="right" width="90" />

<!-- badges: start -->
<!-- badges: end -->

BRREW-ABC (Batched Resilient and Rapid Estimation Workflow through
Approximate Bayesian Computation) : an R package designed to facilitate
inference through a parallelized Approximate Bayesian Computation
Sequential Monte Carlo (ABC SMC) algorithm. This package streamlines the
process of conducting Bayesian inference for complex models by
implementing efficient parallelization techniques.

## Installation

You can install the development version of BRREWABC from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("GaelBn/BRREW-ABC")
```

## Usage

This is a basic example which shows you how to solve a common problem:

``` r
library(BRREWABC)

# model definition
compute_dist = function(x, ss_obs){
    ss_sim = c( x[["alpha"]] + x[["beta"]] + rnorm(1,0,0.1),
           x[["alpha"]] * x[["beta"]] + rnorm(1,0,0.1) ) # a very simple toy model
    dist = sum((ss_sim-ss_obs)^2)
    return(c(dist))
}

MODEL_LIST <- list("m1" = compute_dist)
PRIOR_DIST <- list("m1" = list(c('alpha', 'unif', 0, 4), c('beta', 'unif', 0, 1)))

# create a reference trajectory
sum_stat_obs = c(2.0,0.75)

# run abc smc procedure
res = abcsmc(model_list = MODEL_LIST, prior_dist = PRIOR_DIST,
  ss_obs = sum_stat_obs, max_number_of_gen = 10, nb_acc_prtcl_per_gen = 1000,
  new_threshold_quantile = 0.8, experiment_folderpath = "",
  max_concurrent_jobs = 2, verbose = TRUE)
#> [1] "tmp/currentABCState.RData"
#> Check folder_path for : tmp
#> Folder created successfully.
#> Check folder_path for : res
#> Folder created successfully.
#> Check folder_path for : res/csv
#> Folder created successfully.
#> Check folder_path for : res/figs
#> Folder created successfully.
#> gen 1 
#> threshold: 
#> prtrbtn_krnl_sd: NA NA 
#> -
#> gen 2 
#> threshold: 3.8241 
#> prtrbtn_krnl_sd: 1.155307 0.2888311 
#> -
#> gen 3 
#> threshold: 2.124396 
#> prtrbtn_krnl_sd: 0.9411677 0.262425 
#> -
#> gen 4 
#> threshold: 1.304413 
#> prtrbtn_krnl_sd: 0.7790466 0.2606397 
#> -
#> gen 5 
#> threshold: 0.7881324 
#> prtrbtn_krnl_sd: 0.6282967 0.2600914 
#> -
#> gen 6 
#> threshold: 0.5236896 
#> prtrbtn_krnl_sd: 0.5290049 0.2651203 
#> -
#> gen 7 
#> threshold: 0.3709691 
#> prtrbtn_krnl_sd: 0.4812579 0.2589308 
#> -
#> gen 8 
#> threshold: 0.2567042 
#> prtrbtn_krnl_sd: 0.42443 0.2441114 
#> -
#> gen 9 
#> threshold: 0.1831002 
#> prtrbtn_krnl_sd: 0.3795465 0.2302809 
#> -
#> gen 10 
#> threshold: 0.132646 
#> prtrbtn_krnl_sd: 0.3398539 0.2133927 
#> -
#> Experiment done!

# get results and plots
all_accepted_particles = res$particles
all_thresholds = res$thresholds
plot_abcsmc_res(data = all_accepted_particles, prior = PRIOR_DIST, colorpal = "GnBu")
#> [1] "Plot saved as '.png'."
#> Registered S3 method overwritten by 'GGally':
#>   method from   
#>   +.gg   ggplot2
plot_densityridges(data = all_accepted_particles, prior = PRIOR_DIST, colorpal = "GnBu")
#> [1] "Plot saved as '.png'."
plot_thresholds(data = all_thresholds, nb_threshold = 1, colorpal = "GnBu")
#> [1] "Plot saved as 'png'."
```

<!-- ```{r pairplot, echo = FALSE}
plot_abcsmc_res(data = all_accepted_particles, prior = PRIOR_DIST, colorpal = "GnBu")
``` -->
