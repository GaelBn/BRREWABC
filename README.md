
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
#> Folder already exists.
#> Check folder_path for : res/csv
#> Folder already exists.
#> Check folder_path for : res/figs
#> Folder already exists.
#> gen 1 
#> threshold: 
#> prtrbtn_krnl_sd: NA NA 
#> -
#> gen 2 
#> threshold: 3.998623 
#> prtrbtn_krnl_sd: 1.155235 0.2888307 
#> -
#> gen 3 
#> threshold: 2.153006 
#> prtrbtn_krnl_sd: 0.9469025 0.2584123 
#> -
#> gen 4 
#> threshold: 1.397381 
#> prtrbtn_krnl_sd: 0.7849039 0.2657151 
#> -
#> gen 5 
#> threshold: 0.8519803 
#> prtrbtn_krnl_sd: 0.646988 0.2736002 
#> -
#> gen 6 
#> threshold: 0.5684694 
#> prtrbtn_krnl_sd: 0.5465287 0.2701191 
#> -
#> gen 7 
#> threshold: 0.4040715 
#> prtrbtn_krnl_sd: 0.4819998 0.2591692 
#> -
#> gen 8 
#> threshold: 0.2762117 
#> prtrbtn_krnl_sd: 0.420591 0.2436711 
#> -
#> gen 9 
#> threshold: 0.1980603 
#> prtrbtn_krnl_sd: 0.3758658 0.2273871 
#> -
#> gen 10 
#> threshold: 0.1379844 
#> prtrbtn_krnl_sd: 0.3317047 0.2123022 
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

<div class="figure">

<img src="man/figures/pairplot_m1.png" alt="Pairplot of all iterations" width="100%" />
<p class="caption">
Pairplot of all iterations
</p>

</div>

<div class="figure">

<img src="man/figures/thresholds_dist1.png" alt="Threshold evolution over iterations" width="100%" />
<p class="caption">
Threshold evolution over iterations
</p>

</div>

<div class="figure">

<img src="man/figures/densityridges_m1_alpha.png" alt="Density estimates for alpha" width="100%" />
<p class="caption">
Density estimates for alpha
</p>

</div>

<div class="figure">

<img src="man/figures/densityridges_m1_beta.png" alt="Density estimates for beta" width="100%" />
<p class="caption">
Density estimates for beta
</p>

</div>
