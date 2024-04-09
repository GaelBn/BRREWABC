
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BRREWABC <img src="man/figures/icon.png" align="right" width="90" />

<!-- badges: start -->
<!-- badges: end -->

BRREWABC (Batched Resilient and Rapid Estimation Workflow through
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
devtools::install_github("GaelBn/BRREWABC")
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
  ss_obs = sum_stat_obs, max_number_of_gen = 20, nb_acc_prtcl_per_gen = 3000,
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
#> threshold: 4.028069 
#> prtrbtn_krnl_sd: 1.154902 0.2887234 
#> -
#> gen 3 
#> threshold: 2.144677 
#> prtrbtn_krnl_sd: 0.9314262 0.2606515 
#> -
#> gen 4 
#> threshold: 1.318785 
#> prtrbtn_krnl_sd: 0.7742273 0.2621097 
#> -
#> gen 5 
#> threshold: 0.8125375 
#> prtrbtn_krnl_sd: 0.6405105 0.2649226 
#> -
#> gen 6 
#> threshold: 0.5382799 
#> prtrbtn_krnl_sd: 0.5345319 0.2671884 
#> -
#> gen 7 
#> threshold: 0.3673938 
#> prtrbtn_krnl_sd: 0.4715125 0.2567703 
#> -
#> gen 8 
#> threshold: 0.2547235 
#> prtrbtn_krnl_sd: 0.4156183 0.2435061 
#> -
#> gen 9 
#> threshold: 0.1790777 
#> prtrbtn_krnl_sd: 0.3651152 0.2244952 
#> -
#> gen 10 
#> threshold: 0.1268142 
#> prtrbtn_krnl_sd: 0.3356348 0.2150104 
#> -
#> gen 11 
#> threshold: 0.09292255 
#> prtrbtn_krnl_sd: 0.3091853 0.2053241 
#> -
#> gen 12 
#> threshold: 0.06938362 
#> prtrbtn_krnl_sd: 0.2935447 0.1979522 
#> -
#> gen 13 
#> threshold: 0.05311042 
#> prtrbtn_krnl_sd: 0.2699788 0.1876449 
#> -
#> gen 14 
#> threshold: 0.04051743 
#> prtrbtn_krnl_sd: 0.2553968 0.177729 
#> -
#> gen 15 
#> threshold: 0.0315603 
#> prtrbtn_krnl_sd: 0.2461897 0.1712717 
#> -
#> gen 16 
#> threshold: 0.02477117 
#> prtrbtn_krnl_sd: 0.2389223 0.1683386 
#> -
#> gen 17 
#> threshold: 0.01964872 
#> prtrbtn_krnl_sd: 0.2346167 0.1649224 
#> -
#> gen 18 
#> threshold: 0.01588233 
#> prtrbtn_krnl_sd: 0.2311844 0.1612676 
#> -
#> gen 19 
#> threshold: 0.01261688 
#> prtrbtn_krnl_sd: 0.2175436 0.1535265 
#> -
#> gen 20 
#> threshold: 0.00999965 
#> prtrbtn_krnl_sd: 0.2104832 0.1460712
#> The distance threshold(s) (epsilon(s)) fall(s) below the predetermined min value!
#> [1] 0.00999965
#> Experiment done!

# get results and plots
all_accepted_particles = res$particles
all_thresholds = res$thresholds
plot_abcsmc_res(data = all_accepted_particles, prior = PRIOR_DIST, colorpal = "Blues")
#> [1] "Plot saved as '.png'."
#> Registered S3 method overwritten by 'GGally':
#>   method from   
#>   +.gg   ggplot2
plot_densityridges(data = all_accepted_particles, prior = PRIOR_DIST, colorpal = "Blues")
#> [1] "Plot saved as '.png'."
plot_thresholds(data = all_thresholds, nb_threshold = 1, colorpal = "Blues")
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
