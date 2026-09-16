# Resume an estimation using ABC-SMC

``` r

library(BRREWABC)
```

## Model definition

``` r

compute_dist <- function(x, ss_obs) {
  ss_sim <- c(x[["alpha"]] + x[["beta"]] + rnorm(1, 0, 0.1),
              x[["alpha"]] * x[["beta"]] + rnorm(1, 0, 0.1))
  dist <- sum((ss_sim - ss_obs)^2)
  return(c(dist))
}

model_list <- list("m1" = compute_dist)
```

## Define prior distribution

``` r

prior_dist <- list("m1" = list(c("alpha", "unif", 0, 4),
                               c("beta", "unif", 0, 1)))
```

## Create a reference trajectory

``` r

sum_stat_obs <- c(2.0, 0.75)
```

## Run abc smc procedure

``` r

res <- abcsmc(model_list = model_list,
              prior_dist = prior_dist,
              ss_obs = sum_stat_obs,
              max_number_of_gen = 10,
              nb_acc_prtcl_per_gen = 2000,
              new_threshold_quantile = 0.8,
              experiment_folderpath = "rsmsmpl",
              max_concurrent_jobs = 5,
              verbose = FALSE)
```

## Plot results

``` r

all_accepted_particles <- res$particles
all_thresholds <- res$thresholds
plot_abcsmc_res(data = all_accepted_particles, prior = prior_dist,
                filename = "rsmsmpl/res/figs/rsmsmpl_pairplot_all.png", colorpal = "Greys")
#> [1] "Plot saved as 'png'."
plot_densityridges(data = all_accepted_particles, prior = prior_dist,
                   filename = "rsmsmpl/res/figs/rsmsmpl_densityridges.png", colorpal = "Greys")
#> [1] "Plot saved as 'png'."
plot_thresholds(data = all_thresholds, nb_threshold = 1,
                filename = "rsmsmpl/res/figs/rsmsmpl_thresholds.png", colorpal = "Greys")
#> [1] "Plot saved as 'png'."
```

![Pairplot of all
iterations](../reference/figures/rsmsmpl_pairplot_all.png)

Pairplot of all iterations

![Threshold evolution over
iterations](../reference/figures/rsmsmpl_thresholds.png)

Threshold evolution over iterations

![Density estimates for
alpha](../reference/figures/rsmsmpl_densityridges_alpha.png)

Density estimates for alpha

![Density estimates for
beta](../reference/figures/rsmsmpl_densityridges_beta.png)

Density estimates for beta

## Re-run abc smc procedure from last iteration of previous results

Be sure to update any parameters defining stopping conditions that may
have caused the previous procedure to terminate (here
`max_number_of_gen`).

``` r

res <- abcsmc(model_list = model_list,
              prior_dist = prior_dist,
              ss_obs = sum_stat_obs,
              max_number_of_gen = 20,
              nb_acc_prtcl_per_gen = 2000,
              new_threshold_quantile = 0.8,
              experiment_folderpath = "rsmsmpl",
              max_concurrent_jobs = 5,
              previous_gens = all_accepted_particles,
              previous_epsilons = all_thresholds,
              verbose = FALSE)
#> The distance threshold(s) (epsilon(s)) fall(s) below the predetermined min value!
#> [1] 0.00950524
```

## Plot new results

``` r

all_accepted_particles <- res$particles
all_thresholds <- res$thresholds
plot_abcsmc_res(data = all_accepted_particles, prior = prior_dist,
                filename = "rsmsmpl/res/figs/rsmsmpl_pairplot_all_rsm.png", colorpal = "OrRd")
#> [1] "Number of generations exceed the threshold (15) allowed by ggpairs, it may cause long processing times. You may (re)define the iter argument to choose which generations to plot."
#> [1] "Plot saved as 'png'."
plot_densityridges(data = all_accepted_particles, prior = prior_dist,
                   filename = "rsmsmpl/res/figs/rsmsmpl_densityridges_rsm.png",
                   colorpal = "OrRd")
#> [1] "Plot saved as 'png'."
plot_thresholds(data = all_thresholds, nb_threshold = 1,
                filename = "rsmsmpl/res/figs/rsmsmpl_thresholds_rsm.png", colorpal = "OrRd")
#> [1] "Plot saved as 'png'."
```

![Pairplot of all
iterations](../reference/figures/rsmsmpl_pairplot_all_rsm.png)

Pairplot of all iterations

![Threshold evolution over
iterations](../reference/figures/rsmsmpl_thresholds_rsm.png)

Threshold evolution over iterations

![Density estimates for
alpha](../reference/figures/rsmsmpl_densityridges_rsm_alpha.png)

Density estimates for alpha

![Density estimates for
beta](../reference/figures/rsmsmpl_densityridges_rsm_beta.png)

Density estimates for beta
