# Use of several distances with ABC-SMC

``` r

library(BRREWABC)
```

## Model definition

``` r

compute_dist <- function(x, ss_obs) {
  ss_sim <- c(x[["alpha"]] + x[["beta"]] + rnorm(1, 0, 0.1),
              x[["alpha"]] * x[["beta"]] + rnorm(1, 0, 0.1))
  dist <- (ss_sim - ss_obs)^2 # the function returns a vector of 2 values
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
              max_number_of_gen = 15,
              nb_acc_prtcl_per_gen = 2000,
              new_threshold_quantile = 0.8,
              nb_threshold = 2, # the number of thresholds/distances used has to be indicated
              experiment_folderpath = "svrldist",
              max_concurrent_jobs = 5,
              verbose = FALSE)
#> The distance threshold(s) (epsilon(s)) fall(s) below the predetermined min value!
#> [1] 0.009842929 0.004205746
```

## Plot results

``` r

all_accepted_particles <- res$particles
all_thresholds <- res$thresholds
plot_abcsmc_res(data = all_accepted_particles, prior = prior_dist,
                filename = "svrldist/res/figs/svrldist_pairplot_all.png", colorpal = "YlGnBu")
#> [1] "Plot saved as 'png'."
plot_densityridges(data = all_accepted_particles, prior = prior_dist,
                   filename = "svrldist/res/figs/svrldist_densityridges.png", colorpal = "YlGnBu")
#> [1] "Plot saved as 'png'."
plot_thresholds(data = all_thresholds, nb_threshold = 2, # the number of thresholds/distances used has to be indicated
                filename = "svrldist/res/figs/svrldist_thresholds.png", colorpal = "YlGnBu")
#> [1] "Plot saved as 'png'."
```

![Pairplot of all
iterations](../reference/figures/svrldist_pairplot_all.png)

Pairplot of all iterations

![Threshold evolution over iterations
(1)](../reference/figures/svrldist_thresholds_dist1.png)

Threshold evolution over iterations (1)

![Threshold evolution over iterations
(2)](../reference/figures/svrldist_thresholds_dist2.png)

Threshold evolution over iterations (2)

![Density estimates for
alpha](../reference/figures/svrldist_densityridges_alpha.png)

Density estimates for alpha

![Density estimates for
beta](../reference/figures/svrldist_densityridges_beta.png)

Density estimates for beta
