# Plot abcsmc results : pairplot and model selection

Plot abcsmc results : pairplot and model selection

## Usage

``` r
plot_abcsmc_res(
  data,
  prior,
  filename = "pairplot.png",
  figtitle = "ABC smc results",
  colorpal = "Greys",
  iter = NA
)
```

## Arguments

- data:

  a dataframe containing the estimation results (set of particles
  accepted during iterations)

- prior:

  a list linking model name (character string) to a list describing the
  prior distribution of each parameter estimated (the same as used for
  the abcsmc function)

- filename:

  the file name to be used to save the plots (the extension defines the
  format: pdf or png)

- figtitle:

  the figure title

- colorpal:

  a palette name as used in the RColorBrewer package

- iter:

  if not NA, the id (number) of the iteration to be plotted

## Value

one or several plots in pdf or png format

## Examples

``` r
# see the abcsmc function help for details on how to plot results
```
