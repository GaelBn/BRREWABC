# Plot abcrejection results : pairplot

Plot abcrejection results : pairplot

## Usage

``` r
plot_abcrejection_res(
  data,
  prior,
  thresholds = NA,
  filename = "pairplot.png",
  figtitle = "ABC rejection results",
  colorpal = "Greys"
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

- thresholds:

  if not NA, the threshold(s) (a value or a vector of thresholds) to use
  to select acceptable particles. Should be use only if abcrejection was
  run using thresholds = NA

- filename:

  the file name to be used to save the plots (the extension defines the
  format: pdf or png)

- figtitle:

  the figure title

- colorpal:

  a palette name as used in the RColorBrewer package

## Value

one or several plots in pdf or png format

## Examples

``` r
# see the abcrejection function help for details on how to plot results
```
