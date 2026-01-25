# Plot abcsmc results : densityridges for estimated parameters

Plot abcsmc results : densityridges for estimated parameters

## Usage

``` r
plot_densityridges(
  data,
  prior,
  filename = "densityridges.png",
  figtitle = "",
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
# see the abcsmc function help for details on how to plot results
```
