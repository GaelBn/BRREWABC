# Plot abcsmc results : thresholds over iterations

Plot abcsmc results : thresholds over iterations

## Usage

``` r
plot_thresholds(
  data,
  nb_threshold = 1,
  filename = "thresholds.png",
  figtitle = "",
  colorpal = "Greys"
)
```

## Arguments

- data:

  a dataframe containing the estimation results (thresholds over
  iterations)

- nb_threshold:

  number of thresholds used (see the abcsmc function help for more
  details)

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
