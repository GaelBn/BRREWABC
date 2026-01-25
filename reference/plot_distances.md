# Plot distribution of distances by model

Plot distribution of distances by model

## Usage

``` r
plot_distances(
  data,
  generation = NULL,
  filename = "distances.png",
  figtitle = "",
  colorpal = "Set2",
  bins = 20,
  alpha = 0.6,
  adjust = 1.5,
  width = 8,
  height = 5
)
```

## Arguments

- data:

  a dataframe containing the estimation results (set of particles
  accepted during iterations, including : gen, model, pWeight and one or
  more dist columns)

- generation:

  generation to plot; defaults to the maximum value in gen

- filename:

  the file name to save the plot (extension png or pdf determines
  format)

- figtitle:

  the figure title

- colorpal:

  a palette name as used in the RColorBrewer package

- bins:

  number of histogram bins

- alpha:

  histogram transparency level

- adjust:

  density curve smoothing adjustment

- width:

  plot width in inches (only applies to png/pdf output)

- height:

  plot height in inches (only applies to png/pdf output)

## Value

a dataframe in long format of distances; the plot is saved to file

## Examples

``` r
# plot_distances(data)
```
