# Function to transform coordinates for each mouse and extract layer boundary estimates

This function transforms the raw x,y coordinates of cells into laminar
(y) and columnar (x) coordinates, bins the data, and estimates layer
boundaries for each mouse in the dataset. It generates plots to
visualize the transformation process and can optionally save these
plots. This function was written for "in house" use by Michael Barkasi,
and was not intended for general use. You can try it, but it may not
work for your data and isn't documented in a user-friendly way.

## Usage

``` r
cortical_coordinate_transform(
  count_data,
  total_bins,
  keep_plots = FALSE,
  L1_removed = TRUE,
  nat_left = FALSE,
  nat_right = FALSE,
  leveling_layer = "L4",
  verbose = TRUE
)
```

## Arguments

- count_data:

  A list containing the count data and slice plots for each mouse, as
  produced by the `make_count_data` or `make_count_data_csv` functions

- total_bins:

  Number of bins into which to divide the spatial axes, e.g., 100

- keep_plots:

  Logical, whether to keep the generated plots in the output (default is
  FALSE)

- L1_removed:

  Logical, whether the data in count_data has L1 removed (default is
  TRUE)

- nat_left:

  Logical, setting used based on manual visual checks to ensure
  consistent orientation, causes angle and coordinate sign flips in
  layer leveling (default is FALSE)

- nat_right:

  Logical, setting used based on manual visual checks to ensure
  consistent orientation, causes angle and coordinate sign flips in
  layer leveling (default is FALSE)

- leveling_layer:

  Character string, the layer to use when leveling tissue patch in
  coordinate transformation (default is)

- verbose:

  Logical, whether to print progress messages (default is TRUE)

## Value

A list containing the transformed count data, layer boundary estimates,
and optionally the generated plots
