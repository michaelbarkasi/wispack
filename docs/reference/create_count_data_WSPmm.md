# Function to convert to WSPmm format

This function converts a data frame of MERFISH data produced by the
`cortical_coordinate_transform` function into a count data format
suitable for input into the WSPmm model. The function bins the data
along a specified dimension (either x or y), collapses along the other
dimension, and includes specified genes and fixed effects in the output.
The resulting count data frame contains counts for each gene in each
cell, along with associated metadata such as bin number, context, mouse
identifier, and cell number.

## Usage

``` r
create_count_data_WSPmm(
  df.merfish,
  bin.dim = c("x_bins", "y_bins"),
  gene.list,
  fixed.effect.names,
  context = NULL
)
```

## Arguments

- df.merfish:

  A data frame of MERFISH data produced by the
  `cortical_coordinate_transform` function

- bin.dim:

  A character string specifying the dimension by which to bin the data;
  must be either "x_bins" or "y_bins"

- gene.list:

  A character vector of genes to include in the count data, which will
  form the "species" level of the model

- fixed.effect.names:

  A character vector of column names in `df.merfish` to include as fixed
  effects in the count data

- context:

  A character string specifying the context level for fixed effects; if
  NULL, will use "context", if not a col name of df.merfish, will rep
  the value

## Value

A data frame in count data format suitable for input into the WSPmm
model

## Examples

``` r
# A tiny toy data frame in the wide, one-column-per-gene shape that
# cortical_coordinate_transform() would normally produce.
toy_merfish <- data.frame(
  mouse      = rep(1, 6),
  cell_num   = 1:6,
  x_bins     = c(1, 1, 2, 2, 3, 3),
  geneA      = c(2, 0, 5, 1, 3, 4),
  geneB      = c(1, 2, 0, 3, 2, 1),
  hemisphere = rep(c("left", "right"), 3)
)
create_count_data_WSPmm(
  toy_merfish,
  bin.dim            = "x_bins",
  gene.list          = c("geneA", "geneB"),
  fixed.effect.names = "hemisphere"
)
#>    count bin context  gene mouse cell_num hemisphere
#> 1      2   1 context geneA     1        1       left
#> 2      0   1 context geneA     1        2      right
#> 3      5   2 context geneA     1        3       left
#> 4      1   2 context geneA     1        4      right
#> 5      3   3 context geneA     1        5       left
#> 6      4   3 context geneA     1        6      right
#> 7      1   1 context geneB     1        1       left
#> 8      2   1 context geneB     1        2      right
#> 9      0   2 context geneB     1        3       left
#> 10     3   2 context geneB     1        4      right
#> 11     2   3 context geneB     1        5       left
#> 12     1   3 context geneB     1        6      right
```
