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
create.count.data.WSPmm(
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
