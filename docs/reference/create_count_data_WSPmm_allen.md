# Function to convert Allen-formatted MERFISH data to WSPmm format

This is the
[`create_count_data_WSPmm`](https://michaelbarkasi.github.io/wispack/reference/create_count_data_WSPmm.md)
counterpart for MERFISH data already in Allen Brain Atlas format, i.e.
already in long form with one row per gene/cell pair (columns
`trscrpt_gene_symb`, `trscrpt_ct`, `cluster_id`, `mouse`), rather than
the wide, one-column-per-gene form that `create_count_data_WSPmm`
expects. Because the data is already long, no `gene.list` argument is
needed here: genes are read directly from `trscrpt_gene_symb`. An `age`
fixed effect, if requested via `fixed.effect.names`, is always set to
`"p5x"` for every row, since this function was written for a single age
cohort; include `"age"` in `fixed.effect.names` only when that is the
cohort you are working with.

## Usage

``` r
create_count_data_WSPmm_allen(
  df.merfish,
  bin.dim = c("x_bins", "y_bins"),
  fixed.effect.names,
  context = NULL
)
```

## Arguments

- df.merfish:

  A data frame of Allen-formatted MERFISH data, in long form with one
  row per gene/cell pair

- bin.dim:

  A character string specifying the dimension by which to bin the data;
  must be either "x_bins" or "y_bins"

- fixed.effect.names:

  A character vector of column names in `df.merfish` to include as fixed
  effects in the count data

- context:

  A character string specifying the context level for fixed effects; if
  NULL, will use "context" literally for every row; otherwise must be
  the name of an existing column in `df.merfish`, whose values are used
  as the context level

## Value

A data frame in count data format suitable for input into the WSPmm
model

## Examples

``` r
# A tiny toy data frame in the long, one-row-per-gene/cell-pair Allen format.
toy_allen <- data.frame(
  mouse             = rep(1, 4),
  cluster_id        = 1:4,
  y_bins            = c(1, 1, 2, 2),
  trscrpt_gene_symb = c("geneA", "geneB", "geneA", "geneB"),
  trscrpt_ct        = c(3, 1, 4, 2),
  hemisphere        = c("left", "right", "left", "right")
)
create_count_data_WSPmm_allen(
  toy_allen,
  bin.dim            = "y_bins",
  fixed.effect.names = "hemisphere"
)
#>   count bin context  gene mouse cell_num hemisphere
#> 1     3   1 context geneA     1        1       left
#> 2     1   1 context geneB     1        2      right
#> 3     4   2 context geneA     1        3       left
#> 4     2   2 context geneB     1        4      right
```
