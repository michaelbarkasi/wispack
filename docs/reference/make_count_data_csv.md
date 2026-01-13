# Load and parse data, csv

Loads MERFISH data from a csv file, or multiple files in a directory,
and combines them into a single data frame. This function was written
for "in house" use by Michael Barkasi, and was not intended for general
use. You can try it, but it may not work for your data and isn't
documented in a user-friendly way.

## Usage

``` r
make_count_data_csv(
  data_path,
  remove_L1 = TRUE,
  initialize_zeros = FALSE,
  load_first_file_only = FALSE,
  verbose = TRUE
)
```

## Arguments

- data_path:

  Path to the directory containing the HDF5 files@aliases

- remove_L1:

  Logical, whether to exclude cells in layer 1 (default is TRUE)

- initialize_zeros:

  Logical, whether to initialize new coordinate columns to zero (default
  is FALSE)

- load_first_file_only:

  Logical, whether to load only the first file found (default is FALSE)

- verbose:

  Logical, whether to print progress messages (default is TRUE)

## Value

A list containing the combined count data and a list of slice plots for
each mouse
