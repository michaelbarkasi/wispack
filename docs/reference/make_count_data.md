# Load and parse data, HDF5

Loads MERFISH data from an HDF5 file, or multiple files in a directory,
and combines them into a single data frame. This function was written
for "in house" use by Michael Barkasi, and was not intended for general
use. You can try it, but it may not work for your data and isn't
documented in a user-friendly way.

## Usage

``` r
make_count_data(
  data_path,
  remove_L1 = TRUE,
  ROIname = "Primary auditory area",
  raw = TRUE,
  load_first_file_only = FALSE,
  verbose = TRUE
)
```

## Arguments

- data_path:

  Path to the directory containing the HDF5 files

- remove_L1:

  Logical, whether to exclude cells in layer 1 (default is TRUE)

- ROIname:

  Name of the region of interest (default is "Primary auditory area")

- raw:

  Logical, whether to use raw transcript counts vs normalized ones
  (default is TRUE)

- load_first_file_only:

  Logical, whether to load only the first file found (default is FALSE)

- verbose:

  Logical, whether to print progress messages (default is TRUE)

## Value

A list containing the combined count data and a list of slice plots for
each mouse
