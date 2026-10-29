# Installation instructions

1.  Clone and `cd` into this git repo
2.  Make `build_install.sh` executable by running
    `chmod +x build_install.sh` in bash terminal\
3.  In terminal, run `./build_install.sh` to build and install the
    package\
4.  To ensure a clean start, run:

``` bash
rm -f src/*.o src/*.so
rm -rf wispack.Rcheck
rm -f wispack_*.tar.gz
./build_install.sh
```

Note that additional system dependencies are required *before* the
building process. Below are instructions for Mac and Linux
(Ubuntu/Debian). Installation on Windows should be possible, but is not
currently tested or supported.

## Mac

First, ensure your Mac has the necessary backend dependencies:

1.  Xcode developer tools: In terminal, run `xcode-select --install` and
    follow the instructions
2.  GNU Fortran compiler: Go to <https://mac.r-project.org/tools/>, then
    download and run the relevant gfortran-\[version\]-universal.pkg
3.  Follow the instructions at this same website to install LaTeX for R

Ensure pkg-config and nlopt are installed, e.g., using Homebrew in a
terminal:

``` bash
brew install pkg-config
brew install nlopt
export PKG_CONFIG_PATH="/opt/homebrew/lib/pkgconfig:$PKG_CONFIG_PATH"
```

Install Homebrew if needed:

``` bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
echo >> /Users/[name]/.zprofile
echo 'eval "$(/opt/homebrew/bin/brew shellenv zsh)"' >> /Users/[name]/.zprofile
eval "$(/opt/homebrew/bin/brew shellenv zsh)"
```

Rebuild and cleanly install Rcpp and Stan, in R:
`install.packages(c("Rcpp", "RcppEigen", "RcppParallel", "StanHeaders"))`

## Linux

Need to install the following system dependencies: - `libxml2` for
roxygen2 - `r-base-dev`, `libnlopt-dev` (or `libnlopt-cxx-dev`) for
Rcpp, RcppEigen, and StanHeaders

For PDF documentation generation, `pdflatex` should also be installed

``` bash
sudo apt install \
    libxml2 \
    r-base-dev \
    libnlopt-dev \
    libnlopt-cxx-dev \
    texlive-latex-base \
    texlive-fonts-recommended \
    texlive-fonts-extra \
    texlive-latex-extra
```

## Build documentation

In R, run:

``` r

devtools::document()
pkgdown::build_site()
pkgdown::build_article("<article name>")
```

## Running the test suite

Wispack has a `testthat` suite under `tests/testthat`, split into two
tiers:

- **Fast tier** (runs by default, in well under a minute): synthetic
  data and precomputed fixtures only – no real model fitting against the
  data used in the tutorials. This is what runs automatically with
  `R CMD check` or:

``` r

devtools::test()
```

- **Slow tier** (opt-in, several minutes): full-scale replications of
  the tutorial vignettes’ own fitted-value checks, using the same
  extdata and the same bootstrap/MCMC settings as the tutorials
  themselves. These are gated behind an environment variable so they
  don’t slow down routine testing:

``` r

Sys.setenv(WISPACK_RUN_SLOW_TESTS = "true")
devtools::test()
```

(equivalently,
`WISPACK_RUN_SLOW_TESTS=true Rscript -e 'devtools::test()'` from the
shell).

Both tiers read their expected numeric bounds from a single shared table
(`.wisp_golden`, in `R/golden_checks.R`), which the tutorial vignettes’
own check chunks also read from – so a model-fitting change that shifts
a tutorial’s narrative numbers will show up the same way whether you
notice it by rebuilding the vignette or by running the test suite. If
you deliberately change model-fitting behavior, re-run the slow tier and
update the bounds in `R/golden_checks.R` (and the corresponding prose in
whichever tutorial it affects) in the same commit.

## System requirements

Wispack is computationally intense. Fitting all but the smallest models
can easily use four or five Gb of RAM. We strongly recommend a modern
multi-core Unix system (Apple Silicon or Linux) with at least 24 Gb of
RAM, although realistically a workstation or compute cluster with 100+
Gb of available RAM and 20 or more cores are needed for realistic run
times of complex models.
