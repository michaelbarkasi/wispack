# Build-stamp helper for the vignettes' YAML "date" field #############################################
#
# Several vignettes/articles/*.Rmd files carry a "date: ... (vX.Y)" header recording when the document
# was last confirmed to run against a given package version. Hand-typing that string on every rebuild is
# easy to forget (it's metadata, not content), so this makes it self-updating -- but only when the
# knitr cache for that document was empty at the start of the knit (i.e. every chunk is about to run
# from scratch), not on an incremental re-knit that reuses cached results. An incremental re-knit keeps
# reporting the date/version of the last from-scratch build, since nothing was actually re-verified.
#
# The "last from-scratch build" stamp is stored as a plain text file inside the document's own knitr
# cache directory, so clearing that cache (the usual way to force a real rebuild) also discards the old
# stamp, keeping the two in lockstep without any separate bookkeeping file to forget to delete.

#' Compute (or recall) the "date (vX.Y)" stamp used in a vignette's YAML header.
#'
#' @param doc_name Base name of the .Rmd file, no extension (its knitr cache directory is
#'   "<doc_name>_cache", by knitr's own naming convention).
#' @noRd
.wisp_build_stamp <- function(doc_name) {

  cache_dir  <- paste0(doc_name, "_cache")
  stamp_file <- file.path(cache_dir, "build_stamp.txt")

  cache_was_empty <- !dir.exists(cache_dir) || length(list.files(cache_dir, recursive = TRUE)) == 0

  if (cache_was_empty || !file.exists(stamp_file)) {
    # Either this is a genuine from-empty-cache build, or the cache has leftover content from before
    # this stamping mechanism existed (no recorded stamp to recall) -- either way, nothing trustworthy
    # to recall, so (re)compute and record a fresh stamp.
    stamp <- sprintf(
      "%s (v%s)",
      format(Sys.Date(), "%B %d, %Y"),
      as.character(utils::packageVersion("wispack"))
    )
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
    writeLines(stamp, stamp_file)
  } else {
    stamp <- readLines(stamp_file, warn = FALSE)[1]
  }

  stamp
}
