
# Internal console helpers for the verbose = TRUE progress output printed during model fitting and
# data preprocessing (R/wspmm.R, R/merfish_preprocessing.R). Not exported -- see leading dot and @noRd.

# Aesthetic parameters ####

.snk_section_break_reps <- 120
.snk_simple_break_reps <- 40

# Horizontal rule marks ####

#' @noRd
.snk_horizontal_rule <- function(
    mark = "-",
    reps = .snk_section_break_reps,
    initial_breaks = 1,
    end_breaks = 1
  ) {
    cat(
      paste(rep("\n", initial_breaks)),
      rep(mark, reps),
      paste(rep("\n", end_breaks)),
      sep=""
      )
  }

# Report headers ####

#' @noRd
.snk_report <- function(
    label,
    initial_breaks = 2,
    end_breaks = 0
  ) {
    cat(
      paste(rep("\n", initial_breaks)),
      label, #": ",
      paste(rep("\n", end_breaks)),
      sep=""
      )
  }

#' @noRd
.snk_report_progress <- function(
    label,
    initial_breaks = 1,
    end_breaks = 0
  ) {
    cat(
      paste(rep("\n", initial_breaks)),
      label, "... ",
      paste(rep("\n", end_breaks)),
      sep=""
    )
  }

# Print large objects ####

#' @noRd
.snk_print_var_list <- function(
    main,
    vars,
    vert = TRUE,
    end_breaks = 1
  ) {
    # vars should be a named list
    cat("\n", main, ":", sep="")
    if (vert) {
      for ( v in seq_along(vars) ) {
        if ( class(vars) == "list" ) cat("\n\t", names(vars)[v], ": ", paste(vars[[v]],collapse = ", "), sep="")
        else cat("\n\t", names(vars)[v], ": ", paste(vars[v],collapse = ", "), sep="")
      }
    } else {
      cat("\n")
      print(vars)
    }
    cat(paste(rep("\n", end_breaks)))
  }

#' @noRd
.snk_print_table <- function(
    label,
    tab,
    head = TRUE,
    initial_breaks = 0,
    end_breaks = 1
  ) {
    if (!any(class(tab) %in% c("data.frame", "array", "matrix", "table"))) stop("tab must be a data frame, array, matrix, or table")
    cat(paste(rep("\n", initial_breaks)))
    if (head) cat("\n", label, " (head only):\n", sep="")
    else cat("\n", label, ":\n", sep="")
    cat(rep("-", 30), "\n", sep="")
    if (head) print(head(tab))
    else print(tab)
    cat(rep("-", 4), sep="")
    cat(paste(rep("\n", end_breaks)))
  }

#' @noRd
.snk_print_vec <- function(
    label,
    vec,
    rnd = 3,
    initial_breaks = 1,
    end_breaks = 0
  ) {
    cat(paste(rep("\n", initial_breaks)))
    if (is.numeric(vec)) cat(label, ": ", paste(round(unlist(vec),rnd), collapse = ", "), sep="")
    else if (class(vec) == "difftime") cat(label, ": ", paste(round(unlist(vec),rnd), collapse = ", "), " (", units(vec), ")", sep="")
    else cat(label, ": ", paste(unlist(vec), collapse = ", "), sep="")
    cat(paste(rep("\n", end_breaks)))
  }
