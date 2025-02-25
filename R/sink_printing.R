
# Define environment
sink_env <- new.env( parent = globalenv() )

# Helper functions for printing to output.txt

# Aesthetic parameters ####

snk.step_number <- 0
snk.section_break_reps <- 120
snk.simple_break_reps <- 40
snk.small_break_reps <- 15

# Horizontal rule marks ####

snk.horizontal_rule <- function(
    mark = "-",
    reps = snk.section_break_reps,
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

snk.report <- function(
    label,
    initial_breaks = 2,
    end_breaks = 0
  ) {
    cat(
      paste(rep("\n", initial_breaks)), 
      label, ": ", 
      paste(rep("\n", end_breaks)),
      sep=""
      )
  }

snk.report... <- function(
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

snk.step_heading <- function(
    title, 
    prefix = "Step", 
    reset = FALSE, 
    initial_breaks = 2,
    end_breaks = 0
  ) {
    snk.horizontal_rule("=", initial_breaks = initial_breaks)
    if (!reset) snk.step_number <<- snk.step_number + 1
    else snk.step_number <<- 1
    cat(
      prefix, " ", snk.step_number, ": ", 
      title, 
      paste(rep("\n", end_breaks)), 
      sep="")
  }

snk.summary_heading <- function() {
    cat("Summary:\n", sep="")
    cat(rep("-", 40), "\n", sep="")
  }

snk.print_details <- function(
    label, 
    details
  ) {
    cat("\n\n", label, ":\n", sep="")
    cat(rep("-", 30), "\n", sep="")
    cat(details)
    cat("\n")
    cat(rep("-", 4), sep="")
  }

# Print large objects ####

snk.print_var_charvec <- function(
    main,
    vars
  ) {
    # vars should be a character vector of variable names
    if (class(vars) != "character") stop("vars must be a character vector")
    cat("\n\n", main, ":", sep="")
    for (var in vars) {
      cat("\n\t", var, ": ", paste(round(get(var), 2), collapse = ", "), sep="")
    }
  }

snk.print_var_list <- function(
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

snk.print_function <- function(func) {
  if (class(func) != "character") stop("func must be a character string naming a function")
  cat("\n\n", func, " definition:\n",rep("-", 30),"\n", sep="")
  print(get(func))
  cat(rep("-", 4), sep = "")
}

snk.print_table <- function(
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

snk.print_vec <- function(
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
