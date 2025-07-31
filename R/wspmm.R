
#' @useDynLib wispack, .registration = TRUE
#' @import Rcpp
#' @import methods
#' @import RcppEigen 
#' @import dtwclust
NULL

.onLoad <- function(libname, pkgname) {
  Rcpp::loadModule("wspc", TRUE)
  loadNamespace("colorspace")
}

# Main function for generating WSPmm model #############################################################################

#' Fit wisp to count data
#'
#' This function takes a data frame of wisp variables (as columns) and fits a wisp model to it. Statistical analyses and plots are generated from the fitted model.
#'
#' @param count.data Data.frame, data to be modeled, with columns for model variables (count, bin, parent, child, ran, fixedeffects), or equivalent variables as specified in the \code{variables} argument.
#' @param variables List, names of the columns in \code{count.data} that correspond to the model variables. The list should contain only (but not necessarily all) named elements: \code{count}, \code{bin}, \code{parent}, \code{child}, \code{ran}, and \code{fixedeffects}.
#' @param use.median Logical, if TRUE, the median of the resamples is used as the final parameter estimates; if FALSE, the initial fit by L-BFGS is used.
#' @param MCMC.settings List, settings for the MCMC simulation, including \code{MCMC.burnin}, \code{MCMC.steps}, \code{MCMC.step.size}, \code{MCMC.prior}, and \code{MCMC.neighbor.filter}. Default values are provided.
#' @param bootstraps.num Integer, number of bootstrap resamples to perform. If 0, only MCMC is run.
#' @param converged.resamples.only Logical, if TRUE, only resamples with a converged fit are used for statistical analysis; if FALSE, all resamples are used. Applies only to bootstraps. 
#' @param max.fork Integer, maximum number of parallel processes to use for bootstrapping.
#' @param dim.bounds Numeric vector, block boundaries for plotting in rate-count plots. If empty, the argument is ignored.
#' @param verbose Logical, if TRUE, prints information during the fitting process.
#' @param print.child.summaries Logical, if TRUE, prints summaries of each child level. 
#' @param model.settings List, settings for the C++ model, including \code{buffer_factor}, \code{ctol}, \code{max_penalty_at_distance_factor}, \code{LROcutoff}, \code{LROwindow_factor}, \code{rise_threshold_factor}, \code{max_evals}, \code{rng_seed}, and \code{warp_precision}. Default values are provided.
#' @return List giving the results of the fitted model, including: \code{model.component.list}, \code{count.data.summed}, \code{fitted.parameters}, \code{gamma.disperson}, \code{param.names}, \code{fix}, \code{treatment}, \code{grouping.variables}, \code{param.idx0}, \code{settings}, \code{sample.params}, \code{sample.params.bs}, \code{sample.params.MCMC}, \code{diagnostics.bs}, \code{diagnostics.MCMC}, \code{stats}, and \code{plots}
#' @export
wisp <- function(
    # Data to model
    count.data, 
    # Variable labels
    variables = list(), 
    # Settings used on R side
    use.median = FALSE,
    MCMC.settings = list(),
    bootstraps.num = 0, 
    converged.resamples.only = TRUE,
    max.fork = 1,
    dim.bounds = c(), 
    verbose = TRUE,
    print.child.summaries = TRUE,
    # Setting to pass to C++ model
    model.settings = list()
  ) {
   
    # Make reproducible
    ran.seed <- 1234
    set.seed(ran.seed)
    
    # Data checks and parsing ####
    
    # Helper function 
    check_list <- function(
      user_input,
      default_list
    ) {
      list_name <- deparse(substitute(user_input))
      default_list.names <- names(default_list)
      # ... check if variables is provided 
      if (length(user_input) > 0) {
        # ... check that provided input is a list with valid names
        if (class(user_input) != "list") {
          stop(paste0(list_name, " must be a list"))
        } else if (is.null(names(user_input))) {
          stop(paste0(list_name, " must be a list with named elements"))
        } else if (!all(names(user_input) %in% default_list.names)) {
          stop(paste0(list_name, " can contain only the following elements: ", paste(default_list.names, collapse = ", ")))
        } else {
          # ... check that each element is of the correct type
          for (v in names(user_input)) {
            if (length(default_list[[v]]) > 0) {
              expected_type <- class(default_list[[v]]) 
              if (class(user_input[[v]]) != expected_type) {
                stop(paste0(list_name, "$", v, " must be of type ", expected_type))
              }
            }
          }
        }
      }
      return(default_list.names)
    }
    
    # Parse variable names
    # ... names of columns in count.data giving each model variable
    variables.internal <- list( 
      count = "count",
      bin = "bin", 
      parent = "parent", 
      child = "child",
      ran = "ran",
      fixedeffects = c()
    )
    # ... check that provided variables is a list with valid names
    variables.names <- check_list(variables, variables.internal)
    variables.names <- variables.names[-c(which(variables.names == "fixedeffects"))]  # remove fixedeffects from list of variable names
    # ... load 
    for (v in names(variables)) {
      variables.internal[[v]] <- variables[[v]]
    }
    
    # Relabel and rearrange data columns 
    if (class(count.data) != "data.frame") {
      stop("count.data must be a data frame")
    }
    old_names <- colnames(count.data)
    ordered_cols <- unlist(variables.internal[variables.names])
    if (!all(ordered_cols %in% old_names)) {
      stop(
        paste0(
          "Not all needed variable names found as column names in count.data, missing: ", 
          paste(ordered_cols[!(ordered_cols %in% old_names)], collapse = ", ")))
    }
    if (length(variables.internal$fixedeffects) == 0) {
      # ... extract fixed effect names, if not given
      fe_cols <- !(old_names %in% ordered_cols)
      variables.internal$fixedeffects <- old_names[fe_cols]
    } else {
      # ... if fixed effects name given, load
      fe_cols <- old_names %in% variables.internal$fixedeffects
      if (!all(variables.internal$fixedeffects %in% old_names)) {
        stop(
          paste0(
            "Not all fixed effect names found as column names in count.data, missing: ", 
            paste(variables.internal$fixedeffects[!(variables.internal$fixedeffects %in% old_names)], collapse = ", ")))
      }
    } 
    if (sum(fe_cols) == 0) {
      stop("No fixed effects found in count.data")
    } 
    new_names <- c(variables.names, old_names[fe_cols])
    data <- cbind(count.data[,ordered_cols], count.data[,fe_cols])
    colnames(data) <- new_names
    
    # Parse model settings 
    # ... define default values
    model.settings.internal <- list(
      buffer_factor = 0.05,                       # buffer factor for minimum distance between t-points
      ctol = 1e-6,                                # convergence tolerance
      max_penalty_at_distance_factor = 0.01,      # maximum penalty at distance from structural parameter values
      LROcutoff = 2.0,                            # cutoff for LROcp, a multiple of standard deviation
      LROwindow_factor = 1.25,                    # controls size of window used in LROcp algorithm (window = LROwindow_factor * bin_num * buffer_factor)
      rise_threshold_factor = 0.8,                # amount of detected rise as fraction of total required to end run
      max_evals = 1000,                           # maximum number of evaluations for optimization
      rng_seed = 42,                              # seed for random number generator
      warp_precision = 1e-7                       # decimal precision to retain when selecting really big number as pseudo infinity for unbound warping
    )
    # ... check that provided model.settings is a list with valid names
    model.settings.names <- check_list(model.settings, model.settings.internal)
    # ... check and load values
    for (s in names(model.settings)) {
      ms <- model.settings[[s]]
      if (!(ms > 0)) {
        stop("All model.settings values must be positive numbers")
      } else if (s == "buffer_factor" && !(ms < 1)) {
        stop("model.settings$buffer_factor must be a number between 0 and 1")
      } 
      if (s == "max_evals" || s == "rng_seed") {
        ms <- as.integer(ms)
      }
      if (
        s == "ctol" || 
        s == "max_penalty_at_distance_factor" || 
        s == "rise_threshold_factor" || 
        s == "warp_precision") {
        if (ms > 1) {
          warning(paste0("model.settings$", s, " should be a number less than 1"))
        }
      }
      # ... load value 
      model.settings.internal[[s]] <- ms
    }
    
    # Add inf_warp 
    model.settings.internal$inf_warp <- model.settings.internal$warp_precision / .Machine$double.eps
    
    # Parse MCMC settings
    MCMC.settings.internal <- list(
      MCMC.burnin = 0,
      MCMC.steps = 1e3,
      MCMC.step.size = 1.0,
      MCMC.prior = 1.0,
      MCMC.neighbor.filter = 2
    )
    # ... check that provided MCMC.settings is a list with valid names
    MCMC.settings.names <- check_list(MCMC.settings, MCMC.settings.internal)
    # ... check and load values
    for (s in names(MCMC.settings)) {
      ms <- MCMC.settings[[s]]
      if (!(ms >= 0)) {
        stop("All MCMC.settings values must be >= 0")
      } else if (s == "MCMC.steps" && ms < 100) {
        warning(paste0("MCMC.settings$", s, " should be at least 100"))
      } else if (s == "MCMC.step.size" && ms >= 2.0) {
        warning(paste0("Consider setting MCMC.settings$", s, " below 2.0"))
      } 
      # ... load value 
      MCMC.settings.internal[[s]] <- ms
    }
    
    if (verbose) {
      snk.report("Parsing data and settings for wisp model")
      snk.horizontal_rule(reps = snk.simple_break_reps, end_breaks = 1)
      snk.print_var_list("Model settings", model.settings.internal, vert = TRUE, end_breaks = 1)
      snk.print_var_list("MCMC settings", MCMC.settings.internal, vert = TRUE, end_breaks = 1)
      snk.print_var_list("Variable dictionary", variables.internal, vert = TRUE, end_breaks = 1)
      snk.print_table("Parsed data", data, end_breaks = 0)
    }
    
    # Initialize cpp model ####
    if (verbose) {
      snk.report("Initializing Cpp (wspc) model")
      snk.horizontal_rule(reps = snk.simple_break_reps, end_breaks = 1)
    }
    cpp_model <- new(
      wspc, 
      data, 
      model.settings.internal,
      verbose
    )
    
    # Estimate model parameters with MCMC or bootstrapping ####
    if (verbose) {
      snk.report("Estimating model parameters", initial_breaks = 1)
      snk.horizontal_rule(reps = snk.simple_break_reps, end_breaks = 0)
    }
    
    # Confirm forking is possible
    if (!(Sys.info()["sysname"] == "Darwin" || Sys.info()["sysname"] == "Linux")) {
      if (bootstraps.num > 0 && max.fork > 1) {
        if (verbose) snk.report...("Forking not available on Windows, cannot bootstrap in parallel, setting max fork to 1")
        max.fork <- 1
      }
    } 
    
    # Run MCMC simulation
    if (verbose) snk.report("Running MCMC stimulations", end_breaks = 1)
    start_time_MCMC <- Sys.time()
    MCMC_walk <- cpp_model$MCMC(
      MCMC.settings.internal$MCMC.steps + MCMC.settings.internal$MCMC.burnin, 
      MCMC.settings.internal$MCMC.neighbor.filter,
      MCMC.settings.internal$MCMC.step.size,
      MCMC.settings.internal$MCMC.prior,
      verbose 
    )
    run_time_MCMC <- Sys.time() - start_time_MCMC
    if (verbose) {
      snk.report...("MCMC simulation complete")
      snk.print_vec("MCMC run time (total), minutes", c(as.numeric(run_time_MCMC, units = "mins")))
      snk.print_vec("MCMC run time (per retained step), seconds", c(as.numeric(run_time_MCMC, units = "secs") / (MCMC.settings.internal$MCMC.steps + MCMC.settings.internal$MCMC.burnin)))
      snk.print_vec("MCMC run time (per step), seconds", c((as.numeric(run_time_MCMC, units = "secs") / (MCMC.settings.internal$MCMC.steps + MCMC.settings.internal$MCMC.burnin))/MCMC.settings.internal$MCMC.neighbor.filter), end_breaks = 0)
    }
    
    # Clear out burn-in, if any
    if (MCMC.settings.internal$MCMC.burnin > 0) {
      MCMC_walk <- MCMC_walk[-c(2:(2+MCMC.settings.internal$MCMC.burnin-1)),]
    }
    
    if (bootstraps.num > 0) {
      
      # Run bootstrap fits
      start_time_bs <- Sys.time()
      if (verbose) snk.report("Running bootstrap fits", end_breaks = 1)
      sample_results <- cpp_model$bs_batch(
        bootstraps.num, 
        max.fork,
        verbose
      )
      run_time_bs <- Sys.time() - start_time_bs
      if (verbose) {
        snk.report...("Bootstrap simulation complete")
        snk.print_vec("Bootstrap run time (total), minutes", c(as.numeric(run_time_bs, units = "mins")))
        snk.print_vec("Bootstrap run time (per sample), seconds", c(as.numeric(run_time_bs, units = "secs") / bootstraps.num))
        snk.print_vec("Bootstrap run time (per sample, per thread), seconds", c(as.numeric(run_time_bs, units = "secs") * max.fork / bootstraps.num), end_breaks = 0)
      }
      
    }
    
    # Extract results and diagnostics
    # ... Save MCMC estimates and diagnostics
    n_params <- ncol(MCMC_walk) - 4
    sample.params.MCMC <- MCMC_walk[,1:n_params]
    diagnostics.MCMC <- data.frame(
      pen.neg.value = MCMC_walk[,n_params + 1],
      neg.loglik = MCMC_walk[,n_params + 2], 
      acceptance.ratio = MCMC_walk[,n_params + 3],
      ctr.num = MCMC_walk[,n_params + 4]
    )
    if (bootstraps.num > 0) {
      
      # ... Save bs estimates and diagnostics
      n_params <- ncol(sample_results) - 4
      sample.params.bs <- sample_results[,1:n_params]
      diagnostics.bs <- data.frame( 
        pen.neg.value = sample_results[,n_params + 1],
        neg.loglik = sample_results[,n_params + 2], 
        success.code = sample_results[,n_params + 3],
        num.evals = sample_results[,n_params + 4]
      )
      
    } else {
      sample.params.bs <- NULL
      diagnostics.bs <- NULL
    }
    
    # Set resamples for analysis 
    if (is.null(sample.params.bs)) {
      sample.params <- sample.params.MCMC
    } else {
      sample.params <- sample.params.bs
    }
    
    # Set final fitted parameters
    if (use.median) {
      if (verbose) snk.report...("Setting median parameter samples as final parameters", initial_breaks = 1, end_breaks = 1)
      final_parameters <- apply(sample.params, 2, function(x) median(x, na.rm = TRUE))
    } else {
      if (verbose) snk.report...("Setting full-data fit as parameters", initial_breaks = 2, end_breaks = 1)
      final_parameters <- sample.params[nrow(sample.params),]
    }
    cpp_model$set_parameters(
      final_parameters,
      verbose
    )
    
    # Grab model results and add samples
    results <- cpp_model$results()
    results[["sample.params"]] <- sample.params
    results[["sample.params.bs"]] <- sample.params.bs
    results[["sample.params.MCMC"]] <- sample.params.MCMC
    results[["diagnostics.bs"]] <- diagnostics.bs
    results[["diagnostics.MCMC"]] <- diagnostics.MCMC
    
    # Add variable names 
    results[["variables"]] <- variables
    
    # Run statistical analysis ####
    
    # Initialize shell to hold stats
    stats <- list(
      parameters = data.frame(),
      tps = data.frame(),
      residuals = data.frame(),
      residuals.log = data.frame(),
      variance = data.frame()
    )
    results[["stats"]] <- stats
    
    # Run stats on samples
    if (bootstraps.num == 0) converged.resamples.only <- FALSE
    results$stats$parameters <- sample.stats(
      wisp.results = results,
      conv.resamples.only = converged.resamples.only,
      verbose = verbose
    )
    
    # Analyze residuals 
    residuals <- analyze.residuals(
      wisp.results = results,
      verbose = verbose
    )
    results$stats$residuals <- residuals$stats
    results$stats$residuals.log <- residuals$stats.log
    plots.residuals <- residuals$plots
    
    # Make plots of results ####
    
    # Plot MCMC walks, both parameters and negloglik
    plots.MCMC <- plot.MCMC.walks(
      wisp.results = results,
      verbose = verbose
    )
    
    # Plot normality comparison of MCMC and bootstrap estimates
    plots.MCMC.bs.comparison <- plot.MCMC.bs.comparison(
      wisp.results = results,
      verbose = verbose
    )
    
    # Plot effect parameter distribution
    plots.effect.dist <- plot.effect.dist(
      wisp.results = results,
      verbose = verbose
    )
    
    # Make rate plots 
    plots.ratecount <- plot.ratecount(
      wisp.results = results,
      pred.type = "pred",
      count.type = "count",
      dim.boundaries = dim.bounds,
      verbose = verbose
    )
    
    # Make parameter plots 
    plots.parameters <- plot.parameters(
      wisp.results = results,
      print.plots = FALSE, 
      verbose = verbose
    )
    
    # Gather plots
    plots <- list(
      residuals = plots.residuals,
      ratecount = plots.ratecount,
      parameters = plots.parameters,
      MCMC = plots.MCMC,
      parameter.normality = plots.MCMC.bs.comparison, 
      effect.dist = plots.effect.dist
    )
    results[["plots"]] <- plots
    
    # Print summary plots
    if (print.child.summaries) {
      plot.child.summary(
        wisp.results = results,
        these.parents = NULL,
        these.childs = NULL,
        verbose = TRUE
      )
    }
    
    return(results)
    
  }

# Analysis methods #####################################################################################################

#' Compute p-values using ecdf from parameter resamples
#'
#' This function takes a vector \code{mu.B} of sampled values of a variable X, and a single observation \code{mu.obs} of the same variable, and computes the p-value of \code{mu.obs} using the empirical cumulative distribution function (ecdf) of \code{mu.B}.
#'
#' @param mu.B Numeric vector of sampled values of a variable X, e.g., bootstrapped or MCMC estimates, used to make an empirical cumulative distribution function (ecdf).
#' @param mu.obs Numeric value of the observed variable X, e.g., the mean of \code{mu.B} or an actual observation.
#' @return Numeric p-value.
#' @export
pvalues.samples <- function(
    mu.B,       # vector of bootstrapped or MCMC estimates
    mu.obs      # observed value, either mean of mu.B or actual observation
  ) {
    # Basic idea: Instead of centering data > bootstrapping > estimate parameter, bootstrap > estimate parameter > center data
    Fn <- ecdf(mu.B - mean(mu.B))
    abs.mu.obs <- abs(mu.obs)
    return(
      1 - Fn(abs.mu.obs) + Fn(-abs.mu.obs)
    )
  }

#' Estimate p-values and confidence intervals from resampled parameters
#'
#' Runs statistical analysis on the resampled parameters from the wisp function. It computes p-values and confidence intervals for each parameter, adjusting for multiple comparisons using either the Bonferroni correction or the Holm-Bonferroni method.
#'
#' @param wisp.results List, output of the wisp function.
#' @param alpha Numeric value giving significance level for p-values and confidence intervals. Default is 0.05.
#' @param Bonferroni Logical, if TRUE, uses the Bonferroni correction for multiple comparisons; if FALSE, uses the Holm-Bonferroni method. Default is FALSE.
#' @param conv.resamples.only Logical, if TRUE, only resamples with a converged fit are used for statistical analysis; if FALSE, all resamples are used. Default is TRUE.
#' @param verbose Logical, if TRUE, prints information during the statistical analysis.
#' @return Data frame giving, for each parameter, its name, estimate, confidence interval (CI.low, CI.high), p-value, adjusted p-value (p.value.adj), adjusted alpha (alpha.adj), and significance level (significance).
#' @export
sample.stats <- function(
    wisp.results,
    alpha = 0.05,
    Bonferroni = FALSE,
    conv.resamples.only = TRUE,
    verbose = TRUE
  ) {
    
    # Multiple comparisons correction explanation:
    #   - If the Bonferroni correction is used, alpha and p_values are adjusted  
    #       by the number of tests performed.
    #   - If the Holm-Bonferroni method is used, alpha and p_values are adjusted by 
    #       the number of tests performed minus the rank of the p-value plus 1. 
    #   - In both cases, the CI are presented as "(1-alpha)%", but calculated with adjustment.
    
    # Run stats on simulation results
    if (verbose) {
      snk.report("Running stats on simulation results", initial_breaks = 1)
      snk.horizontal_rule(reps = snk.simple_break_reps)
    }
    
    # Grab simulation results
    n_conv <- sum(wisp.results$diagnostics.bs$success.code == 3)
    n_total <- nrow(wisp.results$diagnostics.bs)
    if (conv.resamples.only) {
      if (n_conv / n_total < 0.9) {
        warning_message <- paste0("Warning: Only ", n_conv, " out of ", n_total, " resamples converged to fit. Results may be unreliable. Using all resamples.")
        if (verbose) snk.report...(
          warning_message
        )
        warning(warning_message)
        conv.resamples.only <- FALSE
      }
      if (verbose) snk.report...("Grabbing sample results, only resamples with converged fit")
      sample_results <- wisp.results$sample.params[wisp.results$diagnostics.bs$success.code == 3,]
    } else {
      if (verbose) snk.report...("Grabbing sample results")
      sample_results <- wisp.results$sample.params
    }
    
    if (length(sample_results) == 0) {
      snk.report...("\nNo samples to analyze") 
    } else {
      
      # Grab parameter values 
      if (verbose) snk.report...("Grabbing parameter values") 
      fitted_params <- wisp.results$fitted.parameters
      n_params <- length(fitted_params)
      
      # Grabbing indexes of parameters to test
      fitted_params_names <- wisp.results$param.names
      test_mask <- grepl("beta_Rt|beta_tslope|beta_tpoint", fitted_params_names)
      
      # Compute 95% confidence intervals
      if (verbose) snk.report...("Computing 95% confidence intervals")
      sample.params_ci <- apply(sample_results, 2, quantile, probs = c(alpha/2, 1 - alpha/2))
      
      # Estimate p_values from samples
      if (verbose) snk.report...("Estimating p-values from resampled parameters")
      p_values <- rep(NA, n_params)
      p_values_adj <- rep(NA, n_params)
      sig_marks <- rep(" ", n_params)
      for (n in seq_along(fitted_params)) {
        if (test_mask[n]) {
          p_values[n] <- pvalues.samples(sample_results[,n], fitted_params[n]) 
        }
        if (is.na(p_values[n]) && test_mask[n]) stop("Problem with p-value calculation, NA")
      }
      
      # Sanity check
      num_of_tests <- sum(test_mask)
      if (num_of_tests != sum(!is.na(p_values))) stop("Problem with p-value calculation")
      
      # Check resample size 
      recommended_resample_size <- num_of_tests / alpha
      if (verbose) {
        snk.report(paste0("Recommended resample size for alpha = ", alpha, ", ", num_of_tests, " tests"))
        snk.print_vec("with bootstrapping/MCMC", c(recommended_resample_size))
        snk.print_vec("Actual resample size", c(nrow(sample_results)))
      }
      
      # Estimate significance, adjusting for multiple comparisons
      if (!Bonferroni) {
        # Use Holm-Bonferroni method to calculate p-values
        p_value_order <- order(p_values, na.last = NA)
        p_values_adj[p_value_order] <- p_values[p_value_order] * (num_of_tests:1)
        alpha_adj <- rep(NA, n_params)
        alpha_adj[p_value_order] <- alpha / (num_of_tests:1)
        # Adjust CI
        for (i in p_value_order) {
          sample.params_ci[,i] <- quantile(
            sample_results[,i], 
            probs = c(alpha_adj[i]/2, 1 - alpha_adj[i]/2),
            na.rm = TRUE
          )
        }
      } else {
        # Use Bonferroni method to calculate p-values
        p_values_adj <- rep(p_values * num_of_tests, n_params) 
        alpha_adj <- rep(alpha / num_of_tests, n_params) 
        # Adjust CI
        sample.params_ci <- apply(sample_results, 2, quantile, probs = c(alpha_adj[1]/2, 1 - alpha_adj[1]/2), na.rm = TRUE)
      }
      sig_marks[!is.na(p_values_adj) & p_values_adj < alpha/50] <- "***"                            # < 0.001
      sig_marks[!is.na(p_values_adj) & p_values_adj < alpha/5 & p_values_adj >= alpha/50] <- "**"   # < 0.01, >= 0.001
      sig_marks[!is.na(p_values_adj) & p_values_adj < alpha & p_values_adj >= alpha/5] <- "*"       # < 0.05, >= 0.01
      sig_marks[!is.na(p_values_adj) & p_values_adj >= alpha] <- "ns"                               # >= 0.05
      
      # Make summary table
      stats.parameters <- data.frame(
        "parameter" = fitted_params_names,
        "estimate" = fitted_params,
        "CI.low" = sample.params_ci[1,], 
        "CI.high" = sample.params_ci[2,], 
        "p.value" = p_values, 
        "p.value.adj" = p_values_adj,
        "alpha.adj" = alpha_adj,
        "significance" = sig_marks
      )
      rownames(stats.parameters) <- NULL
      
      # Print summary table of stat analysis 
      if (verbose) snk.print_table("Stat summary", stats.parameters, head = TRUE, initial_breaks = 1)
      return(stats.parameters)
      
    }
    
  }

#' Analyze residuals from wisp fit
#'
#' This function takes a vector \code{mu.B} of sampled values of a variable X, and a single observation \code{mu.obs} of the same variable, and computes the p-value of \code{mu.obs} using the empirical cumulative distribution function (ecdf) of \code{mu.B}.
#'
#' @param wisp.results List, output of the wisp function.
#' @param verbose Logical, if TRUE, prints information during the statistical analysis.
#' @return Numeric p-value.
#' @export
analyze.residuals <- function(
    wisp.results,
    verbose = TRUE
  ) {
    
    if (verbose) {
      snk.report("Analyzing residuals", initial_breaks = 1)
      snk.horizontal_rule(reps = snk.simple_break_reps)
    }
    
    # Font sizes 
    label_size <- 5.5
    title_size <- 20 
    axis_size <- 12 
    legend_size <- 12
    
    if (verbose) snk.report...("Computing residuals")
    # Compute residuals 
    wisp.results$count.data.summed$residuals <- wisp.results$count.data.summed$pred - wisp.results$count.data.summed$count
    wisp.results$count.data.summed$residuals.log <- wisp.results$count.data.summed$pred.log - wisp.results$count.data.summed$count.log
    # ... positive residuals mean we overshot, negative residuals mean undershooting
    
    qq_convolved <- function(
      gp, 
      resids = wisp.results$count.data.summed$residuals.log,
      parents = wisp.results$count.data.summed$parent,
      child = wisp.results$count.data.summed$child,
      dispersion.matrix = wisp.results$gamma.dispersion
    ) {
      
      clean_parents <- parents[!is.na(resids)]
      clean_child <- child[!is.na(resids)]
      resids <- resids[!is.na(resids)] 
      n <- length(resids)
      sd_resid <- sd(resids)
      mean_resid <- mean(resids)
      dispersion <- rep(NA, n)
      for (j in colnames(dispersion.matrix)) {
        for (i in rownames(dispersion.matrix)) {
          mask <- clean_parents == j & clean_child == i
          dispersion[mask] <- dispersion.matrix[i, j]
        }
      }
      resid_order <- order(resids)
      resids <- resids[resid_order]
      dispersion <- dispersion[resid_order]
      
      qnorm_dgamma_sd <- function(norm_obs, norm_mean, X, gamma_expected, gamma_variance) {
        gamma_shape <- gamma_expected * gamma_expected / gamma_variance
        out <- qnorm(norm_obs, mean = norm_mean, sd = X) * dgamma(x = X, shape = gamma_shape, rate = gamma_shape / gamma_expected)
        return(out)
      }
      
      convolved_quantile <- rep(NA, n)
      for (i in 1:n) {
        this_quantile <- integrate(
          f = qnorm_dgamma_sd, 
          lower = 0,
          upper = Inf, 
          norm_obs = (i - 0.5)/n, 
          norm_mean = 0,
          gamma_expected = sd_resid,
          gamma_variance = dispersion[i]
        )$value
        convolved_quantile[i] <- this_quantile
      }
      
      qq1 <- qqnorm(resids, plot.it = FALSE)
      
      # Put into one data frame with a group column
      df <- rbind(
        data.frame(sample = qq1$y, theoretical = qq1$x, group = "Normal"),
        data.frame(sample = resids, theoretical = convolved_quantile, group = "convolved")
      )
      
      # Create the ggplot
      resid_plot <- ggplot(df, aes(x = theoretical, y = sample, color = group)) +
        geom_point(size = 1.5) +  
        geom_abline(intercept = mean_resid, slope = sd_resid, color = "black", linewidth = 0.8, linetype = "dashed") +  
        labs(
          y = "Ordered Log-residuals",
          x = "Theoretical Quantiles",
          title = paste0("Q-Q Plot of Log-residuals (", gp, ")")
        ) +
        scale_color_manual(
          name = "Quantile Distribution",
          labels = c("Normal", "Gamma Convolved"),
          values = c("Normal" = "steelblue", "convolved" = "red")
        ) +
        theme_minimal() +
        theme(
          legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, size = title_size),
          axis.title = element_text(size = axis_size),
          axis.text = element_text(size = axis_size),
          legend.title = element_text(size = legend_size),
          legend.text = element_text(size = legend_size)
        )
      
      return(resid_plot)
      
    }
    
    if (verbose) snk.report...("Making masks")
    # Grab masks for subsets of interest
    mask_list <- list()
    mask_list[["all"]] <- rep(TRUE, ncol(wisp.results$count.data.summed))
    mask_ctr <- 1
    for (rn_lvl in wisp.results$grouping.variables$ran.lvls) {
      mask_list[[paste0("ran_lvl_", rn_lvl)]] <- wisp.results$count.data.summed$ran == rn_lvl
      mask_ctr <- mask_ctr + 1
    }
    for (trt_lvl in wisp.results$treatment$names) {
      mask_list[[paste0("treatment_", trt_lvl)]] <- wisp.results$count.data.summed$treatment == trt_lvl
    }
    for (gvp_lvl in wisp.results$grouping.variables$parent.lvls) {
      for (gv_lvl in wisp.results$grouping.variables$child.lvls) {
        mask_list[[paste0(gvp_lvl,"_",gv_lvl)]] <- wisp.results$count.data.summed$child == gv_lvl & 
          wisp.results$count.data.summed$parent == gvp_lvl
      }
    }
    
    if (verbose) snk.report...("Making plots and saving stats")
    plots.residuals <- list()
    stats.residuals <- data.frame()
    stats.residuals.log <- data.frame()
    for (gp in names(mask_list)) {
      
      # Basic histogram
      df_wide <- na.omit(wisp.results$count.data.summed[mask_list[[gp]],])
      hist_resids_plot <- ggplot(df_wide, aes(x=residuals.log)) +
        geom_histogram(bins = 50, fill = "steelblue", color = "black", na.rm = TRUE) + 
        labs(x = "Log-residual Value", y = "Frequency", title = paste0("Histogram of Log-residuals (",gp,")")) +
        theme_minimal() + 
        theme(
          legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, size = title_size),
          axis.title = element_text(size = axis_size),
          axis.text = element_text(size = axis_size),
          legend.title = element_text(size = legend_size),
          legend.text = element_text(size = legend_size)
        )
      
      # qq-plot
      qq_resids_plot <- qq_convolved(
        gp = gp,
        resids = df_wide$residuals.log,
        parents = df_wide$parent,
        child = df_wide$child,
        dispersion.matrix = wisp.results$gamma.dispersion
      )
      # qq_resids_plot <- ggplot(df_wide, aes(sample = residuals)) +
      #   stat_qq(color = "steelblue", size = 1) +  # Q-Q points
      #   stat_qq_line(color = "black", linetype = "dashed") +  # Reference line
      #   labs(x = "Theoretical Quantiles", y = "Sample Quantiles", 
      #        title = paste0("Q-Q Plot of Log-residuals (", gp, ")")) +
      #   theme_minimal()
      
      # Save separately then combine
      plots.residuals[[paste0(gp,"_hist")]] <- hist_resids_plot
      plots.residuals[[paste0(gp,"_qq")]] <- qq_resids_plot
      residual_plots <- list(hist_resids_plot, qq_resids_plot)
      residual_plots <- do.call(arrangeGrob, c(residual_plots, ncol = 2))
      plots.residuals[[paste0(gp,"_resid")]] <- residual_plots
      
      # Find and save residual summary stats
      stats.residuals.gp <- data.frame(
        group = gp,
        mean = mean(df_wide$residuals, na.rm = TRUE),
        sd = sd(df_wide$residuals, na.rm = TRUE),
        variance = var(df_wide$residuals, na.rm = TRUE)
      )
      stats.residuals <- rbind(stats.residuals, stats.residuals.gp)
      
      # Find and save residual summary stats for log residuals
      stats.residuals.log.gp <- data.frame(
        group = gp,
        mean = mean(df_wide$residuals.log, na.rm = TRUE),
        sd = sd(df_wide$residuals.log, na.rm = TRUE),
        variance = var(df_wide$residuals.log, na.rm = TRUE)
      )
      stats.residuals.log <- rbind(stats.residuals.log, stats.residuals.log.gp)
      
    }
    
    # Print results 
    if (verbose) {
      snk.print_table(
        "Log-residual summary by grouping variables", 
        stats.residuals.log,
        initial_breaks = 1
      )
    }
    
    return(
      list(
        stats = stats.residuals,
        stats.log = stats.residuals.log,
        plots = plots.residuals
      )
    )
    
  }

# Plotting methods #####################################################################################################

#' Plot fitted model and data
#'
#' This function takes wisp model results (a fitted line) and observed data (counts) and plots them together for visual comparison. It can also include independent block boundaries for comparison if provided.
#'
#' @param wisp.results List, output of the wisp function.
#' @param pred.type Character string, the name of the predicted rate column in the count data (e.g., "pred.log" or "pred").
#' @param count.type Character string, the name of the observed count column in the count data (e.g., "count.log" or "count").
#' @param dim.boundaries Numeric vector, independent block boundaries to plot for comparison. If empty, the argument is ignored.
#' @param print.all Logical, if TRUE, prints all plots; if FALSE, only returns plots in list without printing any. 
#' @param y.lim Numeric vector of length 2, limits for the y-axis of the plots. If NA, defaults to automatic limits.
#' @param count.alpha.none Numeric, transparency for count points when random level is "none". If left NA, defaults to 0.25.
#' @param count.alpha.ran Numeric, transparency for count points when random level is not "none". If left NA, defaults to 0.25.
#' @param pred.alpha.none Numeric, transparency for predicted lines when random level is "none". If left NA, defaults to 1.0.
#' @param pred.alpha.ran Numeric, transparency for predicted lines when random level is not "none". If left NA, defaults to 0.9.
#' @param rans.to.print Character vector, list of random levels to include on each child plot. If c(), all random levels are included.
#' @param childs.to.print Character vector, list of child levels to place on their own plot. If c(), all child levels are plotted individually.
#' @param verbose Logical, if TRUE, prints updates about the plotting process.
#' @return List of ggplot objects for rate-count plots.
#' @export
plot.ratecount <- function(
    wisp.results,
    pred.type = "pred.log",
    count.type = "count.log",
    dim.boundaries = c(),
    print.all = FALSE,
    y.lim = NA,
    count.alpha.none = NA, # These values have defaults which will be used if left NA
    count.alpha.ran = NA,
    pred.alpha.none = NA,
    pred.alpha.ran = NA,
    rans.to.print = c(),
    childs.to.print = c(),
    verbose = TRUE
  ) {
    
    if (verbose) snk.report...("Making rate-count plots")
   
    # Grab data and run checks 
    df <- wisp.results$count.data.summed 
    df$treatment <- as.factor(df$treatment)
    df$treatment <- relevel(df$treatment, ref = "ref")
    y_lab <- paste0("Rate (", pred.type, ")")
    if (sum(colnames(df) == pred.type) == 0) stop("pred.type not found in count.data.summed")
    if (sum(colnames(df) == count.type) == 0) stop("count.type not found in count.data.summed")
    
    # aes (aesthetics) settings 
    make_parent_ref <- FALSE
    if (is.na(count.alpha.ran)) count.alpha.ran <- 0.25
    if (is.na(count.alpha.none)) count.alpha.none <- 0.25
    if (is.na(pred.alpha.ran)) pred.alpha.ran <- 0.9
    if (is.na(pred.alpha.none)) pred.alpha.none <- 1.0
    if (length(rans.to.print) == 0) rans.to.print <- unique(df[,"ran"])
    if (length(childs.to.print) == 0) {
      make_parent_ref <- TRUE
      childs.to.print <- wisp.results$grouping.variables$child.lvls
      } 
    count_size <- 1.5
    ran_size <- 0.75
    ran_linetype <- "longdash"
    col_fixEff_ref <- "black"
    boundary_linewidth <- 1
    boundary_linetype <- "dotdash"
    boundary_color <- "darkgray"
    tpoint_linewidth <- 1
    tpoint_linetype <- "dashed"
    
    # Make color palettes
    num_of_colors <- length(wisp.results$grouping.variables$child.lvls)
    colors_hues <- seq(0,360,length.out = num_of_colors + 1)[2:num_of_colors]
    child_colors <- colorspace::qualitative_hcl(
      n = num_of_colors,
      h = c(colors_hues[1], colors_hues[num_of_colors]),
      c = 80,
      l = 60,
      fixup = TRUE,
      alpha = 1,
      palette = NULL,
      rev = FALSE,
      register = ""
    )
    
    num_of_colors <- length(wisp.results$treatment$names) - 1
    colors_hues <- seq(0,360,length.out = num_of_colors + 1)[2:num_of_colors]
    treatment_colors <- colorspace::qualitative_hcl(
      n = num_of_colors,
      h = c(colors_hues[1], colors_hues[num_of_colors]),
      c = 80,
      l = 60,
      fixup = TRUE,
      alpha = 1,
      palette = NULL,
      rev = FALSE,
      register = ""
    )
    treatment_colors <- c( "black", treatment_colors )
    names(treatment_colors) <- wisp.results$treatment$names
    
    # Create function for plotting reference and random effects
    plot_parent_ref <- function(df, fvp) {
      
      # Reference-level filtering 
      ref_idx_fvp <- df[,"treatment"] == "ref" & df[,"parent"] == fvp
      
      # Initial ggplot with jittered points from df
      plot <- ggplot() +
        geom_jitter(
          data = df[ ref_idx_fvp & df[,"ran"] == "none", ], 
          aes(x = bin, y = .data[[count.type]], color = child), 
          width = 0.5, height = 0, alpha = count.alpha.none, size = count_size, na.rm = TRUE
        ) +  
        geom_line(
          data = df[ ref_idx_fvp & df[,"ran"] == "none", ],  
          aes(x = bin, y = .data[[pred.type]], color = child),
          linewidth = 2, alpha = pred.alpha.none, na.rm = TRUE
        ) +  
        labs(y = y_lab, x = "Bin", color = "fixed GV") +
        scale_colour_manual(values = child_colors ) +  
        theme_minimal() +
        ggtitle(paste0("Ref-Class Dynamics and RE (", pred.type, "), ", fvp))
      
      # Add lines for random effects from df
      for (rl in seq_along(unique(df[,"ran"]))) {
        plot <- plot + 
          geom_jitter(
            data = df[ ref_idx_fvp & df[,"ran"] == rl, ], 
            aes(x = bin, y = .data[[count.type]], color = child), 
            width = 0.5, height = 0, alpha = count.alpha.ran, size = count_size, na.rm = TRUE
          ) + 
          geom_line(
            data = df[ ref_idx_fvp & df[,"ran"] == rl, ],                
            aes(x = bin, y = .data[[pred.type]], color = child), 
            linetype = ran_linetype, linewidth = ran_size, alpha = pred.alpha.ran, na.rm = TRUE
          )  
      }
      
      if (length(dim.boundaries) > 0) {
        plot <- plot + geom_vline(
          xintercept = dim.boundaries, 
          color = boundary_color, 
          linetype = boundary_linetype, 
          linewidth = boundary_linewidth, 
          na.rm = TRUE
        )
      }
      
      return(plot)
      
    }
   
    # Create function for plotting fixed effects
    plot_parent_fixEff <- function(df, fvp, fv) {
      
      # Grab treatment level names
      treatment_levels = unique(df[,"treatment"])
      
      # Grab index mask for this parent-child level combination
      gvf_idx_fvp <- df[,"child"] == fv &                         # only rates for this level of fixed-effects grouping 
        df[,"parent"] == fvp                                      # only rates for this parent fixed grouping variable
      
      # Grab index mask for case without random effects
      gvf_idx_fvp_noRanEff <- gvf_idx_fvp & df[,"ran"] == "none"  # only rates without random effects
      
      # Plot fixed effects for each treatment case / interaction, including baseline
      plot <- ggplot() 
      for (fe in 1:length(treatment_levels)) {
        
        # Grab this fixed-effect level
        fe_name <- treatment_levels[fe]
        
        # plot effect without random effects
        plot <- plot + 
          geom_jitter(
            data = df[gvf_idx_fvp_noRanEff & df[,"treatment"] == treatment_levels[fe], ], 
            aes(x = bin, y = .data[[count.type]], color = treatment), 
            width = 0.5, height = 0, alpha = count.alpha.none, size = count_size, na.rm = TRUE
          ) +  
          geom_line(
            data = df[gvf_idx_fvp_noRanEff & df[,"treatment"] == treatment_levels[fe], ], 
            aes(x = bin, y = .data[[pred.type]], color = treatment), 
            linetype = "solid", linewidth = 1.5, alpha = pred.alpha.none, na.rm = TRUE
          ) +  
          labs(y = y_lab, x = "Bin") +
          theme_minimal() +
          ggtitle(paste0("Fixed Effects for ", fv, " within ", fvp, " (", pred.type, ")")) 
        
        # plot effect with random effects 
        for (rl in rans.to.print) { 
          plot <- plot + 
            geom_jitter(
              data = df[ gvf_idx_fvp & df[,"ran"] == rl & df[,"treatment"] == treatment_levels[fe], ], 
              aes(x = bin, y = .data[[count.type]], color = treatment), 
              width = 0.5, height = 0, alpha = count.alpha.ran, size = count_size, na.rm = TRUE
            ) + 
            geom_line(
              data = df[ gvf_idx_fvp & df[,"ran"] == rl & df[,"treatment"] == treatment_levels[fe], ],                
              aes(x = bin, y = .data[[pred.type]], color = treatment), 
              linetype = ran_linetype, linewidth = ran_size, alpha = pred.alpha.ran, na.rm = TRUE
            )  
        }
        
      }
      plot <- plot + scale_color_manual(values = treatment_colors)
      if (length(y.lim) == 2) plot <- plot + ylim(y.lim)
      
      if (length(dim.boundaries) > 0) {
        plot <- plot + geom_vline(
          xintercept = dim.boundaries, 
          color = boundary_color, 
          linetype = boundary_linetype, 
          linewidth = boundary_linewidth, 
          na.rm = TRUE
        )
      }
      
      return(plot)
      
    }
    
    # Helper function
    sqrt_or_next <- function(x) {
      sqrt_x <- sqrt(x)          # Take the square root
      if (sqrt_x == floor(sqrt_x)) {
        return(sqrt_x)           # Return the square root if it's an integer
      } else {
        return(ceiling(sqrt_x))  # Return the next integer up if not
      }
    }
    
    # Start list to hold plots
    plot_list <- list()
    
    # Generate, save, and print the plots
    for (fvp in wisp.results$grouping.variables$parent.lvls) {
      
      # Make reference class and random-effects plot
      if (make_parent_ref) {
        plot_list[[paste0("plot_",pred.type,"_parent_",fvp)]] <- plot_parent_ref(df,fvp)
        if (print.all) {
          print(plot_list[[paste0("plot_",pred.type,"_parent_",fvp)]])
        }
      }
      
      # Make fixed effects plot, for each child level
      for (fv in childs.to.print) {
        # Make the plot
        plot_list[[paste0("plot_",pred.type,"_parent_",fvp,"_fixEff_",fv)]] <- plot_parent_fixEff(df,fvp,fv)
        if (print.all) {
          print(plot_list[[paste0("plot_",pred.type,"_parent_",fvp,"_fixEff_",fv)]])
        }
      }
      
    }
    
    return(plot_list)
    
  }

#' Plot of wisp parameters
#'
#' Function to make nicely formatted violin (or bar) plots of the fitted wisp parameters, including confidence intervals if stat information is available.
#'
#' @param wisp.results List, output of the wisp function.
#' @param child.lvl Character string, the child level to be plotted. If NULL, all child levels are plotted.
#' @param violin Logical, if TRUE, plots violin plots for each parameter; if FALSE, uses bar plots. 
#' @param print.plots Logical, if TRUE, prints the plots to the console; if FALSE, only returns a list of plots without printing.
#' @param child.classes List, a list of character vectors specifying which child levels to include together in plots. If NULL, all child levels are included in a single plot.
#' @param verbose Logical, if TRUE, prints updates about the plotting process.
#' @return List of ggplot objects for parameter plots.
#' @export
plot.parameters <- function(
    wisp.results,
    child.lvl = NULL, # NULL (plot all) or a single child level to be plotted
    violin = TRUE,
    print.plots = TRUE,
    child.classes = NULL,
    verbose = TRUE
  ) {
    
    # aes (aesthetics) settings
    star_gap_factor <- 0.15
    expand_factor <- 0.1
    vscale <- "area"       # default is "area", can also try "width"
    mean_dash_size <- 8
    sig_marks_unit <- "mm" # default is "mm", can try "pt"
    sig_marks_size <- 8
    
    # Checks
    print_stats <- TRUE
    sample.params <- wisp.results$sample.params
    if (length(sample.params) == 0) {
      print_stats <- FALSE
      if (violin) {
        violin <- FALSE
        if (verbose) snk.report...("sample.params not found, setting violin = FALSE and not printing stats")
      } 
    } else {
      downsample_size <- min(1e2, nrow(sample.params))
    }
    
    if (verbose) snk.report...("Making parameter plots")
    
    # Grab fitted parameters 
    fitted_params <- wisp.results$fitted.parameters
    
    # Grab masks and set legend position
    param_names <- wisp.results$param.names
    point_mask <- grepl("tpoint", param_names)
    rate_mask <- grepl("Rt", param_names)
    slope_mask <- grepl("tslope", param_names)
    pointR_mask <- grepl("wfactor_point", param_names) # for random effect on point
    rateR_mask <- grepl("wfactor_rate", param_names) # for random effect on rate
    slopeR_mask <- grepl("wfactor_slope", param_names) # for random effect on slope
    legpos <- "none"
    
    # Format parameter names for nice printing
    param_names_saved <- param_names
    for (gvp_lvl in as.character(wisp.results$grouping.variables$parent.lvls)) {
      param_names <- gsub(paste0("baseline_",gvp_lvl,"_"), "", param_names)
      param_names <- gsub(paste0("_",gvp_lvl), "", param_names)
    }
    for (gv_lvl in as.character(wisp.results$grouping.variables$child.lvls)) {
      param_names <- gsub(paste0("_",gv_lvl,"_"), "_", param_names)
    }
    param_names <- gsub(paste0("wfactor_point_"), "level_", param_names)
    param_names <- gsub(paste0("wfactor_rate_"), "level_", param_names)
    param_names <- gsub(paste0("wfactor_slope_"), "level_", param_names)
    param_names <- gsub(paste0("_X"), "", param_names)
    param_names <- gsub(paste0("beta_"), "", param_names)
    param_names <- gsub(paste0("Tns/Blk"), "T/B", param_names)
    treatment_names_rev <- as.character(wisp.results$treatment$names)
    treatment_names_rev <- treatment_names_rev[length(treatment_names_rev):1]
    for (treatment_name in treatment_names_rev) {
      param_names <- gsub(paste0(gsub("\\*", "\\\\*", treatment_name), "_"), "", param_names)
    }
    
    # Rename fitted_params vector with nice names
    names(fitted_params) <- param_names
    if (print_stats) {
      sample_ci <- rbind(wisp.results$stats$parameters$CI.low, wisp.results$stats$parameters$CI.high)
      sig_marks <- wisp.results$stats$parameters$significance
    }
    
    # Check child.classes
    if (length(child.classes) != 0) {
      if (class(child.classes) != "list") stop("child.classes must be a list")
      for (cc in 1:length(child.classes)) {
        if (length(child.classes[[cc]]) == 0) stop("child.classes must have at least one element")
        if (class(child.classes[[cc]] != "character")) stop("child.classes must be a list of character vectors")
        if (!all(child.classes[[cc]] %in% as.character(wisp.results$grouping.variables$child.lvls))) {
          stop("child.classes must be a list of character vectors containing only levels of the child grouping variable")
          }
      }
    }
    if (length(child.classes) == 0 || length(child.lvl) != 0) {
      # If classes not provided, plot all children together
      # If a single child is provided for plotting, disregard any classes provided
      child.classes <- list(all = as.character(wisp.results$grouping.variables$child.lvls))
      }
    child_class_names <- names(child.classes)
    
    for (gvp_lvl in as.character(wisp.results$grouping.variables$parent.lvls)) {
      
      for (cc in 1:length(child.classes)) {
        
        baseline_df <- data.frame() 
        ranEff_df <- data.frame() 
        
        parameter_comparison_plots <- list()
        
        # Make baseline and ranEff data frames ####
        
        for (gv_lvl in child.classes[[cc]]) {
          
          if (length(child.lvl) != 0 && gv_lvl != child.lvl) next
          
          # Make data frame for baseline plot
          baseline_mask <- grepl(gvp_lvl, param_names_saved) & grepl(gv_lvl, param_names_saved) & grepl("baseline", param_names_saved) 
          
          # Baseline tpoint parameters
          if (violin) {
            baseline_fitted_tpoint <- c(as.matrix(sample.params[sample(1:nrow(sample.params), downsample_size, replace = FALSE), baseline_mask & point_mask]))
            params <- rep(param_names[baseline_mask & point_mask], each = downsample_size)
            means <- rep(fitted_params[baseline_mask & point_mask], each = downsample_size)
          } else {
            baseline_fitted_tpoint <- fitted_params[baseline_mask & point_mask]
            params <- param_names[baseline_mask & point_mask]
            means <- fitted_params[baseline_mask & point_mask]
          }
          baseline_fitted_tpoint_df <- data.frame(
            means = means,
            value = baseline_fitted_tpoint, 
            parameter = params,
            type = rep("tpoint", length(baseline_fitted_tpoint)),
            child = rep(gv_lvl, length(baseline_fitted_tpoint))
          )
          if (print_stats) {
            if (violin) {
              baseline_fitted_tpoint_df$value_low <- rep(sample_ci[1,baseline_mask & point_mask], each = downsample_size)
              baseline_fitted_tpoint_df$value_high <- rep(sample_ci[2,baseline_mask & point_mask], each = downsample_size)
              baseline_fitted_tpoint_df$sig_marks <- rep(sig_marks[baseline_mask & point_mask], each = downsample_size)
            } else {
              baseline_fitted_tpoint_df$value_low <- sample_ci[1,baseline_mask & point_mask]
              baseline_fitted_tpoint_df$value_high <- sample_ci[2,baseline_mask & point_mask] 
              baseline_fitted_tpoint_df$sig_marks <- sig_marks[baseline_mask & point_mask]
            }
          }
          baseline_df <- rbind(baseline_df, baseline_fitted_tpoint_df) 
          
          # Baseline rate parameters
          if (violin) {
            baseline_fitted_rate <- c(as.matrix(sample.params[sample(1:nrow(sample.params), downsample_size, replace = FALSE), baseline_mask & rate_mask]))
            params <- rep(param_names[baseline_mask & rate_mask], each = downsample_size)
            means <- rep(fitted_params[baseline_mask & rate_mask], each = downsample_size)
          } else {
            baseline_fitted_rate <- fitted_params[baseline_mask & rate_mask]
            params <- param_names[baseline_mask & rate_mask]
            means <- fitted_params[baseline_mask & rate_mask]
          }
          baseline_fitted_rate_df <- data.frame(
            means = means, 
            value = baseline_fitted_rate, 
            parameter = params,
            type = rep("rate", length(baseline_fitted_rate)),
            child = rep(gv_lvl, length(baseline_fitted_rate))
          )
          if (print_stats) {
            if (violin) {
              baseline_fitted_rate_df$value_low <- rep(sample_ci[1,baseline_mask & rate_mask], each = downsample_size)
              baseline_fitted_rate_df$value_high <- rep(sample_ci[2,baseline_mask & rate_mask], each = downsample_size)
              baseline_fitted_rate_df$sig_marks <- rep(sig_marks[baseline_mask & rate_mask], each = downsample_size)
            } else {
              baseline_fitted_rate_df$value_low <- sample_ci[1,baseline_mask & rate_mask]
              baseline_fitted_rate_df$value_high <- sample_ci[2,baseline_mask & rate_mask]
              baseline_fitted_rate_df$sig_marks <- sig_marks[baseline_mask & rate_mask]
            }
          }
          baseline_df <- rbind(baseline_df, baseline_fitted_rate_df) 
          
          # Baseline tslope parameters 
          if (violin) {
            baseline_fitted_slope <- c(as.matrix(sample.params[sample(1:nrow(sample.params), downsample_size, replace = FALSE), baseline_mask & slope_mask]))
            params <- rep(param_names[baseline_mask & slope_mask], each = downsample_size)
            means <- rep(fitted_params[baseline_mask & slope_mask], each = downsample_size)
          } else {
            baseline_fitted_slope <- fitted_params[baseline_mask & slope_mask]
            params <- param_names[baseline_mask & slope_mask]
            means <- fitted_params[baseline_mask & slope_mask]
          }
          baseline_fitted_slope_df <- data.frame(
            means = means, 
            value = baseline_fitted_slope, 
            parameter = params,
            type = rep("tslope", length(baseline_fitted_slope)),
            child = rep(gv_lvl, length(baseline_fitted_slope))
          )
          if (print_stats) {
            if (violin) {
              baseline_fitted_slope_df$value_low <- rep(sample_ci[1,baseline_mask & slope_mask], each = downsample_size)
              baseline_fitted_slope_df$value_high <- rep(sample_ci[2,baseline_mask & slope_mask], each = downsample_size)
              baseline_fitted_slope_df$sig_marks <- rep(sig_marks[baseline_mask & slope_mask], each = downsample_size)
            } else {
              baseline_fitted_slope_df$value_low <- sample_ci[1,baseline_mask & slope_mask]
              baseline_fitted_slope_df$value_high <- sample_ci[2,baseline_mask & slope_mask]
              baseline_fitted_slope_df$sig_marks <- sig_marks[baseline_mask & slope_mask]
            }
          }
          baseline_df <- rbind(baseline_df, baseline_fitted_slope_df) 
          
          # Make data frame for ranEff plot
          ranEff_mask <- grepl(gv_lvl, param_names_saved) & grepl("wfactor", param_names_saved) 
          
          # Random effects on tpoints
          if (violin) {
            ranEff_fitted_point <- c(as.matrix(sample.params[sample(1:nrow(sample.params), downsample_size, replace = FALSE), ranEff_mask & pointR_mask]))
            params <- rep(param_names[ranEff_mask & pointR_mask], each = downsample_size)
            means <- rep(fitted_params[ranEff_mask & pointR_mask], each = downsample_size)
          } else {
            ranEff_fitted_point <- fitted_params[ranEff_mask & pointR_mask]
            params <- param_names[ranEff_mask & pointR_mask]
            means <- fitted_params[ranEff_mask & pointR_mask]
          }
          if (sum(ranEff_fitted_point, na.rm = TRUE) != 0) {
            # If zero, this was a degree-zero child, and we don't want to plot it
            ranEff_fitted_point_df <- data.frame(
              means = means, 
              value = ranEff_fitted_point, 
              parameter = params, 
              type = rep("point", length(ranEff_fitted_point)),
              child = rep(gv_lvl, length(ranEff_fitted_point))
            )
            if (print_stats) {
              if (violin) {
                ranEff_fitted_point_df$value_low <- rep(sample_ci[1,ranEff_mask & pointR_mask], each = downsample_size)
                ranEff_fitted_point_df$value_high <- rep(sample_ci[2,ranEff_mask & pointR_mask], each = downsample_size)
                ranEff_fitted_point_df$sig_marks <- rep(sig_marks[ranEff_mask & pointR_mask], each = downsample_size)
              } else {
                ranEff_fitted_point_df$value_low <- sample_ci[1,ranEff_mask & pointR_mask]
                ranEff_fitted_point_df$value_high <- sample_ci[2,ranEff_mask & pointR_mask]
                ranEff_fitted_point_df$sig_marks <-sig_marks[ranEff_mask & pointR_mask]
              }
            }
            ranEff_df <- rbind(ranEff_df, ranEff_fitted_point_df) 
          }
          
          # Random effects on rates
          if (violin) {
            ranEff_fitted_rate <- c(as.matrix(sample.params[sample(1:nrow(sample.params), downsample_size, replace = FALSE), ranEff_mask & rateR_mask]))
            params <- rep(param_names[ranEff_mask & rateR_mask], each = downsample_size)
            means <- rep(fitted_params[ranEff_mask & rateR_mask], each = downsample_size)
          } else {
            ranEff_fitted_rate <- fitted_params[ranEff_mask & rateR_mask]
            params <- param_names[ranEff_mask & rateR_mask]
            means <- fitted_params[ranEff_mask & rateR_mask]
          }
          ranEff_fitted_rate_df <- data.frame(
            means = means, 
            value = ranEff_fitted_rate, 
            parameter = params, 
            type = rep("rate", length(ranEff_fitted_rate)),
            child = rep(gv_lvl, length(ranEff_fitted_rate))
          )
          if (print_stats) {
            if (violin) {
              ranEff_fitted_rate_df$value_low <- rep(sample_ci[1,ranEff_mask & rateR_mask], each = downsample_size)
              ranEff_fitted_rate_df$value_high <- rep(sample_ci[2,ranEff_mask & rateR_mask], each = downsample_size)
              ranEff_fitted_rate_df$sig_marks <- rep(sig_marks[ranEff_mask & rateR_mask], each = downsample_size)
            } else {
              ranEff_fitted_rate_df$value_low <- sample_ci[1,ranEff_mask & rateR_mask]
              ranEff_fitted_rate_df$value_high <- sample_ci[2,ranEff_mask & rateR_mask]
              ranEff_fitted_rate_df$sig_marks <- sig_marks[ranEff_mask & rateR_mask]
            }
          }
          ranEff_df <- rbind(ranEff_df, ranEff_fitted_rate_df) 
          
          # Random effects on slopes
          if (violin) {
            ranEff_fitted_slope <- c(as.matrix(sample.params[sample(1:nrow(sample.params), downsample_size, replace = FALSE), ranEff_mask & slopeR_mask]))
            params <- rep(param_names[ranEff_mask & slopeR_mask], each = downsample_size)
            means <- rep(fitted_params[ranEff_mask & slopeR_mask], each = downsample_size)
          } else {
            ranEff_fitted_slope <- fitted_params[ranEff_mask & slopeR_mask]
            params <- param_names[ranEff_mask & slopeR_mask]
            means <- fitted_params[ranEff_mask & slopeR_mask]
          }
          ranEff_fitted_slope_df <- data.frame(
            means = means, 
            value = ranEff_fitted_slope, 
            parameter = params, 
            type = rep("slope", length(ranEff_fitted_slope)),
            child = rep(gv_lvl, length(ranEff_fitted_slope))
          )
          if (print_stats) {
            if (violin) {
              ranEff_fitted_slope_df$value_low <- rep(sample_ci[1,ranEff_mask & slopeR_mask], each = downsample_size)
              ranEff_fitted_slope_df$value_high <- rep(sample_ci[2,ranEff_mask & slopeR_mask], each = downsample_size)
              ranEff_fitted_slope_df$sig_marks <- rep(sig_marks[ranEff_mask & slopeR_mask], each = downsample_size)
            } else {
              ranEff_fitted_slope_df$value_low <- sample_ci[1,ranEff_mask & slopeR_mask]
              ranEff_fitted_slope_df$value_high <- sample_ci[2,ranEff_mask & slopeR_mask]
              ranEff_fitted_slope_df$sig_marks <- sig_marks[ranEff_mask & slopeR_mask]
            }
          }
          ranEff_df <- rbind(ranEff_df, ranEff_fitted_slope_df)
          
        }
        
        # Make character columns into factors
        baseline_df$child <- as.factor(baseline_df$child)
        baseline_df$parameter <- as.factor(baseline_df$parameter)
        baseline_df$type <- as.factor(baseline_df$type)
        ranEff_df$child <- as.factor(ranEff_df$child)
        ranEff_df$parameter <- as.factor(ranEff_df$parameter)
        ranEff_df$type <- as.factor(ranEff_df$type)
        
        # Make baseline plots ####
        
        title_bl <- paste0("Baseline parameters, fitted (", gvp_lvl, ")")
        plot_name <- paste0("plot_baseline_", gvp_lvl, "_", child_class_names[cc])
        if (length(child.lvl) != 0) plot_name <- paste0("plot_baseline_", gvp_lvl, "_", child.lvl)
        
        parameter_comparison_plots[[plot_name]] <- ggplot(baseline_df, aes(x = parameter, y = value, fill = type)) 
        if (violin) {
          parameter_comparison_plots[[plot_name]] <- parameter_comparison_plots[[plot_name]] +
            geom_violin(aes(y = value), scale = vscale, na.rm = TRUE) +
            geom_point(aes(y = means), shape = 95, size = mean_dash_size, na.rm = TRUE) 
        } else {
          parameter_comparison_plots[[plot_name]] <- parameter_comparison_plots[[plot_name]] +
            geom_bar(aes(y = value), stat = "identity", na.rm = TRUE)
        }
        
        parameter_comparison_plots[[plot_name]] <- parameter_comparison_plots[[plot_name]] +
          facet_wrap(~ interaction(child,type), scales = "free") +
          theme_minimal() +
          labs(title = title_bl) +
          theme(legend.position = legpos) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1))
        
        if (print_stats) {
          parameter_comparison_plots[[plot_name]] <- parameter_comparison_plots[[plot_name]] + 
            geom_errorbar(
              aes(ymin = value_low, ymax = value_high), 
              width = 0.2, na.rm = TRUE) +
            geom_text(
              aes(x = parameter, 
                  y = value_high + star_gap_factor*(value_high-value_low),
                  label = sig_marks,
                  size = sig_marks_size),
              size.unit = sig_marks_unit, 
              na.rm = TRUE
            ) 
        }
        
        # Make ranEff plots ####
        
        title_ran <- paste0("RanEff parameters, fitted (", gvp_lvl, ")")
        plot_name <- paste0("plot_ranEff_", gvp_lvl, "_", child_class_names[cc])
        if (length(child.lvl) != 0) plot_name <- paste0("plot_ranEff_", gvp_lvl, "_", child.lvl)
        
        parameter_comparison_plots[[plot_name]] <- ggplot(ranEff_df, aes(x = parameter, y = value, fill = type)) 
        if (violin) {
          parameter_comparison_plots[[plot_name]] <- parameter_comparison_plots[[plot_name]] + 
            geom_violin(scale = vscale, na.rm = TRUE) + 
            geom_point(aes(y = means), shape = 95, size = mean_dash_size, na.rm = TRUE)
        } else {
          parameter_comparison_plots[[plot_name]] <- parameter_comparison_plots[[plot_name]] + 
            geom_bar(stat = "identity", na.rm = TRUE)
        }
        
        parameter_comparison_plots[[plot_name]] <- parameter_comparison_plots[[plot_name]] +
          facet_wrap(~ interaction(child, type), scales = "free") +
          theme_minimal() +
          labs(title = title_ran) +
          theme(legend.position = legpos) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
          scale_fill_manual(values = c("point" = "#c356ea", "rate" = "#ffc100", "slope" = "#00a99d")) 
        
        if (print_stats) {
          parameter_comparison_plots[[plot_name]] <- parameter_comparison_plots[[plot_name]] + 
            geom_errorbar(
              aes(ymin = value_low, ymax = value_high),
              width = 0.2, na.rm = TRUE) +
            geom_text(
              aes(x = parameter, 
                  y = value_high + star_gap_factor*(value_high-value_low),
                  label = sig_marks,
                  size = sig_marks_size),
              size.unit = sig_marks_unit, 
              na.rm = TRUE
            ) +
            scale_y_continuous(
              expand = expansion(mult = c(0, expand_factor))
            )
        }
        
        # Make treatment data frames and plots ####
        
        treatment_dfL <- list() 
        for (treatment in as.character(wisp.results$treatment$names)) { 
          if (treatment == "ref") next # ... covered in baseline plots above
          
          # Make treatment data frames ####
          treatment_dfL[[treatment]] <- data.frame()
          for (gv_lvl in child.classes[[cc]]) {
            
            # If only printing one child level, skip the rest
            if (length(child.lvl) != 0 && gv_lvl != child.lvl) next
            
            # Grab index masks
            treatment_mask <- grepl(gvp_lvl, param_names_saved) & 
              grepl(gv_lvl, param_names_saved) & 
              grepl("beta", param_names_saved) & 
              !grepl("beta_shape", param_names_saved) &
              grepl(gsub("\\*", "\\\\*", paste0("_",treatment,"_")), param_names_saved)
            
            # Add rows for rate effects
            if (violin) {
              treatment_fitted <- c(as.matrix(sample.params[sample(1:nrow(sample.params), downsample_size, replace = FALSE), treatment_mask & rate_mask]))
              params <- rep(param_names[treatment_mask & rate_mask], each = downsample_size)
              means <- rep(fitted_params[treatment_mask & rate_mask], each = downsample_size)
            } else {
              treatment_fitted <- fitted_params[treatment_mask & rate_mask]
              params <- param_names[treatment_mask & rate_mask]
              means <- fitted_params[treatment_mask & rate_mask]
            }
            treatment_fitted_df <- data.frame(
              means = means,
              value = treatment_fitted, 
              parameter = params, 
              type = rep("rate", length(treatment_fitted)),
              child = rep(gv_lvl, length(treatment_fitted))
            )
            if (print_stats) {
              if (violin) {
                treatment_fitted_df$value_low <- rep(sample_ci[1,treatment_mask & rate_mask], each = downsample_size)
                treatment_fitted_df$value_high <- rep(sample_ci[2,treatment_mask & rate_mask], each = downsample_size)
                treatment_fitted_df$sig_marks <- rep(sig_marks[treatment_mask & rate_mask], each = downsample_size)
              } else {
                treatment_fitted_df$value_low <- sample_ci[1,treatment_mask & rate_mask]
                treatment_fitted_df$value_high <- sample_ci[2,treatment_mask & rate_mask]
                treatment_fitted_df$sig_marks <- sig_marks[treatment_mask & rate_mask]
              }
            }
            treatment_dfL[[treatment]] <- rbind(treatment_dfL[[treatment]], treatment_fitted_df) 
            
            # Add rows for slope effects
            if (violin) {
              treatment_fitted <- c(as.matrix(sample.params[sample(1:nrow(sample.params), downsample_size, replace = FALSE), treatment_mask & slope_mask]))
              params <- rep(param_names[treatment_mask & slope_mask], each = downsample_size)
              means <- rep(fitted_params[treatment_mask & slope_mask], each = downsample_size)
            } else {
              treatment_fitted <- fitted_params[treatment_mask & slope_mask]
              params <- param_names[treatment_mask & slope_mask]
              means <- fitted_params[treatment_mask & slope_mask]
            }
            treatment_fitted_df <- data.frame(
              means = means, 
              value = treatment_fitted, 
              parameter = params, 
              type = rep("tslope", length(treatment_fitted)),
              child = rep(gv_lvl, length(treatment_fitted))
            )
            if (print_stats) {
              if (violin) {
                treatment_fitted_df$value_low <- rep(sample_ci[1,treatment_mask & slope_mask], each = downsample_size)
                treatment_fitted_df$value_high <- rep(sample_ci[2,treatment_mask & slope_mask], each = downsample_size)
                treatment_fitted_df$sig_marks <- rep(sig_marks[treatment_mask & slope_mask], each = downsample_size)
              } else {
                treatment_fitted_df$value_low <- sample_ci[1,treatment_mask & slope_mask]
                treatment_fitted_df$value_high <- sample_ci[2,treatment_mask & slope_mask]
                treatment_fitted_df$sig_marks <- sig_marks[treatment_mask & slope_mask]
              }
            }
            treatment_dfL[[treatment]] <- rbind(treatment_dfL[[treatment]], treatment_fitted_df) 
            
            # Add rows for point effects
            if (violin) {
              treatment_fitted <- c(as.matrix(sample.params[sample(1:nrow(sample.params), downsample_size, replace = FALSE), treatment_mask & point_mask]))
              params <- rep(param_names[treatment_mask & point_mask], each = downsample_size)
              means <- rep(fitted_params[treatment_mask & point_mask], each = downsample_size)
            } else {
              treatment_fitted <- fitted_params[treatment_mask & point_mask]
              params <- param_names[treatment_mask & point_mask]
              means <- fitted_params[treatment_mask & point_mask]
            }
            treatment_fitted_df <- data.frame(
              means = means, 
              value = treatment_fitted, 
              parameter = params, 
              type = rep("tpoint", length(treatment_fitted)),
              child = rep(gv_lvl, length(treatment_fitted))
            )
            if (print_stats) {
              if (violin) {
                treatment_fitted_df$value_low <- rep(sample_ci[1,treatment_mask & point_mask], each = downsample_size)
                treatment_fitted_df$value_high <- rep(sample_ci[2,treatment_mask & point_mask], each = downsample_size)
                treatment_fitted_df$sig_marks <- rep(sig_marks[treatment_mask & point_mask], each = downsample_size)
              } else {
                treatment_fitted_df$value_low <- sample_ci[1,treatment_mask & point_mask]
                treatment_fitted_df$value_high <- sample_ci[2,treatment_mask & point_mask]
                treatment_fitted_df$sig_marks <- sig_marks[treatment_mask & point_mask]
              }
            }
            treatment_dfL[[treatment]] <- rbind(treatment_dfL[[treatment]], treatment_fitted_df) 
            
          }
          
          # Make treatment plots ####
          title_fe <- paste0("treatment parameters, fitted (", gvp_lvl, ", ", treatment, ")")
          plot_name <- paste0("plot_treatment_",treatment,"_",gvp_lvl, "_", child_class_names[cc])
          if (length(child.lvl) != 0) plot_name <- paste0("plot_treatment_",treatment,"_", gvp_lvl, "_", child.lvl)
          
          treatment_dfL[[treatment]]$child <- as.factor(treatment_dfL[[treatment]]$child)
          treatment_dfL[[treatment]]$parameter <- as.factor(treatment_dfL[[treatment]]$parameter)
          treatment_dfL[[treatment]]$type <- as.factor(treatment_dfL[[treatment]]$type)
          
          parameter_comparison_plots[[plot_name]] <- ggplot(treatment_dfL[[treatment]], aes(x = parameter, y = value, fill = type)) 
          if (violin) {
            parameter_comparison_plots[[plot_name]] <- parameter_comparison_plots[[plot_name]] + 
              geom_violin(scale = vscale, na.rm = TRUE) + 
              geom_point(aes(y = means), shape = 95, size = mean_dash_size, na.rm = TRUE)
          } else {
            parameter_comparison_plots[[plot_name]] <- parameter_comparison_plots[[plot_name]] + 
              geom_bar(stat = "identity", na.rm = TRUE)
          }
          
          parameter_comparison_plots[[plot_name]] <- parameter_comparison_plots[[plot_name]] +
            facet_wrap(~ interaction(child,type), scales = "free") +
            theme_minimal() +
            labs(title = title_fe) +
            theme(legend.position = legpos) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))
          
          if (print_stats) {
            parameter_comparison_plots[[plot_name]] <- parameter_comparison_plots[[plot_name]] +
              geom_errorbar(
                aes(ymin = value_low, ymax = value_high),
                width = 0.2, na.rm = TRUE) +
              geom_text(
                aes(x = parameter, 
                    y = value_high + star_gap_factor*(value_high-value_low),
                    label = sig_marks,
                    size = sig_marks_size),
                size.unit = sig_marks_unit, 
                na.rm = TRUE
              )  +
              scale_y_continuous(
                expand = expansion(mult = c(0, expand_factor))
              )
          }
          
        }
        
        # Print and return raw plots ####
        
        if (print.plots) {
          for (i in seq_along(parameter_comparison_plots)) {
            print(parameter_comparison_plots[[i]])
          }
        }
        
        return(parameter_comparison_plots)
        
      }
      
    }
    
  }

#' Plot parameter distributions from WISP results as histograms
#'
#' Function to make nicely formatted histograms of fitted parameters from WISP results.
#'
#' @param wisp.results List, output of the wisp function.
#' @param verbose Logical, if TRUE, prints information during plotting.
#' @return List of ggplot objects containing histograms of fitted parameters.
#' @export
plot.effect.dist <- function(
    wisp.results,
    verbose = TRUE 
  ) {
    
    if (verbose) snk.report...("Making effect parameter distribution plots")
    
    plots.effects_dist <- list()
    bs_fitted_params <- wisp.results$sample.params
    
    # Warping factor distributions
    if (verbose) snk.report...("Warping factors distributions")
    
    # ... for point
    wfactors_point_mask <- grepl("wfactor_point", wisp.results$param.names)
    wfactors_point <- c(bs_fitted_params[,wfactors_point_mask])
    plot_wfactor_point_effects_dist <- ggplot() +
      geom_histogram(
        data = data.frame(vals = wfactors_point), aes(x = vals, y = after_stat(ndensity)),
        bins = 75, fill = "skyblue", alpha = 0.5, na.rm = TRUE) +
      labs(title = "Distribution of random point effects (warping factors)", x = "Warping factor, point", y = "Count") +
      theme_minimal() 
    
    # ... for rate
    wfactors_rate_mask <- grepl("wfactor_rate", wisp.results$param.names)
    wfactors_rate <- c(bs_fitted_params[,wfactors_rate_mask])
    plot_wfactor_rate_effects_dist <- ggplot() +
      geom_histogram(
        data = data.frame(vals = wfactors_rate), aes(x = vals, y = after_stat(ndensity)),
        bins = 75, fill = "skyblue", alpha = 0.5, na.rm = TRUE) +
      labs(title = "Distribution of random rate effects (warping factors)", x = "Warping factor, rate", y = "Count") +
      theme_minimal() 
    
    # ... for slope 
    wfactors_slope_mask <- grepl("wfactor_slope", wisp.results$param.names)
    wfactors_slope <- c(bs_fitted_params[,wfactors_slope_mask])
    plot_wfactor_slope_effects_dist <- ggplot() +
      geom_histogram(
        data = data.frame(vals = wfactors_slope), aes(x = vals, y = after_stat(ndensity)),
        bins = 75, fill = "skyblue", alpha = 0.5, na.rm = TRUE) +
      labs(title = "Distribution of random slope effects (warping factors)", x = "Warping factor, slope", y = "Count") +
      theme_minimal()
    
    if (verbose) print(plot_wfactor_point_effects_dist)
    plots.effects_dist[["plot_wfactor_point_effects_dist"]] <- plot_wfactor_point_effects_dist
    if (verbose) print(plot_wfactor_rate_effects_dist)
    plots.effects_dist[["plot_wfactor_rate_effects_dist"]] <- plot_wfactor_rate_effects_dist
    if (verbose) print(plot_wfactor_slope_effects_dist)
    plots.effects_dist[["plot_wfactor_slope_effects_dist"]] <- plot_wfactor_slope_effects_dist
    
    # Rate effects distribution
    if (verbose) snk.report...("Rate effects distribution")
    rate_effs_mask <- grepl("Rt", wisp.results$param.names) & grepl("beta", wisp.results$param.names)
    rate_effects <- c(bs_fitted_params[,rate_effs_mask])
    plot_rate_effects_effects_dist <- ggplot() +
      geom_histogram(
        data = data.frame(vals = rate_effects), aes(x = vals, y = after_stat(ndensity)),
        bins = 75, fill = "skyblue", alpha = 0.5, na.rm = TRUE) +
      xlim(min(rate_effects),max(rate_effects)) +
      labs(title = "Distribution of fixed rate effects", x = "Rate effect", y = "Count") +
      theme_minimal() 
    
    if (verbose) print(plot_rate_effects_effects_dist)
    plots.effects_dist[["plot_rate_effects_effects_dist"]] <- plot_rate_effects_effects_dist
    
    # Slope effects distribution
    if (verbose) snk.report...("Slope effects distribution")
    tslope_effs_mask <- grepl("tslope", wisp.results$param.names) & grepl("beta", wisp.results$param.names)
    if (any(tslope_effs_mask)) {
      
      tslope_effects <- c(bs_fitted_params[,tslope_effs_mask])
      plot_slope_effects_effects_dist <- ggplot() +
        geom_histogram(
          data = data.frame(vals = tslope_effects), aes(x = vals, y = after_stat(ndensity)),
          bins = 75, fill = "skyblue", alpha = 0.5, na.rm = TRUE) +
        xlim(min(tslope_effects), max(tslope_effects)) +
        labs(title = "Distribution of fixed slope effects", x = "t-slope effect", y = "Count") +
        theme_minimal() 
      
      if (verbose) print(plot_slope_effects_effects_dist)
      plots.effects_dist[["plot_slope_effects_effects_dist"]] <- plot_slope_effects_effects_dist
      
    }
    
    # Point effects distribution
    if (verbose) snk.report...("tpoint effects distribution")
    tpoint_effs_mask <- grepl("tpoint", wisp.results$param.names) & grepl("beta", wisp.results$param.names)
    if (any(tpoint_effs_mask)) {
      
      tpoint_effects <- c(bs_fitted_params[,tpoint_effs_mask])
      plot_tpoint_effects_effects_dist <- ggplot() +
        geom_histogram(
          data = data.frame(vals = tpoint_effects), aes(x = vals, y = after_stat(ndensity)),
          bins = 75, fill = "skyblue", alpha = 0.5, na.rm = TRUE) +
        labs(title = "Distribution of fixed t-point effects", x = "t-point effect", y = "Count") +
        theme_minimal() 
      
      if (verbose) print(plot_tpoint_effects_effects_dist)
      plots.effects_dist[["plot_tpoint_effects_effects_dist"]] <- plot_tpoint_effects_effects_dist
      
    }
    
    return(plots.effects_dist)
    
  }

#' Print rate-count, residual, and parameter plots for one child level together.
#'
#' Function to summarize all important information for an individual child level on one plot.
#'
#' @param wisp.results List, output of the wisp function.
#' @param these.parents Character vector, optional, specifies which parent levels to summarize. Defaults to all. 
#' @param these.childs Character vector, optional, specifies which child levels to summarize. Defaults to all.
#' @param verbose Logical, if TRUE, prints information during plotting.
#' @return Nothing, plots are printed to the current graphics device.
#' @export
plot.child.summary <- function(
    wisp.results,
    these.parents = NULL,
    these.childs = NULL,
    verbose = TRUE
  ) {
    
    # Check for needed plots
    if (
      length(wisp.results$plots$ratecount) == 0 || 
      length(wisp.results$plots$parameters) == 0 ||
      length(wisp.results$plots$residuals) == 0
      ) {
      stop("No rate or parameter plots found in wisp.results")
    }
    
    # Specify parent and child levels to summarize
    gvp_lvls <- as.character(wisp.results$grouping.variables$parent.lvls)
    if (length(these.parents) != 0) gvp_lvls <- gvp_lvls[gvp_lvls %in% these.parents]
    gv_lvls <- as.character(wisp.results$grouping.variables$child.lvls)
    if (length(these.childs) != 0) gv_lvls <- gv_lvls[gv_lvls %in% these.childs]
    
    for (gvp_lvl in gvp_lvls) {
      
      if (verbose) snk.report(paste0("Printing child summary plots for ", gvp_lvl))
      first_print <- TRUE
      
      for (gv_lvl in gv_lvls) {
        
        if (verbose && first_print) {
          cat(gv_lvl)
          first_print <- FALSE
        } else if (verbose) {
          cat(",", gv_lvl)
        }
        wisp.results$plots$parameters <- plot.parameters(
          wisp.results = wisp.results,
          child.lvl = gv_lvl, 
          print.plots = FALSE, 
          verbose = FALSE 
        )
        
        # Grab / make plots
        p1 <- wisp.results$plots$ratecount
        p2 <- wisp.results$plots$parameters
        p3 <- wisp.results$plots$residuals
        
        # Find plots for this parent
        p_mask1 <- grepl(gvp_lvl, names(p1))
        p_mask2 <- grepl(gvp_lvl, names(p2))
        p_mask3 <- grepl(gvp_lvl, names(p3))
        
        # Find plots for this child
        c_mask1 <- grepl(gv_lvl, names(p1))
        c_mask2 <- grepl(gv_lvl, names(p2)) 
        c_mask3 <- grepl(gv_lvl, names(p3))
        
        # Find residual hist and qq plots 
        histqq <- grepl("hist|qq", names(p3))
        
        # Find treatment plots
        iX_mask <- grepl("treatment", names(p2))
        
        # Subset
        p_rates <- p1[p_mask1 & c_mask1]
        p_treatment <- p2[p_mask2 & c_mask2 & iX_mask]
        p_otherparams <- p2[p_mask2 & c_mask2 & !iX_mask]
        p_residuals <- p3[p_mask3 & c_mask3 & histqq]
        
        # Combine and print
        resids <- do.call(arrangeGrob, c(p_residuals, ncol = length(p_residuals)))
        rates_and_residuals <- arrangeGrob(ggplotGrob(p_rates[[1]]), resids, ncol = 1)
        treatments <- do.call(arrangeGrob, c(as.list(p_treatment), ncol = length(p_treatment)))
        other_params <- do.call(arrangeGrob, c(as.list(p_otherparams), ncol = length(p_otherparams)))
        params <- arrangeGrob(treatments, other_params, ncol = 1)
        p <- arrangeGrob(rates_and_residuals, params, ncol = 2, widths = c(0.4,0.6))
        grid.arrange(p)
        
      }
      
    }
    
  }

#' Plot sampling of random walks and negative log likelihood from MCMC simulations
#'
#' Function to make nicely formatted histograms of fitted parameters from WISP results.
#'
#' @param wisp.results List, output of the wisp function.
#' @param verbose Logical, if TRUE, prints information during plotting.
#' @param low_samples Integer, number of low-value parameters to plot. Defaults to 10.
#' @return List of ggplot objects containing plots of walks for low-value parameters, high-value parameters, and normalized negative log likelihood.
#' @export
plot.MCMC.walks <- function(
    wisp.results,
    verbose = TRUE,
    low_samples = 10
  ) {
    
    # Font sizes 
    label_size <- 5.5
    title_size <- 20 
    axis_size <- 12 
    legend_size <- 10
    
    # Grab sampled parameters
    sampled_params <- wisp.results$sample.params.MCMC
    if (length(sampled_params) == 0 || is.null(sampled_params)) return(NULL)
    if (verbose) snk.report...("Making MCMC walks plots")
    sampled_params <- sampled_params[-c(nrow(sampled_params)),]  # Remove last row, which contains L-BFGS values
    
    # Find high-value parameters
    high_val <- 20 < colMeans(sampled_params)
    
    # Down-sample 
    these_rows <- as.integer(seq(1, nrow(sampled_params), length.out = 1000))
    sampled_params_high <- sampled_params[these_rows,high_val]
    sampled_params_low <- sampled_params[these_rows,!high_val]
    if (length(low_samples) > 1) these_cols <- low_samples
    else these_cols <- sample(1:ncol(sampled_params_low), low_samples, replace = FALSE)
    sampled_params_low <- sampled_params_low[,these_cols]
    sampled_params_low <- sampled_params_low[,order(colMeans(abs(sampled_params_low)), decreasing = TRUE)]
    
    # Make long-format data frame
    walks_low <- data.frame(
      value = c(sampled_params_low),
      param = as.factor(rep(1:ncol(sampled_params_low), each = nrow(sampled_params_low))),
      step = rep(these_rows, ncol(sampled_params_low))
    )
    walks_high <- data.frame(
      value = c(sampled_params_high),
      param = as.factor(rep(1:ncol(sampled_params_high), each = nrow(sampled_params_high))),
      step = rep(these_rows, ncol(sampled_params_high))
    )
    
    # Make colors 
    getPalette <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))  # Set1 has max 9 colors
    colors_low <- getPalette(ncol(sampled_params_low))
    colors_high <- getPalette(ncol(sampled_params_high))
    
    # Make plot
    plot.walks.parameters_low <- ggplot(
      data = walks_low, 
      aes(x = step, y = value, group = param, color = param)
    ) + 
      geom_line() + 
      ggtitle("Parameter Walks from MCMC Estimation, Low-Values") +
      scale_color_manual(values = colors_low) +
      theme_minimal() + 
      theme(
        plot.title = element_text(hjust = 0.5, size = title_size),
        axis.title = element_text(size = axis_size),
        axis.text = element_text(size = axis_size),
        legend.title = element_text(size = legend_size),
        legend.text = element_text(size = legend_size),
        legend.position = "none"
      )
    plot.walks.parameters_high <- ggplot(
      data = walks_high, 
      aes(x = step, y = value, group = param, color = param)
    ) + 
      geom_line() + 
      ggtitle("Parameter Walks from MCMC Estimation, High-Values") +
      scale_color_manual(values = colors_high) +
      theme_minimal() + 
      theme(
        plot.title = element_text(hjust = 0.5, size = title_size),
        axis.title = element_text(size = axis_size),
        axis.text = element_text(size = axis_size),
        legend.title = element_text(size = legend_size),
        legend.text = element_text(size = legend_size),
        legend.position = "none"
      )
    
    # Grab trace of negloglik and normalize relative to L-BFGS value
    nll <- unlist(wisp.results[["diagnostics.MCMC"]][["neg.loglik"]])
    pnll <- unlist(wisp.results[["diagnostics.MCMC"]][["pen.neg.value"]])
    nll <- nll / nll[length(nll)]
    pnll <- pnll / pnll[length(pnll)]
    nll <- nll[-c(length(nll))]  # Remove last value, which is the L-BFGS value
    pnll <- pnll[-c(length(pnll))]  # Remove last value, which is the L-BFGS value
    ymin <- min(c(nll, pnll))
    ymax <- max(c(nll, pnll))
    
    # Manually create long-format data frame
    df_long <- data.frame(
      iteration = c(seq_along(nll), seq_along(pnll)),
      value = c(nll, pnll),
      type = c(rep("neg.loglik", length(nll)), rep("penalized", length(pnll)))
    )
    
    # Plot
    plot.walks.nll <- ggplot(df_long, aes(x = iteration, y = value, color = type)) +
      geom_line() +
      scale_color_manual(values = c("neg.loglik" = "blue", "penalized" = "red")) +
      coord_cartesian(ylim = c(ymin, ymax)) +  
      theme_minimal() +
      labs(
        x = "MCMC step number", 
        y = "Normalized Value", 
        title = "Negative log likelihood over MCMC random walk", 
        color = "Type"
      ) + 
      theme(
        plot.title = element_text(hjust = 0.5, size = title_size),
        axis.title = element_text(size = axis_size),
        axis.text = element_text(size = axis_size),
        legend.title = element_text(size = legend_size),
        legend.text = element_text(size = legend_size)
      )
    
    grid.arrange(
      plot.walks.parameters_low, 
      plot.walks.parameters_high, 
      plot.walks.nll,
      ncol = 1)
    
    return(
      list(
        plot.walks.parameters_low = plot.walks.parameters_low,
        plot.walks.parameters_high = plot.walks.parameters_high,
        plot.walks.nll = plot.walks.nll
      )
    )
    
  }

#' Plot individual components of wisp fit for a single child level
#'
#' Extension of \code{plot.ratecount} which plots the individual pieces of the rate-count plot for a single child separately. Returns individual plots for data points only, fit lines only, and data points plus fit lines for individual random levels.
#'
#' @param wisp.results List, output of the wisp function.
#' @param child Character string, the child level to plot. Must be provided, and only one at a time. 
#' @param log Logical, if TRUE, plots on a log scale. Defaults to FALSE.
#' @param dim.boundaries Numeric vector, independent block boundaries to plot for comparison. If empty, the argument is ignored.
#' @param y.lim Numeric vector of length 2, limits for the y-axis of the plots. If NA, defaults to automatic limits.
#' @return List of ggplot objects containing the decomposed rate-count plots.
#' @export
plot.decomposition <- function(
    wisp.results,
    child, # "Nptxr"
    log = FALSE, 
    dim.boundaries = c(), # colMeans(count_data_WSPmm.y$db)
    y.lim = NULL
  ) {
    
    # Comment for future development: Function assumes there is only one parent level. Needs to handle case of multiple parent levels. 
    
    ran_lvls <- wisp.results$grouping.variables$ran.lvls
    ran_lvls <- ran_lvls[ran_lvls != "none"]
    
    ptype <- "pred"
    ctype <- "count"
    if (log) {
      ptype <- paste0(ptype, ".log")
      ctype <- paste0(ctype, ".log")
    }
    
    input_rows <- 5 + length(ran_lvls)
    input <- data.frame(
      count.alpha.none = as.numeric(rep(NA, input_rows)),
      count.alpha.ran = as.numeric(rep(NA, input_rows)),
      pred.alpha.none = as.numeric(rep(NA, input_rows)),
      pred.alpha.ran = as.numeric(rep(NA, input_rows)),
      rans.to.print = as.character(rep(NA, input_rows)),
      stringsAsFactors = FALSE
    )
    rownames(input) <- c("all", "points.obs", "points.extrap", "lines", "none", as.character(ran_lvls))
    
    for (r in 1:input_rows) {
      if (r == 2) {
        input[r,1:4] <- as.numeric(c(0,1,0,0))
      } else if (r == 3) {
        input[r,1:4] <- as.numeric(c(1,0,0,0))
      } else if (r == 4) {
        input[r,1:4] <- as.numeric(c(0,0,1,1))
      } else if (r == 5) {
        input[r,1:4] <- as.numeric(c(1,0,1,0))
        input[r,5] <- "none"
      } else if (r > 5) {
        input[r,1:4] <- as.numeric(c(0,1,0,1))
        input[r,5] <- ran_lvls[r - 5]
      }
    }
    
    plots.decomp <- list()
    length(plots.decomp) <- input_rows
    names(plots.decomp) <- rownames(input)
    for (r in 1:input_rows) {
      
      plots.decomp[[r]] <- plot.ratecount(
        wisp.results = wisp.results,
        pred.type = ptype,
        count.type = ctype,
        dim.boundaries = dim.boundaries,
        y.lim = y.lim,
        count.alpha.none = input$count.alpha.none[r],
        count.alpha.ran = input$count.alpha.ran[r],
        pred.alpha.none = input$pred.alpha.none[r],
        pred.alpha.ran = input$pred.alpha.ran[r],
        rans.to.print = input$rans.to.print[r],
        childs.to.print = c(child)
      )
      
    }
    
    return(plots.decomp)
    
  }

#' Visually compare normality and autocorrelation of bootstrap and MCMC parameter estimates
#'
#' Function allowing for visual comparison of the parameter estimates from bootstrapping vs MCMC simulation. Shows density for ten representative parameters from bootstrapping and MCMC walk, the distributions of Shapiro-Walk p-values for bootstrapping vs MCMC walks, and the density of autocorrelation results for bootstrapping vs MCMC walks.
#'
#' @param wisp.results List, output of the wisp function.
#' @param verbose Logical, if TRUE, prints information during plotting.
#' @return List of ggplot objects containing plots of representative parameter distributions, Shapiro-Wilk normality test results, and autocorrelation plots for bootstrap and MCMC parameter estimates.
#' @export
plot.MCMC.bs.comparison <- function(
    wisp.results,
    verbose = TRUE
  ) {
    
    # For each parameter in the model, both bootstrapping (bs) and the MCMC simulation will 
    # produce a distribution of parameter estimates. We expect that distribution to be normal,
    # which suggests that the estimates are stable. For each parameter and each method (bs and MCMC), 
    # we can quantify "how normal" the estimate distribution is by computing the Shapiro-Wilk normality 
    # p-value on small resamples of the estimate. Higher p-values mean more like a normal distribution, 
    # lower p-values mean less like a normal distribution. We get a single p-value for each parameter 
    # by taking the mean of the resamples, then, for each of bs and MCMC, we'll have a distribution of p-values, 
    # one p-value per parameter. We can then compare the distributions of p-values for bs and MCMC by 
    # looking at their density plots to see if one systematically produces more normal distributions than the other.
    
    # Ensure reproducibility (sampling via Shaprio test not very stable)
    set.seed(1234)
    
    # Font sizes 
    label_size <- 5.5
    title_size <- 20 
    axis_size <- 12 
    legend_size <- 10
    
    # Grab data
    param_mcmc <- wisp.results[["sample.params.MCMC"]]
    param_bs <- wisp.results[["sample.params.bs"]]
    if (is.null(param_bs) || is.null(param_mcmc) || length(param_bs) == 0 || length(param_mcmc) == 0) {
      return(NULL)
    }
    
    if (verbose) snk.report...("Making MCMC vs bootstrap parameter distribution comparison plots")
    
    # Find mean "normality" (p-value of shaprio test) for each parameter, on each method
    dens_list <- lapply(1:ncol(param_mcmc), function(pnum) {
      mcmc_vals <- param_mcmc[, pnum]
      bs_vals   <- param_bs[, pnum]
      n_resamples <- 100
      mcmc_sw_stat <- rep(NA, n_resamples)
      bs_sw_stat   <- rep(NA, n_resamples)
      
      sample_size <- 35
      throw_bs_warning <- FALSE
      throw_MCMC_warning <- FALSE
      for (i in 1:n_resamples) {
        mcmc_vals_sample <- sample(mcmc_vals, sample_size, replace = TRUE)
        bs_vals_sample   <- sample(bs_vals, sample_size, replace = TRUE)
        if (length(unique(mcmc_vals_sample)) > 1) {
          mcmc_sw_stat[i] <- shapiro.test(mcmc_vals_sample)$p.value
        } else {
          mcmc_sw_stat[i] <- NA
          throw_MCMC_warning <- TRUE
        }
        if (length(unique(bs_vals_sample)) > 1) {
          bs_sw_stat[i]   <- shapiro.test(bs_vals_sample)$p.value
        } else {
          bs_sw_stat[i] <- NA
          throw_bs_warning <- TRUE
        }
      }
      if (throw_bs_warning) {
        warning("Some parameters had only one unique bootstrap resample value, Shapiro test not computed for them.")
      }
      if (throw_MCMC_warning) {
        warning("Some parameters had only one unique MCMC resample value, Shapiro test not computed for them.")
      }
      
      return(
        data.frame(
          sw_stat = c(mean(mcmc_sw_stat), mean(bs_sw_stat)),
          method = as.factor(c("MCMC", "Bootstrap"))
        )
      )
    })
    
    df_dens <- do.call(rbind, dens_list)
    
    # Plot
    plot_comparison_Shaprio <- ggplot(df_dens, aes(x = sw_stat, color = method)) +
      geom_density(linewidth = 1.2) +
      theme_minimal() +
      labs(
        title = "Shaprio-Wilk Normality Test on Resampled Parameters", 
        x = "Shaprio-Wilk Normality Test p-value", 
        y = "Density"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, size = title_size),
        axis.title = element_text(size = axis_size),
        axis.text = element_text(size = axis_size),
        legend.title = element_text(size = legend_size),
        legend.text = element_text(size = legend_size)
      )
    
    # For each method, order parameters by "normality" (p-value of shaprio test)
    mcmc_pvals <- rep(NA, ncol(param_mcmc))
    bs_pvals   <- rep(NA, ncol(param_bs))
    for (i in 1:ncol(param_mcmc)) {
      mcmc_vals <- param_mcmc[, i]
      bs_vals   <- param_bs[, i]
      n_resamples <- 100
      mcmc_sw_stat <- rep(NA, n_resamples)
      bs_sw_stat   <- rep(NA, n_resamples)
      
      sample_size <- 35
      for (i in 1:n_resamples) {
        mcmc_vals_sample <- sample(mcmc_vals, sample_size, replace = TRUE)
        bs_vals_sample   <- sample(bs_vals, sample_size, replace = TRUE)
        mcmc_sw_stat[i] <- shapiro.test(mcmc_vals_sample)$p.value
        bs_sw_stat[i]   <- shapiro.test(bs_vals_sample)$p.value
      }
      mcmc_pvals[i] <- mean(mcmc_sw_stat)
      bs_pvals[i] <- mean(bs_sw_stat)
    }
    mcmc_pvals_order <- order(mcmc_pvals)
    bs_pvals_order <- order(bs_pvals)
    
    # Take 10 equally spaced parameteres from the ordered list
    take <- 10
    mcmc_pvals_order <- mcmc_pvals_order[as.integer(seq(1,length(mcmc_pvals_order), length.out = take))]
    bs_pvals_order <- bs_pvals_order[as.integer(seq(1,length(bs_pvals_order), length.out = take))]
    
    # Precompute density estimates for each param and group
    dens_list <- lapply(1:take, function(pnum) {
      mcmc_vals <- param_mcmc[, mcmc_pvals_order[pnum]]
      bs_vals   <- param_bs[, bs_pvals_order[pnum]]
      
      # Density and centering
      d_mcmc <- density(mcmc_vals)
      d_bs   <- density(bs_vals)
      mcmc_peak <- d_mcmc$x[which.max(d_mcmc$y)]
      bs_peak   <- d_bs$x[which.max(d_bs$y)]
      
      # Return centered density curves
      data.frame(
        x = c(d_mcmc$x - mcmc_peak, d_bs$x - bs_peak),
        y = c(d_mcmc$y, d_bs$y),
        method = rep(c("MCMC", "Bootstrap"), each = length(d_mcmc$x)),
        param = factor(pnum)
      )
    })
    
    df_dens <- do.call(rbind, dens_list)
    df_dens$y <- log(df_dens$y+1)
    
    # Plot precomputed densities
    plot_comparison_density <- ggplot(df_dens, aes(x = x, y = y, color = method, method = interaction(method, param))) +
      geom_line(linewidth = 1) +
      theme_minimal() +
      labs(
        title = "Resampled Parameter Distribution Comparison",
        x = "Parameter Value Centered with Peak Density at Zero", y = "Log Density + 1") +
      theme(
        plot.title = element_text(hjust = 0.5, size = title_size),
        axis.title = element_text(size = axis_size),
        axis.text = element_text(size = axis_size),
        legend.title = element_text(size = legend_size),
        legend.text = element_text(size = legend_size)
      )
    
    # Grab sample results
    sample_results_MCMC <- wisp.results$sample.params.MCMC
    sample_results_bs <- wisp.results$sample.params.bs
    
    # Compute correlations
    n_params <- ncol(sample_results_MCMC)
    if (n_params != ncol(sample_results_bs)) {
      stop("Sample results from MCMC and bootstrap have different number of parameters")
    }
    # ... for MCMC random walk steps
    cor_MCMC <- rep(NA, n_params)
    n_samples_MCMC <- nrow(sample_results_MCMC)
    for (i in seq_len(n_params)) {
      cor_MCMC[i] <- cor(
        x = sample_results_MCMC[c(1:(n_samples_MCMC-1)),i], 
        y = sample_results_MCMC[c(2:n_samples_MCMC),i],
        method = "pearson"
      )
    }
    # ... for bootstraps (control, expect no correlation) 
    cor_bs <- rep(NA, n_params)
    n_samples_bs <- nrow(sample_results_bs)
    for (i in seq_len(n_params)) {
      cor_bs[i] <- cor(
        x = sample_results_bs[c(1:(n_samples_bs-1)),i], 
        y = sample_results_bs[c(2:n_samples_bs),i],
        method = "pearson"
      )
    }
    
    # Make summary table
    den_MCMC <- density(cor_MCMC)
    den_bs <- density(cor_bs)
    sample_correlations <- data.frame(
      density = log(c(den_MCMC$y, den_bs$y) + 1),
      auto_cor = c(den_MCMC$x, den_bs$x),
      method = c(rep("MCMC", length(den_MCMC$x)), rep("Bootstrap", length(den_bs$x)))
    )
    
    # Make summary density plot 
    plot_sample_correlations <- ggplot(sample_correlations, aes(x = auto_cor, y = density, color = method)) +
      geom_line(linewidth = 1.2) +
      theme_minimal() +
      labs(
        title = "Next-step Autocorrelation Between Parameter Resamples", 
        x = "Autocorrelation", 
        y = "Log Density + 1"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, size = title_size),
        axis.title = element_text(size = axis_size),
        axis.text = element_text(size = axis_size),
        legend.title = element_text(size = legend_size),
        legend.text = element_text(size = legend_size)
      )
    
    grid.arrange(
      plot_comparison_density,
      plot_comparison_Shaprio,
      plot_sample_correlations,
      ncol = 1
    )
    
    return(
      list(
        plot_sample_correlations = plot_sample_correlations,
        plot_comparison_Shaprio = plot_comparison_Shaprio,
        plot_comparison_density = plot_comparison_density
      )
    )
    
  }

# Debugging and misc ###################################################################################################

# Functions used (in Cpp) in the initial change-point estimation
find_centroid <- function(
    X_lik,     # array with values as lik ratios, rows as bins, columns as trt x ran interactions
    window     # window size for dtwclust::DBA, constraint on slope
  ) {
    
    # Uses R package from: https://doi.org/10.1016/j.patcog.2010.09.013
    # R package: 
      # citHeader("To cite the R package, use:")
      # citation(auto = meta)
      # 
      # bibentry(
      #   "Article",
      #   header = "If the vignette's content was useful, consider citing the summarized version published in:",
      #   title = "Time-Series Clustering in R Using the dtwclust Package",
      #   author = person("Alexis", "Sardá-Espinosa"),
      #   year = "2019",
      #   journal = "The R Journal",
      #   doi = "10.32614/RJ-2019-023"
      # )
    
    # Take transpose, as dtwclust::DBA expects columns as time points
    X <- t(X_lik)
    # Find centroid of X, i.e., the "average" series
    d <- dtwclust::DBA(X, window.size = window) 
    return(d)
  }

project_cp <- function(
    found_cp, # vector of change points 
    centroid, # centroid from which those cp were estimated
    X_lik     # data from which that centroid was computed
  ) {
    # Function will return the implied change points in the original data
    # ... rows are change points (by deg), columns as trt x ran interactions
    
    # Uses R package dtw, with step pattern derived from: 
    # H. Sakoe and S. Chiba, "Dynamic programming algorithm optimization for spoken word recognition," 
    #   in IEEE Transactions on Acoustics, Speech, and Signal Processing, vol. 26, no. 1, pp. 43-49, February 1978, doi: 10.1109/TASSP.1978.1163055. 
    #   This step pattern chosen to enforce a very shallow slope. 
    
    # Get alignment indices
    # ... take transpose, as dtwclust::DBA expects columns as time points
    X <- t(X_lik)
    # ... initialize list to store alignment indices
    aX <- list()
    # ... for each row of X (i.e., each trt x ran interaction), find the indices to align it to the centroid
    for (i in 1:nrow(X)) {
      Xi <- X[i,]
      Xi[is.na(Xi)] <- 1
      alignment_idx <- dtw::dtw(Xi, centroid, step.pattern = dtw::symmetricP05, distance.only = FALSE)
      aX[[i]] <- cbind(alignment_idx$index1, alignment_idx$index2)
    }
    
    # Project alignment indices to original data
    X_likp <- array(NA, dim = c(length(found_cp), ncol(X_lik)))
    for (c in 1:ncol(X_lik)) {
      
      a_idx <- aX[[c]]
      
      for (cpi in 1:length(found_cp)) {
        cp <- found_cp[cpi]                # index of change point on the centroid
        cp_a_idx <- which(a_idx[,2] == cp) # second column is the centroid alignment indices
        X_likp[cpi, c] <- round(mean(a_idx[cp_a_idx, 1])) # project back, take mean, and round to integer
      }
      
    }
    
    return(X_likp)
    
  }

#' Wisp warping function, implemented in R
#' 
#' This function provides an R implementation of the warping function used to warp model components before input into \code{wisp.sigmoid}.
#' 
#' @param z Numeric vector or matrix, values to warp. Scalar values fine as well. Matrix must be 2D.
#' @param b Numeric, warping bound (must be a single value).
#' @param w Numeric, warping factor (must be a single value).
#' @return Numeric vector or matrix, warped values. Returned object will have the same dimensions as input z.
wisp.warp <- function( 
    z, # value to warp, either a scalar, vector, or 2D array/matrix
    b, # warping bound
    w  # warping factor
  ) {
    
    if (length(b) != 1) stop("b must be a single value for warping bound")
    if (length(w) != 1) stop("w must be a single value for warping factor")
    if (!(length(z) > 0)) stop("z must be a vector or matrix of values to warp, but has length zero")
    if (max(z) > b) stop("z must be less than or equal to b")
    
    if (is.null(dim(z)) || length(dim(z)) == 1) {
      out <- rep(NA, length(z))
      for (i in 1:length(z)) {
        out[i] <- warp_mc_R(z[i], b, w)
      }
    } else if (length(dim(z)) == 2) {
      out <- array(NA, dim = dim(z))
      for (i in 1:(dim(z)[1])) {
        for (j in 1:(dim(z)[2])) {
          out[i,j] <- warp_mc_R(z[i,j], b, w)
        }
      }
    } else {
      stop("z must be a vector or matrix of values to warp, but unrecognized dim attribute detected")
    }
    return(out)
    
  }

#' Wisp sigmoid function, implemented in R
#' 
#' This function provides an R implementation of the wisp sigmoid function. It assumes the parameters Rt, tslope, and tpoint have already been warped by \code{wisp.warp}.
#' 
#' @param x Numeric, input spatial position, either a scalar, vector, or 2D array/matrix
#' @param Rt Numeric vector, rate parameters for the wisp function. Degree of the wisp model will be length of this vector minus 1.
#' @param tslope Numeric vector, slope scalars for the wisp function. Must be one less than the length of Rt.
#' @param tpoint Numeric vector, transition points for the wisp function. Must be one less than the length of Rt.
#' @return Numeric vector or matrix, wisp sigmoid values for given input. Returned object will have the same dimensions as input x.
wisp.sigmoid <- function(
    x, 
    Rt,
    tslope,
    tpoint
  ) {
    
    # Check inputs
    if (length(Rt) != length(tslope) + 1) {stop("Rt must be one element longer than tslope")}
    if (length(tslope) != length(tpoint)) {stop("tslope must be the same length as tpoint")}
    deg <- length(Rt) - 1
    
    if (is.null(dim(x)) || length(dim(x)) == 1) {
      out <- rep(NA, length(x))
      for (i in 1:length(x)) {
        out[i] <- poly_sigmoid_R(x[i], deg, Rt, tslope, tpoint)
      }
    } else if (length(dim(x)) == 2) {
      out <- array(NA, dim = dim(x))
      for (i in 1:(dim(x)[1])) {
        for (j in 1:(dim(x)[2])) {
          out[i,j] <- poly_sigmoid_R(x[i,j], deg, Rt, tslope, tpoint)
        }
      }
    } else {
      stop("x must be a vector or matrix of values to warp, but unrecognized dim attribute detected")
    }
    return(out)
    
  }

#' Make plot demonstrating how the wisp function is warped by the warping factors
#' 
#' This function takes user-provided values for wisp function parameters and produces a panel of plots showing how the wisp model is warped by the warping factors
#' 
#' @param w Numeric, warping factor (must be a single value).
#' @param point_pos Numeric, x coordinate at which to place positive warp segment (must be a single value).
#' @param point_neg Numeric, x coordinate at which to place negative warp segment (must be a single value).
#' @param Rt Numeric vector, rate parameters for the wisp function. Degree of the wisp model will be length of this vector minus 1.
#' @param tslope Numeric vector, slope scalars for the wisp function. Must be one less than the length of Rt.
#' @param tpoint Numeric vector, transition points for the wisp function. Must be one less than the length of Rt.
#' @param w_factors Numeric vector, warping factors for the wisp function. Must be length 3. First element is the warping factor for Rt, second for tslope, and third for tpoint. Note that this is more restrictive than the real wisp model, which not only allows for different warping factors across the model components (Rt, tslope, and tpoint), but also across the different elements within each model component as well.
#' @return Nothing. A ggplot object is printed to console.
#' @export
demo.warp.plots <- function(
    w = 2,                       # warping factor
    point_pos = 60,
    point_neg = 40,
    Rt = c(6, 3, 0.2, 6)*4.65,   # rates for poly-sigmoid
    tslope = c(0.4, 0.75, 1),    # slope scalars for poly-sigmoid
    tpoint = c(15, 38, 80),      # transition points for poly-sigmoid
    w_factors = c(0.6, -0.9, 0.5)      # warping factors for poly sigmoid
  ) {
    
    # Font sizes 
    label_size <- 5.5
    title_size <- 20 
    axis_size <- 12 
    legend_size <- 10
    
    # Data
    n <- 1000
    x <- (1:n)/10
    b <- 100
    y <- x
    y1 <- wisp.warp(x, b, w)
    y2 <- wisp.warp(x, b, -w)
    y_pos <- wisp.warp(point_pos, b, w)
    y_neg <- wisp.warp(point_neg, b, -w)
    
    # Organize into a data frame
    df <- data.frame(
      x = rep(x, 3),
      y = c(y, y1, y2),
      curve = factor(rep(c("w = 0", "w > 0", "w < 0"), each = length(x)))
    )
    df$curve <- relevel(df$curve, ref = "w = 0")
    
    df_segments <- data.frame(
      point_pos = point_pos,
      point_neg = point_neg,
      y_pos = y_pos,
      y_neg = y_neg
    )
    
    # Make the ggplot
    demo_plot_warpfunction <- ggplot(df, aes(x = x, y = y, color = curve)) +
      geom_line(linewidth = 1.25) +
      geom_hline(yintercept = 100, linetype = "dashed", color = "darkgray", linewidth = 1) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray", linewidth = 1) +
      #coord_fixed(ratio = 1) +
      geom_segment(
        data = df_segments,
        aes(x = point_pos, xend = point_pos, y = point_pos, yend = y_pos),
        color = "cyan4", linetype = "dashed", linewidth = 0.75) + 
      geom_segment(
        data = df_segments,
        aes(x = point_neg, xend = point_neg, y = point_neg, yend = y_neg),
        color = "deeppink4", linetype = "dashed", linewidth = 0.75) +
      annotate("text", x = 10, y = 95, label = "upper asymptote", size = label_size, color = "black") +
      annotate("text", x = 90, y = 5, label = "lower asymptote", size = label_size, color = "black") +
      annotate("text", x = point_pos - 10, y = (y_pos + point_pos)/2, label = expression(varphi * "(z)(b - z)"), size = label_size, color = "black") +
      annotate("text", x = point_neg + 10, y = (y_neg + point_neg)/2, label = expression(varphi * "(b - z)z"), size = label_size, color = "black") +
      labs(
        x = "z",
        y = expression(omega * "(z, " * rho * ", b)"),
        title = "WSP Warping Function",
        color = "Direction"
      ) +
      scale_color_manual(
        values = c("black", "deeppink2", "cyan3"),
        labels = c(expression(rho * " = 0"), expression(rho * " < 0"), expression(rho * " > 0"))
      ) +
      theme_minimal(base_size = 16) +
      theme(
        plot.title = element_text(hjust = 0.5, size = title_size),
        axis.title = element_text(size = axis_size),
        axis.text = element_text(size = axis_size),
        legend.title = element_text(size = legend_size),
        legend.text = element_text(size = legend_size)
      )
    
    # Construct poly-sigmoid curve, no warping
    x <- seq(0, 100, length.out = n)
    y <- rep(NA, n)
    deg <- length(Rt)-1
    for (i in 1:n) y[i] <- poly_sigmoid_R(
      x[i], # input 
      deg,    # degree 
      Rt,     # Rates
      tslope, # slope scalars
      tpoint  # inflection points
    )
    
    # Construct poly-sigmoid curve, with warping
    yw <- rep(NA, n)
    Rtw <- wisp.warp(Rt, 100000, w_factors[1])
    tslopew <- wisp.warp(tslope, 100000, w_factors[2])
    tpointw <- wisp.warp(tpoint, 100, w_factors[3])
    for (i in 1:n) yw[i] <- poly_sigmoid_R(
      x[i],    # input 
      deg,     # degree 
      Rtw,     # Rates
      tslopew, # slope scalars
      tpointw  # inflection points
    )
    df <- data.frame(x = c(x, x), y = c(y, yw), type = c(rep("unwarped", n), rep("warped", n)))
    
    # Compute slopes 
    m <- tslope 
    for (i in 1:length(tslope)) {
      m[i] <- m[i] * (Rt[i+1] - Rt[i]) / 4
    }
    mw <- tslopew
    for (i in 1:length(tslopew)) {
      mw[i] <- mw[i] * (Rtw[i+1] - Rtw[i]) / 4
    }
    
    # Make line segments to show slopes
    y0 <- tpoint 
    y0w <- tpointw
    for (i in 1:length(y0)) {
      y0[i] <- poly_sigmoid_R(
        y0[i],  # input 
        deg,    # degree 
        Rt,     # Rates
        tslope, # slope scalars
        tpoint  # inflection points
      )
    }
    for (i in 1:length(y0w)) {
      y0w[i] <- poly_sigmoid_R(
        y0w[i],  # input 
        deg,     # degree 
        Rtw,     # Rates
        tslopew, # slope scalars
        tpointw  # inflection points
      )
    }
    slope_run <- ((max(y)-min(y))/abs(diff(Rt))) * 3 + c(6, 0, 0.5)
    slope_runw <- ((max(yw)-min(yw))/abs(diff(Rtw))) * 3 + c(6, 0, 0.5)
    slope_seg_neg <- y0 - slope_run*m 
    slope_seg_negw <- y0w - slope_runw*mw
    slope_seg_pos <- y0 + slope_run*m 
    slope_seg_posw <- y0w + slope_runw*mw
    
    # Data frame to hold slope segments for plotting
    def_segments_slopes <- data.frame(
      slope_seg_neg_x = tpoint - slope_run,
      slope_seg_pos_x = tpoint + slope_run,
      slope_seg_neg_y = slope_seg_neg,
      slope_seg_pos_y = slope_seg_pos
    )
    def_segments_slopesw <- data.frame(
      slope_seg_neg_x = tpointw - slope_runw,
      slope_seg_pos_x = tpointw + slope_runw,
      slope_seg_neg_y = slope_seg_negw,
      slope_seg_pos_y = slope_seg_posw
    )
    
    # Data frame to hold rate segments for plotting
    def_segments_rates <- data.frame(
      x0 = c(0, tpoint)-max(x)*0.1,
      x1 = c(tpoint, max(x))+max(x)*0.1,
      y0 = Rt,
      y1 = Rt
    )
    def_segments_ratesw <- data.frame(
      x0 = c(0, tpointw)-max(x)*0.1,
      x1 = c(tpointw, max(x))+max(x)*0.1,
      y0 = Rtw,
      y1 = Rtw
    )
    
    # Compute y limits
    ylower <- min(
      min(c(
        def_segments_slopes$slope_seg_neg_y, def_segments_slopes$slope_seg_pos_y,
        def_segments_slopesw$slope_seg_neg_y, def_segments_slopesw$slope_seg_pos_y)),
      -mean(y)*0.05
    ) 
    yupper <- max(
      max(c(
        def_segments_slopes$slope_seg_neg_y, def_segments_slopes$slope_seg_pos_y,
        def_segments_slopesw$slope_seg_neg_y, def_segments_slopesw$slope_seg_pos_y)),
      max(y)*1.2
    )
    
    # Make plot
    demo_plot_warpedsigmoid <- ggplot(df) +
      geom_line(aes(x = x, y = y, color = type), linewidth = 1.5) +  
      geom_vline(xintercept = tpoint, linetype = "dashed", color = "red4", linewidth = 1) +
      geom_vline(xintercept = tpointw, linetype = "dashed", color = "red", linewidth = 1) +
      ylim(ylower, yupper) +
      coord_fixed()  +
      geom_segment(
        data = def_segments_rates,
        aes(x = x0, xend = x1, y = y0, yend = y1),
        color = "darkgray", linetype = "dashed", linewidth = 1) +
      geom_segment(
        data = def_segments_ratesw,
        aes(x = x0, xend = x1, y = y0, yend = y1),
        color = "gray", linetype = "dashed", linewidth = 1) +
      geom_segment(
        data = def_segments_slopes,
        aes(x = slope_seg_neg_x, xend = slope_seg_pos_x, y = slope_seg_neg_y, yend = slope_seg_pos_y),
        color = "blue4", linetype = "dashed", linewidth = 0.75) + 
      geom_segment(
        data = def_segments_slopesw,
        aes(x = slope_seg_neg_x, xend = slope_seg_pos_x, y = slope_seg_neg_y, yend = slope_seg_pos_y),
        color = "blue", linetype = "dashed", linewidth = 0.75) +
      labs(
        x = "x",
        y = expression(Psi * "(x, r, s, p)"),
        title = "WSP Sigmoid Function, Warped",
        color = "Type"
      )  +
      scale_color_manual(
        values = c("black", "orange3")
      ) +
      theme_minimal(base_size = 16) +
      theme(
        plot.title = element_text(hjust = 0.5, size = title_size),
        axis.title = element_text(size = axis_size),
        axis.text = element_text(size = axis_size),
        legend.title = element_text(size = legend_size),
        legend.text = element_text(size = legend_size)
      ) 
    
    grid.arrange(demo_plot_warpfunction, demo_plot_warpedsigmoid, ncol = 1)
    # Save at 1156 x 843
    return(NULL)
    
  }

#' Make plot demonstrating how the wisp function is determined by the Rt, slope, and tpoint parameters
#' 
#' This function takes user-provided values for wisp function parameters and produces a panel of plots showing how the wisp model is determined by the Rt, slope, and tpoint parameters.
#' 
#' @param r Numeric, upper asymptote for logistic (must be a single value).
#' @param s Numeric, slope scalar at inflection point (must be a single value).
#' @param Rt Numeric vector, rate parameters for the wisp function. Degree of the wisp model will be length of this vector minus 1.
#' @param tslope Numeric vector, slope scalars for the wisp function. Must be one less than the length of Rt.
#' @param tpoint Numeric vector, transition points for the wisp function. Must be one less than the length of Rt.
#' @return Nothing. A ggplot object is printed to console.
#' @export
demo.sigmoid.plots <- function(
    r = 4,                       # upper asymptote for logistic
    s = 1,                       # slope scalar at inflection point
    Rt = c(6, 3, 0.2, 6)*4.65,   # rates for poly-sigmoid
    tslope = c(0.4, 0.75, 1),    # slope scalars for poly-sigmoid
    tpoint = c(15, 38, 80)       # transition points for poly-sigmoid
  ) {
    
    # Font sizes 
    label_size <- 5.5
    title_size <- 20 
    axis_size <- 12 
    legend_size <- 10
    
    # Check inputs
    if (length(Rt) != length(tslope) + 1) {
      stop("Rt must be one element longer than tslope")
    }
    if (length(tslope) != length(tpoint)) {
      stop("tslope must be the same length as tpoint")
    }
    if (length(r) != 1) {
      stop("r must be a single value for upper asymptote")
    }
    if (length(s) != 1) {
      stop("s must be a single value for slope scalar at inflection point")
    }
    
    # Plot 1, logistic function ####
    
    # Construct logistic curve
    n <- 1000
    x <- seq(-10, 10, length.out = n)
    y <- rep(NA, n)
    for (i in 1:n) y[i] <- sigmoid_stable_R(-x[i]*s)*r
    df <- data.frame(x = x, y = y)
    
    # Compute slope at inflection point (y0)
    y0 <- sigmoid_stable_R(-0*s)*r
    m <- s*(y0*r - y0*y0)/r
    m <- -m # ... flip because the x axis was flipped
    
    # Make line segment to show slope at inflection point
    slope_run <- r/s
    slope_seg_neg <- y0 - slope_run*m
    slope_seg_pos <- y0 + slope_run*m
    
    # Data frame to hold segment for plotting
    df_segments <- data.frame(
      slope_seg_neg_x = -slope_run,
      slope_seg_pos_x = slope_run,
      slope_seg_neg_y = slope_seg_neg,
      slope_seg_pos_y = slope_seg_pos
    )
    
    # Plot logistic function 
    demo_plot_logistic <- ggplot(df, aes(x = x, y = y)) +
      geom_line(linewidth = 1.5)  +
      coord_fixed() +
      geom_hline(yintercept = r, linetype = "dashed", color = "darkgray", linewidth = 1) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray", linewidth = 1) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "red4", linewidth = 1) +
      geom_segment(
        data = df_segments,
        aes(x = slope_seg_neg_x, xend = slope_seg_pos_x, y = slope_seg_neg_y, yend = slope_seg_pos_y),
        color = "blue4", linetype = "dashed", linewidth = 0.75) + 
      annotate("text", x = -7, y = r+0.75, label = "upper asymptote", size = label_size, color = "black") +
      annotate("text", x = 7, y = -0.75, label = "lower asymptote", size = label_size, color = "black") +
      annotate("text", x = 3, y = r+0.75, label = "inflection point", size = label_size, color = "red4") +
      annotate("text", x = -2, y = r+1, angle = atan(m)*57.3, label = "slope", size = label_size, color = "blue4") +
      labs(
        x = "x",
        y = expression(psi * "(x, r, s)"),
        title = "The Logistic Function"
      )  +
      theme_minimal(base_size = 16) +
      theme(
        plot.title = element_text(hjust = 0.5, size = title_size),
        axis.title = element_text(size = axis_size),
        axis.text = element_text(size = axis_size),
        legend.title = element_text(size = legend_size),
        legend.text = element_text(size = legend_size)
      )
    # saved at 1145 x 647
    
    # Plot 2, poly-sigmoid function ####
    
    # Construct poly-sigmoid curve 
    x <- seq(0, 100, length.out = n)
    y <- rep(NA, n)
    deg <- length(Rt)-1
    for (i in 1:n) y[i] <- poly_sigmoid_R(
      x[i], # input 
      deg,    # degree 
      Rt,     # Rates
      tslope, # slope scalars
      tpoint  # inflection points
    )
    df <- data.frame(x = x, y = y)
    
    # Compute slopes 
    m <- tslope 
    for (i in 1:length(tslope)) {
      m[i] <- m[i] * (Rt[i+1] - Rt[i]) / 4
    }
    
    # Make line segments to show slopes
    y0 <- tpoint 
    for (i in 1:length(y0)) {
      y0[i] <- poly_sigmoid_R(
        y0[i],  # input 
        deg,    # degree 
        Rt,     # Rates
        tslope, # slope scalars
        tpoint  # inflection points
      )
    }
    slope_run <- ((max(y)-min(y))/abs(diff(Rt))) * 3 + c(6, 0, 0.5)
    slope_seg_neg <- y0 - slope_run*m 
    slope_seg_pos <- y0 + slope_run*m 
    
    # Data frame to hold slope segments for plotting
    def_segments_slopes <- data.frame(
      slope_seg_neg_x = tpoint - slope_run,
      slope_seg_pos_x = tpoint + slope_run,
      slope_seg_neg_y = slope_seg_neg,
      slope_seg_pos_y = slope_seg_pos
    )
    
    # Data frame to hold rate segments for plotting
    def_segments_rates <- data.frame(
      x0 = c(0, tpoint)-max(x)*0.1,
      x1 = c(tpoint, max(x))+max(x)*0.1,
      y0 = Rt,
      y1 = Rt
    )
    
    # Compute y limits
    ylower <- min(
      min(c(def_segments_slopes$slope_seg_neg_y, def_segments_slopes$slope_seg_pos_y)),
      -mean(y)*0.05
    ) 
    yupper <- max(
      max(c(def_segments_slopes$slope_seg_neg_y, def_segments_slopes$slope_seg_pos_y)),
      max(y)*1.2
    )
    
    # Compute block midpoints 
    block_midpoints <- c(0, tpoint) 
    block_midpoints <- block_midpoints + diff(c(block_midpoints, max(x)))/2
    
    # Make plot
    demo_plot_sigmoid <- ggplot(df, aes(x = x, y = y)) +
      geom_line(linewidth = 1.5) +  
      geom_vline(xintercept = tpoint, linetype = "dashed", color = "red4", linewidth = 1) +
      ylim(ylower, yupper) +
      coord_fixed() +
      geom_segment(
        data = def_segments_rates,
        aes(x = x0, xend = x1, y = y0, yend = y1),
        color = "darkgray", linetype = "dashed", linewidth = 1) +
      geom_segment(
        data = def_segments_slopes,
        aes(x = slope_seg_neg_x, xend = slope_seg_pos_x, y = slope_seg_neg_y, yend = slope_seg_pos_y),
        color = "blue4", linetype = "dashed", linewidth = 0.75) + 
      annotate("text", x = block_midpoints, y = Rt+mean(Rt)*0.15, label = paste("rate", 1:length(Rt)), size = label_size, color = "black") +
      annotate("text", x = tpoint-max(x)*0.08, y = -max(y)*0.1, label = paste("t-point", 1:length(tpoint)), size = label_size, color = "red4") +
      annotate("text", x = tpoint-max(x)*0.035, y = y0*0.95, angle = atan(m)*57.3, label = paste("slope", 1:length(tslope)), size = label_size, color = "blue4") +
      labs(
        x = "x",
        y = expression(Psi * "(x, r, s, p)"),
        title = "The WSP Sigmoid Function"
      )  +
      theme_minimal(base_size = 16) +
      theme(
        plot.title = element_text(hjust = 0.5, size = title_size),
        axis.title = element_text(size = axis_size),
        axis.text = element_text(size = axis_size),
        legend.title = element_text(size = legend_size),
        legend.text = element_text(size = legend_size)
      )
    
    grid.arrange(demo_plot_logistic, demo_plot_sigmoid, ncol = 1)
    # Saved at 1186 x 1032
    return(NULL)
    
  }
