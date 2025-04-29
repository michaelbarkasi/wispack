
.onLoad <- function(libname, pkgname) {
  Rcpp::loadModule("wspc", TRUE)
  loadNamespace("colorspace")
}

# Main function for generating WSPmm model #############################################################################

# Function for fitting model to raw count data list 
wisp <- function(
    # Data to model
    count.data.raw, 
    # Variable labels
    variables = list( # names of columns in count.data.raw giving each model variable type
      count = "count",
      bin = "bin", 
      parent = "parent", 
      child = "child",
      ran = "ran",
      fixedeffects = c()
    ),
    # Settings used on R side
    use.median = FALSE,
    MCMC.burnin = 0,
    MCMC.steps = 1e4,
    MCMC.step.size = 0.05,
    MCMC.prior = 10.0,
    bootstraps.num = 0, 
    converged.resamples.only = TRUE,
    max.fork = 10,
    null.rate = log(2),
    null.slope = 1,
    dim.bounds = NULL, 
    verbose = TRUE,
    print.child.summaries = TRUE,
    # Setting to pass to C++ model
    model.settings = list(
      buffer_factor = 0.05,                       # buffer factor for penalizing distance from structural parameter values
      ctol = 1e-6,                                # convergence tolerance
      max_penalty_at_distance_factor = 0.01,      # maximum penalty at distance from structural parameter values
      LROcutoff = 2.0,                            # cutoff for LROcp
      LROwindow_factor = 2.0,                     # controls size of window used in LROcp algorithm (window = LROwindow_factor * bin_num * buffer_factor)
      rise_threshold_factor = 0.8,                # amount of detected rise as fraction of total required to end run
      max_evals = 1000,                           # maximum number of evaluations for optimization
      rng_seed = 42,                              # seed for random number generator
      warp_precision = 1e-7                       # decimal precision to retain when selecting really big number as pseudo infinity for unbound warping
    )
  ) {
    
    # Make reproducible
    ran.seed <- 1234
    set.seed(ran.seed)
   
    # Relabel and rearrange data columns 
    old_names <- colnames(count.data.raw)
    required_cols <- c("count", "bin", "parent", "child", "ran")
    ordered_cols <- unlist(variables[required_cols])
    if (length(variables$fixedeffects) == 0) {
      fe_cols <- !(old_names %in% ordered_cols)
    } else {
      fe_cols <- old_names %in% variables$fixedeffects
    }
    new_names <- c(required_cols, old_names[fe_cols])
    data <- cbind(count.data.raw[,ordered_cols], count.data.raw[,fe_cols])
    colnames(data) <- new_names
    
    # Add inf_warp 
    model.settings$inf_warp <- model.settings$warp_precision / .Machine$double.eps
    
    # Initialize cpp model ####
    if (verbose) {
      snk.report("Initializing Cpp (wspc) model")
      snk.horizontal_rule(reps = snk.simple_break_reps, end_breaks = 2)
    }
    cpp_model <- new(
      wspc, 
      data, 
      model.settings,
      verbose
    )
    
    # Estimate model parameters with MCMC or bootstrapping ####
    if (verbose) {
      snk.report("Estimating model parameters", initial_breaks = 1)
      snk.horizontal_rule(reps = snk.small_break_reps, end_breaks = 0)
    }
    
    # Confirm forking is possible
    if (!(Sys.info()["sysname"] == "Darwin" || Sys.info()["sysname"] == "Linux")) {
      if (bootstraps.num > 0) {
        if (verbose) snk.report...("Forking not available on Windows, cannot bootstrap, switching to MCMC")
        bootstraps.num <- 0 
        MCMC.steps <- 1e4
      }
    } else if (bootstraps.num > 0) {
      if (verbose) {
        snk.report...("Forking available and bootstrap requested.")
        snk.report...("If bootstrapping not desired, set bootstraps.num = 0")
      }
    }
    
    if (bootstraps.num == 0) {
      
      # Run MCMC simulation
      if (verbose) snk.report("Running MCMC stimulations (single-threaded)", end_breaks = 1)
      MCMC_walk <- cpp_model$MCMC(
        MCMC.steps + MCMC.burnin, 
        MCMC.step.size,
        MCMC.prior,
        verbose 
      )
      if (MCMC.burnin > 0) {
        MCMC_walk <- MCMC_walk[-c(2:(2+MCMC.burnin-1)),]
      }
      
    } else {
      
      # Run bootstrap fits in parallel with forking
      if (verbose) snk.report("Running bootstrap fits (with forking)", end_breaks = 1)
      sample_results <- cpp_model$bs_batch(
        bootstraps.num, 
        max.fork,
        verbose
      )
      
    }
    
    # Extract bs results and diagnostics
    if (bootstraps.num > 0) {
      n_params <- ncol(sample_results) - 4
      sample.params <- sample_results[,1:n_params]
      bs.diagnostics <- data.frame( 
        pen.neg.value = sample_results[,n_params + 1],
        neg.loglik = sample_results[,n_params + 2], 
        success.code = sample_results[,n_params + 3],
        num.evals = sample_results[,n_params + 4]
      )
      MCMC.diagnostics <- NULL
    } else {
      n_params <- ncol(MCMC_walk) - 4
      sample.params <- MCMC_walk[,1:n_params]
      MCMC.diagnostics <- data.frame(
        pen.neg.value = MCMC_walk[,n_params + 1],
        neg.loglik = MCMC_walk[,n_params + 2], 
        acceptance.ratio = MCMC_walk[,n_params + 3],
        ctr.num = MCMC_walk[,n_params + 4]
      )
      bs.diagnostics <- NULL
    }
    
    # Set final fitted parameters
    if (use.median) {
      if (verbose) snk.report...("Setting median parameter samples as final parameters", initial_breaks = 1, end_breaks = 1)
      final_parameters <- apply(sample.params, 2, function(x) median(x, na.rm = TRUE))
    } else {
      if (verbose) snk.report...("Setting full-data fit as parameters", initial_breaks = 1, end_breaks = 1)
      final_parameters <- sample.params[nrow(sample.params),]
    }
    cpp_model$set_parameters(
      final_parameters,
      verbose
    )
    
    # Grab model results and add samples
    results <- cpp_model$results()
    results[["sample.params"]] <- sample.params
    results[["bs.diagnostics"]] <- bs.diagnostics
    results[["MCMC.diagnostics"]] <- MCMC.diagnostics
    
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
      null.rate = null.rate,
      null.slope = null.slope,
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
    
    # Check tpoint stability
    results$stats$tps <- check.tpoint.stability(
      wisp.results = results,
      verbose = verbose
    )
    
    # Make plots of results ####
    
    # Plot MCMC walks, both parameters and negloglik
    if (bootstraps.num == 0) {
      plots.MCMC <- plot.MCMC.walks(
        wisp.results = results
      )
    } else {
      plots.MCMC <- NULL
    }
    
    # Plot structural parameter distribution
    plots.effect.dist <- plot.effect.dist(
      wisp.results = results,
      verbose = verbose
    )
    
    # Make rate plots 
    if (verbose) snk.report...("Making rate-count plots")
    plots.ratecount <- plot.ratecount(
      wisp.results = results,
      pred.type = "pred",
      count.type = "count",
      dim.boundaries = dim.bounds
    )
    
    # Make parameter plots 
    plots.parameters <- plot.parameters(
      wisp.results = results,
      print.plots = FALSE, 
      verbose = FALSE
    )
    
    # Gather plots
    plots <- list(
      residuals = plots.residuals,
      ratecount = plots.ratecount,
      parameters = plots.parameters,
      MCMC = plots.MCMC,
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

# Helper function for computing p_values from bootstraps or MCMC samples
pvalues.samples <- function(
    mu.B,       # vector of bootstrapped or MCMC estimates
    mu.obs,     # observed value, either mean of mu.B or actual observation
    mu.null = 0 # null hypothesis value, usually 0
  ) {
    # Basic idea: Instead of centering data > bootstrapping > estimate parameter, bootstrap > estimate parameter > center data
    # equivalent to: mean(abs(mu.B - mean(mu.B) + mu.null) >= abs(mu.obs))
    Fn <- ecdf(abs(mu.B - mean(mu.B) + mu.null))
    return(
      1 - Fn(abs(mu.obs))
    )
  }

# Function for running stat analysis of bootstraps
sample.stats <- function(
    wisp.results,
    alpha = 0.05,
    Bonferroni = FALSE,
    conv.resamples.only = TRUE,
    null.rate = log(2),
    null.slope = 1,
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
    if (conv.resamples.only) {
      if (verbose) snk.report...("Grabbing sample results, only resamples with converged fit")
      sample_results <- wisp.results$sample.params[wisp.results$bs.diagnostics$success.code == 3,]
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
      
      # Grab count data 
      log_count <- wisp.results$count.data.summed$count.log
      child_col <- wisp.results$count.data.summed$child
      parent_col <- wisp.results$count.data.summed$parent
      ref_none_mask <- wisp.results$count.data.summed$ran == "none" & wisp.results$count.data.summed$treatment == "ref"
      max_bin <- max(wisp.results$count.data.summed$bin)
      
      # Grabbing indexes of parameters to test
      fitted_params_names <- wisp.results$param.names
      baselineRt_mask <- grepl("baseline", fitted_params_names) & grepl("Rt", fitted_params_names)        # mask for baseline rate values
      baseline_mask <- grepl("baseline", fitted_params_names) & !grepl("tpoint", fitted_params_names)     # mask for baseline rate and tslope values
      tslope_mask <- grepl("tslope", fitted_params_names)                                                 # mask for tslope values
      beta_w_mask <- grepl("beta_Rt|beta_tslope|beta_tpoint|wfactor", fitted_params_names)                # only want to test fixed-effect and wfactor parameters (no baseline parameters)
      empty_w_mask <- colSums(sample_results) == 0                                                        # Can only be zero for point wfactors of degree-zero children (which must be zero)
      test_mask <- beta_w_mask & !empty_w_mask
      
      # Compute 95% confidence intervals
      if (verbose) snk.report...("Computing 95% confidence intervals")
      sample.params_ci <- apply(sample_results, 2, quantile, probs = c(alpha/2, 1 - alpha/2))
      
      # Estimate p_values from bs params
      if (verbose) snk.report...("Estimating p-values from resampled parameters")
      p_values <- rep(NA, n_params)
      p_values_adj <- rep(NA, n_params)
      sig_marks <- rep(" ", n_params)
      for (n in seq_along(fitted_params)) {
        if (test_mask[n]) {
          p_values[n] <- pvalues.samples(sample_results[,n], fitted_params[n]) 
        } else if (baseline_mask[n]) {
          if (tslope_mask[n]) { # testing slope
            # Find rise: 
            parent_name <- sub(".*baseline_(.*?)_tslope.*", "\\1", fitted_params_names[n])
            child_name <- sub(".*tslope_(.*?)_Tns/Blk.*", "\\1", fitted_params_names[n])
            parent_child_Rt_mask <- grepl(parent_name, fitted_params_names) & grepl(child_name, fitted_params_names) & baselineRt_mask
            block_num <- sub(".*_Tns/Blk(.*)", "\\1", fitted_params_names[n])
            block_num_next <- as.character(as.integer(block_num) + 1)
            block_num <- paste0("Tns/Blk", block_num)
            block_num_next <- paste0("Tns/Blk", block_num_next)
            this_Rt_mask <- parent_child_Rt_mask & grepl(block_num, fitted_params_names)
            next_Rt_mask <- parent_child_Rt_mask & grepl(block_num_next, fitted_params_names)
            this_Rt <- fitted_params[this_Rt_mask]
            next_Rt <- fitted_params[next_Rt_mask]
            rise <- next_Rt - this_Rt
            if (is.na(rise) || is.infinite(rise)) stop("Problem with rise calculation")
            # Find mean of log-linked count
            count_mask <- parent_col == parent_name & child_col == child_name & ref_none_mask
            mean_log_count <- mean(log_count[count_mask], na.rm = TRUE)
            rise_run_scalar <- mean_log_count / max_bin
            # Compute p-value
            p_values[n] <- pvalues.samples(sample_results[,n], fitted_params[n], null.slope*rise_run_scalar/rise) 
          } else { # testing "log-linked" rate, i.e., true rate is exp of this number
            p_values[n] <- pvalues.samples(sample_results[,n], fitted_params[n], null.rate) 
          }
        }
        if (is.na(p_values[n]) && (test_mask[n] || baseline_mask[n])) stop("Problem with p-value calculation, NA")
      }
      
      # Sanity check
      num_of_tests <- sum(test_mask | baseline_mask)
      if (num_of_tests != sum(!is.na(p_values))) stop("Problem with p-value calculation")
      
      # Check resample size 
      recommended_resample_size <- num_of_tests / alpha
      if (verbose) {
        snk.report(paste0("Recommended resample size for alpha = ", alpha, ", ", num_of_tests, " tests"))
        snk.print_vec("with bootstrapping", c(recommended_resample_size))
        snk.print_vec("with MCMC", c(recommended_resample_size*10))
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

# Function for checking whether tpoints are stable, i.e., a genuine rate transition
check.tpoint.stability <- function(
    wisp.results,
    alpha = 0.05,
    verbose = TRUE
  ) {
   
    # "TPS" = transition point stability. 
    # ... a stable transition point is one where the slope of the rate is significantly greater than 1
    
    # Check that stats have been run
    if (length(wisp.results$stats$parameters) == 0) {
      snk.report...("Can't make TPS table without running stats")
    } else {
      
      # Grab names and other needed values
      fixed_effect_names <- wisp.results$treatment$names
      fixed_effect_components <- wisp.results$treatment$components
      names(fixed_effect_components) <- fixed_effect_names
      mc_names <- wisp.results$model.component.list
      p_val_adj_names <- wisp.results$stats$parameters$parameter
      parent_names <- as.character(wisp.results$grouping.variables$parent.lvls)
      child_names <- as.character(wisp.results$grouping.variables$child.lvls)
      sample.params <- wisp.results$sample.params
      
      # Grab baseline-tslope and beta masks
      beta_mask <- grepl("beta_", p_val_adj_names)
      tslope_mask <- grepl("tslope", p_val_adj_names)
      tslope_bl_mask <- grepl("baseline", p_val_adj_names) & tslope_mask
      tslope_beta_mask <- beta_mask & tslope_mask
      
      # Make column names 
      col_names <- c(
        wisp.results$variables$parent,
        wisp.results$variables$child,
        "TPS.ref", 
        paste0("TPS.", fixed_effect_names[-c(1)]), 
        as.character(levels(interaction(mc_names, fixed_effect_names[-c(1)])))
      )
      m <- length(as.character(levels(interaction(mc_names, fixed_effect_names[-c(1)])))) + length(fixed_effect_names[-c(1)])
      
      # Make row names
      row_names <- as.character(levels(interaction(child_names, parent_names)))
      n <- length(row_names)
      
      # Initialize data frame
      df <- data.frame(rep("", n), rep("", n), rep(FALSE, n))
      for (c in 1:m) df <- cbind(df, rep(FALSE, n))
      colnames(df) <- col_names
      rownames(df) <- row_names
      
      # Test effects
      for (gvp_lvl in parent_names) {
        for (gv_lvl in child_names) {
          
          # Grab parent and child names
          r <- paste0(gv_lvl,".",gvp_lvl)
          df[r,wisp.results$variables$parent] <- gvp_lvl
          df[r,wisp.results$variables$child] <- gv_lvl
          
          # Make mask for this parent-child combination
          pc_mask <- grepl(gvp_lvl, p_val_adj_names) & grepl(gv_lvl, p_val_adj_names)
          
          # Grab adjusted p-values, if any, for this child under this parent
          baseline_tslope_mask <- tslope_bl_mask & pc_mask
          
          # Check if a stable t-point in reference condition
          any_slopes <- any(baseline_tslope_mask)
          if (any_slopes) {
            if (any(wisp.results$stats$parameters$p.value.adj[baseline_tslope_mask] < alpha)) {
              df[r,"TPS.ref"] <- TRUE
            }
          }
          
          # Check for fixed effects on model components 
          for (fe in fixed_effect_names[-c(1)]) { # skip "ref" condition, as it is the baseline
            for (mc in mc_names) {
              
              # Grab the mask for this fixed effect and model component
              fe_mask <- grepl(paste0("_",fe,"_"), p_val_adj_names, fixed = TRUE)
              fe_mc_mask <- beta_mask & pc_mask & fe_mask & grepl(mc, p_val_adj_names)
              
              # Check if there are any significant p-values
              if (any(fe_mc_mask)) {
                if (any(wisp.results$stats$parameters$p.value.adj[fe_mc_mask] < alpha)) {
                  df[r,paste0(mc,".",fe)] <- TRUE
                }
              } else {
                df[r,paste0(mc,".",fe)] <- NA # means deg 1
              }
              
              # Check if TPS in this treatment condition
              if (mc == "tslope" && any_slopes) {
                
                # Grab list of components 
                fe_components <- fixed_effect_components[[fe]]
                
                # Check implied slope for each transition point
                baseline_tslope_idx <- which(baseline_tslope_mask)
                for (i in seq_along(baseline_tslope_idx)) {
                  
                  # Construct index of bs matrix columns to add
                  bs_col_idx <- c(baseline_tslope_idx[i])
                  for (fec in fe_components) {
                    
                    # Grab mask for this tslope effect for this parent-child pair
                    fec_mask <- grepl(paste0("_",paste0(fec),"_"), p_val_adj_names, fixed = TRUE)
                    fec_tslope_mask <- tslope_beta_mask & pc_mask & fec_mask
                    
                    # Sanity check, both should be the number of transition points
                    if (sum(fec_tslope_mask) != sum(baseline_tslope_mask)) stop("Problem testing for TPS in treatment condition.")
                    
                    # Add new column index
                    fec_tslope_idx <- which(fec_tslope_mask)
                    bs_col_idx <- c(bs_col_idx, fec_tslope_idx[i])
                    
                  }
                  
                  # Compute implied slopes for this effects condition in each bs resample
                  if (length(bs_col_idx) < 2) stop("Problem computing TPS in treatment condition (no effects present, baseline only)")
                  implied_bs_slopes <- rowSums(sample.params[,bs_col_idx])
                  
                  # Compute p_value ("- 1" because checking if significantly greater than 1)
                  p_value <- pvalues.samples(implied_bs_slopes - 0, mean(implied_bs_slopes) - 0) 
                  
                  # Test for significance
                  if (p_value < min(wisp.results$stats$parameters$alpha.adj[bs_col_idx], na.rm = TRUE)) df[r,paste0("TPS.",fe)] <- TRUE
                  
                }
                
              }
              
            }
          }
          
        }
      }
      
      # Print summary table of stat analysis 
      if (verbose) snk.print_table("TPS test results", df, head = TRUE, initial_breaks = 0, end_breaks = 0)
      
      # return
      return(df)
      
    }
    
  }

# Method for analyzing residuals
analyze.residuals <- function(
    wisp.results,
    verbose = TRUE
  ) {
    
    if (verbose) {
      snk.report("Analyzing residuals", initial_breaks = 1)
      snk.horizontal_rule(reps = snk.simple_break_reps)
    }
    
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
      
      clean_parents <- parents[ !is.na(resids)]
      clean_child <- child[ !is.na(resids)]
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
      
      df <- data.frame(convolved_quantile, resids)
      
      # Create the ggplot
      resid_plot <- ggplot(df, aes(x = convolved_quantile, y = resids)) +
        geom_point(color = "steelblue", size = 2) +  
        geom_abline(intercept = mean_resid, slope = sd_resid, color = "black", linewidth = 0.8, linetype = "dashed") +  
        labs(
          x = "Ordered Log-residuals",
          y = "Theoretical gamma-convolved Quantiles",
          title = paste0("Q-Q Plot of Log-residuals (", gp, ")")
        ) +
        theme_minimal()
      
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
        labs(x = "Log-residual Value", y = "Frequency", title = paste0("Histograms of Log-residuals (",gp,")")) +
        theme_minimal()
      
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

# Method for plotting model and data
plot.ratecount <- function(
    wisp.results,
    pred.type = "pred.log",
    count.type = "count.log",
    dim.boundaries = NULL,
    print.all = FALSE,
    y.lim = NA,
    count.alpha.none = NA, # These values have defaults which will be used if left NA
    count.alpha.ran = NA,
    pred.alpha.none = NA,
    pred.alpha.ran = NA,
    rans.to.print = NA,
    childs.to.print = NA
  ) {
   
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
    if (is.na(rans.to.print)) rans.to.print <- unique(df[,"ran"])
    if (is.na(childs.to.print)) {
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
      
      if (!is.null(dim.boundaries)) {
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
      
      if (length(dim.boundaries) != 0) {
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

# Method for printing plot of model parameters
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
    
    if (verbose) {
      snk.report("Plotting model parameters")
      snk.horizontal_rule(reps = snk.simple_break_reps)
    }
    
    # Grab fitted parameters 
    fitted_params <- wisp.results$fitted.parameters
    
    # Grab masks and set legend position
    if (verbose) snk.report...("Grabbing masks and setting legend position") 
    param_names <- wisp.results$param.names
    point_mask <- grepl("tpoint", param_names)
    rate_mask <- grepl("Rt", param_names)
    slope_mask <- grepl("tslope", param_names)
    pointR_mask <- grepl("wfactor_point", param_names) # for random effect on point
    rateR_mask <- grepl("wfactor_rate", param_names) # for random effect on rate
    slopeR_mask <- grepl("wfactor_slope", param_names) # for random effect on slope
    legpos <- "none"
    
    # Format parameter names for nice printing
    if (verbose) snk.report...("Formatting parameter names for nice printing") 
    param_names_saved <- param_names
    for (gvp_lvl in as.character(wisp.results$grouping.variables$parent.lvls)) {
      param_names <- gsub(paste0("baseline_",gvp_lvl,"_"), "", param_names)
      param_names <- gsub(paste0("_",gvp_lvl), "", param_names)
    }
    for (gv_lvl in as.character(wisp.results$grouping.variables$child.lvls)) {
      param_names <- gsub(paste0("_",gv_lvl,"_"), "_", param_names)
    }
    param_names <- gsub(paste0("wfactor_"), "", param_names)
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
    
    if (verbose) snk.report...("Making plots") 
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

# Method for plotting structural parameter distributions
plot.effect.dist <- function(
    wisp.results,
    verbose = TRUE 
  ) {
    
    if (verbose) {
      snk.report("Printing effect distributions")
      snk.horizontal_rule(reps = snk.simple_break_reps)
    }
    
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
      labs(title = "Distribution of random point effects (warping factors)", x = "Warping factor, point", y = "Density") +
      theme_minimal() 
    
    # ... for rate
    wfactors_rate_mask <- grepl("wfactor_rate", wisp.results$param.names)
    wfactors_rate <- c(bs_fitted_params[,wfactors_rate_mask])
    plot_wfactor_rate_effects_dist <- ggplot() +
      geom_histogram(
        data = data.frame(vals = wfactors_rate), aes(x = vals, y = after_stat(ndensity)),
        bins = 75, fill = "skyblue", alpha = 0.5, na.rm = TRUE) +
      labs(title = "Distribution of random rate effects (warping factors)", x = "Warping factor, rate", y = "Density") +
      theme_minimal() 
    
    # ... for slope 
    wfactors_slope_mask <- grepl("wfactor_slope", wisp.results$param.names)
    wfactors_slope <- c(bs_fitted_params[,wfactors_slope_mask])
    plot_wfactor_slope_effects_dist <- ggplot() +
      geom_histogram(
        data = data.frame(vals = wfactors_slope), aes(x = vals, y = after_stat(ndensity)),
        bins = 75, fill = "skyblue", alpha = 0.5, na.rm = TRUE) +
      labs(title = "Distribution of random slope effects (warping factors)", x = "Warping factor, slope", y = "Density") +
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
      labs(title = "Distribution of fixed rate effects", x = "Rate effect", y = "Density") +
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
        labs(title = "Distribution of fixed slope effects", x = "t-slope effect", y = "Density") +
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
        labs(title = "Distribution of fixed t-point effects", x = "t-point effect", y = "Density") +
        theme_minimal() 
      
      if (verbose) print(plot_tpoint_effects_effects_dist)
      plots.effects_dist[["plot_tpoint_effects_effects_dist"]] <- plot_tpoint_effects_effects_dist
      
    }
    
    return(plots.effects_dist)
    
  }

# Method for printing all child plots on one figure
plot.child.summary <- function(
    wisp.results,
    these.parents = NULL,
    these.childs = NULL,
    verbose = TRUE
  ) {
    
    if (verbose) {
      snk.report("Printing child summary plots", initial_breaks = 2)
      snk.horizontal_rule(reps = snk.simple_break_reps, end_breaks = 0)
    }
    
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
      
      if (verbose) snk.report(paste0("Making summary plots for ", gvp_lvl))
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

# Method to plot random walks from MCMC simulations 
plot.MCMC.walks <- function(
    wisp.results
  ) {
    
    # Grab sampled parameters
    sampled_params <- wisp.results$sample.params
    
    # Down-sample 
    these_rows <- as.integer(seq(1, nrow(sampled_params), length.out = 1000))
    sampled_params <- sampled_params[these_rows,]
    
    # Take abs and log
    sampled_params <- log(abs(sampled_params) + 1)
    
    # Make long-format data frame
    walks <- data.frame(
      value = c(sampled_params),
      param = rep(1:ncol(sampled_params), each = nrow(sampled_params)),
      step = rep(1:nrow(sampled_params), ncol(sampled_params))
    )
    
    # Make plot
    plot.walks.parameters <- ggplot(data = walks, 
                   aes(x = step, y = value, group = param, color = param)) + 
      geom_line() + 
      ggtitle("Parameter Walks from MCMC Estimation") +
      theme_minimal()
    print(plot.walks.parameters)
    
    # Grab trace of negloglik and normalize
    nll <- unlist(wisp.results[["MCMC.diagnostics"]][["neg.loglik"]])
    pnll <- unlist(wisp.results[["MCMC.diagnostics"]][["pen.neg.value"]])
    nll <- (nll - nll[1]) / nll[1]
    pnll <- (pnll - pnll[1]) / pnll[1]
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
        )
    print(plot.walks.nll)
    
  }

plot.decomposition <- function(
    wisp.results,
    child, # "Nptxr"
    log = FALSE, 
    dim.boundaries = NULL, # colMeans(count_data_WSPmm.y$db)
    y.lim = NULL
  ) {
    
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

# For trying out warping function in R
WSP.warp <- function( 
    x, # value to warp, either a scalar, vector, or 2D array/matrix
    b, # warping bound
    w  # warping factor
  ) {
    
    if (length(b) != 1) stop("b must be a single value for warping bound")
    if (length(w) != 1) stop("w must be a single value for warping factor")
    if (!(length(x) > 0)) stop("x must be a vector or matrix of values to warp, but has length zero")
    if (max(x) > b) stop("x must be less than or equal to b")
    
    if (is.null(dim(x)) || length(dim(x)) == 1) {
      out <- rep(NA, length(x))
      for (i in 1:length(x)) {
        out[i] <- warp_mc_R(x[i], b, w)
      }
    } else if (length(dim(x)) == 2) {
      out <- array(NA, dim = dim(x))
      for (i in 1:(dim(x)[1])) {
        for (j in 1:(dim(x)[2])) {
          out[i,j] <- warp_mc_R(x[i,j], b, w)
        }
      }
    } else {
      stop("x must be a vector or matrix of values to warp, but unrecognized dim attribute detected")
    }
    return(out)
    
  }

# Pre-made plot to explain warping function
demo_warp <- function(
    w = 2,                       # warping factor
    point_pos = 60,
    point_neg = 40,
    Rt = c(6, 3, 0.2, 6)*4.65,   # rates for poly-sigmoid
    tslope = c(0.4, 0.75, 1),    # slope scalars for poly-sigmoid
    tpoint = c(15, 38, 80),      # transition points for poly-sigmoid
    w_factors = c(0.6, -0.9, 0.5)      # warping factors for poly sigmoid
  ) {
    
    # Data
    n <- 1000
    x <- (1:n)/10
    b <- 100
    y <- x
    y1 <- WSP.warp(x, b, w)
    y2 <- WSP.warp(x, b, -w)
    y_pos <- WSP.warp(point_pos, b, w)
    y_neg <- WSP.warp(point_neg, b, -w)
    
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
      geom_line(linewidth = 1.5) +
      geom_hline(yintercept = 100, linetype = "dashed", color = "darkgray", linewidth = 1) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray", linewidth = 1) +
      #coord_fixed(ratio = 1) +
      geom_segment(
        data = df_segments,
        aes(x = point_pos, xend = point_pos, y = point_pos, yend = y_pos),
        color = "blue4", linetype = "dashed", linewidth = 0.75) + 
      geom_segment(
        data = df_segments,
        aes(x = point_neg, xend = point_neg, y = point_neg, yend = y_neg),
        color = "red4", linetype = "dashed", linewidth = 0.75) +
      annotate("text", x = 10, y = 95, label = "upper asymptote", size = 7.5, color = "black") +
      annotate("text", x = 90, y = 5, label = "lower asymptote", size = 7.5, color = "black") +
      annotate("text", x = point_pos - 10, y = (y_pos + point_pos)/2, label = expression(varphi * "(z)(b - z)"), size = 7.5, color = "black") +
      annotate("text", x = point_neg + 10, y = (y_neg + point_neg)/2, label = expression(varphi * "(b - z)z"), size = 7.5, color = "black") +
      labs(
        x = "z",
        y = expression(omega * "(z, " * rho * ", b)"),,
        title = "WSP Warping Function",
        color = "Direction"
      ) +
      scale_color_manual(
        values = c("black", "red", "blue"),
        labels = c(expression(rho * " = 0"), expression(rho * " < 0"), expression(rho * " > 0"))
      ) +
      theme_minimal(base_size = 16) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 30),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18)
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
    Rtw <- WSP.warp(Rt, 100000, w_factors[1])
    tslopew <- WSP.warp(tslope, 100000, w_factors[2])
    tpointw <- WSP.warp(tpoint, 100, w_factors[3])
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
        plot.title = element_text(hjust = 0.5, size = 30),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18)
      ) 
    
    grid.arrange(demo_plot_warpfunction, demo_plot_warpedsigmoid, ncol = 1)
    
    # Save at 1156 x 843
    return(demo_plot_warpfunction)
    
  }

# Pre-made plot to explain poly-sigmoid function 
demo_sigmoid <- function(
    r = 4,                       # upper asymptote for logistic
    s = 1,                       # slope scalar at inflection point
    Rt = c(6, 3, 0.2, 6)*4.65,   # rates for poly-sigmoid
    tslope = c(0.4, 0.75, 1),    # slope scalars for poly-sigmoid
    tpoint = c(15, 38, 80)       # transition points for poly-sigmoid
  ) {
    
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
      annotate("text", x = -7, y = r+0.75, label = "upper asymptote", size = 7.5, color = "black") +
      annotate("text", x = 7, y = -0.75, label = "lower asymptote", size = 7.5, color = "black") +
      annotate("text", x = 3, y = r+0.75, label = "inflection point", size = 7.5, color = "red4") +
      annotate("text", x = -2, y = r+1, angle = atan(m)*57.3, label = "slope", size = 7.5, color = "blue4") +
      labs(
        x = "x",
        y = expression(psi * "(x, r = 4, s = 1)"),
        title = "The Logistic Function"
      )  +
      theme_minimal(base_size = 16) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 30),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18)
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
        deg,      # degree 
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
      coord_fixed()  +
      geom_segment(
        data = def_segments_rates,
        aes(x = x0, xend = x1, y = y0, yend = y1),
        color = "darkgray", linetype = "dashed", linewidth = 1) +
      geom_segment(
        data = def_segments_slopes,
        aes(x = slope_seg_neg_x, xend = slope_seg_pos_x, y = slope_seg_neg_y, yend = slope_seg_pos_y),
        color = "blue4", linetype = "dashed", linewidth = 0.75) + 
      annotate("text", x = block_midpoints, y = Rt+mean(Rt)*0.15, label = paste("rate", 1:length(Rt)), size = 7.5, color = "black") +
      annotate("text", x = tpoint-max(x)*0.08, y = -max(y)*0.1, label = paste("t-point", 1:length(tpoint)), size = 7.5, color = "red4") +
      annotate("text", x = tpoint-max(x)*0.035, y = y0*0.95, angle = atan(m)*57.3, label = paste("slope", 1:length(tslope)), size = 7.5, color = "blue4") +
      labs(
        x = "x",
        y = expression(Psi * "(x, r, s, p)"),
        title = "The WSP Sigmoid Function"
      )  +
      theme_minimal(base_size = 16) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 30),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18)
      )
    
    grid.arrange(demo_plot_logistic, demo_plot_sigmoid, ncol = 1)
    # Saved at 1186 x 1032
    
  }

# R.sum and M.hat helper functions to enforce block-rate constraints
#  ... These are no longer used as part of the prediction or optimization functions (most of that is 
#      now in cpp only), but these are used in the "make" functions for when generating parameters 
#      and data for simulations. 
#  ... Eventually want to move everything, including these, fully to Rcpp.

R.sum <- function(
    R, 
    m.hat
) { 
  # When the rate of the first block is greater than this sum the rate cannot be negative
  lengthR <- length(R)
  if (lengthR != length(m.hat) + 1) stop("R must have length 1 greater than m.hat")
  r_sum <- sum((R[-c(lengthR)] - R[2:lengthR])/m.hat)
  return(r_sum)
}

M.hat <- function(
    tslopes, 
    tpoints, 
    max.bin.dim
) {
  
  # Original formula, derived analytically (do not delete):
  # m.hat <- vector("list", length(tslopes))
  # for (j in seq_along(tslopes)) m.hat[j] <- min(1+exp(-tslopes[j]*(1:max.bin.dim - tpoints[j]))) 
  
  #  ... However, this is very slow. It can be shown (also analytically) that the above formula is equivalent to: 
  return(1+exp(-tslopes*(max.bin.dim - tpoints)))
  #  ... which is computationally much faster.
  
}



