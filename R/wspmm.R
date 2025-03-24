
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
    use.median = TRUE,
    MCMC.burnin = 0,
    MCMC.steps = 1e4,
    MCMC.step.size = 0.1,
    MCMC.prior = 0.5,
    bootstraps.num = 0, 
    converged.resamples.only = FALSE,
    max.fork = 10,
    dim.bounds = NULL, 
    verbose = TRUE,
    print.child.summaries = TRUE,
    # Setting to pass to C++ model
    model.settings = list(
      struc_values = c(5.0, 5.0, 1.0),            # values of structural parameters to test
      buffer_factor = 0.05,                       # buffer factor for penalizing distance from structural parameter values
      ctol = 1e-6,                                # convergence tolerance
      max_penalty_at_distance_factor = 0.01,      # maximum penalty at distance from structural parameter values
      LROcutoff = 2.0,                            # cutoff for LROcp
      LROwindow_factor = 2.0,                     # controls size of window used in LROcp algorithm (window = LROwindow_factor * bin_num * buffer_factor)
      LROfilter_ws_divisor = 2.0,                 # divisor for filter window size in likelihood ratio outlier detection (bigger is smaller window)
      tslope_initial = 1.0,                       # initial value for tslope
      wf_initial = 0.1,                           # initial value for wfactor
      max_evals = 1000,                           # maximum number of evaluations for optimization
      rng_seed = 42                               # seed for random number generator
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
    
    # Return initial values and compute effect difference ratio 
    fe_diff_ratio <- analyze.diff(
      wisp.results = cpp_model$results(),
      count.type = "count.log"
    )
    cpp_model$import_fe_diff_ratio_Rt(fe_diff_ratio$diff.ratio.c, verbose)
    cpp_model$import_fe_diff_ratio_tpoint(fe_diff_ratio$diff.ratio.tp, verbose)
    
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
      MCMC_walk <- MCMC_walk[-c(2:MCMC.burnin),]
     
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
    } else {
      n_params <- ncol(MCMC_walk)
      sample.params <- MCMC_walk
      bs.diagnostics <- NULL
    }
    
    # Set final fitted parameters
    if (use.median) {
      if (verbose) snk.report...("Setting median parameter samples as final parameters", initial_breaks = 1, end_breaks = 1)
      final_parameters <- apply(sample.params, 2, function(x) median(x, na.rm = TRUE))
    } else {
      if (verbose) snk.report...("Setting full-data fit as parameters", initial_breaks = 1, end_breaks = 1)
      if (bootstraps.num > 0) {
        final_parameters <- sample.params[bootstraps.num + 1,]
      } else {
        final_parameters <- sample.params[1,]
      }
    }
    cpp_model$set_parameters(
      final_parameters,
      verbose
    )
    
    # Grab model results and add samples
    results <- cpp_model$results()
    results[["sample.params"]] <- sample.params
    results[["bs.diagnostics"]] <- bs.diagnostics
    
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
    
    # Check tpoint stability
    results$stats$tps <- check.tpoint.stability(
      wisp.results = results,
      verbose = verbose
    )
    
    # Make plots of results ####
    
    # Plot MCMC walks 
    if (bootstraps.num == 0) {
      plots.MCMC <- plot.MCMC.walks(
        wisp.results = results
      )
    } else {
      plots.MCMC <- NULL
    }
    
    # Plot structural parameter stats
    plots.struc <- plot.struc.stats(
      wisp.results = results,
      verbose = verbose
    )
    
    # Make rate plots 
    if (verbose) snk.report...("Making rate-count plots")
    plots.ratecount <- plot.ratecount(
      wisp.results = results,
      pred.type = "pred.log",
      count.type = "count.log",
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
      struc = plots.struc
    )
    results[["plots"]] <- plots
    
    # Print summary plots
    if (print.child.summaries) {
      plot.child.summary(
        wisp.results = results,
        these.parents = NULL,
        these.childs = NULL,
        log.scale = TRUE,
        verbose = TRUE
      )
    }
    
    return(results)
    
  }

# Analysis methods #####################################################################################################

# Helper function for computing p_values from bootstraps or MCMC samples
pvalues.samples <- function(
    mu.B,    # vector of bootstrapped or MCMC estimates
    mu.obs,  # observed value, either mean of mu.B or actual observation
    lower.tail.only = FALSE
  ) {
    
    if (lower.tail.only) return(mean((mu.B - mean(mu.B)) + mu.obs <= 0))
    else return(mean(abs(mu.B - mean(mu.B)) >= abs(mu.obs)))
    
    # Basic idea (two-tailed): Instead of centering data > bootstrapping > estimate parameter, bootstrap > estimate parameter > center data
    
    # Basic idea (one-tailed): Want to test if mu.obs < 0, based on variation of values in mu.B < mu.obs
    # ... "(mu.B - mean(mu.B))" is < 0 iff the bs is less than the mean
    # ... If "(mu.B - mean(mu.B)) + mu.obs" is < 0 and mu.obs > 0, then "(mu.B - mean(mu.B))" is undersized by amount greater than mu.obs's distance from 0
    # ... If "(mu.B - mean(mu.B)) + mu.obs" is < 0 and mu.obs < 0, then "(mu.B - mean(mu.B))" is not oversized by amount greater than mu.obs's distance from 0
    # ... Idea for p value: Want to quantify the number of bs that are either undersized more than the obs exceeds 0, or not oversized enough to compensate for the distance of mu.obs under 0
    
  }

# Function for running stat analysis of bootstraps
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
      
      # Grabbing indexes of parameters to test
      fitted_params_names <- wisp.results$param.names
      struc_mask <- !grepl("baseline|beta_tslope|beta_tpoint|wfactor", fitted_params_names) 
      baseline_mask <- grepl("baseline", fitted_params_names) & !grepl("tpoint", fitted_params_names)     # mask for baseline rate and tslope values
      tslope_mask <- grepl("tslope", fitted_params_names)                                                 # mask for tslope values
      beta_w_mask <- grepl("beta_Rt|beta_tslope|beta_tpoint|wfactor", fitted_params_names)                # only want to test fixed-effect and wfactor parameters (no baseline or structural parameters)
      empty_w_mask <- colSums(sample_results) == 0                                                        # Can only be zero for point wfactors of degree-zero children (which must be zero)
      test_mask <- beta_w_mask & !empty_w_mask
      
      # Compute 95% confidence intervals
      if (verbose) snk.report...("Computing 95% confidence intervals")
      sample.params_ci <- apply(sample_results, 2, quantile, probs = c(alpha/2, 1 - alpha/2))
      
      # Estimate p_values from bs params
      if (verbose) snk.report...("Estimating p-values from bootstraped parameters")
      p_values <- rep(NA, n_params)
      p_values_adj <- rep(NA, n_params)
      sig_marks <- rep(" ", n_params)
      for (n in seq_along(fitted_params)) {
        if (test_mask[n]) {
          p_values[n] <- pvalues.samples(sample_results[,n], fitted_params[n]) 
        } else if (baseline_mask[n]) {
          if (tslope_mask[n]) { # testing slope
            p_values[n] <- pvalues.samples(sample_results[,n] - 0, fitted_params[n] - 0, lower.tail.only = TRUE) 
            # Basic idea: A slope < 1 is unstably shallow and suggests that there is no transition point here
            #   ... hence, expect to rest the value - 1; however, the model uses the exp of the parameter, so zero becomes our 1. 
          } else { # testing "log-linked" rate, i.e., true rate is exp of this number
            p_values[n] <- pvalues.samples(sample_results[,n], fitted_params[n]) 
            # ^ ... These values can be negative, meaning unlinked (i.e., taking exp) they are less than 1, which we interpret as not expressed
            #   ... As there are multiple cells per bin, perhaps the cutoff shouldn't be 1, but the number of cells (so avg is 1)??
          }
        }
      }
      
      # Sanity check
      num_of_tests <- sum(test_mask | baseline_mask)
      if (num_of_tests != sum(!is.na(p_values))) stop("Problem with p-value calculation")
      
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
                  
                  # Compute p_value ("- 1" and lower-tail only because checking if significantly greater than 1)
                  p_value <- pvalues.samples(implied_bs_slopes - 0, mean(implied_bs_slopes) - 0, lower.tail.only = TRUE) 
                  
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
      qq_resids_plot <- ggplot(df_wide, aes(sample = residuals)) +
        stat_qq(color = "steelblue", size = 1) +  # Q-Q points
        stat_qq_line(color = "black", linetype = "dashed") +  # Reference line
        labs(x = "Theoretical Quantiles", y = "Sample Quantiles", 
             title = paste0("Q-Q Plot of Log-residuals (", gp, ")")) +
        theme_minimal()
      
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

analyze.diff <- function(
    wisp.results,
    count.type     # character, either "count", "pred", "count.log", or "pred.log"
  ) {
    
    # Grab count data frame
    count_data <- wisp.results$count.data.summed
    
    # Grab random effects levels 
    ran_lvls <- wisp.results$grouping.variables$ran.lvls[-c(1)] # remove "none" level
    n_ran_lvls <- length(ran_lvls)
    
    # Get treatment levels instanced by each ran level
    trts <- wisp.results$treatment$names
    trts_comps <- wisp.results$treatment$components
    names(trts_comps) <- trts
    n_trts <- length(trts)
    trt_lvls_in_ran <- list()
    for (m in ran_lvls) {
      mask <- count_data$ran == m & !is.na(count_data$count)
      trt_lvls_in_ran[[m]] <- unique(count_data$treatment[mask])
    }
    
    # Grab fixed effects levels
    fix <- wisp.results$fix
    names(fix$lvls) <- fix$name
    names(fix$ref.lvl) <- fix$name
    names(fix$treat.lvl) <- fix$name
    n_fix <- length(fix$name)
    
    # Make treatment component array
    trt_comp_array <- array(
      as.character(NA),
      dim = c(n_trts, n_fix),
      dimnames = list(trts, fix$name)
    )
    for (t in trts) {
      if (t == "ref") {
        trt_comp_array[t,] <- fix$ref.lvl
      } else {
        for (fe in fix$name) {
          if (any(trts_comps[[t]] %in% fix$lvls[[fe]])) {
            trt_comp_array[t, fe] <- fix$treat.lvl[[fe]]
          } else {
            trt_comp_array[t, fe] <- fix$ref.lvl[fe]
          }
        }
      }
    }
    
    # Identify which treatment pairs differ only on one component 
    trt_oneoff_array <- array(
      FALSE,
      dim = c(n_trts, n_trts),
      dimnames = list(trts, trts)
    )
    for (ti in 1:(n_trts - 1)) {
      for (tj in (ti + 1):n_trts) {
        if (sum(trt_comp_array[ti,] != trt_comp_array[tj,]) == 1) {
          trt_oneoff_array[ti,tj] <- TRUE
          trt_oneoff_array[tj,ti] <- TRUE
        }
      }
    }
    
    # Function to make data matrices 
    # ... will randomly nudge the provided tpoints
    make_data_matrices <- function(
      count_type 
    ) {
      
      # Grab prelim info
      parents <- wisp.results$grouping.variables$parent.lvls
      n_parents <- length(parents)
      children <- wisp.results$grouping.variables$child.lvls
      n_children <- length(children)
      bins <- unique(count_data$bin)
      n_bins <- length(bins)
      param_names <- names(wisp.results$fitted.parameters)
      degs <- array( 
        0,
        dim = c(n_children, n_parents),
        dimnames = list(child = children, parent = parents)
      )
      tpoints <- wisp.results$change.points
      
      # Make full count data matrices
      data0 <- array(
        NA, 
        dim = c(
          n_bins, n_children, n_parents,
          n_ran_lvls,
          n_trts
        ),
        dimnames = list(bin = bins, child = children, parent = parents, ran = ran_lvls, fixedeffect = trts)
      )
      for (m in ran_lvls) {
        for (t in trt_lvls_in_ran[[m]]) {
          for (g in children) {
            for (p in parents) {
              mask <- count_data$ran == m & count_data$treatment == t & count_data$child == g & count_data$parent == p
              if (sum(mask) == n_bins) data0[, g, p, m, t] <- count_data[,count_type][mask]
            }
          }
        }
      }
      
      # Find t-points and degrees for each parent-child pair
      for (p in parents) {
        for (g in children) {
          degs[g, p] <- 0
          tpoints_pg <- tpoints[[p]][[g]]
          if (length(tpoints_pg) > 0) degs[g, p] <- nrow(tpoints_pg)
        }
      }
      max_deg <- max(degs)
      max_blocks <- max_deg + 1
      
      # Make tpoint matrices 
      datat <- array(
        NA,
        dim = c(
          max_deg, n_children, n_parents,
          n_ran_lvls + 1,
          n_trts
        ),
        dimnames = list(tp = 1:max_deg, child = children, parent = parents, ran = c("none", ran_lvls), fixedeffect = trts)
      )
      for (p in parents) {
        for (g in children) {
          for (m in ran_lvls) {
            for (t in trt_lvls_in_ran[[m]]) {
              t_num <- which(trts == t)
              m_num <- which(c("none", ran_lvls) == m)
              tpoints_array_col_num <- (t_num - 1) * (n_ran_lvls + 1) + m_num
              tpoints_pg <- tpoints[[p]][[g]]
              if (length(tpoints_pg) > 0) {
                tpoints_pg_this <- tpoints_pg[, tpoints_array_col_num]
                # Randomly nudge 
                jit <- sample((-3):3, 1)
                tpoints_pg_this <- tpoints_pg_this + jit
                # Ensure in bounds
                if (any(tpoints_pg_this < 3)) tpoints_pg_this[tpoints_pg_this < 3] <- 3
                if (any(tpoints_pg_this > n_bins - 2)) tpoints_pg_this[tpoints_pg_this > n_bins - 2] <- n_bins - 2
                datat[1:degs[g, p], g, p, m, t] <- tpoints_pg_this
              }
            }
          }
        }
      }
      
      # Collapse full count data matrices down to blocks
      data <- array(
        NA, 
        dim = c(
          max_blocks, n_children, n_parents,
          n_ran_lvls,
          n_trts
        ),
        dimnames = list(bin = 1:max_blocks, child = children, parents, ran = ran_lvls, fixedeffect = trts)
      )
      for (m in ran_lvls) {
        for (t in trt_lvls_in_ran[[m]]) {
          if (length(data0[, , , m, t]) > 0) {
            for (g in children) {
              for (p in parents) {
                if (degs[g, p] == 0) {
                  data[1, g, p, m, t] <- mean(data0[, g, p, m, t])
                } else {
                  tpoints_pgtm <- datat[, g, p, m, t]
                  for (i in 1:(degs[g, p] + 1)) {
                    if (i == 1) {
                      row_batch <- 1:(tpoints_pgtm[i] - 1) 
                    } else if (i > degs[g, p]) {
                      row_batch <- tpoints_pgtm[i - 1]:n_bins
                    } else {
                      row_batch <- tpoints_pgtm[i - 1]:(tpoints_pgtm[i] - 1) 
                    }
                    data[i, g, p, m, t] <- mean(data0[row_batch, g, p, m, t])
                  }
                }
              }
            }
          }
        }
      }
      
      return(
        list(
          count_data = data,
          tpoint_data = datat
        )
      )
      
    }
    
    # Function for finding normalized differences between random levels
    find_norm_diffs <- function(
      data,
      random # If not random, then one-off
    ) {
      
      # Initialize array to track which random levels and treatments have been compared
      running_tally <- array(
        as.character(NA),
        dim = c(2, 4)
      )
      
      # Initialize variable to hold normalized differences
      norm_diffs <- c()
      
      # ... then for each treatment level
      for (t1 in trts) {
        for (t2 in trts) {
          
          if (random) {
            if (t1 != t2) next
          } else if (!trt_oneoff_array[t1, t2]) {
            next 
          }
          
          # For each unique random-level pair 
          for (m1 in ran_lvls) {
            if (t1 %in% trt_lvls_in_ran[[m1]]) {
              for (m2 in ran_lvls) {
                
                # Looking at variation between random levels, so don't compare levels to themselves
                if (t2 %in% trt_lvls_in_ran[[m2]] && m1 != m2) {
                  
                  # Grab random levels and treatments for this comparison
                  this_row <- c(m1, m2, t1, t2)
                  
                  # Check if these random levels and treatments have been compared and skip if so
                  if (any(apply(running_tally, 1, function(x) all(this_row %in% x)))) {
                    next
                  } else {
                    
                    # Grab data matrices for each random level and treatment level
                    values1 <- data[, , , m1, t1]
                    values2 <- data[, , , m2, t2]
                    
                    # Find the difference between these two random levels and associated treatments 
                    diffs <- abs(values1 - values2)
                    
                    # If this is a legit comparison (no NA), normalize the differences and add to the list
                    if (sum(!is.na(c(diffs))) > 0) {
                      diffs <- diffs/values1
                      diffs <- na.omit(c(diffs))
                      diffs <- diffs[is.finite(diffs)]
                      norm_diffs <- c(norm_diffs, diffs)
                      running_tally <- rbind(running_tally, this_row)
                    }
                    
                  }
                  
                }
              }
            }
          }
        }
      }
      
      return(norm_diffs)
      
    }
    
    # Run analysis a bunch and take means 
    n_runs <- 100
    diff.ratio.c_i <- rep(0, n_runs)
    diff.ratio.tp_i <- rep(0, n_runs)
    for (i in 1:n_runs) {
      
      # Make data matrices with random nudge of tpoints
      data.both <- make_data_matrices(count.type)
      data.c <- data.both$count_data
      data.tp <- data.both$tpoint_data
      
      # For reach treatment level of each fixed effect, find normalized differences between random levels
      ran.diffs.c <- find_norm_diffs(data.c, TRUE)
      oneoff.diffs.c <- find_norm_diffs(data.c, FALSE)
      diff.ratio.c <- mean(oneoff.diffs.c)/mean(ran.diffs.c)
      
      ran.diffs.tp <- find_norm_diffs(data.tp, TRUE)
      oneoff.diffs.tp <- find_norm_diffs(data.tp, FALSE)
      diff.ratio.tp <- mean(oneoff.diffs.tp)/mean(ran.diffs.tp)
      
      # Mathematically, this must be greater than zero, or model fitting will fail
      # ... from a physical point of view, even if the true effect is zero, it will 
      #      be measured as non-zero due to noise, so we should expect this to be greater than 1. 
      if (diff.ratio.c < 1.05) {diff.ratio.c <- 1.05}
      if (diff.ratio.tp < 1.05) {diff.ratio.tp <- 1.05}
      
      diff.ratio.c_i[i] <- diff.ratio.c
      diff.ratio.tp_i[i] <- diff.ratio.tp
      
    }
    
    # Take means 
    diff.ratio.c <- mean(diff.ratio.c_i)
    diff.ratio.tp <- mean(diff.ratio.tp_i)
    
    return(
      list(
        diff.ratio.c = diff.ratio.c,
        diff.ratio.tp = diff.ratio.tp
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
    param_names <- gsub(paste0("beta_shape"), "bshape", param_names) # protect struc parameters
    param_names <- gsub(paste0("beta_"), "", param_names)
    param_names <- gsub(paste0("bshape"), "beta_shape", param_names) # put back
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
            ranEff_fitted_rate <- c(as.matrix(sample.params[sample(1:nrow(sample.params), downsample_size, replace = FALSE), ranEff_mask & !pointR_mask]))
            params <- rep(param_names[ranEff_mask & !pointR_mask], each = downsample_size)
            means <- rep(fitted_params[ranEff_mask & !pointR_mask], each = downsample_size)
          } else {
            ranEff_fitted_rate <- fitted_params[ranEff_mask & !pointR_mask]
            params <- param_names[ranEff_mask & !pointR_mask]
            means <- fitted_params[ranEff_mask & !pointR_mask]
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
              ranEff_fitted_rate_df$value_low <- rep(sample_ci[1,ranEff_mask & !pointR_mask], each = downsample_size)
              ranEff_fitted_rate_df$value_high <- rep(sample_ci[2,ranEff_mask & !pointR_mask], each = downsample_size)
              ranEff_fitted_rate_df$sig_marks <- rep(sig_marks[ranEff_mask & !pointR_mask], each = downsample_size)
            } else {
              ranEff_fitted_rate_df$value_low <- sample_ci[1,ranEff_mask & !pointR_mask]
              ranEff_fitted_rate_df$value_high <- sample_ci[2,ranEff_mask & !pointR_mask]
              ranEff_fitted_rate_df$sig_marks <- sig_marks[ranEff_mask & !pointR_mask]
            }
          }
          ranEff_df <- rbind(ranEff_df, ranEff_fitted_rate_df) 
          
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
          facet_wrap(~ interaction(child,type), scales = "free") +
          theme_minimal() +
          labs(title = title_ran) +
          theme(legend.position = legpos) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
          scale_fill_manual(values = c("point" = "#c356ea", "rate" = "#ffc100")) 
        
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
        
        # Make structural parameter plots ####
        if (length(child.lvl) == 0 && all(names(parameter_comparison_plots) != "plot_struc")) {
          struc_df <- data.frame(
            value = wisp.results$struc.params, 
            parameter = names(wisp.results$struc.params)
          )
          if (print_stats) {
            struc_df$value_low <- rep(NA, nrow(struc_df))
            struc_df$value_high <- rep(NA, nrow(struc_df))
            for (n in names(wisp.results$struc.params)) {
              idx <- which(param_names == n)
              struc_df$value_low[struc_df$parameter == n] <- sample_ci[1,idx]
              struc_df$value_high[struc_df$parameter == n] <- sample_ci[2,idx]
            }
          }
          parameter_comparison_plots[["plot_struc"]] <- ggplot(struc_df, aes(x = parameter, y = value)) +
            geom_bar(stat = "identity", fill = "cadetblue1", na.rm = TRUE) +
            theme_minimal() +
            labs(title = "Structural parameters") +
            theme(legend.position = "none") +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))
          if (print_stats) {
            parameter_comparison_plots[["plot_struc"]] <- parameter_comparison_plots[["plot_struc"]] +
              geom_errorbar(
                data = struc_df[struc_df$parameter != "buffer_factor", ],
                aes(ymin = value_low, ymax = value_high),
                width = 0.2, 
                na.rm = TRUE) 
            # ... all structural values must be above zero ...
            # geom_text(
            #   data = struc_df[struc_df$parameter != "buffer_factor", ], 
            #   aes(x = parameter, 
            #       y = value_high + star_gap_factor*(value_high-value_low), 
            #       label = sig_marks,
            #       size = sig_marks_size),
            #       size.unit = sig_marks_unit
            #   )  +
            # scale_y_continuous(
            #   expand = expansion(mult = c(0, expand_factor))
            # )
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
plot.struc.stats <- function(
    wisp.results,
    verbose = TRUE 
  ) {
    
    if (verbose) {
      snk.report("Printing structural summary distributions")
      snk.horizontal_rule(reps = snk.simple_break_reps)
    }
    
    plots.struc_stats <- list()
    bs_fitted_params <- wisp.results$sample.params
    
    # Warping factors, beta distributions
    if (verbose) snk.report...("Warping factors, beta distributions")
    shape_point <- wisp.results$struc.params["beta_shape_point"]
    shape_rate <- wisp.results$struc.params["beta_shape_rate"]
    rates <- seq(0,1,0.01)
    data <- data.frame(
      rate = rates*2 - 1,
      probability_point = dbeta(rates, shape_point, shape_point),
      probability_rate = dbeta(rates, shape_rate, shape_rate)
    )
    wfactors_point_mask <- grepl("wfactor_point", wisp.results$param.names)
    wfactors_point <- c(bs_fitted_params[,wfactors_point_mask])
    data$probability_point <- data$probability_point / max(data$probability_point)
    plot_wfactor_point_struc_stats <- ggplot() +
      geom_histogram(
        data = data.frame(vals = wfactors_point), aes(x = vals, y = after_stat(ndensity)),
        bins = 75, fill = "skyblue", alpha = 0.5, na.rm = TRUE) +
      geom_line(aes(x = rate, y = probability_point), data = data, linewidth = 1.2, na.rm = TRUE) +
      labs(title = "Beta distribution model of random point effects", x = "Warping factor, point", y = "Density") +
      theme_minimal() 
    
    wfactors_rate_mask <- grepl("wfactor_rate", wisp.results$param.names)
    wfactors_rate <- c(bs_fitted_params[,wfactors_rate_mask])
    data$probability_rate <- data$probability_rate /max(data$probability_rate)
    plot_wfactor_rate_struc_stats <- ggplot() +
      geom_histogram(
        data = data.frame(vals = wfactors_rate), aes(x = vals, y = after_stat(ndensity)),
        bins = 75, fill = "skyblue", alpha = 0.5, na.rm = TRUE) +
      geom_line(aes(x = rate, y = probability_rate), data = data, linewidth = 1.2, na.rm = TRUE) +
      labs(title = "Beta distribution model of random rate effects", x = "Warping factor, rate", y = "Density") +
      theme_minimal() 
    
    if (verbose) print(plot_wfactor_point_struc_stats)
    plots.struc_stats[["plot_wfactor_point_struc_stats"]] <- plot_wfactor_point_struc_stats
    if (verbose) print(plot_wfactor_rate_struc_stats)
    plots.struc_stats[["plot_wfactor_rate_struc_stats"]] <- plot_wfactor_rate_struc_stats
    
    # Rate effects, Gaussian distribution
    if (verbose) snk.report...("Rate effects, Gaussian distribution")
    Rt_shape <- wisp.results$computed_sd_Rt_effect
    rate_effs_mask <- grepl("Rt", wisp.results$param.names) & grepl("beta", wisp.results$param.names)
    rate_effects <- c(bs_fitted_params[,rate_effs_mask])
    rates <- seq(min(rate_effects), max(rate_effects), length.out = 100)
    data <- data.frame(
      rate = rates,
      probability = dnorm(rates, 0, Rt_shape)
    )
    data$probability <- data$probability / max(data$probability)
    plot_rate_effects_struc_stats <- ggplot() +
      geom_histogram(
        data = data.frame(vals = rate_effects), aes(x = vals, y = after_stat(ndensity)),
        bins = 75, fill = "skyblue", alpha = 0.5, na.rm = TRUE) +
      geom_line(aes(x = rate, y = probability), data = data, linewidth = 1.2, na.rm = TRUE) +
      xlim(min(rate_effects),max(rate_effects)) +
      labs(title = "Gaussian distribution model of fixed rate effects", x = "Rate effect", y = "Density") +
      theme_minimal() 
    
    if (verbose) print(plot_rate_effects_struc_stats)
    plots.struc_stats[["plot_rate_effects_struc_stats"]] <- plot_rate_effects_struc_stats
    
    # Slope effects, Gaussian distribution
    if (verbose) snk.report...("Slope effects, Gaussian distribution")
    tslope_shape <- wisp.results$struc.params["sd_tslope_effect"]
    tslope_effs_mask <- grepl("tslope", wisp.results$param.names) & grepl("beta", wisp.results$param.names)
    if (any(tslope_effs_mask)) {
      
      tslope_effects <- c(bs_fitted_params[,tslope_effs_mask])
      rates <- seq(min(tslope_effects), max(tslope_effects), length.out = 100)
      data <- data.frame(
        rate = rates,
        probability = dnorm(rates, 0, tslope_shape)
      )
      data$probability <- data$probability / max(data$probability)
      plot_slope_effects_struc_stats <- ggplot() +
        geom_histogram(
          data = data.frame(vals = tslope_effects), aes(x = vals, y = after_stat(ndensity)),
          bins = 75, fill = "skyblue", alpha = 0.5, na.rm = TRUE) +
        geom_line(aes(x = rate, y = probability), data = data, linewidth = 1.2, na.rm = TRUE) +
        xlim(min(tslope_effects),max(tslope_effects)) +
        labs(title = "Gaussian distribution model of fixed slope effects", x = "t-slope effect", y = "Density") +
        theme_minimal() 
      
      if (verbose) print(plot_slope_effects_struc_stats)
      plots.struc_stats[["plot_slope_effects_struc_stats"]] <- plot_slope_effects_struc_stats
      
    }
    
    # tpoint effects, Gaussian distribution
    if (verbose) snk.report...("tpoint effects, Gaussian distribution")
    tpoint_shape <- wisp.results$computed_sd_tpoint_effect
    tpoint_effs_mask <- grepl("tpoint", wisp.results$param.names) & grepl("beta", wisp.results$param.names)
    if (any(tpoint_effs_mask)) {
      
      tpoint_effects <- c(bs_fitted_params[,tpoint_effs_mask])
      rates <- seq(min(tpoint_effects), max(tpoint_effects), length.out = 100)
      data <- data.frame(
        rate = rates,
        probability = dnorm(rates, 0, tpoint_shape)
      )
      data$probability <- data$probability / max(data$probability)
      plot_tpoint_effects_struc_stats <- ggplot() +
        geom_histogram(
          data = data.frame(vals = tpoint_effects), aes(x = vals, y = after_stat(ndensity)),
          bins = 75, fill = "skyblue", alpha = 0.5, na.rm = TRUE) +
        geom_line(aes(x = rate, y = probability), data = data, linewidth = 1.2, na.rm = TRUE) +
        labs(title = "Gaussian distribution model of fixed t-point effects", x = "t-point effect", y = "Density") +
        theme_minimal() 
      
      if (verbose) print(plot_tpoint_effects_struc_stats)
      plots.struc_stats[["plot_tpoint_effects_struc_stats"]] <- plot_tpoint_effects_struc_stats
      
    }
    
    return(plots.struc_stats)
    
  }

# Method for printing all child plots on one figure
plot.child.summary <- function(
    wisp.results,
    these.parents = NULL,
    these.childs = NULL,
    log.scale = TRUE,
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
    
    # Take abs and log
    sampled_params <- log(abs(sampled_params) + 1)
    
    # Make long-format data frame
    walks <- data.frame(
      value = c(sampled_params),
      param = rep(1:ncol(sampled_params), each = nrow(sampled_params)),
      step = rep(1:nrow(sampled_params), ncol(sampled_params))
    )
    
    plot.walks <- ggplot(data = walks, 
                   aes(x = step, y = value, group = param, color = param)) + 
      geom_line() + 
      ggtitle("Parameter Walks from MCMC Estimation") +
      theme_minimal()
    print(plot.walks)
    
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



