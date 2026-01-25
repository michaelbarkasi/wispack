# Functions for making attractor-based simulations for benchmarking ####################################################

# Helper function to make replicates
make_replicate <- function(
    data, 
    rate_scalar = 0.5, # number between 0 and 1 
    affine_transform = NULL,
    spatial_scalar = 0.05
  ) {
    data_transformed <- data
    # If none provided, make affine transformation to shift cells around
    if (is.null(affine_transform)) {
      Ax <- matrix(c(1,0,rnorm(1,0,spatial_scalar),1),2,2)  # shear x
      Ay <- matrix(c(1,rnorm(1,0,spatial_scalar),0,1),2,2)  # shear y
      As <- matrix(c(rnorm(1,0,spatial_scalar)+1, 0, 0, rnorm(1,0,spatial_scalar)+1),2,2) # scale
      A <- As %*% Ay %*% Ax
    } else {
      A <- affine_transform
    }
    # Make scales to keep cells roughly in bounds
    x_mid <- diff(range(data$coord_x))/2 + min(data$coord_x)
    distx <- (data$coord_x - x_mid)^2
    distx <- distx/max(distx)
    y_mid <- diff(range(data$coord_y))/2 + min(data$coord_y)
    disty <- (data$coord_y - y_mid)^2
    disty <- disty/max(disty)
    # Applied scaled affine transformation
    data_transformed$coord_x = distx*data$coord_x + (1-distx) * (A[1,1]*data$coord_x + A[1,2]*data$coord_y)
    data_transformed$coord_y = disty*data$coord_y + (1-disty) * (A[2,1]*data$coord_x + A[2,2]*data$coord_y)
    # Poisson resampling of genes with scaled rate
    data_transformed$count <- rpois(nrow(data), data_transformed$count * (2*rate_scalar))
    return(data_transformed)
  }

# Helper function to simulate non-SVGs
smooth_data <- function(
    data
  ) {
    idx_shuffled <- sample(nrow(data), nrow(data), replace = FALSE)
    data$coord_x[idx_shuffled] <- data$coord_x
    data$coord_y[idx_shuffled] <- data$coord_y
    return(data)
  }

# Helper function to simulate SV from an attractor point and spatial effect
induce_SV <- function(
    data_mixed, 
    gene, 
    attractor,    # number between 0 and 1
    effect = 0.5  # number between 0 and 1, to simulate fixed effects on rate
  ) {
    # Get attractor coordinates
    n_rows <- nrow(data_mixed)
    attractor_idx <- as.integer(n_rows * attractor)
    attractor_idx <- min(max(1, attractor_idx), n_rows)
    coord_attractor <- c(data_mixed$coord_x[attractor_idx], data_mixed$coord_y[attractor_idx])
    # Find distances from cells (i.e., transcripts) to attractor points
    gene_mask <- data_mixed$gene == gene
    x_diff <- data_mixed$coord_x[gene_mask] - coord_attractor[1]
    y_diff <- data_mixed$coord_y[gene_mask] - coord_attractor[2]
    d <- sqrt(x_diff^2 + y_diff^2)
    # Normalize distance
    d_norm <- 1 - d/max(d)
    # Add noise to normalization by using it as mean for beta distribution
    shape1 <- 1
    shape2 <- (shape1 - (shape1 * d_norm)) / d_norm
    n_points <- sum(gene_mask)
    d_norm_noise <- rbeta(n_points, shape1, shape2)
    # Scale differences with noisy normalized distance
    x_diff_norm <- x_diff * d_norm_noise
    y_diff_norm <- y_diff * d_norm_noise
    # Apply attraction
    data_attracted <- data_mixed
    data_attracted$coord_x[gene_mask] <- data_mixed$coord_x[gene_mask] - x_diff_norm 
    data_attracted$coord_y[gene_mask] <- data_mixed$coord_y[gene_mask] - y_diff_norm 
    # Apply effect
    data_attracted$count[gene_mask] <- rpois(sum(gene_mask), data_attracted$count[gene_mask] * (2*effect)^2)
    return(data_attracted)
  }

# Helper function to bin data
bin_data <- function(
    data,
    n_bins = 100
  ) { 
    # Bin coordinates
    data$bin_x <- cut(
      data$coord_x,
      breaks = seq(min(data$coord_x), max(data$coord_x), length.out = n_bins + 1),
      include.lowest = TRUE,
      labels = FALSE
    )
    data$bin_y <- cut(
      data$coord_y,
      breaks = seq(min(data$coord_y), max(data$coord_y), length.out = n_bins + 1),
      include.lowest = TRUE,
      labels = FALSE
    )
    # Return binned data
    return(data)
  }

#' Function to extract ground-truth from attractor-based simulation made with the function \code{attractor_simulation}
#' 
#' This function extracts the ground-truth values for functional spatial effect values, replicate random effect values, and indicators for spatially variable genes (SVGs) and genes with functional spatial effects (FSEs) from a simulation object created by the \code{attractor_simulation} function.
#' 
#' @param sim A list object returned by the \code{attractor_simulation} function, containing the simulation data and parameters.
#' @return A list containing the following components: \code{fse_true}, \code{ran_true}, \code{SVGs}, and \code{FSEs}.
#' @export
attractor_simulation_ground_truth <- function(
    sim
  ) {
    # Extract ground-truth values for functional spatial effect values
    fse_true <- sim$effect
    # Extract ground-truth values for replicate random effect values
    ran_true <- sim$replicate_rate_scalars
    # Extract ground-truth values for whether a gene is SV
    SVGs <- rep(FALSE, length(sim$genes))
    names(SVGs) <- sim$genes
    SVGs[sim$SVGs] <- TRUE
    # Extract ground-truth values for whether a gene has FSE
    FSEs <- rep(FALSE, length(sim$genes))
    names(FSEs) <- sim$genes
    FSEs[sim$FSEs] <- TRUE 
    # Return as list
    return(
      list(
        fse_true = fse_true,
        ran_true = ran_true,
        SVGs = SVGs,
        FSEs = FSEs
      )
    )
  }

#' Function to simulate data using attractor-based spatial variability
#' 
#' This function takes a seed dataset and simulates spatially variable genes (SVGs) with functional spatial effects (FSEs) based on attractor points. It generates multiple replicates for reference and treatment conditions, applying random spatial transformations and rate scalars.
#' 
#' @param seed_data A data frame containing the seed dataset with at least the columns: \code{gene}, \code{coord_x}, \code{coord_y}, and \code{count}.
#' @param n_bins An integer specifying the number of spatial bins to divide the data coordinates into. Default is 100.
#' @param n_replicates An integer specifying the number of replicates to generate for each treatment condition. Default is 4.
#' @param replicate_spatial_scalar A numeric value controlling the amount of spatial variation introduced in each replicate for each simulation. Default is 0.05.
#' @param min_effect_size A numeric value specifying the minimum effect size for functional spatial effects (FSEs) in each simulation. Default is 0.05. Max positive effect size is always 4, i.e., a 4x change in rate. There is no max to the min effect size, i.e., an effect can drop the rate to zero.
#' @param print_plots A logical value indicating whether to return plots of the simulation process for each gene. Default is FALSE.
#' @return A list containing the following components: \code{genes}, \code{SVGs_idx}, \code{SVGs}, \code{FSEs_idx}, \code{FSEs}, \code{attractor}, \code{effect}, \code{replicate_rate_scalars}, \code{data}, and \code{plots}.
#' @export
attractor_simulation <- function(
    seed_data, 
    n_bins = 100,
    n_replicates = 4,
    replicate_spatial_scalar = 0.05, 
    min_effect_size = 0.05,
    print_plots = FALSE
  ) {
    
    # Check column names 
    if (!all(c("gene", "coord_x", "coord_y", "count") %in% colnames(seed_data))) {
      stop("seed_data must contain columns: gene, coord_x, coord_y, count")
    }
    
    # Get number and list of genes
    genes <- sort(unique(seed_data$gene))
    n_genes <- length(genes)
    
    # Randomly select genes to be spatially variable 
    # ... will select at least one
    SVGs <- sort(sample(n_genes, sample(n_genes, 1)))
    n_SVGs <- length(SVGs)
    
    # Randomly select SVGs to have FSEs
    # ... may select none
    FSEs <- c()
    n_FSEs <- sample(c(0:n_SVGs), 1)
    if (n_FSEs > 0) FSEs <- sort(sample(SVGs, n_FSEs))
    
    # Select attractors (some may not be used)
    attractor <- runif(n_genes)
    
    # Select effects
    # ... default is 0.5 * 2 = 1, i.e. no effect
    effect <- rep(0.5, n_genes)
    # ... apply effects to FSE genes, ensuring they are at least 5% different from no-effect
    if (n_FSEs > 0) {
      FSeffects <- runif(n_FSEs)
      min_pos_effect <- sqrt(1 + min_effect_size)/2
      min_neg_effect <- sqrt(1 - min_effect_size)/2
      FSeffects[FSeffects >= 0.5 & FSeffects < min_pos_effect] <- min_pos_effect
      FSeffects[FSeffects < 0.5 & FSeffects > min_neg_effect] <- min_neg_effect
      effect[FSEs] <- FSeffects
    }
    
    # Smooth seed data (no spatial variation)
    seed_data_smoothed <- smooth_data(seed_data)
    
    # Initialize variables for reference and treatment data
    sim_data_ref <- seed_data_smoothed
    sim_data_trt <- seed_data_smoothed
    
    # Smooth expression distribution, apply attractors, apply effects
    all_plots <- list()
    for (g in c(1:n_genes)) {
      
      # Print seed data and smoothed gene
      if (print_plots) {
        mask <- seed_data_smoothed$gene == genes[g]
        plt_seed <- ggplot(seed_data[mask,], aes(x = coord_x, y = coord_y, size = log(count + 1))) +
          geom_point(alpha = 0.5) + 
          ggtitle(paste0(genes[g], ", seed data")) +
          theme_minimal() +
          guides(size = "none") +
          theme(legend.position = "bottom")
        plt_smoothed <- ggplot(seed_data_smoothed[mask,], aes(x = coord_x, y = coord_y, size = log(count + 1))) +
          geom_point(alpha = 0.5) + 
          ggtitle(paste0(genes[g], ", smoothed")) +
          theme_minimal() +
          guides(size = "none") +
          theme(legend.position = "bottom")
        plt <- list(plt_seed, plt_smoothed)
      }
      if (any(g == SVGs)) {
        # Use attractor to induce spatial variability
        sim_data_ref <- induce_SV(sim_data_ref, genes[g], attractor[g])
        # Apply effects
        sim_data_trt <- induce_SV(sim_data_trt, genes[g], attractor[g], effect[g])
        # Print spatially variable gene
        if (print_plots) {
          FSE_title <- "w/ FSE"
          if (!(g %in% FSEs)) FSE_title <- "w/o FSE"
          mask <- sim_data_ref$gene == genes[g]
          plt_FSE_ref <- ggplot(sim_data_ref[mask,], aes(x = coord_x, y = coord_y, size = log(count + 1))) +
            geom_point(alpha = 0.5) + 
            ggtitle(paste0(genes[g], " ref, SV ", FSE_title)) +
            theme_minimal() +
            guides(size = "none") +
            theme(legend.position = "bottom")
          plt_FSE_trt <- ggplot(sim_data_trt[mask,], aes(x = coord_x, y = coord_y, size = log(count + 1))) +
            geom_point(alpha = 0.5) + 
            ggtitle(paste0(genes[g], " trt, SV ", FSE_title)) +
            theme_minimal() +
            guides(size = "none") +
            theme(legend.position = "bottom")
          plt <- c(plt, list(plt_FSE_ref, plt_FSE_trt))
        }
      }
      
      # Print
      if (print_plots) all_plots[[genes[g]]] <- plt
    }
    
    # Select replicate rate scalars 
    replicate_rate_scalars <- runif(n_replicates)
    
    # Make replicates and bin data
    rep_names <- paste0("replicate", c(1:n_replicates))
    nrow_rep <- nrow(sim_data_ref)
    ncol_rep <- ncol(sim_data_ref) + 2
    sim_data_ref_reps <- as.data.frame(matrix(NA, nrow = n_replicates * nrow_rep, ncol = ncol_rep))
    sim_data_trt_reps <- as.data.frame(matrix(NA, nrow = n_replicates * nrow_rep, ncol = ncol_rep))
    for (r in c(1:n_replicates)) {
      # Make affine transform for replicate 
      Ax <- matrix(c(1,0,rnorm(1,0,replicate_spatial_scalar),1),2,2)  # shear x
      Ay <- matrix(c(1,rnorm(1,0,replicate_spatial_scalar),0,1),2,2)  # shear y
      As <- matrix(c(rnorm(1,0,replicate_spatial_scalar)+1, 0, 0, rnorm(1,0,replicate_spatial_scalar)+1),2,2) # scale
      A <- As %*% Ay %*% Ax
      # Set index
      idx_start <- (r - 1) * nrow_rep + 1
      idx_end <- r * nrow_rep
      idx <- c(idx_start:idx_end)
      # Make replicate and bin data
      sim_data_ref_reps[idx,] <- bin_data(make_replicate(data = sim_data_ref, rate_scalar = replicate_rate_scalars[r], affine_transform = A), n_bins)
      sim_data_trt_reps[idx,] <- bin_data(make_replicate(data = sim_data_trt, rate_scalar = replicate_rate_scalars[r], affine_transform = A), n_bins)
    }
    colnames(sim_data_ref_reps) <- c(colnames(sim_data_ref), "bin_x", "bin_y") 
    colnames(sim_data_trt_reps) <- c(colnames(sim_data_trt), "bin_x", "bin_y") 
    # Create labels for replicates and fixed effect conditions
    sim_data_ref_reps$replicate <- rep(rep_names, each = nrow_rep)
    sim_data_trt_reps$replicate <- rep(rep_names, each = nrow_rep)
    sim_data_ref_reps$fixedeffect <- "ref"
    sim_data_trt_reps$fixedeffect <- "trt"
    
    # Combine binned data
    sim_data <- rbind(sim_data_ref_reps, sim_data_trt_reps)
    
    # Save effects, adjusting so that no-effect is zero
    effect <- effect - 0.5
    names(effect) <- genes
    replicate_rate_scalars <- replicate_rate_scalars - 0.5
    names(replicate_rate_scalars) <- rep_names
    
    # Combine simulation components into a list and return
    sim <- list(
      genes = genes,
      SVGs_idx = SVGs,
      SVGs = genes[SVGs],
      FSEs_idx = FSEs, 
      FSEs = genes[FSEs],
      attractor = attractor, 
      effect = effect, 
      replicate_rate_scalars = replicate_rate_scalars,
      data = sim_data,
      plots = all_plots
    )
    return(sim)
    
  }

# Functions for modeling attractor simulations with wisp ###############################################################

# Rate effect is consistent across all of space, so, average over all blocks for each gene
extract_rate_effect_attractor_simulation_wisp <- function(model) {
  genes <- model[["grouping.variables"]][["species.lvls"]]
  rate_effect <- rep(NA, length(genes))
  names(rate_effect) <- genes
  for (g in genes) {
    mask <- grepl(paste0("beta_Rt_context_", g, "_trt_X"), names(model[["fitted.parameters"]]))
    rate_effect[g] <- mean(model[["fitted.parameters"]][mask])
  }
  return(rate_effect)
}

# Function to extract mean p-value across blocks for each gene 
extract_pvalue_attractor_simulation_wisp <- function(model) {
  genes <- model[["grouping.variables"]][["species.lvls"]]
  p_values <- rep(NA, length(genes))
  names(p_values) <- genes
  for (g in genes) {
    mask <- grepl(paste0("beta_Rt_context_", g, "_trt_X"), model[["stats"]][["parameters"]]$parameter)
    p_values[g] <- mean(model[["stats"]][["parameters"]]$p.value.adj[mask]) / sum(mask) 
    # ^ ... Values are adjusted for multiple tests, but we know that there is only a single degree of freedom in the rate effect per gene
  }
  return(p_values)
}

# Replicate random noise is consistent across genes, so, average over all genes for each replicate
extract_random_effect_attractor_simulation_wisp <- function(model) {
  replicates <- model[["grouping.variables"]][["ran.lvls"]]
  random_effect <- rep(NA, length(replicates))
  names(random_effect) <- replicates
  for (r in replicates) {
    mask <- grepl(paste0("wfactor_rate_", r), names(model[["fitted.parameters"]]))
    random_effect[r] <- mean(model[["fitted.parameters"]][mask])
  }
  random_effect <- random_effect[-c(1)]
  return(random_effect)
} 

# Extract model results from wisp
extract_wisp_attractor_simulation_results <- function(model) {
  fse_est <- extract_rate_effect_attractor_simulation_wisp(model)
  ran_est <- extract_random_effect_attractor_simulation_wisp(model)
  fse_pvalues <- extract_pvalue_attractor_simulation_wisp(model)
  return(
    list(
      fse_est = fse_est,
      ran_est = ran_est,
      fse_pvalues = fse_pvalues
    )
  )
}

#' Function to model simulation with wisp
#' 
#' This function fits a wisp model to an attractor-based simulation dataset and extracts the estimated functional spatial effects, replicate random effects, and p-values for functional spatial effects. It then compiles these results along with the ground-truth values from the simulation for comparison.
#' 
#' @param sim A list object returned by the \code{attractor_simulation} function, containing the simulation data and parameters.
#' @param sim_num An integer specifying the simulation number for setting the random seed and labeling results.
#' @return A data frame containing the estimated values (\code{est}), ground-truth values (\code{true}), parameter types (\code{param}), parameter-associated gene or replicate (\code{id}), modeling method used (\code{method}), and simulation number (\code{sim}). 
model_attractor_simulation_wisp <- function(
    sim, 
    sim_num,
    bs_num = 1e3,
    max_fork = 1
  ) {
    
    # Fit model
    model <- wisp(
      count.data = sim$data,
      variables = list(
        bin = "bin_x", 
        count = "count",
        species = "gene",
        ran = "replicate",
        fixedeffects = c("fixedeffect")
      ),
      fit_only = FALSE,
      MCMC.settings = list(MCMC.steps = 0, MCMC.burnin = 0),
      bootstraps.num = bs_num,
      max.fork = max_fork,
      plot.settings = list(print.plots = FALSE),
      verbose = FALSE,
      ran.seed = sim_num
    )
    
    # Extract model results for comparing to ground truth
    MR <- extract_wisp_attractor_simulation_results(model)
    
    # Extract ground-truth from simulation
    GT <- attractor_simulation_ground_truth(sim)
    
    # Compile results in named vector and return
    return(
      data.frame(
        est = c(MR$fse_est, MR$ran_est, MR$fse_pvalues),
        true = c(GT$fse_true, GT$ran_true, GT$FSEs),
        param = c(
          rep("rate_effect", length(MR$fse_est)), 
          rep("random_effect", length(MR$ran_est)), 
          rep("FSE", length(MR$fse_pvalues))
        ),
        id = c(
          names(MR$fse_est), 
          names(MR$ran_est), 
          names(MR$fse_pvalues)
        ),
        method = "wisp",
        sim = sim_num
      )
    )
    
  }

# Functions to run and analyze benchmarks with attractor simulations ###################################################

#' Function to run benchmarks on attractor-based simulations
#' 
#' This function benchmarks a provided list of modeling methods on attractor-based simulation datasets. It simulates data multiple times, fits each modeling method to the simulated data, and records the results and computation times for each step.
#' 
#' @param seed_data A data frame containing the seed dataset with columns: \code{gene}, \code{coord_x}, \code{coord_y}, and \code{count}. Will be given to \code{attractor_simulation} to generate simulations.
#' @param n_sims An integer specifying the number of simulations to run. Default is 100. 
#' @param n_bins An integer specifying the number of spatial bins to divide the data coordinates into for each simulation. Default is 100.
#' @param n_replicates An integer specifying the number of replicates to generate for each treatment condition in each simulation. Default is 4.
#' @param replicate_spatial_scalar A numeric value controlling the amount of spatial variation introduced in each replicate for each simulation. Default is 0.05.
#' @param min_effect_size A numeric value specifying the minimum effect size for functional spatial effects (FSEs) in each simulation. Default is 0.05. Max positive effect size is always 4, i.e., a 4x change in rate. There is no max to the min effect size, i.e., an effect can drop the rate to zero.
#' @param modeling_functions A named list of modeling functions to benchmark. Each function should take at least the arguments \code{sim} (the simulation object produced by \code{attractor_simulation}) and \code{sim_num} (the simulation number). Default is a list containing only the \code{model_attractor_simulation_wisp} function. Functions provided in this list should return a dataframe with columns \code{est}, \code{true}, \code{param}, \code{id}, \code{method}, and \code{sim}. The \code{true} column should contain, for each simulation and each applicable parameter, the ground-truth value returned by the \code{attractor_simulation_ground_truth} function. The \code{est} column should contain the corresponding estimated value from the modeling function. The \code{param} column should contain the name of the parameter (one of "rate_effect", "random_effect", "FSE", or "SVG"). The \code{id} column should contain the name of the gene or replicate associated with the parameter. The \code{method} column should contain the name of the modeling method used. The \code{sim} column should contain the simulation number.
#' @param modeling_function_args A named list of lists, where each sub-list contains additional arguments to pass to the corresponding modeling function in \code{modeling_functions}. Default is a list containing arguments for the \code{model_attractor_simulation_wisp} function.
#' @return A list containing two components: \code{results}, a data frame compiling the results from all simulations and modeling methods; and \code{times}, a data frame recording the computation times for data simulation and each modeling method for each simulation.
#' @export
run_attractor_sim_benchmarks <- function(
    seed_data, 
    n_sims = 100,
    n_bins = 100,
    n_replicates = 4,
    replicate_spatial_scalar = 0.05,
    min_effect_size = 0.05,
    modeling_functions = list(wisp = model_attractor_simulation_wisp),
    modeling_function_args = list(wisp = list(bs_num = 1e3, max_fork = 1))
  ) {
    
    results <- data.frame()
    times <- matrix(NA, nrow = n_sims, ncol = length(modeling_functions) + 1)
    colnames(times) <- c("data_simulation", names(modeling_functions))
    s <- 1
    error_counter <- 0
    while (s <= n_sims) {
      
      cat("\n\n----------------------")
      cat("\nRunning simulation ", s, "/", n_sims)
      
      # Simulate data and extract count data for modeling
      t_stim_start <- Sys.time()
      sim <- attractor_simulation(seed_data)
      dsim <- Sys.time() - t_stim_start
      units(dsim) <- "secs"
      times[s, "data_simulation"] <- dsim
      cat("\n\tdata sim time:", round(dsim,3), "s") 
      
      # Model with each modeling function
      error_flag <- FALSE
      results_s <- data.frame()
      for (mf in names(modeling_functions)) {
        cat("\n\nModeling with ", mf, "...")
        t_model_start <- Sys.time()
        model_func <- modeling_functions[[mf]]
        model_results <- tryCatch(
          {
            if (mf %in% names(modeling_function_args)) {
              model_args <- modeling_function_args[[mf]]
              do.call(model_func, c(list(sim = sim, sim_num = s), model_args))
            } else {
              model_func(sim = sim, sim_num = s)
            }
          },
          error = function(e) {
            error_flag <<- TRUE
            NULL
          }
        )
        dmodel <- Sys.time() - t_model_start
        units(dmodel) <- "secs"
        times[s, mf] <- dmodel
        cat("\n\t", mf, "time:", round(dmodel,3), "s")
        results_s <- rbind(results_s, model_results)
      }
      
      if (error_flag) {
        error_counter <- error_counter + 1
        cat("\n\nError encountered during modeling, number:", error_counter, "Rerunning simulation", s, "\n")
      } else {
        results <- rbind(results, results_s)
        s <- s + 1
      }
      
    }
    
    return(
      list(
        results = results,
        times = as.data.frame(times)
      )
    )
    
  }

#' Function to analyze benchmark results from attractor-based simulations
#' 
#' This function analyzes the results from benchmarking modeling methods on attractor-based simulation datasets. It computes the performance metrics: correlation coefficients for rate and random effect parameters, false positive rates (FPR), false discovery rates (FDR), and power for spatially variable gene (SVG) and functional spatial effect (FSE) parameters.
#' 
#' @param results The output of the \code{run_attractor_sim_benchmarks} function, either in its native form as a list, or converted into a data frame. 
#' @param sig_thresh A named list specifying the significance threshold for each modeling method when computing FPR, FDR, and power. Default is a list with a single entry for the "wisp" method set to 0.05. If a value is not specified for a method found in results, will use 0.05.
#' @return A data frame summarizing the performance metrics for each modeling method.
#' @export
analyze_attractor_sim_benchmarks <- function(
    results,
    sig_thresh = list(wisp = 0.05)
  ) {
   
    # Check input
    if (class(results) == "data.frame") {
      res_cols <- grepl("results.", colnames(results))
      res_col_names <- gsub("results.", "", colnames(results)[res_cols])
      time_cols <- grepl("times.", colnames(results))
      time_col_names <- gsub("times.", "", colnames(results)[time_cols])
      if (all(res_cols | time_cols)) {
        times <- results[,time_cols]
        colnames(times) <- time_col_names
        results <- results[,res_cols]
        colnames(results) <- res_col_names
      } else {
        stop("results data.frame not in recognized format")
      }
    } else if (class(results) == "list") {
      if (length(results) == 2 && all(names(results) == c("results","times"))) {
        times <- results$times
        results <- results$results
      } else {
        stop("results list not in recognized format")
      }
    } else {
      stop("results not in recognized format")
    }
    
    # Initialize summary dataframe
    methods <- sort(unique(results$method))
    metrics <- c(
      "cor_rate", 
      "cor_ran", 
      "FPR_svg", 
      "FDR_svg", 
      "power_svg", 
      "FPR_fse", 
      "FDR_fse", 
      "power_fse",
      "mean_sec_to_model"
    )
    summary <- as.data.frame(matrix(NA, nrow = length(metrics), ncol = length(methods)))
    colnames(summary) <- methods
    rownames(summary) <- metrics
   
    for (m in methods) {
      
      sigt <- 0.05
      if (m %in% names(sig_thresh)) sigt <- sig_thresh[[m]]
      
      m_mask <- results$method == m
      
      # Compute correlation for rate effect parameter, true vs est
      if (any("rate_effect" == results$param[m_mask])) {
        r_m_mask <- m_mask & results$param == "rate_effect"
        summary["cor_rate", m] <- cor(results[r_m_mask,"true"], results[r_m_mask,"est"], method = "pearson")
      }
      # Compute correlation for random effect parameter, true vs est
      if (any("random_effect" == results$param[m_mask])) {
        r_m_mask <- m_mask & results$param == "random_effect"
        summary["cor_ran", m] <- cor(results[r_m_mask,"true"], results[r_m_mask,"est"], method = "pearson")
      }
      # Compute FPR, FDR, power for SVG parameter
      if (any("SVG" == results$param[m_mask])) {
        r_m_mask <- m_mask & results$param == "SVG"
        FP <- sum(results[r_m_mask & results$true == FALSE, "est"] < sigt)
        TN <- sum(results[r_m_mask & results$true == FALSE, "est"] >= sigt)
        TP <- sum(results[r_m_mask & results$true == TRUE, "est"] < sigt)
        FN <- sum(results[r_m_mask & results$true == TRUE, "est"] >= sigt)
        summary["FPR_svg", m] <- FP / (FP + TN)
        summary["FDR_svg", m] <- FP / (FP + TP)
        summary["power_svg", m] <- TP / (TP + FN)
      }
      # Compute FPR, FDR, power for FSE parameter
      if (any("FSE" == results$param[m_mask])) {
        r_m_mask <- m_mask & results$param == "FSE"
        FP <- sum(results[r_m_mask & results$true == FALSE, "est"] < sigt)
        TN <- sum(results[r_m_mask & results$true == FALSE, "est"] >= sigt)
        TP <- sum(results[r_m_mask & results$true == TRUE, "est"] < sigt)
        FN <- sum(results[r_m_mask & results$true == TRUE, "est"] >= sigt)
        summary["FPR_fse", m] <- FP / (FP + TN)
        summary["FDR_fse", m] <- FP / (FP + TP)
        summary["power_fse", m] <- TP / (TP + FN)
      }
      # Compute mean modeling time
      summary["mean_sec_to_model", m] <- mean(times[,m], na.rm = TRUE)
      
    }
    
    return(summary)
  }
