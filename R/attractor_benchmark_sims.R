# Functions for making attractor-based simulations for benchmarking

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
    data_transformed$coord_x = distx*data$coord_x + (1-distx)* (A[1,1]*data$coord_x + A[1,2]*data$coord_y)
    data_transformed$coord_y = disty*data$coord_y + (1-disty)* (A[2,1]*data$coord_x + A[2,2]*data$coord_y)
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
#' @param seed_data A data frame containing the seed dataset with columns: \code{gene}, \code{coord_x}, \code{coord_y}, and \code{count}.
#' @param n_bins An integer specifying the number of spatial bins to divide the data coordinates into. Default is 100.
#' @param n_replicates An integer specifying the number of replicates to generate for each treatment condition. Default is 4.
#' @param replicate_spatial_scalar A numeric value controlling the amount of spatial variation introduced in each replicate. Default is 0.05.
#' @param min_effect_size A numeric value specifying the minimum effect size for functional spatial effects (FSEs). Default is 0.05. Max positive effect size is always 4, i.e., a 4x change in rate. There is no max to the min effect size, i.e., an effect can drop the rate to zero.
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
        plt_seed <- ggplot(seed_data[mask,], aes(x = coord_x, y = coord_y, color = gene, size = log(count + 1))) +
          geom_point(alpha = 0.5) + 
          ggtitle("Seed data with transcript counts per cell") +
          theme_minimal() +
          guides(size = "none") +
          theme(legend.position = "bottom")
        plt_smoothed <- ggplot(seed_data_smoothed[mask,], aes(x = coord_x, y = coord_y, color = gene, size = log(count + 1))) +
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
          FSE_title <- "with FSE"
          if (!(g %in% FSEs)) FSE_title <- "without FSE"
          mask <- sim_data_ref$gene == genes[g]
          sim_data_plt <- rbind(
            sim_data_ref[mask,],
            sim_data_trt[mask,]
          )
          sim_data_plt$fixedeffect <- c(rep("ref", sum(mask)), rep("trt", sum(mask)))
          plt_FSE <- ggplot(sim_data_plt, aes(x = coord_x, y = coord_y, color = fixedeffect, size = log(count + 1))) +
            geom_point(alpha = 0.5) + 
            ggtitle(paste0(genes[g], ", spatially variable ", FSE_title)) +
            theme_minimal() +
            guides(size = "none") +
            theme(legend.position = "bottom")
          plt <- c(plt, list(plt_FSE))
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


