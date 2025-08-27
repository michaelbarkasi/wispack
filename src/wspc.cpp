
// wspc.cpp
#include "wspc.h"

// WSPmm class and methods *********************************************************************************************

/*
 * Object class to hold and fit Warped Sigmoidal Poisson-Process Mixed-Effect Model (WSPmm, aka "WiSP") model. 
 */

// Constructor
wspc::wspc(
    Rcpp::DataFrame count_data,
    Rcpp::List settings,
    bool verbose
  ) { 
    
    /*
     * Input data should be an R data frame with the columns specified below, and rows 
     *  as token count observations. For example, a spatial axis of measure could be 
     *  divided into 100 bins, each containing some number of cells with various transcript 
     *  counts. Each cell and its transcript count would be a different token observation 
     *  for its bin. 
     */
    
    // Extract settings
    buffer_factor = sdouble((double)settings["buffer_factor"]);
    ctol = (double)settings["ctol"];
    max_penalty_at_distance_factor = sdouble((double)settings["max_penalty_at_distance_factor"]);
    LROcutoff = (double)settings["LROcutoff"];
    LROwindow_factor = (double)settings["LROwindow_factor"];
    rise_threshold_factor = (double)settings["rise_threshold_factor"];
    max_evals = (int)settings["max_evals"];
    rng_seed = (unsigned int)settings["rng_seed"];
    warp_precision = (sdouble)settings["warp_precision"];
    inf_warp = (sdouble)settings["inf_warp"];
    model_settings = Rcpp::clone(settings);
    
    // Report warp_inf 
    const sdouble eps_ = std::numeric_limits<double>::epsilon(); // machine epsilon
    if (verbose) {
      vprint_header("Infinity handling");
      vprint("machine epsilon: (eps_): ", eps_);
      vprint("pseudo-infinity (inf_): ", inf_);
      vprint("warp_precision: ", warp_precision);
      vprint("implied pseudo-infinity for unbounded warp (inf_warp): ", inf_warp);
    }
    
    // Check structure of input data
    CharacterVector col_names = count_data.names();
    CharacterVector required_cols = CharacterVector::create("count", "bin", "parent", "child", "ran");
    int n_cols = col_names.size();
    int r_cols = required_cols.size();
    for (int i = 0; i < r_cols; i++) {
      if (col_names[i] != required_cols[i]) {
        Rcpp::stop("Input data is missing required column (or out of order): " + required_cols[i]);
      }
    }
    vprint("Data structure check passed", false);
    
    // Save tokenized count column before collapsing to sums 
    count_tokenized = to_sVec(Rcpp::as<NumericVector>(count_data["count"]));
    vprint("Saving tokenized count", false);
    
    // Find max bins and set warp bounds
    vprint_header("Extracting variables and data information", verbose);
    bin_num = smax(to_sVec(Rcpp::as<NumericVector>(count_data["bin"])));
    warp_bounds.resize(3); 
    for (int i = 0; i < 3; i++) {warp_bounds[i] = inf_warp;}  // initialize to infinity (no bounds)
    warp_bounds_idx.names() = CharacterVector::create("Rt", "tslope", "tpoint");
    int warp_bound_tpoint_idx = warp_bounds_idx["tpoint"]; 
    warp_bounds[warp_bound_tpoint_idx] = bin_num;  // set tpoint bound to max bin
    vprint("Found max bin: " + std::to_string(bin_num.val()), verbose);
    
    // Extract fixed effects 
    int n_fix = n_cols - r_cols;
    fix_names = CharacterVector(n_fix);                                 // names of fixed effect variables 
    fix_ref = CharacterVector(n_fix);                                   // reference level for each fixed effect
    fix_lvls.resize(n_fix);                                             // levels for each fixed effect
    fix_trt.resize(n_fix);                                              // treatments for each fixed effect
    for (int i = 0; i < n_fix; i++) {
      fix_names[i] = col_names[i + r_cols];
      CharacterVector lvls = Rcpp::sort_unique(Rcpp::as<CharacterVector>(count_data[i + r_cols]));
      if (lvls.size() < 2) {
        Rcpp::stop("Fixed effect " + fix_names[i] + " has less than 2 levels.");
      }
      fix_lvls[i] = lvls;
      fix_ref[i] = lvls[0];                                             // assume first level is reference
      fix_trt[i] = lvls[Rcpp::Range(1, lvls.size() - 1)]; 
    }
    vprint("Fixed effects:", verbose);
    vprintV(fix_names, verbose);
    vprint("Ref levels:", verbose);
    vprintV(fix_ref, verbose);
    
    // Create all possible treatment combinations 
    CharacterVector ref_lvl = {"ref"};                                  // define a baseline (i.e., reference) level
    treatment_components = make_treatments(fix_trt);                    // make components of all possible treatment levels
    treatment_components.insert(treatment_components.begin(), ref_lvl); // add "ref" to represent reference level
    treatment_num = treatment_components.size();                        // grab number of treatments
    treatment_lvls = CharacterVector(treatment_num);                    // resize variable to hold treatment level names
    for (int t = 0; t < treatment_num; t++) {                           // collapse components into level names
      treatment_lvls[t] = Rcpp::collapse(treatment_components[t]);
    }
    vprint("Created treatment levels:", verbose);
    vprintV(treatment_lvls, verbose);
    
    // Create matrix to translate between treatment levels and fixed-effects columns
    // ... for making summed-count fixed-effect column 
    // ... each cell should contain the value of the treatment column n_fix, for the given treatment level
    CharacterMatrix effects_rows(treatment_num, n_fix);
    for (int tr = 0; tr < treatment_num; tr++) {
      CharacterVector trt_given = treatment_components[tr];
      for (int f = 0; f < n_fix; f++) {
        CharacterVector lvls = fix_lvls[f];
        effects_rows(tr, f) = lvls[0]; 
        // ^ ... Assume it's the reference level
        for (String l : lvls) {
          // If treatment level found, replace
          if (any_true(eq_left_broadcast(trt_given, l))) {effects_rows(tr, f) = l;}
        }
      }
    }
    vprint("Created treatment-to-fix translation matrix", false);
    
    // Pre-compute weight-matrix rows 
    // ... for making weights matrix
    weight_rows.resize(treatment_num, treatment_num);
    weight_rows.setOnes();
    for (int tr = 0; tr < treatment_num; tr++) {
      CharacterVector trt_given = treatment_components[tr];
      for (int tc = 0; tc < treatment_num; tc++) {
        CharacterVector trt_testing = treatment_components[tc];
        for (String trt : trt_testing) {
          if(!any_true(eq_left_broadcast(trt_given, trt)) && trt != "ref") {weight_rows(tr, tc) = 0.0;}
          // ^ ... ref level must always have weight 1
        }
      }
    }
    vprint("Pre-computed weight matrix rows", false);
    
    // Extract grouping variables 
    parent_lvls = Rcpp::sort_unique(Rcpp::as<CharacterVector>(count_data["parent"]));
    child_lvls = Rcpp::sort_unique(Rcpp::as<CharacterVector>(count_data["child"]));
    ran_lvls = Rcpp::sort_unique(Rcpp::as<CharacterVector>(count_data["ran"]));
    // ... add "none" to represent no random effect (reference level)
    ran_lvls.push_front("none"); 
    // ... print extracted grouping variables
    vprint("Parent grouping variables:", verbose);
    vprintV(parent_lvls, verbose);
    vprint("Child grouping variables:", verbose);
    vprintV(child_lvls, verbose);
    vprint("Random-effect grouping variables:", verbose);
    vprintV(ran_lvls, verbose);
    
    // Temporarily extract tokenized-count columns 
    IntegerVector binT = Rcpp::as<IntegerVector>(count_data["bin"]); 
    CharacterVector parentT = Rcpp::as<CharacterVector>(count_data["parent"]);
    CharacterVector childT = Rcpp::as<CharacterVector>(count_data["child"]);
    CharacterVector ranT = Rcpp::as<CharacterVector>(count_data["ran"]);
    vprint("Extracted tokenized count columns", false);
    
    // Create summed count data, size constants
    int idx = 0;
    int n_ran = ran_lvls.size();
    int n_child = child_lvls.size();
    int n_parent = parent_lvls.size();
    int bin_num_i = (int)bin_num.val();
    n_count_rows = bin_num_i * n_parent * n_child * n_ran * treatment_num;
    count_row_nums = Rcpp::seq(0, n_count_rows - 1);
    vprint("Total rows in summed count data table: " + std::to_string(n_count_rows), verbose);
    
    // Create summed count data rows, initializations
    count.resize(n_count_rows);
    bin.resize(n_count_rows);
    parent = CharacterVector(n_count_rows);
    child = CharacterVector(n_count_rows);
    ran = CharacterVector(n_count_rows);
    treatment = CharacterVector(n_count_rows);
    weights.resize(n_count_rows, treatment_num);
    vprint("Initialized columns for summed count data", false);
    
    // Initiate count indexes 
    int idx_mcu = 0;
    idx_mc_unique = IntegerVector(n_count_rows/bin_num_i);
    token_pool.resize(n_count_rows);
    count_not_na_mask = LogicalVector(n_count_rows);
    count_not_na_mask.fill(false);
    vprint("Number of rows with unique model components: " + std::to_string(idx_mc_unique.size()), verbose);
    
    // Pre-compute bin masks
    vprint("Pre-computing bin masks", false);
    LogicalMatrix bin_masks(binT.size(), bin_num_i);
    for (int b = 0; b < bin_num_i; b++) { 
      bin_masks.column(b) = eq_left_broadcast(binT, b + 1);
    }
    
    // Create summed count data columns and weight matrix
    vprint_header("Creating summed-count data columns", verbose); 
    LogicalVector nan_mask = !Rcpp::is_na(to_NumVec(count_tokenized));
    for (int r = 0; r < n_ran; r++) {
      LogicalVector ran_mask = eq_left_broadcast(ranT, ran_lvls[r]) & nan_mask;
      
      for (int p = 0; p < n_parent; p++) {
        LogicalVector parent_mask = ran_mask & eq_left_broadcast(parentT, parent_lvls[p]);
        
        for (int c = 0; c < n_child; c++) {
          LogicalVector child_mask = parent_mask & eq_left_broadcast(childT, child_lvls[c]);
          
          for (int t = 0; t < treatment_num; t++) {
            LogicalVector treatment_mask = Rcpp::clone(child_mask);
            for (int f = 0; f < n_fix; f++) {
              treatment_mask = treatment_mask & eq_left_broadcast(Rcpp::as<CharacterVector>(count_data[f + r_cols]), effects_rows(t, f));
            }
            
            // Save mc-unique index
            idx_mc_unique[idx_mcu] = idx; 
            idx_mcu++;
            
            for (int b = 0; b < bin_num_i; b++) {
              
              // Build columns 
              count(idx) = 0.0;
              bin(idx) = static_cast<sdouble>(b + 1.0);
              parent(idx) = parent_lvls[p];
              child(idx) = child_lvls[c];
              ran(idx) = ran_lvls[r];
              treatment(idx) = treatment_lvls[t];
              weights.row(idx) = weight_rows.row(t);
              // ^ ... weights (and weight_rows) is a matrix saying whether the effect from a given treatment level should apply when computing the effect of another treatment level.
              
              // Find token pool
              LogicalVector token_mask = treatment_mask & bin_masks.column(b);
              
              // Sum count 
              IntegerVector token_pool_idx = Rwhich(token_mask);
              if (token_pool_idx.size() > 0 || ran_lvls[r] == "none") {
                token_pool[idx] = token_pool_idx; 
                // ^ ... save for bootstrap resampling
                for (int rw : token_pool_idx) {
                  count(idx) += count_tokenized[rw];
                }
                count_not_na_mask(idx) = true;
              } else {
                count(idx) = stan::math::NOT_A_NUMBER;
              }
              
              // Increment index
              idx++; 
              
            }
            
          }
        }
      }
      
      vprint("Random level " + std::to_string(r) + ", " + std::to_string(r + 1) + "/" + std::to_string(n_ran) + " complete", verbose);
      
    }
    
    // Extract idx from count_not_na_mask
    count_not_na_idx = Rwhich(count_not_na_mask);
    vprint("Extracted non-NA indexes", false);
    
    // Make extrapolation pool and extrapolate "none" rows
    vprint_header("Making extrapolation pool", verbose);
    extrapolation_pool.resize(count.size());
    extrapolation_pool = make_extrapolation_pool(bin, count, parent, child, ran, treatment, verbose); 
    
    // Extrapolate "none" rows
    vprint_header("Making initial parameter estimates", verbose);
    count = extrapolate_none(count, ran, extrapolation_pool, true);
    vprint("Extrapolated 'none' rows", verbose);
    
    // Take log of observed counts 
    count_log.resize(n_count_rows); 
    for (int r = 0; r < n_count_rows; r++) {
      count_log(r) = slog(count(r) + 1.0);
    }
    vprint("Took log of observed counts", verbose);
    
    // Compute running and filter window sizes for LRO change-point detection
    ws = static_cast<int>(std::round(LROwindow_factor * (double)bin_num_i * buffer_factor.val()));
    
    // Compute tpoint buffer
    tpoint_buffer = bin_num * buffer_factor;
    double tpoint_buffer_d = tpoint_buffer.val();
    int tpoint_buffer_int = (int)tpoint_buffer_d;
    if (tpoint_buffer_d != tpoint_buffer_int) {tpoint_buffer_int++;} 
    tpoint_buffer = (sdouble)tpoint_buffer_int;
    
    // Estimate gamma dispersion of raw counts
    List gamma_ests = compute_gamma_dispersion(
      count,
      count_not_na_mask,
      parent,
      child,
      parent_lvls,
      child_lvls
    );
    gamma_dispersion = to_sMat(Rcpp::as<NumericMatrix>(gamma_ests["gamma_dispersion"]));
    gd_child_idx = Rcpp::as<IntegerVector>(gamma_ests["gd_child_idx"]);
    gd_parent_idx = Rcpp::as<IntegerVector>(gamma_ests["gd_parent_idx"]);
    vprint("Estimated gamma dispersion of raw counts", verbose);
    
    // Estimate degree of each parent-child combination at baseline using LRO change-point detection 
    List cp_ests = estimate_change_points(
      bin,
      count_log,
      count_not_na_mask, 
      tpoint_buffer.val(),
      ws,
      bin_num_i,
      LROcutoff,
      parent,
      child,
      ran,
      treatment,
      parent_lvls,
      child_lvls,
      ran_lvls,  
      treatment_lvls,
      treatment_components
    );
    degMat = Rcpp::as<IntegerMatrix>(cp_ests["degMat"]);                              // matrix of degrees of each parent-child combination
    found_cp_list = Rcpp::as<List>(cp_ests["found_cp_list"]);                         // list of found change points for each parent-child combination
    found_cp_trt_list = Rcpp::as<List>(cp_ests["found_cp_trt_list"]);                 // list of found change points for each parent-child combination, averaged across treatments
    vprint("Estimated change points", verbose);
    
    // Find average log counts for each parent-child combination
    count_log_avg_mat_list = find_count_log_means(
      bin,
      count_log,
      count_not_na_mask,
      bin_num_i, 
      parent,
      child,
      treatment,
      parent_lvls,
      child_lvls,
      treatment_lvls,
      treatment_components
    );
    vprint("Found average log counts for each parent-child combination", verbose);
    
    // Find initial parameter estimates for fixed-effect treatments 
    List init_params = estimate_initial_parameters(
      count, 
      count_not_na_mask,
      treatment_num,
      rise_threshold_factor,
      min_initialization_slope,
      weight_rows, 
      mc_list,
      parent,
      child,
      parent_lvls,
      child_lvls,
      degMat,
      found_cp_list, 
      found_cp_trt_list,
      count_log_avg_mat_list 
    );
    List ref_values = Rcpp::as<List>(init_params["ref_values"]);       // reference values for each parent-child combination
    List RtEffs = Rcpp::as<List>(init_params["RtEffs"]);               // rate effects for each parent-child combination
    List tpointEffs = Rcpp::as<List>(init_params["tpointEffs"]);       // tpoint effects for each parent-child combination
    List tslopeEffs = Rcpp::as<List>(init_params["tslopeEffs"]);       // tslope effects for each parent-child combination
    vprint("Estimated initial parameters for fixed-effect treatments", verbose);
    
    // Build default fixed-effects matrices in shell
    List beta = build_beta_shell(mc_list, treatment_lvls, parent_lvls, child_lvls, ref_values, RtEffs, tpointEffs, tslopeEffs, degMat);
    vprint("Built initial beta (ref and fixed-effects) matrices", verbose); 
    
    // Initialize random effect warping factors 
    List wfactors = make_initial_random_effects(
      wfactors_names,
      ran_lvls.size(),
      child_lvls.size()
    );
    vprint("Initialized random effect warping factors", verbose);
    
    // Make and map parameter vector
    List params = make_parameter_vector(
      beta, wfactors,
      parent_lvls, child_lvls, ran_lvls,
      mc_list, 
      treatment_lvls,
      degMat
    );
    vprint("Made and mapped parameter vector", verbose);
    
    // Extract parameter vector information
    fitted_parameters = params["param_vec"];
    param_names = params["param_names"];
    param_wfactor_rate_idx = params["param_wfactor_rate_idx"];
    param_wfactor_point_idx = params["param_wfactor_point_idx"]; 
    param_wfactor_slope_idx = params["param_wfactor_slope_idx"];
    param_beta_Rt_idx = params["param_beta_Rt_idx"];
    param_beta_tslope_idx = params["param_beta_tslope_idx"];
    param_beta_tpoint_idx = params["param_beta_tpoint_idx"];
    param_baseline_idx = params["param_baseline_idx"];
    beta_idx = params["beta_idx"];
    wfactor_idx = params["wfactor_idx"];
    if (verbose) {vprint("Number of parameters: ", (int)fitted_parameters.size());}
    
    // Construct grouping variable ids as indexes for warping factor matrices, by count row
    gv_ran_idx = IntegerVector(n_count_rows);
    gv_fix_idx = IntegerVector(n_count_rows);
    for (int r = 0; r < n_count_rows; r++) {
      CharacterVector::iterator it_ran = std::find(ran_lvls.begin(), ran_lvls.end(), ran[r]);
      CharacterVector::iterator it_fix = std::find(child_lvls.begin(), child_lvls.end(), child[r]);
      gv_ran_idx[r] = std::distance(ran_lvls.begin(), it_ran);
      gv_fix_idx[r] = std::distance(child_lvls.begin(), it_fix);
    }
    vprint("Constructed grouping variable IDs", verbose);
    
    // Compute size of the parameter boundary vector 
    for (int r : idx_mc_unique) {
      // Grab degree for this row
      int c_num = Rwhich(eq_left_broadcast(child_lvls, child[r]))[0];
      int p_num = Rwhich(eq_left_broadcast(parent_lvls, parent[r]))[0];
      int deg = degMat(c_num, p_num);
      if (deg > 0){
        // Add slots for the boundary distance at each tpoint
        boundary_vec_size += deg + 1;
        // Add slots for the tslope boundary distance at each tslope
        boundary_vec_size += deg;
        // Add slots for the boundary distance at each block rate
        boundary_vec_size += deg + 1;
      } else {
        // Add one slot for the single block rate value
        boundary_vec_size++;
      } 
    }
    vprint("Computed size of boundary vector: " + std::to_string(boundary_vec_size), verbose);
    
    // Initialize list to hold results from model fit
    optim_results = List::create(
      _["fitted_parameters"] = NumericVector(), 
      _["penalized_neg_loglik"] = NA_REAL,
      _["neg_loglik"] = NA_REAL, 
      _["success_code"] = NA_INTEGER,
      _["num_evals"] = NA_INTEGER,
      _["bs_times"] = NumericVector()
    );
    
  }

// Destructor
wspc::~wspc() {}

// R copy 
wspc wspc::clone() const {
  wspc this_copy = wspc(*this);
  return this_copy;
}

// Clear stan
void wspc::clear_stan_mem() {
    
    // Temporarily save stan variable values 
    double dbin_num = bin_num.val();
    double dbuffer_factor = buffer_factor.val();
    double dtpoint_buffer = tpoint_buffer.val();
    double dwarp_precision = warp_precision.val();
    double dinf_warp = inf_warp.val();
    double dmax_penalty_at_distance_factor = max_penalty_at_distance_factor.val();
    dVec dbin = to_dVec(bin);
    dVec dcount = to_dVec(count);
    dVec dcount_log = to_dVec(count_log);
    dVec dcount_tokenized = to_dVec(count_tokenized);
    dVec dbp_coefs = to_dVec(bp_coefs);
    dVec dwarp_bounds = to_dVec(warp_bounds);
    NumericMatrix Numweights = to_NumMat(weights);
    NumericMatrix Numweight_rows = to_NumMat(weight_rows);
    NumericMatrix Numgamma_dispersion = to_NumMat(gamma_dispersion);
    
    // Recover memory from stan
    stan::math::recover_memory();
    
    // Re-assign stan variables
    bin_num = (sdouble)dbin_num;
    buffer_factor = (sdouble)dbuffer_factor;
    tpoint_buffer = (sdouble)dtpoint_buffer;
    warp_precision = (sdouble)dwarp_precision;
    inf_warp = (sdouble)dinf_warp;
    max_penalty_at_distance_factor = (sdouble)dmax_penalty_at_distance_factor;
    bin = to_sVec(dbin);
    count = to_sVec(dcount);
    count_log = to_sVec(dcount_log);
    count_tokenized = to_sVec(dcount_tokenized);
    bp_coefs = to_sVec(dbp_coefs);
    warp_bounds = to_sVec(dwarp_bounds);
    weights = to_sMat(Numweights);
    weight_rows = to_sMat(Numweight_rows);
    gamma_dispersion = to_sMat(Numgamma_dispersion);
    
  }

/*
 * *************************************************************************
 * Computing predicted values from parameters
 */

// Compute model component values for rows of summed count data
sVec wspc::compute_warped_mc(
    const String& mc,          // Model component for which to compute values
    const int& r,              // Row of summed count data for which to compute model component 
    const sVec& parameters,    // Parameters to use in computation
    const sdouble& wf          // Warping factor to apply 
  ) const {
    
    // Extract the parameter vector indexes for the current rate row, beta matrices
    List beta_idx_mc = Rcpp::as<List>(beta_idx[mc]);
    List beta_idx_mc_prt = Rcpp::as<List>(beta_idx_mc[(String)parent[r]]);
    IntegerVector betas_mc_idx = beta_idx_mc_prt[Rcpp::as<std::string>(child[r])];
    
    // Grab degree
    int c_num = Rwhich(eq_left_broadcast(child_lvls, child[r]))[0];
    int p_num = Rwhich(eq_left_broadcast(parent_lvls, parent[r]))[0];
    int deg = degMat(c_num, p_num);
    int block_num = deg; 
    if (mc == "Rt") {block_num++;}
    
    // Re-roll beta matrices
    int idx_r = 0;
    sMat betas_mc_r(treatment_num, block_num);
    if (block_num > 0) {
      for (int t = 0; t < block_num; t++) {
        for (int i = 0; i < treatment_num; i++) {
          betas_mc_r(i, t) = parameters(betas_mc_idx[idx_r]);
          idx_r++;
        }
      } 
    }
    
    // Extract weight matrix row
    sVec weight_row = weights.row(r).transpose();
    
    // Compute unwarped model component value 
    sVec mc_vec(block_num);
    mc_vec.setZero();
    for (int t = 0; t < treatment_num; t++) {
      sdouble wi = weight_row(t);
      sVec betai = betas_mc_r.row(t);
      for (int bt = 0; bt < block_num; bt++) {
        mc_vec(bt) += wi*betai(bt);
      }
    }
    
    // Apply warping factor 
    int warp_bound_idx = warp_bounds_idx[mc];
    for (int bt = 0; bt < block_num; bt++) {
      mc_vec(bt) = warp_mc(mc_vec(bt), warp_bounds[warp_bound_idx], wf);
    }
    
    // Send out
    return mc_vec;
    
  }

// Predict log of rates
sVec wspc::predict_rates(
    const sVec& parameters,
    const bool& all_rows 
  ) const {
    
    /*
     * If "all_rows" is true, will compute all summed count rows, even with a count value of NA.
     * 
     * For compatibility with NLopt optimization signatures, this function does 
     *   not use the model's currently set parameters to make predictions, but instead 
     *   uses the parameters which have been passed to it. It also does not set the model's
     *   parameters to the ones passed to it. 
     *   
     * Computing all rows is significantly slower and not necessary for calculating likelihood;
     *   only needed when actually making predictions with the model. 
     */
    
    // Initialize variable to hold predicted rates 
    sVec predicted_rates(n_count_rows);
    iVec mc_unique_rows = Rcpp::as<iVec>(idx_mc_unique);
    
    // Initialize variables to hold model components
    sVec Rt; 
    sVec tslope;
    sVec tpoint;
    
    // Grab warping indices and initiate variables to hold them
    NumericMatrix wfactor_idx_point = wfactor_idx["point"];
    NumericMatrix wfactor_idx_rate = wfactor_idx["rate"];
    NumericMatrix wfactor_idx_slope = wfactor_idx["slope"];
    int f_pw_idx;
    int f_rw_idx;
    int f_sw_idx;
    sdouble f_pw;
    sdouble f_rw;
    sdouble f_sw;
    
    // Compute predicted rate for rows of the summed count data
    for (int r = 0; r < n_count_rows; r++) {
      
      // Update predicted model components if r begins a new batch of unique values 
      int cnt = std::count(mc_unique_rows.begin(), mc_unique_rows.end(), r);
      
      // Unless all_rows, only compute values for non-NA rows or rows that begin a new batch
      if (count_not_na_mask[r] || cnt > 0 || all_rows) {
        
        // Grab warping factors
        if (gv_ran_idx[r] == 0) {
          f_pw = 0.0;
          f_rw = 0.0; 
          f_sw = 0.0;
        } else {
          f_pw_idx = wfactor_idx_point(gv_ran_idx[r], gv_fix_idx[r]);
          f_rw_idx = wfactor_idx_rate(gv_ran_idx[r], gv_fix_idx[r]);
          f_sw_idx = wfactor_idx_slope(gv_ran_idx[r], gv_fix_idx[r]);
          f_pw = parameters(f_pw_idx);
          f_rw = parameters(f_rw_idx); 
          f_sw = parameters(f_sw_idx);
        } 
        
        // Only update predicted model components if r begins a new batch of unique values 
        int cnt = std::count(mc_unique_rows.begin(), mc_unique_rows.end(), r);
        if (cnt > 0) { 
          
          // Compute warped model components for this row r
          Rt = compute_warped_mc("Rt", r, parameters, f_rw);        
          tslope = compute_warped_mc("tslope", r, parameters, f_sw); 
          tpoint = compute_warped_mc("tpoint", r, parameters, f_pw);
          
        } 
        
        // Compute the predicted rate
        int c_num = Rwhich(eq_left_broadcast(child_lvls, child[r]))[0];
        int p_num = Rwhich(eq_left_broadcast(parent_lvls, parent[r]))[0];
        
        // Compute the actual poly-sigmoid!! 
        predicted_rates(r) = poly_sigmoid(
          bin(r),
          degMat(c_num, p_num),
          Rt,
          tslope,
          tpoint
        ); 
        
      }
      
    }
    
    return predicted_rates;
    
  }

/*
 * *************************************************************************
 * Computing objective function (i.e., fitting model and parameter boundary distances)
 */

// Compute neg_loglik of observations under the given parameters
sdouble wspc::neg_loglik(
    const sVec& parameters
  ) const {
   
    // Initialize variable to hold (what will become) nll
    sdouble log_lik = 0.0;
    
    // Predict rates under these parameters
    sVec predicted_rates_log_var = predict_rates(
      parameters, // model parameters for generating prediction
      false       // compute all summed count rows, even with a count value of NA?
      );
    
    // Compute the log-likelihood of the count data, assuming a Poisson distribution with Gamma kernel for over-dispersion
    for (int r : count_not_na_idx) {
      
      if (std::isinf(predicted_rates_log_var(r)) || predicted_rates_log_var(r) < rt_lower_bound || std::isnan(predicted_rates_log_var(r))) {
        return sdouble(inf_);
      } else {
        
        // Find gamma variance for this row
        // ... grab parent and child index numbers
        int n_c = gd_child_idx[(String)child[r]];
        int n_p = gd_parent_idx[(String)parent[r]];
        // ... pull predicted rate from log space
        sdouble pred_rate_var = sexp(predicted_rates_log_var(r)) - 1.0;
        // ... estimate variance of rate outside log space
        sdouble gamma_variance = pred_rate_var + gamma_dispersion(n_c, n_p) * pred_rate_var * pred_rate_var;
        // ... estimate the corresponding variance of the rate back in log space 
        gamma_variance = delta_var_est(gamma_variance, pred_rate_var);
        
        // Pull up if only approximately zero 
        if (predicted_rates_log_var(r) < 0.0) {predicted_rates_log_var(r) = 0.0;}
        
        // Analytic solution to the log of the integral from 0 to positive infinity of the Poisson times Gamma densities
        if (gamma_variance == 0) { 
          // if no over-dispersion, just use Poisson
          log_lik += stan::math::poisson_lpmf(count_log(r), predicted_rates_log_var(r));
        } else if (count_log(r) > 0.0 && predicted_rates_log_var(r) > 0.0) { 
          // otherwise, use Poisson-Gamma integral
          log_lik += slog(poisson_gamma_integral(count_log(r), predicted_rates_log_var(r), gamma_variance));
        }
        
      }
      
    }
    
    // Take negative
    sdouble negloglik = -log_lik;
    
    // Check for infinities (zero likelihood)
    if (std::isinf(negloglik) || negloglik > sdouble(inf_)) {
      negloglik = sdouble(inf_);
    }
    
    return negloglik;
    
  } 

// Compute boundary distances
sVec wspc::boundary_dist(
    const sVec& parameters
  ) const {
    
    // Initialize vector to hold boundary distances
    sVec boundary_dist_vec(boundary_vec_size);
    int ctr = 0;
    
    // Grab warping indices and initiate variables to hold them
    NumericMatrix wfactor_idx_point = wfactor_idx["point"];
    NumericMatrix wfactor_idx_rate = wfactor_idx["rate"];
    NumericMatrix wfactor_idx_slope = wfactor_idx["slope"];
    int f_pw_idx;
    int f_rw_idx;
    int f_sw_idx;
    sdouble f_pw;
    sdouble f_rw;
    sdouble f_sw;
    
    // Compute the boundary distance, for ...
    // ... transition slopes (must be positive)
    // ... transition points (enforces tpoint buffer)
    // ... block rate values (must be positive)
    for (int r : idx_mc_unique) {
      
      // Grab warping factors
      if (gv_ran_idx[r] == 0) {
        f_pw = 0.0;
        f_rw = 0.0; 
        f_sw = 0.0;
      } else {
        f_pw_idx = wfactor_idx_point(gv_ran_idx[r], gv_fix_idx[r]);
        f_rw_idx = wfactor_idx_rate(gv_ran_idx[r], gv_fix_idx[r]);
        f_sw_idx = wfactor_idx_slope(gv_ran_idx[r], gv_fix_idx[r]);
        f_pw = parameters(f_pw_idx);
        f_rw = parameters(f_rw_idx); 
        f_sw = parameters(f_sw_idx);
      } 
      
      // Compute Rt for this row r
      // ... note: Not computing rate for all rows, just ones with unique model components
      // ... this reduces number of rows to compute by a factor of bin_num, which is significant
      sVec Rt = compute_warped_mc("Rt", r, parameters, f_rw);
      
      // Grab degree for this row
      int deg = Rt.size() - 1;
      
      // Compute t-point and R_sum boundary distances
      if (deg > 0){
        
        // Compute tslope and tpoint for this row r
        sVec tslope = compute_warped_mc("tslope", r, parameters, f_sw); 
        sVec tpoint = compute_warped_mc("tpoint", r, parameters, f_pw);
        
        // Transition slopes most be positive
        for (int d = 0; d < deg; d++) {
          boundary_dist_vec(ctr) = tslope(d); 
          ctr++;
        }
        
        // Find tpoint boundary distances
        // WARNING: this code is duplicated in test_tpoint
        sVec tpoint_ext(deg + 2);
        for (int bt = 0; bt < tpoint_ext.size(); bt++) {
          if (bt == 0) {
            tpoint_ext(bt) = 0.0;
          } else if (bt <= deg) {
            tpoint_ext(bt) = tpoint(bt - 1);
          } else { 
            tpoint_ext(bt) = bin_num;
          } 
        }
        
        // Transition points must be in increasing order 
        // ... and can't be closer than the buffer
        // ... and first point must be > buffer, 
        // ... and last point must be < bin_num - buffer
        for (int d = 0; d < deg + 1; d++) {
          sdouble buffer_dist = (tpoint_ext(d + 1) - tpoint_ext(d)) - tpoint_buffer; 
          boundary_dist_vec(ctr) = buffer_dist;
          ctr++;
        }
        
        // All block rate values must be positive
        for (int d = 0; d < deg + 1; d++) {
          boundary_dist_vec(ctr) = Rt(d) - rt_lower_bound;
          ctr++;
        }
        
      } else {
        
        // ... trivial to check if rate (Rt) is positive
        boundary_dist_vec(ctr) = Rt(0) - rt_lower_bound;
        ctr++;
        
      }
      
    } 
    
    // Check for nan
    for (int i = 0; i < boundary_vec_size; i++) {
      if (std::isnan(boundary_dist_vec(i))) {
        boundary_dist_vec(i) = 0.0;
      }
    }
    
    return boundary_dist_vec;
    
  }

// Test for tpoints below the buffer
bool wspc::test_tpoints(
    const sVec& parameters,
    const bool& verbose
  ) const {
    
    // Initialize vector to hold tpoints and unique model-component rows
    sVec tpoint;
    iVec mc_unique_rows = Rcpp::as<iVec>(idx_mc_unique);
    
    // Grab point-warping indices and initiate variables to hold them
    NumericMatrix wfactor_idx_point = wfactor_idx["point"];
    int f_pw_idx;
    sdouble f_pw;
    
    // Compute predicted rate for rows of the summed count data
    for (int r = 0; r < n_count_rows; r++) {
      
      // Only update predicted model components if r begins a new batch of unique values 
      int cnt = std::count(mc_unique_rows.begin(), mc_unique_rows.end(), r);
      if (cnt > 0) { 
        
        // Grab warping factors
        if (gv_ran_idx[r] == 0) {
          f_pw = 0.0;
        } else {
          f_pw_idx = wfactor_idx_point(gv_ran_idx[r], gv_fix_idx[r]);
          f_pw = parameters(f_pw_idx);
        }
        
        // Find degree
        int c_num = Rwhich(eq_left_broadcast(child_lvls, child[r]))[0];
        int p_num = Rwhich(eq_left_broadcast(parent_lvls, parent[r]))[0];
        int deg = degMat(c_num, p_num);
        if (deg > 0) {
          
          // Compute tpoints for this row r
          tpoint = compute_warped_mc("tpoint", r, parameters, f_pw);
          
          // Find tpoint boundary distances
          // WARNING: this code is duplicated from boundary_dist
          sVec tpoint_ext(deg + 2);
          for (int bt = 0; bt < tpoint_ext.size(); bt++) {
            if (bt == 0) {
              tpoint_ext(bt) = 0.0;
            } else if (bt <= deg) {
              tpoint_ext(bt) = tpoint(bt - 1);
            } else { 
              tpoint_ext(bt) = bin_num;
            } 
          }
          // Transition points must be in increasing order 
          // ... and can't be closer than the buffer
          // ... and first point must be > buffer, 
          // ... and last point must be < bin_num - buffer
          for (int d = 0; d < deg + 1; d++) {
            sdouble buffer_dist = (tpoint_ext(d + 1) - tpoint_ext(d)) - tpoint_buffer;
            if (buffer_dist < 0.0) {
              if (verbose) {
                vprint(
                  "Found tpoint below buffer, dist: " + std::to_string(buffer_dist.val()) +
                    ", deg: " + std::to_string(d) +
                    ", row: " + std::to_string(r) +
                    ", tpoint_ext(d + 1): " + std::to_string(tpoint_ext(d + 1).val()) +
                    ", tpoint_ext(d): " + std::to_string(tpoint_ext(d).val()),
                    true);
              }
              return false;
            }
          }
          
        }
        
      }  
      
    }  
    
    return true;
    
  } 

// Compute min boundary penalty
sdouble wspc::min_boundary_dist(
    const sVec& parameters
  ) const {
   
    // Compute boundary_dist and take min
    sVec bd = boundary_dist(parameters);
    sdouble bd_min = smin(bd);
    
    return bd_min;
    
  }

// Wrap neg_min_boundary_dist in form needed for NLopt constraint function
double wspc::min_boundary_dist_NLopt(
    const dVec& x,
    dVec& grad,
    void* data
  ) { 
    
    // Grab model
    wspc* model = static_cast<wspc*>(data);
    
    // Convert dVec to Eigen with stan
    sVec parameters_var = to_sVec(x);
    
    // Compute max_boundary_dist
    sdouble fx = model->min_boundary_dist(parameters_var);
    
    // Compute gradient if needed
    if (!grad.empty()) {
      Eigen::VectorXd grad_eigen = model->grad_min_boundary_dist(parameters_var);
      grad.assign(grad_eigen.data(), grad_eigen.data() + grad_eigen.size());
    }  
    
    // Return the value of the neg_loglik
    return fx.val(); 
    
  } 

// Compute neg_loglik plus boundary penalty (main objective function) 
sdouble wspc::bounded_nll(
    const sVec& parameters
  ) const {
    
    // Compute weighted negative log-likelihood
    sdouble bnll = neg_loglik(parameters);
    
    // Compute boundary distance
    sVec bd = boundary_dist(parameters);
    
    // Add boundary penalty
    for (int i = 0; i < boundary_vec_size; i++) {
      bnll += boundary_penalty_transform(bd(i), bp_coefs(i));
    }
    
    /*
     * Idea of boundary penalty transform: When "far" from boundary, the total penalty will be at most 
     *  some specified fraction (e.g., 0.1) of the magnitude of the neg_loglik, but if any one component of 
     *  the boundary distance approaches zero, the penalty smoothly goes to infinity. 
     */
    
    return bnll;
    
  }

// Wrap bounded_nll in form needed for NLopt objective function
double wspc::bounded_nll_NLopt(
    const dVec& x, 
    dVec& grad, 
    void* data
  ) {
    
    // Grab model
    wspc* model = static_cast<wspc*>(data);
    
    // Convert dVec to Eigen with stan
    sVec parameters_var = to_sVec(x);
    
    // Compute bounded_nll
    sdouble fx = model->bounded_nll(parameters_var);
    
    // Compute gradient if needed
    if (!grad.empty()) {
      Eigen::VectorXd grad_eigen = model->grad_bounded_nll(parameters_var);
      grad.assign(grad_eigen.data(), grad_eigen.data() + grad_eigen.size());
    } 
    
    // Return the value of the neg_loglik
    return fx.val(); 
    
  }

/*
 * *************************************************************************
 * Computing gradients with stan reverse-mode autodiff
 */

// Compute the gradient of the bounded_nll function
// ... this is the gradient function used in model optimization
Eigen::VectorXd wspc::grad_bounded_nll(
    const sVec& p_
  ) const { 
   
    // Create nested autodiff context
    stan::math::nested_rev_autodiff nested;
    
    // Make copy
    sVec p = p_;
    
    // Initialize bounded_nll variable
    sdouble bnll = bounded_nll(p);
    
    // Initialize variable to hold gradient
    Eigen::VectorXd gr_bnll(p.size());
    
    // Compute bounded_nll and its gradient
    stan::math::grad(bnll, p, gr_bnll);
    
    // Return bounded_nll gradient
    return gr_bnll;
    
  }

// Compute the gradient of the max_boundary_dist function
// ... this is the gradient function used in the search for feasible parameters
Eigen::VectorXd wspc::grad_min_boundary_dist(
    const sVec& p_
  ) const { 
    
    // Create nested autodiff context
    stan::math::nested_rev_autodiff nested;
    
    // Make copy
    sVec p = p_;
    
    // Initialize min_boundary_dist variable
    sdouble mbd = min_boundary_dist(p);
    
    // Initialize variable to hold gradient
    Eigen::VectorXd gr_mbd(p.size());
    
    // Compute min_boundary_dist and its gradient
    stan::math::grad(mbd, p, gr_mbd);
    
    // Return min_boundary_dist gradient
    return gr_mbd;
    
  }

/*
 * *************************************************************************
 * Bootstrapping and MCMC model fitting, for statistical testing
 */

// Fit model using NLopt
void wspc::fit(
    const bool set_params, 
    const bool verbose
  ) { 
    
    // Set boundary-penalty coefficients 
    sVec initial_params_var = to_sVec(fitted_parameters);
    sdouble max_penalty_at_distance = neg_loglik(initial_params_var) * max_penalty_at_distance_factor;
    sdouble coefs_square = static_cast<double>(boundary_vec_size)/max_penalty_at_distance;
    sdouble coefs = ssqrt(coefs_square);
    bp_coefs = boundary_dist(initial_params_var);
    for (int i = 0; i < boundary_vec_size - 1; i++) {
      bp_coefs(i) = coefs/bp_coefs(i);
    } 
    
    // Grab initial parameters
    dVec x = to_dVec(fitted_parameters); 
    size_t n = x.size();
    
    // Set up NLopt optimizer
    nlopt::opt opt(nlopt::LD_LBFGS, n); // Might try? LD_LBFGS, LN_SBPLX, LN_COBYLA, GN_DIRECT
    opt.set_min_objective(wspc::bounded_nll_NLopt, this);
    opt.set_ftol_rel(ctol);             // Stop when iteration changes objective fn value by less than this fraction 
    opt.set_maxeval(max_evals);         // Maximum number of evaluations to try
    
    // Find and print initial neg_loglik and total objective values
    if (verbose) {
      sdouble initial_nll = neg_loglik(initial_params_var);
      vprint("Initial neg_loglik: ", initial_nll);
      sdouble initial_obj = bounded_nll(initial_params_var);
      vprint("Initial neg_loglik with penalty: ", initial_obj);
      sdouble mb_dist = min_boundary_dist(initial_params_var);
      vprint("Initial min boundary distance: ", mb_dist);
    } 
    
    // Fit model
    int success_code = 0;
    double min_fx;
    try {
      nlopt::result sc = opt.optimize(x, min_fx);
      success_code = static_cast<int>(sc);
    } catch (std::exception& e) { 
      if (verbose) {
        Rcpp::Rcout << "Optimization failed: " << e.what() << std::endl;
      } 
      success_code = 0;
    } 
    
    // Find final neg_loglik for bs diagnostics
    sVec parameters_final = to_sVec(x);
    sdouble final_nll = neg_loglik(parameters_final);
    
    // Print final neg_loglik, total objective, and min boundary distance values
    if (verbose) {
      vprint("Final neg_loglik: ", final_nll);
      vprint("Final neg_loglik with penalty: ", min_fx);
      sdouble mb_dist = min_boundary_dist(parameters_final);
      vprint("Final min boundary distance: ", mb_dist);
      int num_evals = opt.get_numevals();
      vprint("Number of evaluations: ", num_evals);
      if (success_code == 0) {
        vprint("Warning: optimization did not converge.");
      }
    } 
    
    if (set_params) { // See above for explanation
      set_parameters(to_NumVec(x), verbose);
    } 
    
    // Save optimization results
    //fitted_parameters = to_NumVec(x);
    optim_results["fitted_parameters"] = to_NumVec(x); // for predicting log-linked count
    optim_results["penalized_neg_loglik"] = min_fx;
    optim_results["neg_loglik"] = final_nll.val();
    optim_results["success_code"] = success_code;
    optim_results["num_evals"] = opt.get_numevals();
    
  } 

// Fit model to bootstrap resample
dVec wspc::bs_fit(
    int bs_num,               // A unique number for this resample
    bool clear_stan           // Recover stan memory at end?
  ) {
    
    // Set random number generator with unique seed based on resample ID
    unsigned int seed = rng_seed + bs_num;
    pcg32 rng(seed);
    
    // Resample (with replacement), re-sum token pool, take logs, recompute gamma_dispersion
    for (int r : count_not_na_idx) {
      
      // ... but only for actual observations of random grouping variable
      if (ran[r] != "none") {
        
        // ... redraw randomly (with replacement) from the token pool of r and re-sum into r's count 
        IntegerVector token_pool_r = token_pool[r];
        int resample_sz = token_pool_r.size();
        if (resample_sz < 1) {
          // ... ensure new point is viable
          Rcpp::Rcout << "Error: resample size < 1, for row " << r << std::endl;
          Rcpp::stop("Error: resample size < 1");
        }
        
        count(r) = 0.0;
        for (int i = 0; i < resample_sz; i++) {
          // ... randomly select integer between 0 and resample_sz
          int resample_idx = rng(resample_sz); 
          count(r) += count_tokenized[token_pool_r[resample_idx]];
        }
        count_log(r) = slog(count(r) + 1.0);
        
      }
    }
    
    // Extrapolate none's and take their logs
    count = extrapolate_none(count, ran, extrapolation_pool, true);
    iVec r_rows = Rcpp::as<iVec>(count_row_nums[eq_left_broadcast(ran,"none")]);
    for (int r : r_rows) {
      count_log(r) = slog(count(r) + 1.0);
    }
    
    // Estimate gamma dispersion of these new raw re-sampled counts
    List gamma_ests = compute_gamma_dispersion(
      count,
      count_not_na_mask,
      parent,
      child,
      parent_lvls,
      child_lvls
    );
    gamma_dispersion = to_sMat(Rcpp::as<NumericMatrix>(gamma_ests["gamma_dispersion"]));
    gd_child_idx = Rcpp::as<IntegerVector>(gamma_ests["gd_child_idx"]);
    gd_parent_idx = Rcpp::as<IntegerVector>(gamma_ests["gd_parent_idx"]);
    
    // Find average these new re-sampled log counts for each parent-child combination
    count_log_avg_mat_list = find_count_log_means(
      bin,
      count_log,
      count_not_na_mask,
      (int)bin_num.val(), 
      parent,
      child,
      treatment,
      parent_lvls,
      child_lvls,
      treatment_lvls,
      treatment_components
    );
    
    // Find initial parameter estimates for fixed-effect treatments, for new re-sampled data
    List init_params = estimate_initial_parameters(
      count, 
      count_not_na_mask,
      treatment_num,
      rise_threshold_factor,
      min_initialization_slope,
      weight_rows, 
      mc_list,
      parent,
      child,
      parent_lvls,
      child_lvls,
      degMat,
      found_cp_list, 
      found_cp_trt_list,
      count_log_avg_mat_list 
    );
    List ref_values = Rcpp::as<List>(init_params["ref_values"]);       // reference values for each parent-child combination
    List RtEffs = Rcpp::as<List>(init_params["RtEffs"]);               // rate effects for each parent-child combination
    List tpointEffs = Rcpp::as<List>(init_params["tpointEffs"]);       // tpoint effects for each parent-child combination
    List tslopeEffs = Rcpp::as<List>(init_params["tslopeEffs"]);       // tslope effects for each parent-child combination
    
    // Build default fixed-effects matrices in shell
    List beta = build_beta_shell(mc_list, treatment_lvls, parent_lvls, child_lvls, ref_values, RtEffs, tpointEffs, tslopeEffs, degMat);
    
    // Initialize new random effect warping factors 
    List wfactors = make_initial_random_effects(
      wfactors_names,
      ran_lvls.size(),
      child_lvls.size()
    );
    
    // Make and map parameter vector
    List params = make_parameter_vector(
      beta, wfactors,
      parent_lvls, child_lvls, ran_lvls,
      mc_list, 
      treatment_lvls,
      degMat
    );
    
    // Extract new parameter vector information
    fitted_parameters = params["param_vec"];
    // ... rest should be the same as before, but extracting for sanity
    param_names = params["param_names"];
    param_wfactor_rate_idx = params["param_wfactor_rate_idx"];
    param_wfactor_point_idx = params["param_wfactor_point_idx"]; 
    param_wfactor_slope_idx = params["param_wfactor_slope_idx"];
    param_beta_Rt_idx = params["param_beta_Rt_idx"];
    param_beta_tslope_idx = params["param_beta_tslope_idx"];
    param_beta_tpoint_idx = params["param_beta_tpoint_idx"];
    param_baseline_idx = params["param_baseline_idx"];
    beta_idx = params["beta_idx"];
    wfactor_idx = params["wfactor_idx"];
    
    // Check feasibility 
    List feasibility_results = check_parameter_feasibility(to_sVec(fitted_parameters), false); 
    fitted_parameters = Rcpp::as<NumericVector>(feasibility_results["parameters"]);
    
    // Fit model
    fit(
      false, // don't bother setting parameters
      false  // don't print anything
    );
    
    // Prepare and return results
    dVec fitted_parameters = to_dVec(Rcpp::as<NumericVector>(optim_results["fitted_parameters"]));
    fitted_parameters.reserve(fitted_parameters.size() + 4);
    fitted_parameters.push_back(Rcpp::as<double>(optim_results["penalized_neg_loglik"]));
    fitted_parameters.push_back(Rcpp::as<double>(optim_results["neg_loglik"]));
    fitted_parameters.push_back(Rcpp::as<double>(optim_results["success_code"]));
    fitted_parameters.push_back(Rcpp::as<double>(optim_results["num_evals"]));
    
    // Clear stan memory
    if (clear_stan) {
      clear_stan_mem();
    }
    
    return fitted_parameters;
    
  }

// Fork bootstraps (parallel processing)
Rcpp::NumericMatrix wspc::bs_batch(
    int bs_num_max,           // Number of bootstraps to perform
    int max_fork,             // Maximum number of forked processes per batch
    bool verbose
  ) {
    
    // Reset initial parameters to check their feasibility 
    set_parameters(fitted_parameters, verbose);
    
    // Fit full model 
    vprint("Performing initial fit of full data", verbose);
    fit(
      false, // don't bother setting parameters
      false  // don't print anything
    );
    double pnll = optim_results["penalized_neg_loglik"]; 
    if (verbose) {vprint("Penalized neg_loglik: ", pnll);}
    
    // Initiate variables to hold results
    const int c_num = fitted_parameters.size() + 4;
    const int r_num = bs_num_max + 1;
    NumericMatrix bs_results(r_num, c_num);
    
    // Save results from initial full fit
    NumericVector these_results = optim_results["fitted_parameters"];
    dVec full_results = to_dVec(these_results);
    full_results.reserve(full_results.size() + 4);
    full_results.push_back(Rcpp::as<double>(optim_results["penalized_neg_loglik"]));
    full_results.push_back(Rcpp::as<double>(optim_results["neg_loglik"]));
    full_results.push_back(Rcpp::as<double>(optim_results["success_code"]));
    full_results.push_back(Rcpp::as<double>(optim_results["num_evals"]));
    
    // Save full fit results in last row of results matrix
    bs_results.row(bs_num_max) = to_NumVec(full_results);
    
    // Perform bootstrap fits in batches
    if (max_fork < 1) {max_fork = 1;} 
    const int batch_num = std::round(bs_num_max / max_fork);
    int tracker_steps = 10;
    IntegerVector tracker = iseq(batch_num/tracker_steps, batch_num, tracker_steps);
    for (int b = 0; b < batch_num; b++) {
      
      // Initiate timer and grab start time
      Timer batch_timer;
      batch_timer.step("start");  
      
      if (max_fork > 1) {
        // Run in parallel with forking
        
        // Pipes for inter-process communication
        iVec pids(bs_num_max);
        std::vector<std::array<int, 2>> pipes(bs_num_max); 
        
        // Initialize pipes and fork processes
        for (int i = 0; i < max_fork; i++) {
          
          int this_row = b * max_fork + i;
          pipe(pipes[i].data()); // Create a pipe
          pid_t pid = fork();
          
          if (pid == 0) { // Child process
            
            // Close read end
            close(pipes[i][0]); 
            
            // Fit bootstrap
            dVec result = bs_fit(this_row, false); 
            
            // Send result
            write(pipes[i][1], result.data(), sizeof(double) * c_num);
            
            // Close write end
            close(pipes[i][1]); 
            
            // Exit child process
            _exit(0); 
            
            // implicitly recovers stan memory by killing the process which initiated it
            
          } else if (pid > 0) { // Parent process
            
            // Grab child pid
            pids[i] = pid;
            
            // Close write end
            close(pipes[i][1]); 
            
          } else {
            Rcpp::stop("Fork failed!");
          }
          
        }
        
        // Fetch results from pipes
        for (int i = 0; i < max_fork; i++) {
          
          // Wait for child process
          waitpid(pids[i], NULL, 0); 
          
          // Create a temporary buffer to hold the row data
          dVec temp_row(c_num); 
          
          // Read the row from the pipe into the buffer
          read(pipes[i][0], temp_row.data(), sizeof(double) * c_num);
          
          // Copy the buffer contents into the corresponding row of the matrix
          int this_row = b * max_fork + i;
          bs_results.row(this_row) = to_NumVec(temp_row);
          
          // Close read end
          close(pipes[i][0]); 
          
        } 
        
      } else {
        // Run in serial
        
        // Fit bootstrap
        dVec result = bs_fit(b, true); 
        
        // Copy into the corresponding row of the matrix
        bs_results.row(b) = to_NumVec(result);
        
      }
      
      batch_timer.step("end");
      NumericVector batch_times(batch_timer);
      double batch_time = 1e-9 * (batch_times[1] - batch_times[0])/max_fork;
      if (verbose) {
        // Tracker
        if (any_true(eq_left_broadcast(tracker, b)) || b == batch_num || b == 1) {
          Rcpp::Rcout << "Batch: " << b << "/" << batch_num << ", " << batch_time << " sec/bs" << std::endl;
        }
      }
      
    }
    
    vprint("All complete!", verbose); 
    
    // Clear stan memory
    if (bs_num_max == 0) {
      clear_stan_mem();
    }
    
    return bs_results;
    
  }

// Markov-chain Monte Carlo (MCMC) simulation
Rcpp::NumericMatrix wspc::MCMC(
    int n_steps,              // Number of steps to take in random walk
    int neighbor_filter,      // Keep only ever neighbor_filter step
    double step_size,         // Step size for random walk
    double prior_sd,          // standard deviation to use in prior
    bool verbose
  ) {
    
    // Reset initial parameters to check their feasibility 
    set_parameters(fitted_parameters, verbose);
    
    // Fit full model 
    vprint("Performing initial fit of full data", verbose);
    fit(
      false,  // don't set fitted parameters
      false   // don't print anything
    );
    double pnll = optim_results["penalized_neg_loglik"]; 
    if (verbose) {vprint("Penalized neg_loglik: ", pnll);}
    
    // Initiate variables to hold results
    const int n_params = fitted_parameters.size();
    const int c_num = n_params + 4;
    const int r_num = n_steps + 1;
    NumericMatrix RMW_steps(r_num, c_num);
    
    
    // Make baseline parameter mask 
    LogicalVector baseline_mask(n_params);
    baseline_mask.fill(false);
    for (int i : param_baseline_idx) {
      baseline_mask(i) = true;
    }
    
    // Make beta_Rt parameter mask
    LogicalVector beta_Rt_mask(n_params);
    beta_Rt_mask.fill(false);
    for (int i : param_beta_Rt_idx) {
      beta_Rt_mask(i) = true;
    } 
    
    // Save results from initial full fit
    NumericVector these_results = optim_results["fitted_parameters"];
    dVec full_results = to_dVec(these_results);
    full_results.reserve(full_results.size() + 4);
    full_results.push_back(Rcpp::as<double>(optim_results["penalized_neg_loglik"]));
    full_results.push_back(Rcpp::as<double>(optim_results["neg_loglik"]));
    full_results.push_back(double(1.0)); // acceptance ratio
    full_results.push_back(double(0.0)); // ctr number
    
    // Save full fit results in last row of results matrix
    RMW_steps.row(n_steps) = to_NumVec(full_results);
    
    // Run MCMC simulation
    int step = 0;
    int ctr = 0;
    int inf_loop_ctr = 0;
    int last_viable_step = 0;
    double acceptance_rate = 1.0;
    int tracker_steps = 10;
    if (tracker_steps > n_steps/2) {tracker_steps = (int)n_steps/2;}
    IntegerVector tracker = iseq(n_steps/10, n_steps, tracker_steps);
    bool printed_step1 = false;
    // Grab current point (model parameters) in random walk
    NumericVector params_current = fitted_parameters;
    NumericVector last_viable_parameters = params_current;
    
    // Initialize neighbor counter
    int neighbor_counter = 0;
    
    // Set random number generator with seed
    unsigned int seed = rng_seed + n_steps;
    pcg32 walk_rng(seed);
    
    // Take steps, Random-walk Metropolois algorithm
    while (step < n_steps) {
      
      // Initialize priors 
      double prior_current = 0.0;
      double prior_next = 0.0;
      
      // Generate random step
      // ... initiate vector to hold new parameters
      NumericVector params_next(n_params);
      // ... test boundary distance
      sVec bd_current_vec = boundary_dist(to_sVec(params_current));
      sdouble bd_current_min = smin(bd_current_vec); 
      sdouble bd_current_transformed = 0.0;
      // ... if above boundary, scale step to boundary distance
      if (bd_current_min > 0.0) {
        
        // ... transform boundary penalty
        for (int i = 0; i < boundary_vec_size; i++) {
          bd_current_transformed += boundary_penalty_transform(bd_current_vec(i), bp_coefs(i));
        }
        bd_current_transformed += 1.0;
        
        // ... for each parameter
        for (int i = 0; i < n_params; i++) {
          
          // ... calculate step size
          double normalized_step_size = step_size * std::log10(std::abs(params_current(i)) + 1.0);
          double bounded_step_size = normalized_step_size / bd_current_transformed.val();
          if (bounded_step_size == 0.0) {
            // ... presumably this case means current parameter is extremely close to zero or very close to boundary
            params_next(i) = pcg_rnorm(params_current(i), step_size, walk_rng);
          } else {
            // ... take next step
            params_next(i) = pcg_rnorm(params_current(i), bounded_step_size, walk_rng);
          }
          
          // Check that a good baseline parameter hasn't been sent below zero
          if (baseline_mask(i) || beta_Rt_mask(i)) {
            if (params_next(i) < 0.0) {
              params_next(i) = 0.0;
            }
          }
          
          // While looping, compute priors for this random step
          std::string this_param = Rcpp::as<std::string>(param_names[i]);
          if (
              !(pattern_match("baseline", this_param) && pattern_match("tpoint", this_param))
              // ... baseline tpoint values are uniformly distributed
          ) {
            prior_current += log_dNorm(params_current(i), 0.0, prior_sd);
            prior_next += log_dNorm(params_next(i), 0.0, prior_sd);
          }
          
        }
        
      } else {
        
        params_current = last_viable_parameters;
        step = last_viable_step;
        clear_stan_mem(); // Clear stan memory
        continue;
        
      }
      
      // Compute likelihoods for this random step
      double loglik_current = -bounded_nll(to_sVec(params_current)).val();
      double loglik_next = -bounded_nll(to_sVec(params_next)).val();
      // Calculate posteriors and acceptance probability 
      double acceptance = std::exp(
        (loglik_next + prior_next) - 
        (loglik_current + prior_current) 
      );
      if (acceptance > 1.0) {acceptance = 1.0;}
      if (std::isnan(acceptance)) {acceptance = 0.0;}
      
      // Accept or reject the proposed step
      double ran_draw = R::runif(0.0, 1.0); 
      if (ran_draw < acceptance) {
        // Set new parameters as current 
        params_current = params_next;
        // Advance neighbor counter 
        neighbor_counter++;
        if (neighbor_counter == neighbor_filter) {
          // Tracker
          if (any_true(eq_left_broadcast(tracker, step))) {
            int this_step_batch = Rwhich(eq_left_broadcast(tracker, step))[0];
            Rcpp::Rcout << "step: " << (this_step_batch + 1) * (n_steps/tracker.size()) << "/" << n_steps << std::endl;
            tracker(this_step_batch) -= 1; // ensure report is only printed once
          } else if (step == 0 && !printed_step1) {
            Rcpp::Rcout << "step: " << 1 << "/" << n_steps << std::endl;
            printed_step1 = true; // ensure report is only printed once
          }
          // Save new parameters and results
          dVec full_results = to_dVec(params_next);
          full_results.reserve(full_results.size() + 4);
          full_results.push_back(-loglik_next);
          full_results.push_back(-loglik_next - (bd_current_transformed.val() - 1.0));
          full_results.push_back(acceptance); // acceptance ratio
          full_results.push_back(double(ctr + 0.0)); // ctr number
          RMW_steps.row(step) = to_NumVec(full_results);
          // Advance step
          step++;
          // Reset neighbor counter
          neighbor_counter = 0;
        }
      } 
      
      // Advance counters
      ctr++; 
      inf_loop_ctr++;
      
      // Check for infinite loop
      acceptance_rate = (static_cast<double>(step)/static_cast<double>(ctr))*static_cast<double>(neighbor_filter);
      if (inf_loop_ctr > 200) {
        if (acceptance_rate < 0.01) {
          // ... go back to last viable step
          step = last_viable_step;
          params_current = last_viable_parameters;
          vprint("Infinite loop suspected, resetting parameters to last viable set. (Try lowering MCMC.step.size)", verbose);
        } else {
          // ... update last viable step
          last_viable_step = step;
        }
        inf_loop_ctr = 0;
      }
      
      // Clear stan memory
      clear_stan_mem();
      
    }
    
    // Report acceptance rate
    vprint("All complete!", verbose); 
    if (verbose) {vprint("Acceptance rate (aim for 0.2-0.3): ", acceptance_rate);}
    
    return RMW_steps;
    
  }

// Resample (demo, not used in package)
NumericMatrix wspc::resample(
    int n_resample            // total number of resamples to draw
  ) {
    
    NumericMatrix resamples(n_count_rows, n_resample);
    IntegerVector NA_idx = Rwhich(!count_not_na_mask);
    
    for (int j = 0; j < n_resample; j++) {
      
      // Set random number generator with unique seed based on resample ID
      unsigned int seed = rng_seed + j;
      pcg32 rng(seed);
      
      sVec count_new(n_count_rows);
      count_new.setZero();
      
      // Resample (with replacement), re-sum token pool, and take logs
      for (int r : count_not_na_idx) {
        if (ran[r] != "none") {
          count_new(r) = 0.0;
          IntegerVector token_pool_r = token_pool[r];
          int resample_sz = token_pool_r.size();
          for (int i = 0; i < resample_sz; i++) {
            int resample_idx = rng(resample_sz); // randomly select integer between 0 and resample_sz
            count_new(r) += count_tokenized[token_pool_r[resample_idx]];
          }
        }
      }
      
      // Extrapolate none's 
      count_new = extrapolate_none(count_new, ran, extrapolation_pool, true);
      
      for (int r : NA_idx) {
        count_new(r) = stan::math::NOT_A_NUMBER;
      }
      resamples.column(j) = to_NumVec(count_new);
      
    }
    
    return resamples;
    
  }

/*
 * *************************************************************************
 * Setting parameters
 */

// Set the model with the given parameters
void wspc::set_parameters(
    const NumericVector& parameters,
    const bool verbose
  ) { 
    
    /*
     * This function: 
     *  - Checks feasibility of provided parameters
     *  - Predicts (and saves) rates based on provided parameters
     *  - Saves parameters in fitted_parameters vector
     */
    
    // Ensure provided parameters are feasible (i.e., don't predict negative rates) 
    sVec parameters_var = to_sVec(parameters);
    List feasibility_results = check_parameter_feasibility(parameters_var, verbose);
    bool feasible = feasibility_results["feasible"];
    
    // Stop if not feasible 
    if (!feasible) {
      Rcpp::stop("Negative or nan rates predicted!");
    }
    
    // Convert back to doubles, remove from log space, and save predicted values
    predicted_rates_log = Rcpp::as<NumericVector>(feasibility_results["predicted_rates_log"]);
    predicted_rates = NumericVector(n_count_rows); 
    for (int i = 0; i < n_count_rows; i++) {
      predicted_rates[i] = std::exp(predicted_rates_log(i)) - 1.0;
    }
    
    // Update saved parameter vector 
    fitted_parameters = Rcpp::as<NumericVector>(feasibility_results["parameters"]);
    optim_results["fitted_parameters"] = Rcpp::as<NumericVector>(feasibility_results["parameters"]);
    
  } 

// Check and correct parameter feasibility
Rcpp::List wspc::check_parameter_feasibility(
    const sVec& parameters_var,
    const bool verbose
  ) { 
    
    vprint("Checking feasibility of provided parameters", verbose); 
    
    // Initialize vectors to return 
    sVec feasible_parameters_var = parameters_var; 
    sVec predicted_rates_log_var;
    
    // Test if any tpoints are below the buffer 
    bool feasible = test_tpoints(parameters_var, verbose);
    if (verbose) {
      if (feasible) {
        vprint("... no tpoints below buffer");
      } else {
        vprint("... tpoints found below buffer");
      }
    }
    if (feasible) {
      
      // Predict rates 
      predicted_rates_log_var = predict_rates(
        parameters_var, // model parameters for generating prediction
        true            // compute all summed count rows, even with a count value of NA?
      );
      
      // Test if provided parameters produce any nan rates
      for (int i = 0; i < n_count_rows; i++) {
        if (std::isnan(predicted_rates_log_var(i))) {
          feasible = false;
        }
      }
      
      if (verbose) {
        if (feasible) {
          vprint("... no NAN rates predicted");
        } else {
          vprint("... NAN rates predicted");
        }
      }
      
      // Test if provided parameters produce any negative rates
      for (int i = 0; i < n_count_rows; i++) {
        if (!std::isnan(predicted_rates_log_var(i)) && predicted_rates_log_var(i) < rt_lower_bound) {
          feasible = false;
        }
      }
      
      if (verbose) {
        if (feasible) {
          vprint("... no negative rates predicted");
        } else { 
          vprint("... negative rates predicted");
        }
      } 
      
    }
    
    if (verbose) {
      if (feasible) {
        vprint("Provided parameters are feasible");
      } else {
        vprint("Provided parameters not feasible, searching nearby");
      }
    }
    
    // Find and report initial distance to boundary
    sdouble initial_dist = min_boundary_dist(parameters_var); 
    if (verbose) {
      vprint("Initial boundary distance (want >0): ", initial_dist);
    }
    
    // If not feasible, attempt to find a feasible starting point
    if (!feasible) {
      
      // Variables for optimization
      dVec x = to_dVec(to_NumVec(parameters_var));
      size_t n = x.size();
      
      // Set up NLopt optimizer
      nlopt::opt opt(nlopt::LD_LBFGS, n);
      opt.set_max_objective(wspc::min_boundary_dist_NLopt, this);
      opt.set_ftol_rel(ctol);               // stop when iteration changes objective fn value by less than this fraction 
      opt.set_maxeval(max_evals);           // Maximum number of evaluations to try 
      opt.set_stopval(0.01);                // ensure boundary distance is at least this much above zero (and then stop)
      double max_fx;
      int success_code = 0;
      
      // Optimize
      try {
        nlopt::result sc = opt.optimize(x, max_fx);
        success_code = static_cast<int>(sc);
      } catch (std::exception& e) {
        if (verbose) {
          Rcpp::Rcout << "Optimization failed: " << e.what() << std::endl;
        }
        success_code = 0;
      }
      
      // Find and report final distance to boundary
      if (verbose) {
        vprint("Numer of evals: ", (int)opt.get_numevals());
        vprint("Success code: ", (int)success_code);
        vprint("Final boundary distance: ", (double)max_fx);
      }
      
      // Final check of boundary distance
      if (max_fx <= 0) {
        success_code = 0;
      }
      
      // Check for success
      if (success_code == 0) {
        vprint("Could not find a nearby feasible parameters, returning provided ones", verbose); 
      } else { // found a feasible starting point, save
        feasible = true;
        vprint("Nearby feasible parameters found!", verbose); 
        feasible_parameters_var = to_sVec(x);
        // Recompute predicted rates 
        predicted_rates_log_var = predict_rates(
          feasible_parameters_var, // model parameters for generating prediction
          true                     // compute all summed count rows, even with a count value of NA?
        );
      }
      
    }
    
    // Return feasible parameters
    return List::create(
      _["feasible"] = feasible,
      _["predicted_rates_log"] = to_NumVec(predicted_rates_log_var),
      _["parameters"] = to_NumVec(feasible_parameters_var)
    );
    
  } 

/*
 * *************************************************************************
 * Export data to R
 */

Rcpp::List wspc::results() {
    
    NumericVector predicted_rates_out(n_count_rows);
    NumericVector predicted_rates_log_out(n_count_rows);
    if (predicted_rates.size() == n_count_rows) {
      // conditional to prevent trying to access an empty vector
      predicted_rates_out = predicted_rates;
      predicted_rates_log_out = predicted_rates_log;
    }
    
    // Create summed count data frame
    DataFrame count_data_summed = DataFrame::create(
      _["row"] = count_row_nums,
      _["bin"] = to_NumVec(bin),
      _["count"] = to_NumVec(count),
      _["pred"] = predicted_rates_out,
      _["count.log"] = to_NumVec(count_log),
      _["pred.log"] = predicted_rates_log_out,
      _["parent"] = parent,
      _["child"] = child,
      _["ran"] = ran,
      _["treatment"] = treatment
    );
    
    // Add parameter names to fitted parameter vector 
    fitted_parameters.names() = param_names;
    
    // Collect fixed-effect names and levels 
    fix_ref.names() = fix_names;
    List fixed_effects = List::create(
      _["name"] = fix_names,
      _["lvls"] = fix_lvls,
      _["treat.lvl"] = fix_trt,
      _["ref.lvl"] = fix_ref
    );
    
    // Pack up treatment name information
    List treat = List::create(
      _["names"] = treatment_lvls, 
      _["components"] = treatment_components
    );
    
    // Put grouping variable information into list
    List grouping_variables = List::create(
      _["parent.lvls"] = parent_lvls,
      _["child.lvls"] = child_lvls,
      _["ran.lvls"] = ran_lvls
    );
    
    // Pack up parameter indexes into list
    List param_idx = List::create(
      _["beta"] = beta_idx,
      _["w.factor"] = wfactor_idx
    );
    
    // Collect token pool
    List token_pool_list(n_count_rows); 
    for (int i = 0; i < n_count_rows; i++) {
      if (token_pool[i].size() > 0) {
        token_pool_list[i] = (IntegerVector)token_pool[i];
      } 
    }
    
    // Reformat gamma dispersion parameters 
    NumericMatrix g_dispersion = to_NumMat(gamma_dispersion);
    rownames(g_dispersion) = child_lvls;
    colnames(g_dispersion) = parent_lvls;
    
    // Make final list to return 
    List results_list = List::create(
      _["model.component.list"] = mc_list,
      _["count.data.summed"] = count_data_summed,
      _["fitted.parameters"] = fitted_parameters,
      _["gamma.dispersion"] = g_dispersion,
      _["param.names"] = param_names,
      _["fix"] = fixed_effects,
      _["treatment"] = treat,
      _["grouping.variables"] = grouping_variables,
      _["param.idx0"] = param_idx, // "0" to indicate this goes out w/ C++ zero-based indexing
      _["token.pool"] = token_pool_list,
      _["settings"] = model_settings
    );
    
    return results_list;
    
  }

/*
 * *************************************************************************
 * Testing and debugging in R
 */

// Wrap neg_loglik in form needed for R
double wspc::neg_loglik_debug(
    const dVec& x
  ) {
    
    // Convert to sVec
    sVec parameters_var = to_sVec(x);
    
    // Compute neg_loglik
    double negloglik = neg_loglik(parameters_var).val();
    
    // Clear stan memory
    clear_stan_mem(); 
    
    // Return as double
    return negloglik;
    
  }

// Wrap bounded_nll in form needed for R
double wspc::bounded_nll_debug(
    const dVec& x
  ) { 
    
    /*
     * This function is just for testing and debugging, e.g., 
     *  comparing stan grad to finite difference in R. 
     */
    
    // Convert dVec to Eigen with stan
    sVec parameters_var = to_sVec(x);
    
    // Compute bounded_nll
    double fx = bounded_nll(parameters_var).val();
    
    // Clear stan memory
    clear_stan_mem(); 
    
    // Return the value of the bounded_nll
    return fx; 
    
  } 

// Wrap grad_bounded_nll_debug in form needed for R
Rcpp::NumericVector wspc::grad_bounded_nll_debug(
    const dVec& x 
  ) const { 
    
    /*
     * This function is just for testing and debugging, e.g., 
     *  comparing stan grad to finite difference in R. 
     */
    
    // Convert dVec to Eigen with stan
    sVec parameters_var = to_sVec(x);
    
    // Compute grdient of bounded_nll
    Eigen::VectorXd grad_fx = grad_bounded_nll(parameters_var);
    
    // Cast to NumericVector
    NumericVector grad_fx_R(grad_fx.size());
    for (int i = 0; i < grad_fx.size(); i++) {
      grad_fx_R[i] = grad_fx(i);
    } 
    
    // Return the value of the grad of bounded_nll
    return grad_fx_R; 
    
  }  

// Export the class constructor and select fields and methods to R
RCPP_EXPOSED_CLASS(wspc)
RCPP_MODULE(wspc) {
    class_<wspc>("wspc")
    .constructor<DataFrame, List, bool>()  
    .field("optim_results", &wspc::optim_results)
    .field("fitted_parameters", &wspc::fitted_parameters)
    .method("set_parameters", &wspc::set_parameters)
    .method("neg_loglik_debug", &wspc::neg_loglik_debug)
    .method("bounded_nll_debug", &wspc::bounded_nll_debug)
    .method("grad_bounded_nll_debug", &wspc::grad_bounded_nll_debug)
    .method("fit", &wspc::fit)
    .method("bs_fit", &wspc::bs_fit)
    .method("bs_batch", &wspc::bs_batch)
    .method("MCMC", &wspc::MCMC)
    .method("resample", &wspc::resample)
    .method("results", &wspc::results);
  }
