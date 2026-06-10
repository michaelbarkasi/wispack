
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
    round_none = (bool)settings["round_none"];
    trtKO = Rcpp::as<CharacterVector>(settings["trtKO"]);
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
    
    // Remove any KO columns 
    CharacterVector col_names = count_data.names();
    for (String this_trtKO : trtKO) {
      IntegerVector trkKO_idx = grep_cpp(col_names, this_trtKO);
      int cnt = 0;
      for (int idx : trkKO_idx) {
        col_names.erase(col_names.begin() + (idx - cnt));  // ... remove from col_names, adjusting for previous removals
      }
    }
    
    // Check structure of input data
    CharacterVector required_cols = CharacterVector::create("count", "bin", "context", "species", "ran");
    col_names = col_names[!col_names.isin(trtKO)];
    int n_cols = col_names.size();
    int r_cols = required_cols.size();
    for (int i = 0; i < r_cols; i++) {
      if (col_names[i] != required_cols[i]) {
        Rcpp::stop("Input data is missing required column (or out of order): " + required_cols[i]);
      }
    }
    
    // Save tokenized count column before collapsing to sums 
    count_tokenized = to_sVec(Rcpp::as<NumericVector>(count_data["count"]));
    
    // Find max bins and set warp bounds
    vprint_header("Extracting variables and data information", verbose);
    bin_num = (sdouble)settings["max_bin"];
    if (bin_num == 0.0) {bin_num = smax(to_sVec(Rcpp::as<NumericVector>(count_data["bin"])));}
    warp_bounds.resize(3); 
    for (int i = 0; i < 3; i++) {warp_bounds[i] = inf_warp;}            // initialize to infinity (no bounds)
    warp_bounds[2] = bin_num;                                           // set tpoint bound to max bin
    vprint("Max bin: " + std::to_string(bin_num.val()), verbose);
    
    // Extract fixed effects 
    int n_fix = n_cols - r_cols;                                        // number of fixed effect variables, assume all non-required columns are fixed effects
    fix_names = CharacterVector(n_fix);                                 // names of fixed effect variables 
    fix_ref = CharacterVector(n_fix);                                   // reference level for each fixed effect
    fix_lvls.resize(n_fix);                                             // levels for each fixed effect, fix_lvls is a std vector holding an Rcpp CharacterVector
    fix_trt.resize(n_fix);                                              // treatments for each fixed effect, fix_trt is a std vector holding an Rcpp CharacterVector
    LogicalVector is_time = {false};                                    // logical vector to track which fixed effect levels are in "timeseries"
    CharacterVector is_time_name = {"ref"};
    CharacterVector time_series_names = {"no time series"};
    for (int i = 0; i < n_fix; i++) {
      fix_names[i] = col_names[i + r_cols];
      
      // Grab column of fixed-effect variable
      SEXP col = count_data[fix_names[i]];
      CharacterVector lvls;
      
      // Sort levels
      if (Rf_isNumeric(col)) {
        // Sort numerically, then convert to character
        NumericVector tmp = Rcpp::sort_unique(Rcpp::as<NumericVector>(col));
        lvls = Rcpp::as<CharacterVector>(tmp);
      } else {
        // Sort lexicographically
        lvls = Rcpp::sort_unique(Rcpp::as<CharacterVector>(col));
      }
      
      // Check number of levels
      if (lvls.size() < 2) {
        Rcpp::stop("Fixed effect " + fix_names[i] + " has less than 2 levels.");
      } else if (lvls.size() > 2 && fix_names[i] != "timeseries") {
        Rcpp::stop("Fixed effect " + fix_names[i] + " has more than 2 levels. Only binary fixed effects are currently supported, except for 'timeseries'.");
      }
      fix_lvls[i] = lvls;
      fix_ref[i] = lvls[0];                                             // assume first level is reference
      fix_trt[i] = lvls[Rcpp::Range(1, lvls.size() - 1)];               // assume all other levels are treatment levels
      
      // Check for time series
      if (fix_names[i] == "timeseries") {                               // make it easy to access time-series element rank from the element name
        timeseries_rank = seq(1, lvls.size());
        timeseries_rank.names() = lvls;
        time_series_names = lvls;
        for (String l : lvls) {
          is_time.push_back(true);
          is_time_name.push_back(l);
        }
      } else {
        for (String l : lvls) {
          is_time.push_back(false);
          is_time_name.push_back(l);
        }
      }
    }
    is_time.names() = is_time_name;
    
    // Report
    if (fix_names.size() == 0) {
      vprint("No fixed effects detected.", verbose);
    } else {
      vprint("Fixed effects:", verbose);
      vprintV(fix_names, verbose);
      vprint("Ref levels:", verbose);
      vprintV(fix_ref, verbose);
    }
    vprint("No time series detected.", verbose && !any_true(is_time));
    vprint("Time series detected:", verbose && any_true(is_time));
    vprintV(time_series_names, verbose && any_true(is_time));
    
    // Create all possible treatment combinations 
    CharacterVector ref_lvl = {"ref"};                                  // define a baseline (i.e., reference) level
    treatment_components = make_treatments(fix_trt);                    // make components of all possible treatment levels
    treatment_components.insert(treatment_components.begin(), ref_lvl); // add "ref" to represent reference level
    treatment_num = treatment_components.size();                        // grab number of treatments
    treatment_lvls = CharacterVector(treatment_num);                    // resize variable to hold treatment level names
    for (int t = 0; t < treatment_num; t++) {                           // collapse components into level names
      treatment_lvls[t] = Rcpp::collapse(treatment_components[t]);
    }
    
    // Check for any treatment-level knockouts and remove
    for (String this_trtKO : trtKO) {
      IntegerVector trkKO_idx = grep_cpp(treatment_lvls, this_trtKO);
      int cnt = 0;
      for (int idx : trkKO_idx) {
        treatment_lvls.erase(treatment_lvls.begin() + (idx - cnt));  // ... remove from treatment_lvls
        treatment_components.erase(treatment_components.begin() + (idx - cnt++));  // ... remove from treatment_components, adjusting for previous removals
        treatment_num--;  // ... decrease
      }
    }
    vprint("Created treatment levels:", verbose);
    vprintV(treatment_lvls, verbose);
    
    // Create matrix to translate between treatment levels and fixed-effects columns
    // ... for making summed-count fixed-effect column 
    // ... each cell should contain the value of the treatment column n_fix, for the given treatment level
    CharacterMatrix effects_rows(treatment_num, n_fix);
    for (int tr = 0; tr < treatment_num; tr++) {
      CharacterVector trt_components = treatment_components[tr];
      for (int f = 0; f < n_fix; f++) {
        CharacterVector lvls = fix_lvls[f];
        effects_rows(tr, f) = lvls[0]; 
        // ^ ... Assume it's the reference level
        for (String l : lvls) {
          // If treatment level found, replace
          if (any_true(eq_left_broadcast(trt_components, l))) {effects_rows(tr, f) = l;}
        }
      }
    }
    
    // When predicting parameter values for treatment level tr_predict, should treatment level tr_input's effect be applied?
    // ... i.e., what's the "weight" of each treatment level (columns of weight matrix) relative to the others (rows of weight matrix)? 
    weight_rows.resize(treatment_num, treatment_num);
    weight_rows.setOnes(); // ... assume "yes"
    for (int tr_predict = 0; tr_predict < treatment_num; tr_predict++) {
      // Grab components of treatment level tr_predict
      CharacterVector tr_predict_components = treatment_components[tr_predict];
      for (int tr_input = 0; tr_input < treatment_num; tr_input++) {
        // Grab components of treatment level tr_input
        CharacterVector tr_input_components = treatment_components[tr_input];
        // ... if tr_input is the base reference level, tr_input must be applied
        if (tr_input_components.size() == 1 && tr_input_components[0] == "ref") {continue;} 
        // ... else, must check whether all trc_i of tr_input are in tr_predict
        LogicalVector input_comp_included;
        for (String trc_i : tr_input_components) {
          // if trc_i is not in tr_predict ...
          if(!any_true(eq_left_broadcast(tr_predict_components, trc_i))) {
            // ... check whether trc_i is a time point from "timeseries" ...
            if (is_time[trc_i]) {
              // ... then check if any component trc_p of tr_predict ...
              bool pred_has_time = false;
              for (String trc_p : tr_predict_components) {
                // ... is a time point in "timeseries" ...
                if (is_time[trc_p]) { 
                  pred_has_time = true;
                  // ... that is greater than or equal to trc_i
                  int i_rank = timeseries_rank[trc_i];
                  int p_rank = timeseries_rank[trc_p];
                  if (i_rank <= p_rank) {
                    // ... if so, trc_i is included in tr_predict
                    input_comp_included.push_back(true);
                  } else {
                    // ... else, trc_i is not included in tr_predict
                    input_comp_included.push_back(false);
                  }
                  // ... if so, tr_input applies and should be left at 1.0
                } 
              }
              if (!pred_has_time) {
                // trc_i is not in tr_predict
                input_comp_included.push_back(false);
              }
            } else {
              // trc_i is not in tr_predict
              input_comp_included.push_back(false);
            }
          } else {
            // trc_i is in tr_predict
            input_comp_included.push_back(true);
          }
        }
        // Sanity check 
        if (input_comp_included.size() != tr_input_components.size()) {
          Rcpp::stop("Error in constructing weight matrix.");
        }
        if (!all_true(input_comp_included)) {weight_rows(tr_predict, tr_input) = 0.0;}
      }
    }
    vprint("Constructed weight_row matrix", verbose);
    
    // Extract grouping variables 
    context_lvls = Rcpp::sort_unique(Rcpp::as<CharacterVector>(count_data["context"]));
    species_lvls = Rcpp::sort_unique(Rcpp::as<CharacterVector>(count_data["species"]));
    ran_lvls = Rcpp::sort_unique(Rcpp::as<CharacterVector>(count_data["ran"]));
    // ... add "none" to represent no random effect (reference level)
    if (ran_lvls.size() > 1) {ran_lvls.push_front("none");} 
    // ... print extracted grouping variables
    vprint("Context grouping levels:", verbose);
    vprintV(context_lvls, verbose);
    vprint("Species grouping levels:", verbose);
    vprintV(species_lvls, verbose);
    vprint("Random-effect grouping levels:", verbose);
    vprintV(ran_lvls, verbose);
    
    // Temporarily extract tokenized-count columns 
    IntegerVector binT = Rcpp::as<IntegerVector>(count_data["bin"]); 
    CharacterVector contextT = Rcpp::as<CharacterVector>(count_data["context"]);
    CharacterVector speciesT = Rcpp::as<CharacterVector>(count_data["species"]);
    CharacterVector ranT = Rcpp::as<CharacterVector>(count_data["ran"]);
    
    // Create summed count data, size constants
    int idx = 0;
    int n_ran = ran_lvls.size();
    int n_species = species_lvls.size();
    int n_context = context_lvls.size();
    int bin_num_i = (int)bin_num.val();
    n_count_rows = bin_num_i * n_context * n_species * n_ran * treatment_num;
    count_row_nums = Rcpp::seq(0, n_count_rows - 1);
    vprint("Total rows in summed count data table: " + std::to_string(n_count_rows), verbose);
    
    // Create summed count data rows, initializations
    count.resize(n_count_rows);
    bin.resize(n_count_rows);
    context = CharacterVector(n_count_rows);
    species = CharacterVector(n_count_rows);
    context_num.reserve(n_count_rows);
    species_num.reserve(n_count_rows);
    ran = CharacterVector(n_count_rows);
    treatment = CharacterVector(n_count_rows);
    weights.resize(n_count_rows, treatment_num);
    
    // Initiate count indexes 
    int idx_mcu = 0;
    idx_mc_unique = IntegerVector(n_count_rows/bin_num_i);
    token_pool.resize(n_count_rows);
    count_not_na_mask = LogicalVector(n_count_rows);
    count_not_na_mask.fill(false);
    vprint("Number of rows with unique model components: " + std::to_string(idx_mc_unique.size()), verbose);
    
    // Pre-compute bin masks
    LogicalMatrix bin_masks(binT.size(), bin_num_i);
    for (int b = 0; b < bin_num_i; b++) { 
      bin_masks.column(b) = eq_left_broadcast(binT, b + 1);
    }
    
    // Create summed count data columns and weight matrix
    vprint_header("Creating summed-count data columns and weight matrix", verbose); 
    LogicalVector nan_mask = !Rcpp::is_na(to_NumVec(count_tokenized));
    for (int r = 0; r < n_ran; r++) {
      LogicalVector ran_mask = eq_left_broadcast(ranT, ran_lvls[r]) & nan_mask;
      
      for (int c = 0; c < n_context; c++) {
        LogicalVector context_mask = ran_mask & eq_left_broadcast(contextT, context_lvls[c]);
        
        for (int s = 0; s < n_species; s++) {
          LogicalVector species_mask = context_mask & eq_left_broadcast(speciesT, species_lvls[s]);
          
          for (int t = 0; t < treatment_num; t++) {
            LogicalVector treatment_mask = Rcpp::clone(species_mask);
            for (int f = 0; f < n_fix; f++) {
              treatment_mask = treatment_mask & eq_left_broadcast(Rcpp::as<CharacterVector>(count_data[fix_names[f]]), effects_rows(t, f));
            }
            
            // Save mc-unique index
            idx_mc_unique[idx_mcu] = idx; 
            idx_mcu++;
            
            for (int b = 0; b < bin_num_i; b++) {
              
              // Build columns 
              count(idx) = 0.0;
              bin(idx) = static_cast<sdouble>(b + 1.0);
              context(idx) = context_lvls[c];
              species(idx) = species_lvls[s];
              context_num.push_back(c);
              species_num.push_back(s);
              ran(idx) = ran_lvls[r];
              treatment(idx) = treatment_lvls[t];
              weights.row(idx) = weight_rows.row(t);
              // ^ ... weights (and weight_rows) is a matrix specifying whether the effect from a given treatment level should apply when computing the effect of another treatment level.
              // note: a complex treatment level like "right, KO, P12" is giving an interaction effect, not the effect of any one of its components. 
              
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
    
    // Make extrapolation pool
    if (n_ran > 1) {vprint_header("Making extrapolation pool", verbose);}
    make_extrapolation_pool(verbose); 
    
    // Wrap up initialization
    vprint_header("Wrapping up initialization", verbose);
    
    // Extrapolate "none" rows
    if (ran_lvls.size() > 1) {
      count = extrapolate_none(count, ran, extrapolation_pool, true, round_none);
      vprint("Extrapolated 'none' rows", verbose);
    }
    
    // Take log of observed counts 
    count_log.resize(n_count_rows); 
    for (int r = 0; r < n_count_rows; r++) {
      count_log(r) = slog(count(r) + 1.0);
    } 
    vprint("Took log of observed counts", verbose);
    
    // Compute tpoint buffer
    tpoint_buffer = bin_num * buffer_factor;
    double tpoint_buffer_d = tpoint_buffer.val();
    int tpoint_buffer_int = (int)tpoint_buffer_d;
    if (tpoint_buffer_d != tpoint_buffer_int) {tpoint_buffer_int++;} 
    tpoint_buffer = (sdouble)tpoint_buffer_int;
    
    // Estimate gamma dispersion of raw counts
    compute_gamma_dispersion();
    vprint("Estimated gamma dispersion of raw counts", verbose);
    
    // Find average log counts for each context-species combination
    find_count_log_means();
    vprint("Found average log counts for each context-species combination", verbose);
    
    // Construct grouping variable ids as indexes for warping factor matrices, by count row
    gv_ran_idx.reserve(n_count_rows);
    gv_sps_idx.reserve(n_count_rows);
    for (int r = 0; r < n_count_rows; r++) {
      CharacterVector::iterator it_ran = std::find(ran_lvls.begin(), ran_lvls.end(), ran[r]);
      CharacterVector::iterator it_fix = std::find(species_lvls.begin(), species_lvls.end(), species[r]);
      gv_ran_idx[r] = std::distance(ran_lvls.begin(), it_ran);
      gv_sps_idx[r] = std::distance(species_lvls.begin(), it_fix);
    }
    vprint("Constructed grouping variable IDs", verbose);
    
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
 * Initial setup
 */

// Computes gamma_dispersion matrix and related vectors
void wspc::compute_gamma_dispersion() {
    
    // Grab number of context and species levels
    int n_species = species_lvls.size();
    int n_context = context_lvls.size();
    
    // Initialize gamma dispersion matrix
    gamma_dispersion = sMat(n_species, n_context);
    // ... initialized with all zeros
    gd_species_idx = IntegerVector(n_species);
    gd_context_idx = IntegerVector(n_context);
    gd_species_idx.names() = species_lvls;
    gd_context_idx.names() = context_lvls;
    
    // Loop through context levels
    for (int c = 0; c < n_context; c++) {
      
      LogicalVector context_mask = eq_left_broadcast(context, context_lvls[c]);
      gd_context_idx[(String)context_lvls[c]] = c;
      
      // Loop through species levels
      for (int s = 0; s < n_species; s++) { 
        
        // Estimate dispersion of raw count (not log)
        LogicalVector species_mask = eq_left_broadcast(species, species_lvls[s]);
        LogicalVector cs_mask = context_mask & species_mask & count_not_na_mask;
        sVec count_cs_masked = masked_vec(count, cs_mask); 
        sdouble count_cs_mean = vmean(count_cs_masked);
        sdouble count_cs_var = vvar(count_cs_masked);
        gamma_dispersion(s, c) = gamma_dispersion_formula(count_cs_mean, count_cs_var);
        gd_species_idx[(String)species_lvls[s]] = s;
        
      }
      
    }
    
  }

void wspc::make_extrapolation_pool(bool verbose) {
    
    // Construct nan and none masks
    LogicalVector nan_mask = !Rcpp::is_na(to_NumVec(count));
    LogicalVector none_mask = !eq_left_broadcast(ran, "none");
    
    // Find indexes of none rows
    IntegerVector r_idx = Rwhich(!none_mask);
    
    // Resize vector to hold extrapolation pools
    extrapolation_pool = std::vector<IntegerVector>(count.size());
    
    // Convert bin vector to numeric and get max bin number
    NumericVector bin_NumVec = to_NumVec(bin);
    int max_bin = Rcpp::max(bin_NumVec);
    
    // Pre-make masks for bin, context, species, and treatment levels
    LogicalMatrix bin_masks(bin_NumVec.size(), max_bin);
    LogicalMatrix context_masks(context.size(), context_lvls.size());
    LogicalMatrix species_masks(species.size(), species_lvls.size());
    LogicalMatrix treatment_masks(treatment.size(), treatment_lvls.size());
    colnames(context_masks) = context_lvls;
    colnames(species_masks) = species_lvls;
    colnames(treatment_masks) = treatment_lvls;
    for (int i = 0; i < bin_masks.ncol(); i++) {bin_masks.column(i) = eq_left_broadcast(bin_NumVec, i + 1);}
    for (int i = 0; i < context_lvls.size(); i++) {context_masks.column(i) = eq_left_broadcast(context, context_lvls[i]);}
    for (int i = 0; i < species_lvls.size(); i++) {species_masks.column(i) = eq_left_broadcast(species, species_lvls[i]);}
    for (int i = 0; i < treatment_lvls.size(); i++) {treatment_masks.column(i) = eq_left_broadcast(treatment, treatment_lvls[i]);}
    
    // Loop through none rows and find their interpolation pools
    int n_rows = r_idx.size();
    IntegerVector tracker = iseq((int)(n_rows/5 - 1), n_rows - 1, 5); 
    for (int ri = 0; ri < n_rows; ri++) {
      int r = r_idx[ri];
      
      // Find all rows with the same fixed effects and bin
      LogicalVector mask = bin_masks.column(bin(r).val() - 1)
        & context_masks.column(Rwhich(eq_left_broadcast(context_lvls, context(r)))[0])
        & species_masks.column(Rwhich(eq_left_broadcast(species_lvls, species(r)))[0])
        & treatment_masks.column(Rwhich(eq_left_broadcast(treatment_lvls, treatment(r)))[0])
        & nan_mask
        & none_mask;
     
      // Extract pool 
      extrapolation_pool[r] = Rwhich(mask);
      
      // Track progress
      if (any_true(eq_left_broadcast(tracker, ri)) || n_rows <= 5) {
        vprint("row: " + std::to_string(ri + 1) + "/" + std::to_string(n_rows), verbose);
      } 
      
    }
  }

// Function to take means of count_log 
void wspc::find_count_log_means() {
    
    // Convert bin num to integer
    int bin_num_i = (int)bin_num.val();
    
    // Grab number of context and species levels
    const int n_species = species_lvls.size();
    const int n_context = context_lvls.size();
    
    // Initialize lists to hold count_log_avg_mat for each context-species combination 
    count_log_avg_mat_list = List(n_context);
    count_log_avg_mat_list.names() = context_lvls;
    
    // Loop through context levels
    for (int c = 0; c < n_context; c++) {
      
      // Make context mask
      LogicalVector context_mask = eq_left_broadcast(context, context_lvls[c]);
      
      // Set up lists for count_log_avg_mat for this context
      count_log_avg_mat_list[(String)context_lvls[c]] = List(n_species);
      name_proxylist(count_log_avg_mat_list[(String)context_lvls[c]], species_lvls);
      
      // Loop through species levels
      for (int s = 0; s < n_species; s++) { 
        
        // Make species mask and context-species mask
        LogicalVector species_mask = eq_left_broadcast(species, species_lvls[s]);
        LogicalVector cs_mask = context_mask & species_mask & count_not_na_mask;
        
        // Find mean of count_log for each species-context pair
        NumericMatrix count_log_avg_mat(bin_num_i, treatment_num);
        for (int t = 0; t < treatment_num; t++) {
          String trt = treatment_lvls[t];
          
          // Make mask for treatment rows of this context-species pair
          LogicalVector trt_mask = eq_left_broadcast(treatment, trt);
          LogicalVector cs_trt_mask = cs_mask & trt_mask;
          
          // Grab average counts for this treatment level (used later to set initial values)
          dVec count_log_avg(bin_num_i); 
          sVec count_trt_masked = masked_vec(count_log, cs_trt_mask);
          sVec bin_trt_masked = masked_vec(bin, cs_trt_mask);
          for (int b = 0; b < bin_num_i; b++) {
            LogicalVector mask_b = eq_left_broadcast(to_NumVec(bin_trt_masked), (double)b + 1.0);
            sVec count_b = masked_vec(count_trt_masked, mask_b);
            count_log_avg[b] = vmean(count_b).val();
          }
          count_log_avg_mat.column(t) = to_NumVec(count_log_avg);
          
        }
       
        // Assign results matrix to the context-species list
        assign_proxylist(count_log_avg_mat_list[(String)context_lvls[c]], (String)species_lvls[s], count_log_avg_mat);
        
      }
    }
    
  }

// Function to estimate change points
void wspc::estimate_change_points() {
   
    // Convert bin number 
    int bin_num_i = (int)bin_num.val();
    
    // Grab number of context and species levels
    const int n_ran = ran_lvls.size();
    const int n_species = species_lvls.size();
    const int n_context = context_lvls.size();
    const int n_ran_trt = n_ran * treatment_num;
    
    // Initialize matrix to hold degrees of each context-species combination
    degMat = IntegerMatrix(n_species, n_context);
    rownames(degMat) = species_lvls;
    colnames(degMat) = context_lvls;
    
    // Initialize lists to hold results matrices for each context-species combination 
    // ... for found_cp
    found_cp_list = List(n_context);
    found_cp_list.names() = context_lvls;
    // ... for found_cp_trt
    found_cp_trt_list = List(n_context);
    found_cp_trt_list.names() = context_lvls;
    
    // Loop through context levels
    for (int c = 0; c < n_context; c++) {
      
      // Make context mask
      LogicalVector context_mask = eq_left_broadcast(context, context_lvls[c]);
      
      // Set up lists for results matrices for this context
      // ... for found_cp
      found_cp_list[(String)context_lvls[c]] = List(n_species);
      name_proxylist(found_cp_list[(String)context_lvls[c]], species_lvls);
      // ... for found_cp_trt
      found_cp_trt_list[(String)context_lvls[c]] = List(n_species);
      name_proxylist(found_cp_trt_list[(String)context_lvls[c]], species_lvls);
      
      // Loop through species levels
      for (int s = 0; s < n_species; s++) { 
        
        // Make species mask and context-species mask
        LogicalVector species_mask = eq_left_broadcast(species, species_lvls[s]);
        LogicalVector cs_mask = context_mask & species_mask & count_not_na_mask;
        
        // Initialize count array for change-point detection
        sMat count_masked_array(bin_num_i, n_ran_trt);
        count_masked_array.setZero();
        LogicalVector good_col(n_ran_trt);
        // Collect count values for each treatment-ran level interaction of this species-context pair
        for (int t = 0; t < treatment_num; t++) {
          String trt = treatment_lvls[t];
          
          // Make mask for treatment rows of this context-species pair
          LogicalVector trt_mask = eq_left_broadcast(treatment, trt);
          LogicalVector cs_trt_mask = cs_mask & trt_mask;
          
          // Collect count values for each ran level and this treatment trt
          for (int r = 0; r < n_ran; r++) {
            // Make mask for ran level rows of this treatment (of this context-species pair)
            LogicalVector ran_mask = eq_left_broadcast(ran, ran_lvls[r]);
            LogicalVector mask = cs_trt_mask & ran_mask;
            
            // Make masked copies of count_log and bin
            sVec count_masked = masked_vec(count_log, mask);  
            sVec bin_masked = masked_vec(bin, mask);
            if (count_masked.size() == bin_num_i && bin_masked.size() == bin_num_i) {
              
              // Ensure count_masked is in correct order
              // ... should be, but sanity check 
              for (int b = 0; b < bin_num_i; b++) {
                if (bin_masked[b] != (b + 1.0)) {
                  stop("Count or bin vectors not in correct order.");
                }
              }
              
              // Set this column and mark good
              count_masked_array.col(t*n_ran + r) = count_masked;
              good_col(t*n_ran + r) = true;
              
            } else {
              good_col(t*n_ran + r) = false;
            }
          }
          
        }
        
        // Extract good column numbers 
        IntegerVector good_col_idx = Rwhich(good_col);
        sMat count_masked_array_good = count_masked_array(Eigen::all, to_iVec(good_col_idx));
        
        // Estimate change points from masked count series
        IntegerMatrix found_cp_good = LROcp_array(
          count_masked_array_good,    // 2D matrix of points to test for change points (columns as series, rows as bins)
          ws,                         // running window size 
          LROcutoff,                  // points more than this times sd considered outliers
          tpoint_buffer.val()         // Minimum distance between two change points
        );
        
        // Estimate degree of this context-species pair 
        int deg = found_cp_good.rows();
        degMat(s, c) = deg;
        
        // Fill columns into the found_cp matrix
        IntegerMatrix found_cp(deg, n_ran_trt);
        // ^ ... Rcpp should initialize these matrices with all zeros
        if (deg > 0) {
          for (int si = 0; si < good_col_idx.size(); si++) {
            // ... grab change points
            found_cp.column(good_col_idx[si]) = found_cp_good.column(si);
          }
        }
        
        // Extract treatment means for each change point
        NumericMatrix found_cp_trt(deg, treatment_num);
        for (int t = 0; t < treatment_num; t++) {
          for (int d = 0; d < deg; d++) {
            found_cp_trt(d, t) = 0.0;
            int n_ran_hit = 0;
            for (int r = 0; r < n_ran; r++) {
              if (good_col(t*n_ran + r)) { // ... ensure there is data here
                found_cp_trt(d, t) += (double)found_cp(d, t*n_ran + r);
                n_ran_hit++;
              }
            }
            found_cp_trt(d, t) = found_cp_trt(d, t) / (double)n_ran_hit;
          }
        }
        
        // Assign results matrices to the context-species list
        assign_proxylist(found_cp_list[(String)context_lvls[c]], (String)species_lvls[s], found_cp);
        assign_proxylist(found_cp_trt_list[(String)context_lvls[c]], (String)species_lvls[s], found_cp_trt);
        
      }
    }
    
  }

// Estimate change points and initial parameters for model fitting
void wspc::LRO_initial_param_ests(
    bool verbose,
    double LROwf,
    double LROco
  ) {
   
    if (LROwf == 0.0) {LROwf = LROwindow_factor;}
    if (LROco != 0.0) {LROcutoff = LROco;}
    
    // Compute running and filter window sizes for LRO change-point detection
    ws = static_cast<int>(std::round(LROwf * bin_num.val() * buffer_factor.val()));
    
    // Estimate degree of each context-species combination at baseline using LRO change-point detection 
    estimate_change_points();
    vprint("Estimated change points", verbose);
    
    // Find initial parameter estimates
    estimate_initial_parameters();
    vprint("Estimated initial parameters", verbose);
    vprint("Number of parameters: " + std::to_string((int)fitted_parameters.size()), verbose);
    
    // Compute size of the parameter boundary vector 
    boundary_vec_size = 0;
    for (int r : idx_mc_unique) {
      // Grab degree for this row
      int deg = degMat(species_num[r], context_num[r]);
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
    vprint("Size of boundary vector: " + std::to_string(boundary_vec_size), verbose);
    
    // Check for any treatment-level knockouts and set effects to zero in param vector
    for (String this_trtKO : trtKO) {fitted_parameters[grep_cpp(param_names, "_" + ((std::string)this_trtKO) + "_")] = 0.0;}
    
    // Check feasibility 
    List feasibility_results = check_parameter_feasibility(to_sVec(fitted_parameters), false); 
    fitted_parameters = Rcpp::as<NumericVector>(feasibility_results["parameters"]);
    
  }

// Search for best LRO change-point detection settings
NumericMatrix wspc::LRO_grid_search(bool verbose) {
    
    if (verbose) {
      vprint_header("Performing grid search to find optimal LRO parameters");
    }
    dVec LROwf_range = {1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 3.00};
    dVec LROco_range = LROwf_range;
    NumericMatrix results(LROwf_range.size() * LROco_range.size(), 8);
    int n = count_not_na_idx.size();
    
    for (int i = 0; i < LROwf_range.size(); ++i) {
      for (int j = 0; j < LROco_range.size(); ++j) {
        if (verbose) {
          vprint("LRO window factor: " + std::to_string(LROwf_range[i]) + ", LRO cutoff: " + std::to_string(LROco_range[j]));
        }
        LRO_initial_param_ests(false, LROwf_range[i], LROco_range[j]);
        fit(false, false);
        int k = fitted_parameters.size();
        int idx = i * LROco_range.size() + j;
        double nll = optim_results["neg_loglik"];
        double bnll = optim_results["penalized_neg_loglik"];
        results(idx, 0) = LROwf_range[i];
        results(idx, 1) = LROco_range[j];
        results(idx, 2) = nll;
        results(idx, 3) = bnll;
        results(idx, 4) = (double)optim_results["success_code"];
        results(idx, 5) = (double)k;
        results(idx, 6) = 2.0 * ((double)k + nll);
        results(idx, 7) = 2.0 * (std::log((double)n) * (double)k + nll);
        clear_stan_mem();
      }
    }
    
    colnames(results) = CharacterVector::create(
      "LROwindow_factor", 
      "LROcutoff", 
      "neg_loglik", 
      "penalized_neg_loglik",
      "success_code", 
      "n_params",
      "AIC", 
      "BIC"
      );
    
    return results;
    
  }   

/*
 * *************************************************************************
 * Computing predicted values from parameters
 */

// Compute model component values for rows of summed count data
sVec wspc::compute_warped_mc(
    int mc,                    // Model component number for which to compute values, 0 = Rt, 1 = tslope, 2 = tpoint
    int r,                     // Row of summed count data for which to compute model component 
    const sVec& parameters,    // Parameters to use in computation
    const sdouble& wf          // Warping factor to apply 
  ) const {
    
    // Grab degree
    int bktp = degMat(species_num[r], context_num[r]);
    if (mc == 0) {bktp++;} // ... 0 is Rt
    
    // Extract weight matrix row
    sVec weight_row = weights.row(r).transpose();
    
    // Compute unwarped model component value  
    sVec mc_vec(bktp);
    mc_vec.setZero();
    for (int bt = 0; bt < bktp; bt++) {
      for (int t = 0; t < treatment_num; t++) {
        mc_vec(bt) += weight_row(t) * parameters(beta_idx[mc][context_num[r]][species_num[r]][bt*treatment_num + t]);
      }
      // ... and apply warping factor
      mc_vec(bt) = warp_mc(mc_vec(bt), warp_bounds[mc], wf);
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
        int wfi = gv_sps_idx[r] * (ran_lvls.size() - 1) + gv_ran_idx[r];
        if (gv_ran_idx[r] == 0) {
          f_pw = 0.0;
          f_rw = 0.0; 
          f_sw = 0.0;
        } else {
          f_pw = parameters(wfactor_idx[0][wfi]);
          f_rw = parameters(wfactor_idx[1][wfi]); 
          f_sw = parameters(wfactor_idx[2][wfi]);
        } 
        
        // Only update predicted model components if r begins a new batch of unique values 
        int cnt = std::count(mc_unique_rows.begin(), mc_unique_rows.end(), r);
        if (cnt > 0) { 
          
          // Compute warped model components for this row r
          Rt = compute_warped_mc("Rt", r, parameters, f_rw);        
          tslope = compute_warped_mc("tslope", r, parameters, f_sw); 
          tpoint = compute_warped_mc("tpoint", r, parameters, f_pw);
          
        } 
        
        // Compute the actual poly-sigmoid!! 
        predicted_rates(r) = poly_sigmoid(
          bin(r),
          degMat(species_num[r], context_num[r]);,
          Rt,
          tslope,
          tpoint
        ); 
        
      }
      
    }
    
    return predicted_rates;
    
  }

// Predict log of rates, R wrapper 
NumericVector wspc::predict_rates_R(
    const NumericVector& parameters_R,
    const bool& all_rows 
  ) {
    
    // Convert parameters to sVec
    sVec parameters = to_sVec(parameters_R);
    
    // Compute predicted rates
    sVec predicted_rates = predict_rates(
      parameters,
      all_rows
    );
    
    // Convert to NumericVector 
    NumericVector predicted_rates_R = to_NumVec(predicted_rates);
    
    // Clear memory and return
    clear_stan_mem();
    return predicted_rates_R;
    
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
        // ... grab context and species index numbers
        int n_c = gd_species_idx[(String)species[r]];
        int n_p = gd_context_idx[(String)context[r]];
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
    sdouble f_pw;
    sdouble f_rw;
    sdouble f_sw;
    
    // Compute the boundary distance, for ...
    // ... transition slopes (must be positive)
    // ... transition points (enforces tpoint buffer)
    // ... block rate values (must be positive)
    for (int r : idx_mc_unique) {
      
      // Grab warping factors
      int wfi = gv_sps_idx[r] * (ran_lvls.size() - 1) + gv_ran_idx[r];
      if (gv_ran_idx[r] == 0) {
        f_pw = 0.0;
        f_rw = 0.0; 
        f_sw = 0.0;
      } else {
        f_pw = parameters(wfactor_idx[0][wfi]);
        f_rw = parameters(wfactor_idx[1][wfi]); 
        f_sw = parameters(wfactor_idx[2][wfi]);
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
          int wfi = gv_sps_idx[r] * (ran_lvls.size() - 1) + gv_ran_idx[r];
          f_pw = parameters(wfactor_idx_point[wfi]);
        }
        
        // Find degree
        int deg = degMat(species_num[r], context_num[r]);
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
          // ... and first point must be >= buffer, 
          // ... and last point must be <= bin_num - buffer
          for (int d = 0; d < deg + 1; d++) {
            sdouble buffer_dist = (tpoint_ext(d + 1) - tpoint_ext(d)) - tpoint_buffer;
            if (buffer_dist < 0.0) {
              if (verbose) {
                vprint("Found tpoint below buffer, distance: " + std::to_string(buffer_dist.val()), true);
                vprint("   deg: " + std::to_string(d) + ", row: " + std::to_string(r), true);
                vprint("   tpoint_ext(d + 1): " + std::to_string(tpoint_ext(d + 1).val()) + ", tpoint_ext(d): " + std::to_string(tpoint_ext(d).val()), true);
                vprint("   species: " + std::string(species[r]) + ", context: " + std::string(context[r]) + ", ran level: " + std::string(ran[r]) + ", treatment: " + std::string(treatment[r]), true);
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
    
    // Compute min_boundary_dist
    sdouble fx = model->min_boundary_dist(parameters_var);
    if (!std::isfinite(fx.val())) {
      Rcpp::Rcout << "Non-finite objective from min_boundary_dist_NLopt" << std::endl;
    }
    
    // Compute gradient if needed
    if (!grad.empty()) {
      Eigen::VectorXd grad_eigen = model->grad_min_boundary_dist(parameters_var);
      grad.assign(grad_eigen.data(), grad_eigen.data() + grad_eigen.size());
    }  
    for (size_t i = 0; i < grad.size(); ++i) {
      if (!std::isfinite(grad[i])) {
        grad[i] = 0.0;
      }
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

// Compute the gradient of the min_boundary_dist function
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
    vprint("Setting boundary penalty coefficients", verbose);
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
    
    // Find final neg_loglik for diagnostics
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
    optim_results["fitted_parameters"] = to_NumVec(x); // for predicting log-linked count
    optim_results["penalized_neg_loglik"] = min_fx;
    optim_results["neg_loglik"] = final_nll.val();
    optim_results["success_code"] = success_code;
    optim_results["num_evals"] = opt.get_numevals();
    
  } 

// Resample count for bootstrapping
void wspc::resample(
    unsigned int seed
  ) {
   
    // Initialize random number generator
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
        // ... re-sum into r's count
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
    count = extrapolate_none(count, ran, extrapolation_pool, true, round_none);
    iVec r_rows = Rcpp::as<iVec>(count_row_nums[eq_left_broadcast(ran,"none")]);
    for (int r : r_rows) {count_log(r) = slog(count(r) + 1.0);}
    
  }

// Fit model to bootstrap resample
dVec wspc::bs_fit(
    int bs_num,               // A unique number for this resample
    bool clear_stan           // Recover stan memory at end?
  ) {
    
    // Resample counts 
    resample(rng_seed + bs_num);
    
    // Estimate gamma dispersion of these new raw re-sampled counts
    compute_gamma_dispersion();
    
    // Find average these new re-sampled log counts for each context-species combination
    find_count_log_means();
    
    // Find initial parameter estimates for fixed-effect treatments, for new re-sampled data
    estimate_initial_parameters();
    
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
    const int batch_num = std::round((double)bs_num_max / (double)max_fork);
    int tracker_steps = 10;
    IntegerVector tracker = iseq(batch_num/tracker_steps, batch_num, tracker_steps);
    for (int b = 0; b < batch_num; b++) {
      
      // Initiate timer and grab start time
      Timer batch_timer;
      batch_timer.step("start");  
      
      if (max_fork > 1) {
        // Run in parallel with forking
        int batch_size = std::min(max_fork, bs_num_max - b * max_fork);
        
        // Pipes for inter-process communication
        iVec pids(batch_size);
        std::vector<std::array<int, 2>> pipes(batch_size); 
        
        // Initialize pipes and fork processes
        for (int i = 0; i < batch_size; i++) {
          
          int this_row = b * max_fork + i;
          pipe(pipes[i].data()); // Create a pipe
          pid_t pid = fork();
          
          if (pid == 0) { // child process
            
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
            
          } else if (pid > 0) { // parent process
            
            // Grab child pid
            pids[i] = pid;
            
            // Close write end
            close(pipes[i][1]); 
            
          } else {
            Rcpp::stop("Fork failed!");
          }
          
        }
        
        // Fetch results from pipes
        for (int i = 0; i < batch_size; i++) {
          
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
    bool start_from_fit,      // Start from parameters found with gradient descent? 
    bool verbose
  ) {
    
    // Reset initial parameters to check their feasibility 
    set_parameters(fitted_parameters, verbose);
    
    // Fit full model 
    vprint("Performing initial fit of full data", verbose);
    fit(
      start_from_fit,   // set fitted parameters?
      false             // don't print anything
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
          int param_oom = static_cast<int>(std::floor(std::log10(std::abs(params_current(i)))));
          double normalized_step_size = step_size * std::pow(10, static_cast<double>(param_oom + 1.0));
          double bounded_step_size = normalized_step_size / bd_current_transformed.val();
          if (bounded_step_size == 0.0) {
            // ... presumably this case means current parameter is extremely close to zero or very close to boundary
            params_next(i) = pcg_rnorm(params_current(i), step_size/10, walk_rng);
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
            vprint("step: " + std::to_string((this_step_batch + 1) * (n_steps/tracker.size())) + "/" + std::to_string(n_steps), verbose);
            tracker(this_step_batch) -= 1; // ensure report is only printed once
          } else if (step == 0 && !printed_step1) {
            vprint("step: 1/" + std::to_string(n_steps), verbose);
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
      if (inf_loop_ctr > 10000) {
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

/*
 * *************************************************************************
 * Setting parameters
 */

// Function for converting structured encoding of parameters into a vector
void wspc::estimate_initial_parameters() {
    
    /*
     * Initializations
     */
    
    // Initialize vector to hold all parameter values (faster to use std and wrap at end)
    dVec param_vector(0);
    
    // Initialize List to hold names
    List param_names_list;
    
    int idx = 0;
    int n_context = context_lvls.size();
    int n_species = species_lvls.size();
    int n_ran = ran_lvls.size();
   
    /*
     * Beta parameters (ref level and fixed effects)
     */
    
    // Loop through model components 
    for (String mc : mc_list) {
      
      // Loop through context values
      std::vector<std::vector<std::vector<int>>> beta_idx_mc;
      for (int c = 0; c < n_context; c++) {
        String cxt = context_lvls[c];
        
        // Make context mask
        LogicalVector context_mask = eq_left_broadcast(context, cxt);
        
        // Extract found_cp, found_cp_trt, and count_log_avg_mat for this context level
        List found_cp_list_c = found_cp_list[cxt];
        List found_cp_trt_list_c = found_cp_trt_list[cxt];
        List count_log_avg_mat_list_c = count_log_avg_mat_list[cxt];
        
        // Loop through species values
        std::vector<std::vector<int>> beta_idx_mc_cxt;
        for (int s = 0; s < n_species; s++) {
          String sps = species_lvls[s];
          
          // Extract deg and block num 
          int deg = degMat(s, c);
          int n_blocks = deg + 1;
          int bktp = deg;
          if (mc == "Rt") {bktp = b_blocks;}
          
          // Initialize sMat to hold Rt values 
          sMat RtVals(treatment_num, n_blocks);
          
          NumericMatrix Effs(treatment_num, bktp);
          NumericMatrix run_estimates(treatment_num, deg);
          
          // Extract found_cp, found_cp_trt, and count_log_avg_mat for this context-species pair
          IntegerMatrix found_cp = found_cp_list_c[sps];
          NumericMatrix found_cp_trt = found_cp_trt_list_c[sps];
          NumericMatrix count_log_avg_mat = count_log_avg_mat_list_c[sps];
          
          // Make species mask and context-species mask
          LogicalVector species_mask = eq_left_broadcast(species, species_lvls[s]);
          LogicalVector cs_mask = context_mask & species_mask & count_not_na_mask;
          
          if (deg > 0 || mc == "Rt") { 
            
            // Compute Effs for this model component, context, and species
            if (mc == "Rt") {
              // Estimate Rt values
              for (int i = 0; i < treatment_num; i++) {
                NumericVector count_log_avg = count_log_avg_mat.column(i);
                NumericVector found_cp_trt_i_num = found_cp_trt.column(i);
                IntegerVector found_cp_trt_i(deg);
                for (int d = 0; d < deg; d++) {
                  found_cp_trt_i[d] = std::round(found_cp_trt_i_num[d]);
                }
                dVec Rt_est(n_blocks); 
                if (n_blocks == 1) { // case when deg = 0
                  Rt_est[0] = vmean(count_log_avg);
                  RtVals(i, 0) = (sdouble)Rt_est[0];
                } else {
                  // Estimate rate values as mean of count values in each block 
                  std::vector<dVec> est_rate_runs = est_bkRates_tRuns(n_blocks, count_log_avg, found_cp_trt_i, rise_threshold_factor);
                  Rt_est = est_rate_runs[0];
                  run_estimates.row(i) = to_NumVec(est_rate_runs[1]);
                  RtVals.row(i) = to_sVec(Rt_est); 
                }
              }
              // Estimate Rt effect values 
              for (int bk = 0; bk < n_blocks; bk++) {Effs.column(bk) = to_NumVec(weight_rows.fullPivLu().solve(RtVals.col(bk)));}
            } else if (mc == "tpoint") {
              // Estimate fixed effects from treatments for each tpoint
              sMat found_cp_trt_transposed = to_sMat(found_cp_trt).transpose(); 
              for (int d = 0; d < deg; d++) {Effs.column(d) = to_NumVec(weight_rows.fullPivLu().solve(found_cp_trt_transposed.col(d)));}
            } else if (mc == "tslope") {
              // Estimate tslopes for each treatment
              NumericMatrix found_slope_trt(deg, treatment_num);
              for (int i = 0; i < treatment_num; i++) {
                for (int d = 0; d < deg; d++) {
                  found_slope_trt(d, i) = 4.0/run_estimates(i, d); 
                  if (found_slope_trt(d, i) < min_initialization_slope) {found_slope_trt(d, i) = min_initialization_slope;}
                  // ^ ... keep from initializing too close to zero
                }
              }
              // Estimate fixed effects from treatments for each tslope
              sMat found_slope_trt_transposed = to_sMat(found_slope_trt).transpose();
              for (int d = 0; d < deg; d++) {Effs.column(d) = to_NumVec(weight_rows.fullPivLu().solve(found_slope_trt_transposed.col(d)));}
              
            }
            
            // For each block/transition point t and treatment i
            std::vector<int> beta_idx_mc_cxt_sps;
            for (int t = 0; t < bktp; t++) {
              String t_name = "Tns/Blk" + std::to_string(t + 1);
              for (int i = 0; i < treatment_num; i++) {
                
                // Add name
                CharacterVector param_name;
                if (i == 0) {param_name = CharacterVector::create("baseline", cxt, mc, sps, t_name);} 
                else {param_name = CharacterVector::create("beta", mc, cxt, sps, treatment_lvls[i], "X", t_name);}
                param_names_list.push_back(param_name);
               
                // Collect indices 
                if (i > 0) {
                  if (mc == "Rt") {
                    param_beta_Rt_idx.push_back(idx);
                  } else if (mc == "tslope") {   
                    param_beta_tslope_idx.push_back(idx);
                  } else if (mc == "tpoint") {   
                    param_beta_tpoint_idx.push_back(idx);
                  }   
                } else {
                  param_baseline_idx.push_back(idx);
                  if (mc == "Rt") {
                    param_baseline_Rt_idx.push_back(idx);
                  } else if (mc == "tslope") {   
                    param_baseline_tslope_idx.push_back(idx);
                  } else if (mc == "tpoint") {   
                    param_baseline_tpoint_idx.push_back(idx);
                  }
                }
                beta_idx_mc_cxt_sps.push_back(idx);
                idx++;  
                
                // Make parameter vector
                param_vector.push_back(Effs(i,t));
                
              } 
            }
          } 
          beta_idx_mc_cxt.pushback(beta_idx_mc_cxt_sps);
        } 
        beta_idx_mc.push_back(beta_idx_mc_cxt);
      }
      beta_idx.push_back(beta_idx_mc);
    }
    
    /*
     * Random-effect parameters
     */
    
    // Unroll the warping factors
    for (int s = 0; s < n_species; s++) {
      // Intentionally skip the first level, which is "none"
      for (int r = 1; r < n_ran; r++) {
        
        // Map
        String c_name = species_lvls[s];
        String r_name = ran_lvls[r];
        // ... Point warp
        CharacterVector param_name_point = CharacterVector::create("wfactor", "point", r_name, "X", c_name);
        param_names_list.push_back(param_name_point);
        param_wfactor_point_idx.push_back(idx);
        wfactor_idx_point.push_back(idx);
        idx++;
        // ... Rate warp
        CharacterVector param_name_rate = CharacterVector::create("wfactor", "rate", r_name, "X", c_name);
        param_names_list.push_back(param_name_rate);
        param_wfactor_rate_idx.push_back(idx);
        wfactor_idx_rate.push_back(idx); 
        idx++;
        // ... Slope warp
        CharacterVector param_name_slope = CharacterVector::create("wfactor", "slope", r_name, "X", c_name);
        param_names_list.push_back(param_name_slope);
        param_wfactor_slope_idx.push_back(idx);
        wfactor_idx_slope.push_back(idx);
        idx++;
       
        // Make
        // ... Point warp
        param_vector.push_back(R::runif(-0.1, 0.1));
        // ... Rate warp
        param_vector.push_back(R::runif(-0.1, 0.1));
        // ... Slope warp
        param_vector.push_back(R::runif(-0.1, 0.1));
        
      }
    } 
    wfactor_idx = {wfactor_idx_point, wfactor_idx_rate, wfactor_idx_slope};
    
    // Reformat parameter names as single strings
    int n_params = param_names_list.size();
    param_names = CharacterVector(n_params);
    for (int n = 0; n < n_params; n++) {
      CharacterVector name_comps = Rcpp::as<CharacterVector>(param_names_list[n]);
      int m = name_comps.size();
      if (m == 0) {Rcpp::stop("Empty parameter name!");}
      String name = name_comps[0]; 
      if (m > 1) {
        for (int i = 1; i < m; i++) {
          name += "_" + name_comps[i];
        }
      }
      param_names[n] = name;
    } 
    
    // Wrap parameter vector
    fitted_parameters = Rcpp::wrap(param_vector);
    
  }

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
      Rcpp::stop("Negative or nan rates predicted, or other boundary condition violated!");
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
    
    if (feasible) {
      if (verbose) {vprint("... no tpoints below buffer");}
      
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
    
    // Report feasibility of provided parameters
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
      vprint("Initial boundary distance (want > 0): ", initial_dist);
    }
    
    // If not feasible, attempt to find a feasible starting point
    if (!feasible) {
      
      // Variables for optimization
      dVec x = to_dVec(parameters_var);
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
      int n_tries = 0;
      int max_tries = 6;
      while (n_tries < max_tries && success_code == 0) {
        // Manually reduce t-point effects 
        for (int i : param_beta_tpoint_idx) {x[i] *= 0.5;}
        for (int i : param_wfactor_point_idx) {x[i] *= 0.5;}
        // Manually check baseline t-points 
        for (int i = 1; i < param_baseline_tpoint_idx.size(); i++) {
          int idx0 = param_baseline_tpoint_idx[i-1];
          int idx1 = param_baseline_tpoint_idx[i];
          if (x[idx1] - x[idx0] < tpoint_buffer) {
            x[idx0] -= tpoint_buffer.val()/4.0;
            x[idx1] += tpoint_buffer.val()/4.0;
          }
        }
        // Run nlopt 
        try {
          nlopt::result sc = opt.optimize(x, max_fx);
          success_code = static_cast<int>(sc);
        } catch (std::exception& e) {
          if (verbose && n_tries == max_tries - 1) {
            Rcpp::Rcout << "Optimization failed: " << e.what() << std::endl;
          }
          success_code = 0;
        }
        n_tries++;
      }
      
      // Find and report final distance to boundary
      if (verbose) {
        vprint("Number of tries: ", n_tries);
        vprint("Number of evals (last try): ", (int)opt.get_numevals());
        vprint("Success code (last try): ", (int)success_code);
        vprint("Final boundary distance (last try): ", (double)max_fx);
      }
      
      // Final check of boundary distance
      if (max_fx <= 0) {success_code = 0;}
      
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
      _["bin"] = to_NumVec(bin),
      _["count"] = to_NumVec(count),
      _["pred"] = predicted_rates_out,
      _["count.log"] = to_NumVec(count_log),
      _["pred.log"] = predicted_rates_log_out,
      _["context"] = context,
      _["species"] = species,
      _["ran"] = ran,
      _["treatment"] = treatment
    );
    
    // Add parameter names to fitted parameter vector 
    fitted_parameters.names() = param_names;
    
    // Collect fixed-effect names and levels 
    List fix_lvls_list = List(fix_lvls.size());
    for (int i = 0; i < fix_lvls.size(); i++) {
      fix_lvls_list[i] = (CharacterVector)fix_lvls[i];
    }
    List fix_trt_list = List(fix_trt.size());
    for (int i = 0; i < fix_trt.size(); i++) {
      fix_trt_list[i] = (CharacterVector)fix_trt[i];
    }
    fix_lvls_list.names() = fix_names;
    fix_trt_list.names() = fix_names;
    fix_ref.names() = fix_names;
    List fixed_effects = List::create(
      _["names"] = fix_names,
      _["lvls"] = fix_lvls_list,
      _["treat.lvl"] = fix_trt_list,
      _["ref.lvl"] = fix_ref
    );
    
    // Pack up treatment name information
    List treatment_components_list = List(treatment_components.size());
    for (int i = 0; i < treatment_components.size(); i++) {
      treatment_components_list[i] = (CharacterVector)treatment_components[i];
    }
    treatment_components_list.names() = treatment_lvls;
    List treat = List::create(
      _["names"] = treatment_lvls, 
      _["components"] = treatment_components_list
    );
    
    // Put grouping variable information into list
    List grouping_variables = List::create(
      _["context.lvls"] = context_lvls,
      _["species.lvls"] = species_lvls,
      _["ran.lvls"] = ran_lvls
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
    rownames(g_dispersion) = species_lvls;
    colnames(g_dispersion) = context_lvls;
    
    // Convert weight matrix 
    NumericMatrix weight_matrix = to_NumMat(weight_rows);
    rownames(weight_matrix) = treatment_lvls;
    colnames(weight_matrix) = treatment_lvls;
    
    // Make final list to return 
    List results_list = List::create(
      _["model.component.list"] = mc_list,
      _["count.data.summed"] = count_data_summed,
      _["fitted.parameters"] = fitted_parameters,
      _["gamma.dispersion"] = g_dispersion,
      _["param.names"] = param_names,
      _["fix"] = fixed_effects,
      _["treatment"] = treat,
      _["weight.matrix"] = weight_matrix,
      _["grouping.variables"] = grouping_variables,
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
    .method("LRO_initial_param_ests", &wspc::LRO_initial_param_ests)
    .method("LRO_grid_search", &wspc::LRO_grid_search)
    .method("set_parameters", &wspc::set_parameters)
    .method("neg_loglik_debug", &wspc::neg_loglik_debug)
    .method("bounded_nll_debug", &wspc::bounded_nll_debug)
    .method("grad_bounded_nll_debug", &wspc::grad_bounded_nll_debug)
    .method("predict_rates_R", &wspc::predict_rates_R)
    .method("fit", &wspc::fit)
    .method("bs_fit", &wspc::bs_fit)
    .method("bs_batch", &wspc::bs_batch)
    .method("MCMC", &wspc::MCMC)
    .method("resample", &wspc::resample)
    .method("clear_stan_mem", &wspc::clear_stan_mem)
    .method("results", &wspc::results);
  }
