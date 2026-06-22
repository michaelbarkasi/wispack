
// wspc.cpp
#include "wspc.h"

// WSPmm class and methods *********************************************************************************************

/*
 * Object class to hold and fit Warped Sigmoidal Poisson-Process Mixed-Effect Model (WSPmm, aka "WiSP") model. 
 */

// Load settings
void wspc::init_settings(
    const Rcpp::List& settings
  ) {
    cp.buffer_factor = Rcpp::as<double>(settings["buffer_factor"]);
    mf.ctol = Rcpp::as<double>(settings["ctol"]);
    mf.max_penalty_at_distance_factor = Rcpp::as<double>(settings["max_penalty_at_distance_factor"]);
    cp.LRO_cost = Rcpp::as<Rcpp::String>(settings["LRO_cost"]);
    cp.LROcutoff = Rcpp::as<double>(settings["LROcutoff"]);
    cp.LROwindow_factor = Rcpp::as<double>(settings["LROwindow_factor"]);
    cp.rise_threshold_factor = Rcpp::as<double>(settings["rise_threshold_factor"]);
    mf.max_evals = Rcpp::as<int>(settings["max_evals"]);
    mf.rng_seed = Rcpp::as<unsigned int>(settings["rng_seed"]);
    mf.warp_precision = Rcpp::as<double>(settings["warp_precision"]);
    mf.inf_warp = Rcpp::as<double>(settings["inf_warp"]);
    round_none = Rcpp::as<bool>(settings["round_none"]);
    ev.trt_KO = Rcpp::as<CharacterVector>(settings["trtKO"]);
    n_.bin = Rcpp::as<int>(settings["max_bin"]);
    model_settings = Rcpp::clone(settings);
  }

// Get grouping variables
void wspc::init_gv(
    const Rcpp::DataFrame& count_data,
    bool verbose
  ) {
    // Extract grouping variables 
    grouping_lvls.context = Rcpp::sort_unique(Rcpp::as<CharacterVector>(count_data["context"]));
    grouping_lvls.species = Rcpp::sort_unique(Rcpp::as<CharacterVector>(count_data["species"]));
    grouping_lvls.ran = Rcpp::sort_unique(Rcpp::as<CharacterVector>(count_data["ran"]));
    // ... add "none" to represent no random effect (reference level)
    if (grouping_lvls.ran.size() > 1) {grouping_lvls.ran.push_front("none");} 
    // ... print extracted grouping variables
    vprint("Context grouping levels:", verbose);
    vprintV(grouping_lvls.context, verbose);
    vprint("Species grouping levels:", verbose);
    vprintV(grouping_lvls.species, verbose);
    vprint("Random-effect grouping levels:", verbose);
    vprintV(grouping_lvls.ran, verbose);
    // ... grab factor sizes
    n_.ran = grouping_lvls.ran.size();
    n_.context = grouping_lvls.context.size();
    n_.species = grouping_lvls.species.size();
    n_.count_rows = n_.bin * n_.context * n_.species * n_.ran * n_.treatments;
  }

// Extract fixed effects 
void wspc::extract_fixeff(
    const Rcpp::DataFrame& count_data,
    bool verbose
  ) {
    
    // Get column names and remove any KO columns 
    CharacterVector col_names = count_data.names();
    for (String trt_KO : ev.trt_KO) {
      int cnt = 0;
      for (int idx : grep_cpp(col_names, trt_KO)) {col_names.erase(col_names.begin() + (idx - cnt++));}  // ... remove from col_names, adjusting for previous removals
    }
   
    // Check structure of column names
    CharacterVector required_cols = CharacterVector::create("count", "bin", "context", "species", "ran");
    int n_cols = col_names.size();
    int r_cols = required_cols.size();
    for (int i = 0; i < r_cols; i++) {
      if (col_names[i] != required_cols[i]) {Rcpp::stop("Input data is missing required column (or out of order): " + required_cols[i]);}
    }
   
    // Loop through fixed effects
    n_.fix = n_cols - r_cols;               // number of fixed effect variables, assume all non-required columns are fixed effects
    ev.fix_names = CharacterVector(n_.fix);        // names of fixed effect variables 
    ev.fix_ref = CharacterVector(n_.fix);          // reference level for each fixed effect
    ev.fix_lvls.resize(n_.fix);                    // levels for each fixed effect, ev.fix_lvls is a std vector holding an Rcpp CharacterVector
    ev.fix_trt.resize(n_.fix);                     // treatments for each fixed effect, ev.fix_trt is a std vector holding an Rcpp CharacterVector
    CharacterVector is_time_name = {"ref"};
    CharacterVector time_series_names = {"no time series"};
    for (int i = 0; i < n_.fix; i++) {
      ev.fix_names[i] = col_names[i + r_cols];
      
      // Grab column of fixed-effect variable
      SEXP col = count_data[(String)ev.fix_names[i]];
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
        Rcpp::stop("Fixed effect " + ev.fix_names[i] + " has less than 2 levels.");
      } else if (lvls.size() > 2 && ev.fix_names[i] != "timeseries") {
        Rcpp::stop("Fixed effect " + ev.fix_names[i] + " has more than 2 levels. Only binary fixed effects are currently supported, except for 'timeseries'.");
      }
      ev.fix_lvls[i] = lvls;
      ev.fix_ref[i] = lvls[0];                                   // assume first level is reference
      ev.fix_trt[i] = lvls[Rcpp::Range(1, lvls.size() - 1)];     // assume all other levels are treatment levels
      
      // Check for time series
      if (ev.fix_names[i] == "timeseries") {                     // make it easy to access time-series element rank from the element name
        ev.timeseries_rank = seq(1, lvls.size());
        ev.timeseries_rank.names() = lvls;
        time_series_names = lvls;
        for (String l : lvls) {
          ev.is_time.push_back(true);
          is_time_name.push_back(l);
        }
      } else {
        for (String l : lvls) {
          ev.is_time.push_back(false);
          is_time_name.push_back(l);
        }
      }
    }
    ev.is_time.names() = is_time_name;
    
    // Report
    if (ev.fix_names.size() == 0) {
      vprint("No fixed effects detected.", verbose);
    } else {
      vprint("Fixed effects:", verbose);
      vprintV(ev.fix_names, verbose);
      vprint("Ref levels:", verbose);
      vprintV(ev.fix_ref, verbose);
    }
    vprint("No time series detected.", verbose && !any_true(ev.is_time));
    vprint("Time series detected:", verbose && any_true(ev.is_time));
    vprintV(time_series_names, verbose && any_true(ev.is_time));
    
  }

// Make and set treatment levels 
CharacterMatrix wspc::set_treatment_levels(
    bool verbose
  ) {
    
    // Create all possible treatment combinations 
    ev.trt_components = make_treatments(ev.fix_trt);                    // make components of all possible treatment levels
    ev.trt_components.insert(ev.trt_components.begin(), {"ref"}); // add "ref" to represent reference level
    n_.treatments = ev.trt_components.size();                         // grab number of treatments
    ev.trt_lvls = CharacterVector(n_.treatments);                     // resize variable to hold treatment level names
    for (int t = 0; t < n_.treatments; t++) {                            // collapse components into level names
      ev.trt_lvls[t] = Rcpp::collapse(ev.trt_components[t]);
    }
    
    // Check for any treatment-level knockouts and remove
    for (String trt_KO : ev.trt_KO) {
      int cnt = 0;
      for (int idx : grep_cpp(ev.trt_lvls, trt_KO)) {
        ev.trt_lvls.erase(ev.trt_lvls.begin() + (idx - cnt));  // ... remove from ev.trt_lvls
        ev.trt_components.erase(ev.trt_components.begin() + (idx - cnt++));  // ... remove from ev.trt_components, adjusting for previous removals
        n_.treatments--;  // ... decrease
      }
    }
    vprint("Created treatment levels:", verbose);
    vprintV(ev.trt_lvls, verbose);
    
    // Create matrix to translate between treatment levels and fixed-effects columns
    // ... for making summed-count fixed-effect column 
    // ... each cell should contain the value of the treatment column n_.fix, for the given treatment level
    CharacterMatrix effects_rows(n_.treatments, n_.fix);
    for (int tr = 0; tr < n_.treatments; tr++) {
      CharacterVector trt_components = ev.trt_components[tr];
      for (int f = 0; f < n_.fix; f++) {
        CharacterVector lvls = ev.fix_lvls[f];
        effects_rows(tr, f) = lvls[0]; 
        // ^ ... Assume it's the reference level
        for (String l : lvls) {
          // If treatment level found, replace
          if (any_true(eq_left_broadcast(trt_components, l))) {effects_rows(tr, f) = l;}
        }
      }
    }
    
    return effects_rows;
    
  }

// Function to make weight rows matrix
// When predicting parameter values for treatment level tr_predict, should treatment level tr_input's effect be applied?
// ... i.e., what's the "weight" of each treatment level (columns of weight matrix) relative to the others (rows of weight matrix)? 
void wspc::make_weight_rows_matrix() {
    ev.weight_rows.resize(n_.treatments, n_.treatments);
    ev.weight_rows.setOnes(); // ... assume "yes"
    for (int tr_predict = 0; tr_predict < n_.treatments; tr_predict++) {
      // Grab components of treatment level tr_predict
      CharacterVector tr_predict_components = ev.trt_components[tr_predict];
      for (int tr_input = 0; tr_input < n_.treatments; tr_input++) {
        // Grab components of treatment level tr_input
        CharacterVector tr_input_components = ev.trt_components[tr_input];
        // ... if tr_input is the base reference level, tr_input must be applied
        if (tr_input_components.size() == 1 && tr_input_components[0] == "ref") {continue;} 
        // ... else, must check whether all trc_i of tr_input are in tr_predict
        for (int i = 0; i < tr_input_components.size(); i++) {
          String trc_i = tr_input_components[i];
          // if trc_i is not in tr_predict ...
          if(!any_true(eq_left_broadcast(tr_predict_components, trc_i))) {
            // ... check whether trc_i is a time point from "timeseries" ...
            if (ev.is_time[trc_i]) {
              // ... then check if any component trc_p of tr_predict ...
              bool pred_has_time = false;
              for (String trc_p : tr_predict_components) {
                // ... is a time point in "timeseries" ...
                if (ev.is_time[trc_p]) { 
                  pred_has_time = true;
                  // ... that is less than trc_p
                  int i_rank = ev.timeseries_rank[trc_i];
                  int p_rank = ev.timeseries_rank[trc_p];
                  if (i_rank > p_rank) {
                    // ... if so, trc_i is not included in tr_predict
                    ev.weight_rows(tr_predict, tr_input) = 0.0;
                    continue;
                  }
                  // ... if so, tr_input applies and should be left at 1.0
                } 
              }
              if (!pred_has_time) {
                // trc_i is not in tr_predict
                ev.weight_rows(tr_predict, tr_input) = 0.0;
                continue;
              }
            } else {
              // trc_i is not in tr_predict
              ev.weight_rows(tr_predict, tr_input) = 0.0;
              continue;
            }
          } 
        }
      }
    }
  }


// Build summed count data struct from tokenized count data
void wspc::init_summed_count(
    const Rcpp::DataFrame& count_data,
    const CharacterMatrix& effects_rows,
    bool verbose
  ) {
    
    // Allocate columns
    vprint("Total rows in summed count data table: " + std::to_string(n_.count_rows), verbose);
    sc.count.resize(n_.count_rows);
    sc.bin.resize(n_.count_rows);
    sc.ran.resize(n_.count_rows);
    sc.cxt.resize(n_.count_rows);
    sc.sps.resize(n_.count_rows);
    sc.trt.resize(n_.count_rows);
    sc.pr.resize(n_.count_rows);
    sc.pr_log.resize(n_.count_rows); 
    sc.weights.resize(n_.count_rows, n_.treatments);
    
    // Save tokenized count column before collapsing to sums 
    sc.tokens = to_dVec(Rcpp::as<NumericVector>(count_data["count"]));
    
    // Allocate index tracking vectors
    int idx_mcu = 0;
    sc.idx_mc_unique     = IntegerVector(n_.count_rows / n_.bin);
    sc.token_pool        = std::vector<iVec>(n_.count_rows);
    sc.not_na_mask = LogicalVector(n_.count_rows);
    sc.not_na_mask.fill(false);
    vprint("Number of rows with unique model components: " + std::to_string(sc.idx_mc_unique.size()), verbose);
    
    vprint_header("Wrapping up initialization", verbose);
    
    // Extract tokenized-count columns from input data frame
    IntegerVector    binT     = Rcpp::as<IntegerVector>(count_data["bin"]);
    CharacterVector  contextT = Rcpp::as<CharacterVector>(count_data["context"]);
    CharacterVector  speciesT = Rcpp::as<CharacterVector>(count_data["species"]);
    CharacterVector  ranT     = Rcpp::as<CharacterVector>(count_data["ran"]);
    
    // Pre-extract fixed-effect columns (avoids repeated DataFrame lookup in the hot loop)
    std::vector<CharacterVector> fixT(n_.fix);
    for (int f = 0; f < n_.fix; f++) {
      fixT[f] = Rcpp::as<CharacterVector>(count_data[(String)ev.fix_names[f]]);
    }
    
    // Build hash map: (ran, context, species, fix_vals..., bin) -> token row indices
    // Only non-NA rows are included.
    int N_tok = (int)sc.tokens.size();
    std::unordered_map<std::string, iVec> token_map;
    token_map.reserve(N_tok);
    for (int row = 0; row < N_tok; row++) {
      if (std::isnan(sc.tokens[row])) { continue; }
      std::string key;
      key += std::string(ranT[row])     + '\0';
      key += std::string(contextT[row]) + '\0';
      key += std::string(speciesT[row]) + '\0';
      for (int f = 0; f < n_.fix; f++) { key += std::string(fixT[f][row]) + '\0'; }
      key += std::to_string(binT[row]);
      token_map[key].push_back(row);
    }
    
    // Fill summed count columns and weight matrix
    vprint("Creating summed-count data columns and weight matrix", verbose);
    int idx = 0;
    for (int r = 0; r < n_.ran; r++) {
      bool is_none = (grouping_lvls.ran[r] == "none");
      for (int c = 0; c < n_.context; c++) {
        for (int s = 0; s < n_.species; s++) {
          for (int t = 0; t < n_.treatments; t++) {
            
            // Build the (ran, context, species, fix_vals) key prefix once per (r,c,s,t)
            std::string base_key;
            base_key += std::string(grouping_lvls.ran[r])     + '\0';
            base_key += std::string(grouping_lvls.context[c]) + '\0';
            base_key += std::string(grouping_lvls.species[s]) + '\0';
            for (int f = 0; f < n_.fix; f++) {base_key += std::string(effects_rows(t, f)) + '\0';}
            
            sc.idx_mc_unique[idx_mcu++] = idx;
            
            for (int b = 0; b < n_.bin; b++) {
              sc.count[idx]       = 0.0;
              sc.bin[idx]         = static_cast<double>(b + 1.0);
              sc.ran[idx]         = r;
              sc.cxt[idx]         = c;
              sc.sps[idx]         = s;
              sc.trt[idx]         = t;
              sc.weights.row(idx) = ev.weight_rows.row(t);
              
              if (is_none) {
                // "none" ran level: no real tokens — count stays 0, pool stays empty
                sc.not_na_mask(idx) = true;
              } else {
                std::string key = base_key + std::to_string(b + 1);
                auto it = token_map.find(key);
                if (it != token_map.end()) {
                  sc.token_pool[idx] = it->second;
                  for (int rw : it->second) { sc.count[idx] += sc.tokens[rw]; }
                  sc.not_na_mask(idx) = true;
                } else {
                  sc.count[idx] = std::numeric_limits<double>::quiet_NaN();
                }
              }
              idx++;
            }
          }
        }
      }
    }
    
    // Take log of observed counts 
    sc.count_log.resize(n_.count_rows); 
    for (int r = 0; r < n_.count_rows; r++) {sc.count_log[r] = std::log(sc.count[r] + 1.0);} 
    
    // Derive index vectors from masks
    sc.not_na_idx     = Rwhich(sc.not_na_mask);
    sc.not_none_mask  = !eq_left_broadcast(sc.ran, 0);
    sc.r_idx          = Rwhich(!sc.not_none_mask);
    
  }

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
    
    // Extract settings and grouping variables 
    init_settings(settings);
    init_gv(count_data, verbose); 
    
    // Report warp_inf 
    const double eps_ = std::numeric_limits<double>::epsilon(); // machine epsilon
    if (verbose) {
      vprint_header("Infinity handling");
      vprint("machine epsilon: (eps_): ", eps_);
      vprint("pseudo-infinity (inf_): ", inf_);
      vprint("warp_precision: ", mf.warp_precision);
      vprint("implied pseudo-infinity for unbounded warp (inf_warp): ", mf.inf_warp);
    }
    
    // Find max bins and set warp bounds
    vprint_header("Extracting variables and data information", verbose);
    if (n_.bin == 0) {n_.bin = (int)Rcpp::max(Rcpp::as<NumericVector>(count_data["bin"]));}
    vprint("Max bin: " + std::to_string(n_.bin), verbose);
    mf.warp_bounds.resize(3);
    mf.warp_bounds << mf.inf_warp, mf.inf_warp, (double)n_.bin;
    
    // Compute tpoint buffer
    cp.tpoint_buffer = (double)n_.bin * cp.buffer_factor;
    int tpoint_buffer_int = (int)cp.tpoint_buffer;
    if (cp.tpoint_buffer != tpoint_buffer_int) {tpoint_buffer_int++;} 
    cp.tpoint_buffer = (double)tpoint_buffer_int;
    
    // Extract fixed effects 
    extract_fixeff(count_data, verbose);
    
    // Find and set treatment levels 
    CharacterMatrix effects_rows = set_treatment_levels(verbose);
    
    // Recalculate n_.count_rows now that n_.treatments is known.
    // (n_.count_rows was first computed in init_gv() before set_treatment_levels()
    // ran, so n_.treatments was still 0 at that point, giving n_.count_rows = 0.)
    n_.count_rows = n_.bin * n_.context * n_.species * n_.ran * n_.treatments;
    
    // Make weight rows matrix
    make_weight_rows_matrix();
    
    // Build summed count data
    init_summed_count(count_data, effects_rows, verbose);
    
    // Make extrapolation pool
    make_extrapolation_pool(); 
    
    // Extrapolate "none" rows
    extrapolate_none();
    
    // Estimate gamma dispersion of raw counts
    compute_gamma_dispersion();
    
    // Find average log counts for each context-species combination
    find_count_log_means();
    
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
    // Recover memory from stan
    stan::math::recover_memory();
  }

/*
 * *************************************************************************
 * Initial setup
 */

// Get row indices of count data from the given factor levels
iVec wspc::fetch_count_idx(
    iVec I // indexes of each factor
  ) {
    const int stack_depth = 5;
    if (stack_depth != I.size()) {Rcpp::stop("Invalid stack size.");}
    const int n[5] = {n_.ran, n_.context, n_.species, n_.treatments, n_.bin};
    iVec idx = {0};
    for (int i = 0; i < stack_depth; ++i) {
      int remainder = 1;
      for (int l = stack_depth - 1; l > i; --l) {remainder *= n[l];}
      if (I[i] < 0) {
        // index all levels of this factor 
        iVec idxj;
        for (int k = 0; k < idx.size(); ++k) {
          for (int j = 0; j < n[i]; ++j) {
            idxj.push_back(idx[k] + j * remainder);
          }
        }
        idx = idxj;
      } else {
        // index specific level of this factor
        for (int k = 0; k < idx.size(); ++k) {idx[k] += remainder * I[i];}
      }
    }
    return idx;
  }

// Computes gamma_dispersion matrix and related vectors
void wspc::compute_gamma_dispersion() {
    // Initialize gamma dispersion matrix
    mf.gamma_dispersion = std::vector<dVec>(n_.species, dVec(n_.context));
    // Loop through context and species levels
    for (int c = 0; c < n_.context; c++) {
      for (int s = 0; s < n_.species; s++) { 
        // Estimate dispersion of raw count (not log)
        dVec count_cs_masked; 
        for (int r : fetch_count_idx({-1, c, s, -1, -1})) {
          if (sc.not_na_mask[r]) {count_cs_masked.push_back(sc.count[r]);}
        }
        mf.gamma_dispersion[s][c] = gamma_dispersion_formula(vmean(count_cs_masked), vvar(count_cs_masked));
      }
    }
  }

// Function to pre-make extrapolation pools for ran-level "none"
void wspc::make_extrapolation_pool() {
    sc.extrap_pool = std::vector<iVec>(n_.count_rows);
    // Loop through none rows and find their interpolation pools
    int n_rows = sc.r_idx.size();
    for (int ri = 0; ri < n_rows; ri++) {
      int r = sc.r_idx[ri];
      iVec idx = fetch_count_idx({-1, sc.cxt[r], sc.sps[r], sc.trt[r], (int)sc.bin[r] - 1});
      for (int rii : idx) {
        if (sc.not_na_mask[rii] && sc.not_none_mask[rii]) {sc.extrap_pool[r].push_back(rii);}
      }
    }
  }

// Function to extrapolate "none" rows
void wspc::extrapolate_none() {
    for (int r : sc.r_idx) {
      int extrap_sz = sc.extrap_pool[r].size();
      if (extrap_sz > 0) {
        sc.count[r] = 0.0; // Find mean of the log values
        for (int i : sc.extrap_pool[r]) {sc.count[r] += std::log(sc.count[i] + 1.0);}
        sc.count[r] /= (double)extrap_sz;
        sc.count[r] = std::exp(sc.count[r]) - 1.0; // Convert back from log space
        if (round_none) {sc.count[r] = std::round(sc.count[r]);}
      }
    }
    for (int r : sc.r_idx) {sc.count_log[r] = std::log(sc.count[r] + 1.0);}
  }

// Function to take means of count_log 
void wspc::find_count_log_means() {
    // Initialize lists to hold sc.count_log_avg for each context-species combination 
    sc.count_log_avg = std::vector<std::vector<NumericMatrix>>(n_.context, std::vector<NumericMatrix>(n_.species));
    // Loop through context and species levels
    for (int c = 0; c < n_.context; c++) {
      for (int s = 0; s < n_.species; s++) { 
        // Find mean of count_log for each species-context pair
        sc.count_log_avg[c][s] = NumericMatrix(n_.bin, n_.treatments);
        for (int t = 0; t < n_.treatments; t++) {
          dVec count_log_avg(n_.bin); 
          for (int b = 0; b < n_.bin; b++) {
            dVec count_b; 
            for (int i : fetch_count_idx({-1, c, s, t, b - 1})) {
              if (sc.not_na_mask[i]) {count_b.push_back(sc.count_log[i]);}
            }
            count_log_avg[b] = vmean(count_b);
          }
          sc.count_log_avg[c][s].column(t) = to_NumVec(count_log_avg);
        }
      }
    }
  }

// Function to estimate change points
void wspc::estimate_change_points() {
   
    // Grab number of context and species levels
    const int n_ran_trt = n_.ran * n_.treatments;
    
    // Initialize matrix to hold degrees of each context-species combination
    cp.deg = IntegerMatrix(n_.species, n_.context);
    rownames(cp.deg) = grouping_lvls.species;
    colnames(cp.deg) = grouping_lvls.context;
    
    // Reserve space to hold results matrices for each context-species combination 
    cp.found = std::vector<std::vector<IntegerMatrix>>(n_.context, std::vector<IntegerMatrix>(n_.species));
    cp.found_trt = std::vector<std::vector<NumericMatrix>>(n_.context, std::vector<NumericMatrix>(n_.species));
    
    // Loop through context and species levels
    for (int c = 0; c < n_.context; c++) {
      for (int s = 0; s < n_.species; s++) {
        eMat count_masked_array(n_.bin, n_ran_trt);
        count_masked_array.setZero();
        LogicalVector good_col(n_ran_trt);
        for (int t = 0; t < n_.treatments; t++) {
          for (int r = 0; r < n_.ran; r++) {
            dVec count_log_masked; 
            for (int i : fetch_count_idx({r, c, s, t, -1})) {
              if (sc.not_na_mask[i]) {count_log_masked.push_back(sc.count_log[i]);}
            }
            // Check and save column
            if (count_log_masked.size() == n_.bin) {
              count_masked_array.col(t*n_.ran + r) = to_eVec(count_log_masked);
              good_col(t*n_.ran + r) = true;
            } else {
              good_col(t*n_.ran + r) = false;
            }
          }
          
        }
        
        // Extract good column numbers 
        IntegerVector good_col_idx = Rwhich(good_col);
        eMat count_masked_array_good = count_masked_array(Eigen::all, to_iVec(good_col_idx));
        
        // Estimate change points from masked count series
        IntegerMatrix found_cp_good = LROcp_array(
          count_masked_array_good,       // 2D matrix of points to test for change points (columns as series, rows as bins)
          cp.ws,                         // running window size 
          cp.LROcutoff,                  // points more than this times sd considered outliers
          cp.tpoint_buffer               // Minimum distance between two change points
        );
        
        // Estimate degree of this context-species pair 
        cp.deg(s, c) = found_cp_good.rows();
        
        // Fill columns into the cp.found matrix
        cp.found[c][s] = IntegerMatrix(cp.deg(s, c), n_ran_trt);
        // ^ ... Rcpp should initialize these matrices with all zeros
        if (cp.deg(s, c) > 0) {
          for (int si = 0; si < good_col_idx.size(); si++) {
            // ... grab change points
            cp.found[c][s].column(good_col_idx[si]) = found_cp_good.column(si);
          }
        }
        
        // Extract treatment means for each change point
        cp.found_trt[c][s] = NumericMatrix(cp.deg(s, c), n_.treatments);
        for (int t = 0; t < n_.treatments; t++) {
          for (int d = 0; d < cp.deg(s, c); d++) {
            double temp = 0.0;
            int n_ran_hit = 0;
            for (int r = 0; r < n_.ran; r++) {
              if (good_col(t*n_.ran + r)) { // ... ensure there is data here
                temp += (double)cp.found[c][s](d, t*n_.ran + r);
                n_ran_hit++;
              }
            }
            cp.found_trt[c][s](d, t) = temp / (double)n_ran_hit;
          }
        }
        
      }
    }
    
    // Compute implied size of the parameter boundary vector 
    mf.boundary_vec_size = 0;
    for (int r : sc.idx_mc_unique) {
      // Grab degree for this row
      if (cp.deg(sc.sps[r], sc.cxt[r]) > 0){
        // Add slots for the boundary distance at each tpoint
        mf.boundary_vec_size += cp.deg(sc.sps[r], sc.cxt[r]) + 1;
        // Add slots for the tslope boundary distance at each tslope
        mf.boundary_vec_size += cp.deg(sc.sps[r], sc.cxt[r]);
        // Add slots for the boundary distance at each block rate
        mf.boundary_vec_size += cp.deg(sc.sps[r], sc.cxt[r]) + 1;
      } else {
        // Add one slot for the single block rate value
        mf.boundary_vec_size++;
      } 
    }
    
  }

// Estimate change points and initial parameters for model fitting
void wspc::LRO_initial_param_ests(
    bool verbose,
    double LROwf,
    double LROco
  ) {
   
    if (LROwf != 0.0) {cp.LROwindow_factor = LROwf;}
    if (LROco != 0.0) {cp.LROcutoff = LROco;} // ... LROcutoff used in estimate_change_points, so set before calling
    
    // Compute running and filter window sizes for LRO change-point detection
    cp.ws = static_cast<int>(std::round(cp.LROwindow_factor * (double)n_.bin * cp.buffer_factor));
    
    // Estimate degree of each context-species combination at baseline using LRO change-point detection 
    estimate_change_points();
    vprint("Estimated change points", verbose);
    
    // Find initial parameter estimates
    estimate_initial_parameters();
    vprint("Estimated initial parameters", verbose);
    vprint("Number of parameters: " + std::to_string((int)mf.params.size()), verbose);
    vprint("Size of boundary vector: " + std::to_string(mf.boundary_vec_size), verbose);
    
  }

// Search for best LRO change-point detection settings
void wspc::LRO_grid_search(bool verbose) {
    
    // Check if a grid search can be done
    if (!(cp.LRO_cost == "AIC" || cp.LRO_cost == "BIC" || cp.LRO_cost == "NLL")) {return;}
    
    // Set up search grid
    vprint("Performing grid search (resolution 0.25) for optimal LRO window factor and cutoff", verbose);
    dVec LROwf_range = {1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 3.00};
    dVec LROco_range = LROwf_range;
    cp.LRO_grid_search_results = NumericMatrix(LROwf_range.size() * LROco_range.size(), 8);
    int n = sc.not_na_idx.size();
    
    // Run LRO estimation for each initial grid point
    for (int i = 0; i < LROwf_range.size(); ++i) {
      vprint("LRO window factor: " + std::to_string(LROwf_range[i]), verbose);
      for (int j = 0; j < LROco_range.size(); ++j) {
        LRO_initial_param_ests(false, LROwf_range[i], LROco_range[j]);
        fit(false);
        int k = mf.params.size();
        int idx = i * LROco_range.size() + j;
        cp.LRO_grid_search_results(idx, 0) = LROwf_range[i];
        cp.LRO_grid_search_results(idx, 1) = LROco_range[j];
        cp.LRO_grid_search_results(idx, 2) = mf.nll;
        cp.LRO_grid_search_results(idx, 3) = mf.bnll;
        cp.LRO_grid_search_results(idx, 4) = (double)mf.success_code;
        cp.LRO_grid_search_results(idx, 5) = (double)k;
        cp.LRO_grid_search_results(idx, 6) = 2.0 * ((double)k + mf.nll);
        cp.LRO_grid_search_results(idx, 7) = 2.0 * (std::log((double)n) * (double)k + mf.nll);
        clear_stan_mem();
      }
    }
    
    // Label results matrix
    colnames(cp.LRO_grid_search_results) = CharacterVector::create(
      "LROwindow_factor", 
      "LROcutoff", 
      "neg_loglik", 
      "penalized_neg_loglik",
      "success_code", 
      "n_params",
      "AIC", 
      "BIC"
      );
    
    // Find lowest cost window_factor and cutoff
    double cost = inf_; 
    for (int i = 0; i < LROwf_range.size() * LROco_range.size(); i++) {
      double this_cost = inf_; 
      if (cp.LRO_cost == "AIC") {
        this_cost = cp.LRO_grid_search_results(i,6);
      } else if (cp.LRO_cost == "BIC") {
        this_cost = cp.LRO_grid_search_results(i,7);
      } else if (cp.LRO_cost == "NLL") {
        this_cost = cp.LRO_grid_search_results(i,2);
      }
      if (this_cost < cost) {
        cost = this_cost;
        cp.LROwindow_factor = cp.LRO_grid_search_results(i,0);
        cp.LROcutoff = cp.LRO_grid_search_results(i,1);
      } 
    }
    
    vprint("Optimal LRO window factor: " + std::to_string(cp.LROwindow_factor) + ", optimal LRO cutoff: " + std::to_string(cp.LROcutoff), verbose);
    
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
    int bktp = cp.deg(sc.sps[r], sc.cxt[r]);
    if (mc == 0) {bktp++;} // ... 0 is Rt
    // Extract weight matrix row
    eVec weight_row = sc.weights.row(r).transpose();
    // Compute unwarped model component value  
    sVec mc_vec(bktp);
    mc_vec.setZero();
    for (int bt = 0; bt < bktp; bt++) {
      for (int t = 0; t < n_.treatments; t++) {
        mc_vec(bt) += weight_row(t) * parameters(beta_idx[mc][sc.cxt[r]][sc.sps[r]][bt*n_.treatments + t]);
      }
      // ... and apply warping factor
      mc_vec(bt) = warp_mc(mc_vec(bt), mf.warp_bounds(mc), wf);
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
    sVec out(n_.count_rows);
    iVec mc_unique_rows = Rcpp::as<iVec>(sc.idx_mc_unique);
    
    // Initialize variables to hold model components
    sVec Rt; sVec tslope; sVec tpoint;
    
    // Grab warping indices and initiate variables to hold them
    sdouble f_pw; sdouble f_rw; sdouble f_sw;
    
    // Compute predicted rate for rows of the summed count data
    for (int r = 0; r < n_.count_rows; r++) {
      
      // Update predicted model components if r begins a new batch of unique values 
      int cnt = std::count(mc_unique_rows.begin(), mc_unique_rows.end(), r);
      
      // Unless all_rows, only compute values for non-NA rows or rows that begin a new batch
      if (sc.not_na_mask[r] || cnt > 0 || all_rows) {
        
        // Grab warping factors
        // WARNING: this code is duplicated in boundary_dist
        int wfi = sc.sps[r] * (grouping_lvls.ran.size() - 1) + sc.ran[r] - 1;
        if (sc.ran[r] == 0) {
          f_pw = 0.0; f_rw = 0.0; f_sw = 0.0;
        } else {
          f_pw = parameters(wfactor_idx[0][wfi]);
          f_rw = parameters(wfactor_idx[1][wfi]); 
          f_sw = parameters(wfactor_idx[2][wfi]);
        } 
        
        // Only update predicted model components if r begins a new batch of unique values 
        if (cnt > 0) { 
          // Compute warped model components for this row r
          Rt = compute_warped_mc(0, r, parameters, f_rw);        
          tslope = compute_warped_mc(1, r, parameters, f_sw); 
          tpoint = compute_warped_mc(2, r, parameters, f_pw);
        } 
        
        // Compute the poly-sigmoid
        out(r) = poly_sigmoid(
          sc.bin[r],
          cp.deg(sc.sps[r], sc.cxt[r]),
          Rt,
          tslope,
          tpoint
        ); 
        
      }
      
    }
    
    return out;
    
  }

// Predict log of rates, R wrapper 
NumericVector wspc::predict_rates_R(
    const NumericVector& parameters_R,
    const bool& all_rows 
  ) {
    // Convert parameters to sVec
    sVec parameters = to_sVec(parameters_R);
    // Compute predicted rates
    sVec out = predict_rates(
      parameters,
      all_rows
    );
    // Convert to NumericVector 
    NumericVector out_R = to_NumVec(out);
    // Clear memory and return
    clear_stan_mem();
    return out_R;
  }

/*
 * *************************************************************************
 * Computing objective function (i.e., fitting model and parameter boundary distances)
 */

// Compute nll of observations under the given parameters
sdouble wspc::compute_nll(
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
    for (int r : sc.not_na_idx) {
      
      if (std::isinf(predicted_rates_log_var(r)) || predicted_rates_log_var(r) < rt_lower_bound || std::isnan(predicted_rates_log_var(r))) {
        return sdouble(inf_);
      } else {
        
        // Find gamma variance for this row
        // ... pull predicted rate from log space
        sdouble pred_rate_var = sexp(predicted_rates_log_var(r)) - 1.0;
        // ... estimate variance of rate outside log space
        sdouble gamma_variance = pred_rate_var + mf.gamma_dispersion[sc.sps[r]][sc.cxt[r]] * pred_rate_var * pred_rate_var;
        // ... estimate the corresponding variance of the rate back in log space 
        gamma_variance = delta_var_est(gamma_variance, pred_rate_var);
        
        // Pull up if only approximately zero 
        if (predicted_rates_log_var(r) < 0.0) {predicted_rates_log_var(r) = 0.0;}
        
        // Analytic solution to the log of the integral from 0 to positive infinity of the Poisson times Gamma densities
        if (gamma_variance == 0) { 
          // if no over-dispersion, just use Poisson (via the existing log_dPois helper)
          log_lik += log_dPois(sc.count_log[r], predicted_rates_log_var(r));
        } else if (sc.count_log[r] > 0.0 && predicted_rates_log_var(r) > 0.0) { 
          // otherwise, use Poisson-Gamma integral
          log_lik += slog(poisson_gamma_integral(sc.count_log[r], predicted_rates_log_var(r), gamma_variance));
        }
        
      }
      
    }
    
    // Take negative
    sdouble negloglik = -log_lik;
    
    // Check for infinities (zero likelihood)
    if (std::isinf(negloglik) || negloglik > sdouble(inf_)) {negloglik = sdouble(inf_);}
    
    return negloglik;
    
  } 

// Compute boundary distances
sVec wspc::boundary_dist(
    const sVec& parameters
  ) const {
    
    // Initialize vector to hold boundary distances
    sVec boundary_dist_vec(mf.boundary_vec_size);
    int ctr = 0;
    
    // Grab warping indices and initiate variables to hold them
    sdouble f_pw; sdouble f_rw; sdouble f_sw;
    
    // Compute the boundary distance, for ...
    // ... transition slopes (must be positive)
    // ... transition points (enforces tpoint buffer)
    // ... block rate values (must be positive)
    for (int r : sc.idx_mc_unique) {
      
      // Grab warping factors
      // WARNING: this code is duplicated from predict_rates
      int wfi = sc.sps[r] * (grouping_lvls.ran.size() - 1) + sc.ran[r] - 1;
      if (sc.ran[r] == 0) {
        f_pw = 0.0; f_rw = 0.0; f_sw = 0.0;
      } else {
        f_pw = parameters(wfactor_idx[0][wfi]);
        f_rw = parameters(wfactor_idx[1][wfi]); 
        f_sw = parameters(wfactor_idx[2][wfi]);
      } 
      
      // Compute Rt for this row r
      // ... note: Not computing rate for all rows, just ones with unique model components
      // ... this reduces number of rows to compute by a factor of n_.bin, which is significant
      sVec Rt = compute_warped_mc(0, r, parameters, f_rw);
      
      // Grab degree for this row
      int deg = Rt.size() - 1;
      
      // Compute t-point and R_sum boundary distances
      if (deg > 0){
        
        // Compute tslope and tpoint for this row r
        sVec tslope = compute_warped_mc(1, r, parameters, f_sw); 
        sVec tpoint = compute_warped_mc(2, r, parameters, f_pw);
        
        // Transition slopes most be positive
        for (int d = 0; d < deg; d++) {boundary_dist_vec(ctr++) = tslope(d);}
        
        // Find tpoint boundary distances
        sVec tpoint_ext(deg + 2);
        for (int bt = 0; bt < tpoint_ext.size(); bt++) {
          if (bt == 0) {tpoint_ext(bt) = 0.0;} 
          else if (bt <= deg) {tpoint_ext(bt) = tpoint(bt - 1);} 
          else {tpoint_ext(bt) = (sdouble)n_.bin;} 
        }
        
        // Transition points must be in increasing order 
        // ... and can't be closer than the buffer
        // ... and first point must be > buffer, 
        // ... and last point must be < n_.bin - buffer
        for (int d = 0; d < deg + 1; d++) {boundary_dist_vec(ctr++) = (tpoint_ext(d + 1) - tpoint_ext(d)) - cp.tpoint_buffer;}
        
        // All block rate values must be positive
        for (int d = 0; d < deg + 1; d++) {boundary_dist_vec(ctr++) = Rt(d) - rt_lower_bound;}
        
      } else {
        // ... trivial to check if rate (Rt) is positive
        boundary_dist_vec(ctr++) = Rt(0) - rt_lower_bound;
      }
      
    } 
    
    // Check for nan
    for (int i = 0; i < mf.boundary_vec_size; i++) {
      if (std::isnan(boundary_dist_vec(i))) {
        boundary_dist_vec(i) = 0.0;
      }
    }
    
    return boundary_dist_vec;
    
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
    if (!std::isfinite(fx.val())) {fx = 0.0;}
    
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
    
    // Return the value
    return fx.val(); 
    
  } 

// Compute nll plus boundary penalty (main objective function) 
sdouble wspc::compute_bnll(
    const sVec& parameters
  ) const {
    
    // Compute weighted negative log-likelihood
    sdouble bnll = compute_nll(parameters);
    // Compute boundary distance
    sVec bd = boundary_dist(parameters);
    // Add boundary penalty
    for (int i = 0; i < mf.boundary_vec_size; i++) {bnll += boundary_penalty_transform(bd(i), mf.bp_coefs(i));}
    
    /*
     * Idea of boundary penalty transform: When "far" from boundary, the total penalty will be at most 
     *  some specified fraction (e.g., 0.1) of the magnitude of the nll, but if any one component of 
     *  the boundary distance approaches zero, the penalty smoothly goes to infinity. 
     */
    
    return bnll;
    
  }

// Wrap compute_bnll in form needed for NLopt objective function
double wspc::compute_bnll_NLopt(
    const dVec& x, 
    dVec& grad, 
    void* data
  ) {
    
    // Grab model
    wspc* model = static_cast<wspc*>(data);
    // Convert dVec to Eigen with stan
    sVec parameters_var = to_sVec(x);
    // Compute bounded_nll
    sdouble fx = model->compute_bnll(parameters_var);
    
    // Compute gradient if needed
    if (!grad.empty()) {
      Eigen::VectorXd grad_eigen = model->grad_compute_bnll(parameters_var);
      grad.assign(grad_eigen.data(), grad_eigen.data() + grad_eigen.size());
    } 
    
    // Return the value of the nll
    return fx.val(); 
    
  }

/*
 * *************************************************************************
 * Computing gradients with stan reverse-mode autodiff
 */

// Compute the gradient of the compute_bnll function
// ... this is the gradient function used in model optimization
Eigen::VectorXd wspc::grad_compute_bnll(
    const sVec& p_
  ) const { 
    // Create nested autodiff context
    stan::math::nested_rev_autodiff nested;
    // Make copy to create var nodes for p
    sVec p = p_;
    // Initialize bnll variable
    sdouble bnll = compute_bnll(p);
    // Initialize variable to hold gradient
    Eigen::VectorXd gr_bnll(p.size());
    // Compute bnll and its gradient
    stan::math::grad(bnll, p, gr_bnll);
    // Return bnll gradient
    return gr_bnll;
  }

// Compute the gradient of the min_boundary_dist function
// ... this is the gradient function used in the search for feasible parameters
Eigen::VectorXd wspc::grad_min_boundary_dist(
    const sVec& p_
  ) const { 
    // Create nested autodiff context
    stan::math::nested_rev_autodiff nested;
    // Make copy to create var nodes for p
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
void wspc::fit(const bool verbose) { 
    
    // Set boundary-penalty coefficients 
    vprint("Setting boundary penalty coefficients", verbose);
    sVec initial_params_var = to_sVec(mf.params);
    double max_penalty_at_distance = compute_nll(initial_params_var).val() * mf.max_penalty_at_distance_factor;
    double coefs = std::sqrt(static_cast<double>(mf.boundary_vec_size)/max_penalty_at_distance);
    mf.bp_coefs = to_eVec(boundary_dist(initial_params_var));
    for (int i = 0; i < mf.boundary_vec_size - 1; i++) {mf.bp_coefs(i) = coefs/mf.bp_coefs(i);} 
    
    // Grab initial parameters
    dVec x = to_dVec(mf.params); 
    size_t n = x.size();
    
    // Set up NLopt optimizer
    nlopt::opt opt(nlopt::LD_LBFGS, n); // Might try? LD_LBFGS, LN_SBPLX, LN_COBYLA, GN_DIRECT
    opt.set_min_objective(wspc::compute_bnll_NLopt, this);
    opt.set_ftol_rel(mf.ctol);             // Stop when iteration changes objective fn value by less than this fraction 
    opt.set_maxeval(mf.max_evals);         // Maximum number of evaluations to try
    
    // Find and print initial nll and total objective values
    if (verbose) {
      sdouble initial_nll = compute_nll(initial_params_var);
      vprint("Initial nll: ", initial_nll);
      sdouble initial_obj = compute_bnll(initial_params_var);
      vprint("Initial nll with penalty: ", initial_obj);
      sdouble mb_dist = min_boundary_dist(initial_params_var);
      vprint("Initial min boundary distance: ", mb_dist);
    } 
    
    // Fit model
    int sc = 0;
    double min_fx;
    try {
      nlopt::result sc_ = opt.optimize(x, min_fx);
      sc = static_cast<int>(sc_);
    } catch (std::exception& e) { 
      if (verbose) {
        Rcpp::Rcout << "Optimization failed: " << e.what() << std::endl;
      } 
      sc = 0;
    } 
    
    // Find final nll for diagnostics
    sVec parameters_final = to_sVec(x);
    sdouble final_nll = compute_nll(parameters_final);
    
    // Print final nll, total objective, and min boundary distance values
    if (verbose) {
      vprint("Final nll: ", final_nll);
      vprint("Final nll with penalty: ", min_fx);
      sdouble mb_dist = min_boundary_dist(parameters_final);
      vprint("Final min boundary distance: ", mb_dist);
      int num_evals = opt.get_numevals();
      vprint("Number of evaluations: ", num_evals);
      if (sc == 0) {
        vprint("Warning: optimization did not converge.");
      }
    } 
    
    // Save optimization results
    mf.params = to_NumVec(x);
    mf.bnll = min_fx;
    mf.nll = final_nll.val();
    mf.success_code = sc;
    mf.num_evals = opt.get_numevals();
    
    // Check feasibility 
    check_parameter_feasibility(to_sVec(mf.params)); 
    
  } 

// Resample count (with replacement) for bootstrapping
void wspc::resample(
    unsigned int seed
  ) {
   
    // Initialize random number generator
    std::mt19937 rng(seed);
    
    // Resample and re-sum token pool
    for (int r : sc.not_na_idx) {
      // ... but only for actual observations of random grouping variable
      if (sc.ran[r] > 0) {
        // ... redraw randomly (with replacement) from the token pool of r and re-sum into r's count 
        int resample_sz = sc.token_pool[r].size();
        if (resample_sz < 1) {
          // ... ensure new point is viable
          Rcpp::Rcout << "Error: resample size < 1, for row " << r << std::endl;
          Rcpp::stop("Error: resample size < 1");
        }
        // ... re-sum into r's count
        sc.count[r] = 0.0;
        for (int i = 0; i < resample_sz; i++) {
          // ... randomly select integer between 0 and resample_sz
          int resample_idx = rUnifi(0, resample_sz - 1, rng);
          sc.count[r] += sc.tokens[sc.token_pool[r][resample_idx]];
        }
        sc.count_log[r] = std::log(sc.count[r] + 1.0);
      }
    }
    
    // Extrapolate none's and take their logs
    extrapolate_none();
    
    // Estimate gamma dispersion of these new raw re-sampled counts
    compute_gamma_dispersion();
    
    // Find average these new re-sampled log counts for each context-species combination
    find_count_log_means();
    
  }

// Fit model to bootstrap resample
dVec wspc::bs_fit(
    int bs_num,               // A unique number for this resample
    bool clear_stan           // Recover stan memory at end?
  ) {
    
    // Resample counts 
    resample(mf.rng_seed + bs_num);
    // Find initial parameter estimates, for new re-sampled data
    estimate_initial_parameters();
    // Fit model
    fit(false);
    
    // Prepare and return results
    dVec res = to_dVec(mf.params);
    res.reserve(mf.params.size() + 4);
    res.push_back(mf.bnll);
    res.push_back(mf.nll);
    res.push_back((double)mf.success_code);
    res.push_back((double)mf.num_evals);
    
    // Clear stan memory
    if (clear_stan) {clear_stan_mem();}
    
    return res;
    
  }

// Fork bootstraps (parallel processing)
Rcpp::NumericMatrix wspc::bs_batch(
    int bs_num_max,           // Number of bootstraps to perform
    int max_fork,             // Maximum number of forked processes per batch
    bool use_median,
    bool verbose
  ) {
    
    // Fit full model 
    vprint("Performing initial fit of full data", verbose);
    fit(false);
    if (verbose) {vprint("Penalized nll: ", mf.bnll);}
    
    // Initiate variables to hold results
    const int c_num = mf.params.size() + 4;
    const int r_num = bs_num_max + 1;
    NumericMatrix bs_results(r_num, c_num);
    
    // Save results from initial full fit
    dVec res = to_dVec(mf.params);
    res.reserve(res.size() + 4);
    res.push_back(mf.bnll);
    res.push_back(mf.nll);
    res.push_back((double)mf.success_code);
    res.push_back((double)mf.num_evals);
    
    // Save full fit results in last row of results matrix
    bs_results.row(bs_num_max) = to_NumVec(res);
    
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
        // Run in serial ... Fit bootstrap and copy into the corresponding row of the matrix
        bs_results.row(b) = to_NumVec(bs_fit(b, true));
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
    
    // Use median instead of initial full-fit
    if (use_median) {
      // Compute column-wise median of all rows except the last
      int r_num1 = r_num - 1;
      for (int j = 0; j < c_num; j++) {
        NumericVector col_vals(r_num1); // r_num elements (rows 0 to r_num-1)
        for (int i = 0; i < r_num1; i++) {col_vals[i] = bs_results(i, j);}
        std::sort(col_vals.begin(), col_vals.end());
        double median;
        if (r_num1 % 2 == 0) {median = (col_vals[r_num1 / 2 - 1] + col_vals[r_num1 / 2]) / 2.0;} 
        else {median = col_vals[r_num1 / 2];}
        bs_results(r_num1, j) = median;
        if (j < mf.params.size()) {mf.params[j] = median;}
      }
    }
    
    // Clear stan memory and return
    vprint("All complete!", verbose); 
    if (bs_num_max == 0) {clear_stan_mem();}
    return bs_results;
    
  }

// Markov-chain Monte Carlo (MCMC) simulation
Rcpp::NumericMatrix wspc::MCMC(
    int n_burnin,             // Number of initial steps to take (to find optimal parameters)
    int n_steps,              // Number of steps to take in random walk (post-burnin)
    int neighbor_filter,      // Keep only ever neighbor_filter step
    double step_size,         // Step size for random walk
    double prior_sd,          // standard deviation to use in prior
    bool start_from_fit,      // Start from parameters found with gradient descent? 
    bool use_median,
    bool verbose
  ) {
    
    // Initiate variables to hold results
    const int n_params = mf.params.size();
    const int c_num = n_params + 4;
    const int r_num = n_steps + 1;
    NumericMatrix RMW_steps(r_num, c_num);
    n_steps += n_burnin;
    
    // Fit full model 
    if (start_from_fit) {
      vprint("Starting MCMC walk from L-BFGS optimized parameters", verbose);
      fit(false);
      if (verbose) {vprint("Penalized nll: ", mf.bnll);}
      // ... save results from initial full fit
      dVec res = to_dVec(mf.params);
      res.reserve(res.size() + 4);
      res.push_back(mf.bnll); 
      res.push_back(mf.nll);
      res.push_back(double(1.0)); // acceptance ratio
      res.push_back(double(0.0)); // ctr number
      // ... save full fit results in last row of results matrix
      RMW_steps.row(r_num - 1) = to_NumVec(res);
    } else {
      // Set boundary-penalty coefficients 
      vprint("Setting boundary penalty coefficients", verbose);
      sVec initial_params_var = to_sVec(mf.params);
      double max_penalty_at_distance = compute_nll(initial_params_var).val() * mf.max_penalty_at_distance_factor;
      double coefs = std::sqrt(static_cast<double>(mf.boundary_vec_size)/max_penalty_at_distance);
      mf.bp_coefs = to_eVec(boundary_dist(initial_params_var));
      for (int i = 0; i < mf.boundary_vec_size - 1; i++) {mf.bp_coefs(i) = coefs/mf.bp_coefs(i);} 
    }
    
    // Make baseline parameter mask 
    LogicalVector baseline_mask(n_params);
    baseline_mask.fill(false);
    for (int i : param_idx.baseline) {baseline_mask(i) = true;}
    
    // Make beta_Rt parameter mask
    LogicalVector beta_Rt_mask(n_params);
    beta_Rt_mask.fill(false);
    for (int i : param_idx.beta_Rt) {beta_Rt_mask(i) = true;} 
    
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
    NumericVector params_initial = mf.params;
    NumericVector params_current = mf.params;
    NumericVector last_viable_parameters = params_current;
    
    // Initialize neighbor counter
    int neighbor_counter = 0;
    
    // Set random number generator with seed
    unsigned int seed = mf.rng_seed + n_steps;
    std::mt19937 walk_rng(seed);
    
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
        for (int i = 0; i < mf.boundary_vec_size; i++) {bd_current_transformed += boundary_penalty_transform(bd_current_vec(i), mf.bp_coefs(i));}
        bd_current_transformed += 1.0;
        
        // ... for each parameter
        for (int i = 0; i < n_params; i++) {
          
          // ... calculate step size
          int param_oom = static_cast<int>(std::floor(std::log10(std::abs(params_current(i)))));
          double normalized_step_size = step_size * std::pow(10, static_cast<double>(param_oom + 1.0));
          double bounded_step_size = normalized_step_size / bd_current_transformed.val();
          if (bounded_step_size == 0.0) {
            // ... presumably this case means current parameter is extremely close to zero or very close to boundary
            params_next(i) = rNorm(params_current(i), step_size/10, walk_rng);
          } else {
            // ... take next step
            params_next(i) = rNorm(params_current(i), bounded_step_size, walk_rng);
          }
          
          // Check that a good baseline parameter hasn't been sent below zero
          if (baseline_mask(i) || beta_Rt_mask(i)) {
            if (params_next(i) < 0.0) {params_next(i) = 0.0;}
          }
          
          // While looping, compute priors for this random step
          std::string this_param = Rcpp::as<std::string>(mf.param_names[i]);
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
      double loglik_current = -compute_bnll(to_sVec(params_current)).val();
      double loglik_next = -compute_bnll(to_sVec(params_next)).val();
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
          if (step >= n_burnin) {
            dVec full_results = to_dVec(params_next);
            full_results.reserve(full_results.size() + 4);
            full_results.push_back(-loglik_next);
            full_results.push_back(-loglik_next - (bd_current_transformed.val() - 1.0));
            full_results.push_back(acceptance); // acceptance ratio
            full_results.push_back(double(ctr + 0.0)); // ctr number
            RMW_steps.row(step - n_burnin) = to_NumVec(full_results);
          }
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
    
    if (use_median) {
      // Compute column-wise median of all rows except the last
      int r_num1 = r_num - 1;
      for (int j = 0; j < c_num; j++) {
        NumericVector col_vals(r_num1); // r_num elements (rows 0 to r_num-1)
        for (int i = 0; i < r_num1; i++) {col_vals[i] = RMW_steps(i, j);}
        std::sort(col_vals.begin(), col_vals.end());
        double median;
        if (r_num1 % 2 == 0) {median = (col_vals[r_num1 / 2 - 1] + col_vals[r_num1 / 2]) / 2.0;} 
        else {median = col_vals[r_num1 / 2];}
        if (j < n_params) {mf.params[j] = median;}
      }
      check_parameter_feasibility(to_sVec(mf.params)); 
      for (int j = 0; j < n_params; j++) {RMW_steps(r_num1, j) = mf.params[j];}
    } else {
      // Reset to initial params 
      mf.params = params_initial;
      // Fit with L-BFGS
      fit(false);
      if (verbose) {vprint("Penalized nll: ", mf.bnll);}
      // ... save results from initial full fit
      dVec res = to_dVec(mf.params);
      res.reserve(res.size() + 4);
      res.push_back(mf.bnll);
      res.push_back(mf.nll);
      res.push_back(double(1.0)); // acceptance ratio
      res.push_back(double(0.0)); // ctr number
      // ... save full fit results in last row of results matrix
      RMW_steps.row(r_num - 1) = to_NumVec(res);
    }
   
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
   
    /*
     * Beta parameters (ref level and fixed effects)
     */
    
    std::vector<std::vector<NumericMatrix>> run_estimates;
    
    // Clear vectors for pushbacks
    beta_idx.clear();
    param_idx = ParamIdx{};
    
    // Loop through model components 
    for (String mc : mf.mc_list) {
      
      // Loop through context values
      std::vector<std::vector<iVec>> beta_idx_mc;
      for (int c = 0; c < n_.context; c++) {
        String cxt = grouping_lvls.context[c];
        
        // Initialize run-estimate vector 
        std::vector<NumericMatrix> run_estimates_cxt;
        
        // Loop through species values
        std::vector<iVec> beta_idx_mc_cxt;
        for (int s = 0; s < n_.species; s++) {
          String sps = grouping_lvls.species[s];
          
          // Extract deg and block num 
          int deg = cp.deg(s, c);
          int n_blocks = deg + 1;
          int bktp = deg;
          if (mc == "Rt") {bktp = n_blocks;}
          
          // Initialize eMat to hold Rt values 
          eMat RtVals(n_.treatments, n_blocks);
          
          NumericMatrix Effs(n_.treatments, bktp);
          NumericMatrix run_estimates_sps(n_.treatments, deg);
          
          iVec beta_idx_mc_cxt_sps;
          if (deg > 0 || mc == "Rt") {
            
            // Compute Effs for this model component, context, and species
            if (mc == "Rt") {
              // Estimate Rt values
              for (int i = 0; i < n_.treatments; i++) {
                NumericVector count_log_avg = sc.count_log_avg[c][s].column(i);
                NumericVector found_cp_trt_i_num = cp.found_trt[c][s].column(i);
                IntegerVector found_cp_trt_i(deg);
                for (int d = 0; d < deg; d++) {
                  found_cp_trt_i[d] = std::round(found_cp_trt_i_num[d]);
                }
                dVec Rt_est(n_blocks); 
                if (n_blocks == 1) { // case when deg = 0
                  Rt_est[0] = vmean(count_log_avg);
                  RtVals(i, 0) = (double)Rt_est[0];
                } else {
                  // Estimate rate values as mean of count values in each block 
                  std::vector<dVec> est_rate_runs = est_bkRates_tRuns(n_blocks, count_log_avg, found_cp_trt_i, cp.rise_threshold_factor);
                  Rt_est = est_rate_runs[0];
                  run_estimates_sps.row(i) = to_NumVec(est_rate_runs[1]);
                  RtVals.row(i) = to_eVec(Rt_est); 
                }
              }
              // Estimate Rt effect values 
              for (int bk = 0; bk < n_blocks; bk++) {Effs.column(bk) = to_NumVec(ev.weight_rows.fullPivLu().solve(RtVals.col(bk)));}
            } else if (mc == "tpoint") {
              // Estimate fixed effects from treatments for each tpoint
              eMat found_cp_trt_transposed = to_eMat(cp.found_trt[c][s]).transpose(); 
              for (int d = 0; d < deg; d++) {
                Effs.column(d) = to_NumVec(ev.weight_rows.fullPivLu().solve(found_cp_trt_transposed.col(d)));
              }
            } else if (mc == "tslope") {
              // Estimate tslopes for each treatment
              NumericMatrix found_slope_trt(deg, n_.treatments);
              for (int i = 0; i < n_.treatments; i++) {
                for (int d = 0; d < deg; d++) {
                  found_slope_trt(d, i) = 4.0/run_estimates[c][s](i, d); 
                  if (found_slope_trt(d, i) < cp.min_initialization_slope) {found_slope_trt(d, i) = cp.min_initialization_slope;}
                  // ^ ... keep from initializing too close to zero
                }
              }
              // Estimate fixed effects from treatments for each tslope
              eMat found_slope_trt_transposed = to_eMat(found_slope_trt).transpose();
              for (int d = 0; d < deg; d++) {Effs.column(d) = to_NumVec(ev.weight_rows.fullPivLu().solve(found_slope_trt_transposed.col(d)));}
              
            }
            
            // For each block/transition point t and treatment i
            for (int t = 0; t < bktp; t++) {
              String t_name = "Tns/Blk" + std::to_string(t + 1);
              for (int i = 0; i < n_.treatments; i++) {
                
                // Add name
                CharacterVector param_name;
                if (i == 0) {param_name = CharacterVector::create("baseline", cxt, mc, sps, t_name);} 
                else {param_name = CharacterVector::create("beta", mc, cxt, sps, ev.trt_lvls[i], "X", t_name);}
                param_names_list.push_back(param_name);
               
                // Collect indices 
                if (i > 0) {
                  if (mc == "Rt") {
                    param_idx.beta_Rt.push_back(idx);
                  } else if (mc == "tslope") {   
                    param_idx.beta_tslope.push_back(idx);
                  } else if (mc == "tpoint") {   
                    param_idx.beta_tpoint.push_back(idx);
                  }   
                } else {
                  param_idx.baseline.push_back(idx);
                  if (mc == "Rt") {
                    param_idx.baseline_Rt.push_back(idx);
                  } else if (mc == "tslope") {   
                    param_idx.baseline_tslope.push_back(idx);
                  } else if (mc == "tpoint") {   
                    param_idx.baseline_tpoint.push_back(idx);
                  }
                }
                beta_idx_mc_cxt_sps.push_back(idx);
                idx++;  
                
                // Make parameter vector
                param_vector.push_back(Effs(i,t));
                
              } 
            }
          } 
          if (mc == "Rt") {run_estimates_cxt.push_back(run_estimates_sps);}
          beta_idx_mc_cxt.push_back(beta_idx_mc_cxt_sps);
        }
        if (mc == "Rt") {run_estimates.push_back(run_estimates_cxt);}
        beta_idx_mc.push_back(beta_idx_mc_cxt);
      }
      beta_idx.push_back(beta_idx_mc);
    }
    
    /*
     * Random-effect parameters
     */
    
    // Unroll the warping factors
    iVec wfactor_idx_point;
    iVec wfactor_idx_rate;
    iVec wfactor_idx_slope;
    for (int s = 0; s < n_.species; s++) {
      // Intentionally skip the first level, which is "none"
      for (int r = 1; r < n_.ran; r++) {
        
        String sps_name = grouping_lvls.species[s];
        String ran_name = grouping_lvls.ran[r];
        // ... Point warp
        CharacterVector param_name_point = CharacterVector::create("wfactor", "point", ran_name, "X", sps_name);
        param_names_list.push_back(param_name_point);
        param_idx.wfactor_point.push_back(idx);
        wfactor_idx_point.push_back(idx);
        param_vector.push_back(R::runif(-0.1, 0.1));
        idx++;
        // ... Rate warp
        CharacterVector param_name_rate = CharacterVector::create("wfactor", "rate", ran_name, "X", sps_name);
        param_names_list.push_back(param_name_rate);
        param_idx.wfactor_rate.push_back(idx);
        wfactor_idx_rate.push_back(idx); 
        param_vector.push_back(R::runif(-0.1, 0.1));
        idx++;
        // ... Slope warp
        CharacterVector param_name_slope = CharacterVector::create("wfactor", "slope", ran_name, "X", sps_name);
        param_names_list.push_back(param_name_slope);
        param_idx.wfactor_slope.push_back(idx);
        wfactor_idx_slope.push_back(idx);
        param_vector.push_back(R::runif(-0.1, 0.1));
        idx++;
        
      }
    } 
    wfactor_idx = {wfactor_idx_point, wfactor_idx_rate, wfactor_idx_slope};
    
    // Reformat parameter names as single strings
    int n_params = param_names_list.size();
    mf.param_names = CharacterVector(n_params);
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
      mf.param_names[n] = name;
    } 
    
    // Wrap parameter vector
    mf.params = Rcpp::wrap(param_vector);
    
    // Check feasibility 
    check_parameter_feasibility(to_sVec(mf.params)); 
    
  }

// Check and correct parameter feasibility
void wspc::check_parameter_feasibility(const sVec& parameters_var) { 
    
    // Initialize vectors to return 
    sVec feasible_parameters_var = parameters_var; 
    sVec predicted_rates_log_var;
    // Compute boundary_dist and take min
    sVec bd = boundary_dist(parameters_var);
    sdouble bd_min = smin(bd);
    bool feasible = bd_min > 0.0;
    // Test if provided parameters produce any nan rates
    if (feasible) {
      predicted_rates_log_var = predict_rates(
        parameters_var, // model parameters for generating prediction
        true            // compute all summed count rows, even with a count value of NA?
      );
      for (int i = 0; i < n_.count_rows; i++) {
        if (std::isnan(predicted_rates_log_var(i))) {
          feasible = false;
          break;
        }
      }
    }
    // Test if provided parameters produce any negative rates
    if (feasible) {
      for (int i = 0; i < n_.count_rows; i++) {
        if (predicted_rates_log_var(i) < rt_lower_bound) {
          feasible = false;
          break;
        }
      }
    }
    
    // If not feasible, attempt to find a feasible starting point
    if (!feasible) {
      
      // Variables for optimization
      dVec x = to_dVec(parameters_var);
      size_t n = x.size();
      
      // Set up NLopt optimizer
      nlopt::opt opt(nlopt::LD_LBFGS, n);
      opt.set_max_objective(wspc::min_boundary_dist_NLopt, this);
      opt.set_ftol_rel(mf.ctol);               // stop when iteration changes objective fn value by less than this fraction 
      opt.set_maxeval(mf.max_evals);           // Maximum number of evaluations to try 
      opt.set_stopval(0.01);                // ensure boundary distance is at least this much above zero (and then stop)
      double max_fx;
      int sc = 0;
      
      // Optimize
      int n_tries = 0;
      int max_tries = 6;
      while (n_tries < max_tries && sc == 0) {
        // Manually reduce t-point effects 
        for (int i : param_idx.beta_tpoint) {x[i] *= 0.5;}
        for (int i : param_idx.wfactor_point) {x[i] *= 0.5;}
        // Manually check baseline t-points 
        for (int i = 1; i < param_idx.baseline_tpoint.size(); i++) {
          int idx0 = param_idx.baseline_tpoint[i-1];
          int idx1 = param_idx.baseline_tpoint[i];
          if (x[idx1] - x[idx0] < cp.tpoint_buffer) {
            x[idx0] -= cp.tpoint_buffer/4.0;
            x[idx1] += cp.tpoint_buffer/4.0;
          }
        }
        // Run nlopt 
        try {
          nlopt::result sc_ = opt.optimize(x, max_fx);
          sc = static_cast<int>(sc_);
        } catch (std::exception& e) {
          sc = 0;
        }
        n_tries++;
      }
      
      // Final check of boundary distance
      if (max_fx <= 0) {sc = 0;}
      
      // Check for success
      if (sc == 0) {
        Rcpp::warning("Could not find nearby feasible parameters (boundary violation or nan rates), returning provided ones");
      } else { // found a feasible starting point, save
        feasible = true;
        feasible_parameters_var = to_sVec(x);
        // Recompute predicted rates 
        predicted_rates_log_var = predict_rates(
          feasible_parameters_var, // model parameters for generating prediction
          true                     // compute all summed count rows, even with a count value of NA?
        );
      }
      
    }
    
    // Save feasible parameters and predicted rates 
    if (feasible) {
      mf.params = to_NumVec(feasible_parameters_var);
      for (int i = 0; i < n_.count_rows; i++) {
        sc.pr_log[i] = predicted_rates_log_var(i).val();
        sc.pr[i] = std::exp(sc.pr_log[i]) - 1.0;
      }
    } 
    
  } 

/*
 * *************************************************************************
 * Export data to R
 */

Rcpp::List wspc::results() {
    
    CharacterVector context(n_.count_rows);
    CharacterVector species(n_.count_rows);
    CharacterVector ran(n_.count_rows);
    CharacterVector treatment(n_.count_rows);
    for (int i = 0; i < n_.count_rows; i++) {
      context[i] = grouping_lvls.context[sc.cxt[i]];
      species[i] = grouping_lvls.species[sc.sps[i]];
      ran[i] = grouping_lvls.ran[sc.ran[i]];
      treatment[i] = ev.trt_lvls[sc.trt[i]];
    }
    // Create summed count data frame
    DataFrame count_data_summed = DataFrame::create(
      _["bin"] = Rcpp::wrap(sc.bin),
      _["count"] = Rcpp::wrap(sc.count),
      _["pred"] = sc.pr,
      _["count.log"] = Rcpp::wrap(sc.count_log),
      _["pred.log"] = sc.pr_log,
      _["context"] = context,
      _["species"] = species,
      _["ran"] = ran,
      _["treatment"] = treatment
    );
    
    // Add parameter names to fitted parameter vector 
    mf.params.names() = mf.param_names;
    
    // Collect fixed-effect names and levels 
    List fix_lvls_list = List(ev.fix_lvls.size());
    for (int i = 0; i < ev.fix_lvls.size(); i++) {
      fix_lvls_list[i] = (CharacterVector)ev.fix_lvls[i];
    }
    List fix_trt_list = List(ev.fix_trt.size());
    for (int i = 0; i < ev.fix_trt.size(); i++) {
      fix_trt_list[i] = (CharacterVector)ev.fix_trt[i];
    }
    fix_lvls_list.names() = ev.fix_names;
    fix_trt_list.names() = ev.fix_names;
    ev.fix_ref.names() = ev.fix_names;
    List fixed_effects = List::create(
      _["names"] = ev.fix_names,
      _["lvls"] = fix_lvls_list,
      _["treat.lvl"] = fix_trt_list,
      _["ref.lvl"] = ev.fix_ref
    );
    
    // Pack up treatment name information
    List treatment_components_list = List(ev.trt_components.size());
    for (int i = 0; i < ev.trt_components.size(); i++) {
      treatment_components_list[i] = (CharacterVector)ev.trt_components[i];
    }
    treatment_components_list.names() = ev.trt_lvls;
    List treat = List::create(
      _["names"] = ev.trt_lvls, 
      _["components"] = treatment_components_list
    );
    
    // Put grouping variable information into list
    List grouping_variables = List::create(
      _["context.lvls"] = grouping_lvls.context,
      _["species.lvls"] = grouping_lvls.species,
      _["ran.lvls"] = grouping_lvls.ran
    );
    
    // Pack up parameter indexes into list
    List paramidx = List::create(
      _["beta"] = Rcpp::wrap(beta_idx),
      _["w.factor"] = Rcpp::wrap(wfactor_idx)
    );
    
    // Reformat gamma dispersion parameters 
    NumericMatrix g_dispersion = to_NumMat(mf.gamma_dispersion);
    rownames(g_dispersion) = grouping_lvls.species;
    colnames(g_dispersion) = grouping_lvls.context;
    
    // Convert weight matrix 
    NumericMatrix weight_matrix = to_NumMat(ev.weight_rows);
    rownames(weight_matrix) = ev.trt_lvls;
    colnames(weight_matrix) = ev.trt_lvls;
    
    // Make final list to return 
    List results_list = List::create(
      _["model.component.list"] = mf.mc_list,
      _["count.data.summed"] = count_data_summed,
      _["LRO.grid.search.results"] = cp.LRO_grid_search_results,
      _["fitted.parameters"] = mf.params,
      _["gamma.dispersion"] = g_dispersion,
      _["param.names"] = mf.param_names,
      _["fix"] = fixed_effects,
      _["treatment"] = treat,
      _["weight.matrix"] = weight_matrix,
      _["grouping.variables"] = grouping_variables,
      _["param.idx0"] = paramidx, // "0" to indicate this goes out w/ C++ zero-based indexing
      _["token.pool"] = Rcpp::wrap(sc.token_pool),
      _["settings"] = model_settings
    );
    
    return results_list;
    
  }

/*
 * *************************************************************************
 * Testing and debugging in R
 */

// Wrap compute_nll in form needed for R
double wspc::compute_nll_debug(
    const dVec& x
  ) {
    
    // Convert to sVec
    sVec parameters_var = to_sVec(x);
    
    // Compute nll
    double negloglik = compute_nll(parameters_var).val();
    
    // Clear stan memory
    clear_stan_mem(); 
    
    // Return as double
    return negloglik;
    
  }

// Wrap compute_bnll in form needed for R
double wspc::compute_bnll_debug(
    const dVec& x
  ) { 
    
    /*
     * This function is just for testing and debugging, e.g., 
     *  comparing stan grad to finite difference in R. 
     */
    
    // Convert dVec to Eigen with stan
    sVec parameters_var = to_sVec(x);
    
    // Compute bounded_nll
    double fx = compute_bnll(parameters_var).val();
    
    // Clear stan memory
    clear_stan_mem(); 
    
    // Return the value of compute_bnll
    return fx; 
    
  } 

// Wrap grad_compute_bnll_debug in form needed for R
Rcpp::NumericVector wspc::grad_compute_bnll_debug(
    const dVec& x 
  ) const { 
    
    /*
     * This function is just for testing and debugging, e.g., 
     *  comparing stan grad to finite difference in R. 
     */
    
    // Convert dVec to Eigen with stan
    sVec parameters_var = to_sVec(x);
    
    // Compute grdient of compute_bnll
    Eigen::VectorXd grad_fx = grad_compute_bnll(parameters_var);
    
    // Cast to NumericVector
    NumericVector grad_fx_R(grad_fx.size());
    for (int i = 0; i < grad_fx.size(); i++) {
      grad_fx_R[i] = grad_fx(i);
    } 
    
    // Return the value of the grad of compute_bnll
    return grad_fx_R; 
    
  }  

// Export the class constructor and select fields and methods to R
RCPP_EXPOSED_CLASS(wspc)
RCPP_MODULE(wspc) {
    class_<wspc>("wspc")
    .constructor<DataFrame, List, bool>()  
    .method("LRO_initial_param_ests", &wspc::LRO_initial_param_ests)
    .method("LRO_grid_search", &wspc::LRO_grid_search)
    .method("compute_nll_debug", &wspc::compute_nll_debug)
    .method("compute_bnll_debug", &wspc::compute_bnll_debug)
    .method("grad_compute_bnll_debug", &wspc::grad_compute_bnll_debug)
    .method("predict_rates_R", &wspc::predict_rates_R)
    .method("fit", &wspc::fit)
    .method("bs_batch", &wspc::bs_batch)
    .method("MCMC", &wspc::MCMC)
    .method("clear_stan_mem", &wspc::clear_stan_mem)
    .method("results", &wspc::results);
  }
