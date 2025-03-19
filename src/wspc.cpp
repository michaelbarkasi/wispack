
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
    NumericVector struc_values_settings = settings["struc_values"];
    double buffer_factor_settings = settings["buffer_factor"];
    double ctol_settings = settings["ctol"];
    double max_penalty_at_distance_factor_settings = settings["max_penalty_at_distance_factor"];
    double LROcutoff_settings = settings["LROcutoff"];
    double LROwindow_factor_settings = settings["LROwindow_factor"];
    double tslope_initial_settings = settings["tslope_initial"];
    double wf_initial_settings = settings["wf_initial"];
    int max_evals_settings = settings["max_evals"];
    int initial_fits_settings = settings["initial_fits"];
    unsigned int rng_seed_settings = settings["rng_seed"];
    for (int i = 0; i < struc_values.size(); i++) {struc_values[i] = struc_values_settings[i];}
    buffer_factor = sdouble(buffer_factor_settings);
    ctol = ctol_settings;
    max_penalty_at_distance_factor = sdouble(max_penalty_at_distance_factor_settings);
    LROcutoff = LROcutoff_settings;
    LROwindow_factor = LROwindow_factor_settings;
    tslope_initial = tslope_initial_settings;
    wf_initial = wf_initial_settings;
    max_evals = max_evals_settings;
    initial_fits = initial_fits_settings;
    rng_seed = rng_seed_settings;
    model_settings = Rcpp::clone(settings);
    
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
    if (tslope_initial < 1.0) {
      Rcpp::stop("tslope_initial must be greater than or equal to 1.0.");
    }
    vprint("Data structure check passed", verbose);
    
    // Save tokenized count column before collapsing to sums 
    count_tokenized = to_sVec(Rcpp::as<NumericVector>(count_data["count"]));
    vprint("Saving tokenized count", verbose);
    
    // Find max bins
    bin_num = smax(to_sVec(Rcpp::as<NumericVector>(count_data["bin"])));
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
    vprint("Extracted fixed effects:", verbose);
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
    vprint("Created treatment levels", verbose);
    vprintV(treatment_lvls, verbose);
    
    // Create matrix to translate between treatment levels and fixed-effects columns
    // ... for making summed-count fixed-effect column 
    // ... each cell should contain the value of the treatment column n_fix, for the given treatment level
    CharacterMatrix effects_rows(treatment_num, n_fix);
    for (int tr = 0; tr < treatment_num; tr++) {
      CharacterVector trt_given = treatment_components[tr];
      for (int f = 0; f < n_fix; f++) {
        CharacterVector lvls = fix_lvls[f];
        effects_rows(tr, f) = lvls[0]; // Assume it's the reference level
        for (String l : lvls) {
          // If treatment level found, replace
          if (any_true(eq_left_broadcast(trt_given, l))) {effects_rows(tr, f) = l;}
        }
      }
    }
    vprint("Created treatment-to-fix translation matrix", verbose);
    
    // Pre-compute weight-matrix rows 
    // ... for making weights matrix
    sMat weight_rows(treatment_num, treatment_num);
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
    vprint("Pre-computed weight matrix rows", verbose);
    
    // Extract grouping variables 
    parent_lvls = Rcpp::sort_unique(Rcpp::as<CharacterVector>(count_data["parent"]));
    child_lvls = Rcpp::sort_unique(Rcpp::as<CharacterVector>(count_data["child"]));
    ran_lvls = Rcpp::sort_unique(Rcpp::as<CharacterVector>(count_data["ran"]));
    ran_lvls.push_front("none");                                        // add "none" to represent no random effect (reference level)
    // ... print extracted grouping variables
    vprint("Extracted parent grouping variables:", verbose);
    vprintV(parent_lvls, verbose);
    vprint("Extracted child grouping variables:", verbose);
    vprintV(child_lvls, verbose);
    vprint("Extracted random-effect grouping variables:", verbose);
    vprintV(ran_lvls, verbose);
    
    // Temporarily extract tokenized-count columns 
    IntegerVector binT = Rcpp::as<IntegerVector>(count_data["bin"]); 
    CharacterVector parentT = Rcpp::as<CharacterVector>(count_data["parent"]);
    CharacterVector childT = Rcpp::as<CharacterVector>(count_data["child"]);
    CharacterVector ranT = Rcpp::as<CharacterVector>(count_data["ran"]);
    vprint("Extracted tokenized count columns", verbose);
    
    // Create summed count data, size constants
    int idx = 0;
    int n_ran = ran_lvls.size();
    int n_child = child_lvls.size();
    int n_parent = parent_lvls.size();
    int bin_num_i = (int)bin_num.val();
    n_count_rows = bin_num_i * n_parent * n_child * n_ran * treatment_num;
    count_row_nums = Rcpp::seq(0, n_count_rows - 1);
    observed_mean_ran_eff.resize(n_ran - 1); 
    observed_mean_ran_eff.setZero(); 
    vprint("Grabbed size constants for summed count data, total rows: " + std::to_string(n_count_rows), verbose);
    
    // Create summed count data rows, initializations
    count.resize(n_count_rows);
    bin.resize(n_count_rows);
    parent = CharacterVector(n_count_rows);
    child = CharacterVector(n_count_rows);
    ran = CharacterVector(n_count_rows);
    treatment = CharacterVector(n_count_rows);
    weights.resize(n_count_rows, treatment_num);
    vprint("Initialized columns for summed count data", verbose);
    
    // Initiate count indexes 
    int idx_mcu = 0;
    idx_mc_unique = IntegerVector(n_count_rows/bin_num_i);
    token_pool.resize(n_count_rows);
    count_not_na_mask = LogicalVector(n_count_rows);
    count_not_na_mask.fill(false);
    vprint("Initialized count indexes, number of rows with unique model components: " + std::to_string(idx_mc_unique.size()), verbose);
    
    // Pre-compute bin masks
    vprint("Pre-computing bin masks", verbose);
    LogicalMatrix bin_masks(binT.size(), bin_num_i);
    for (int b = 0; b < bin_num_i; b++) { 
      bin_masks.column(b) = eq_left_broadcast(binT, b + 1);
    }
    
    // Create summed count data columns and weight matrix
    vprint("Creating summed-count data columns ...", verbose); 
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
    vprint("Extracted non-NA indexes", verbose);
    
    // Make extrapolation pool and extrapolate "none" rows
    vprint("Making extrapolation pool ...", verbose);
    extrapolation_pool.resize(count.size());
    extrapolation_pool = make_extrapolation_pool(bin, count, parent, child, ran, treatment, verbose); 
    
    // Extrapolate "none" rows
    count = extrapolate_none(count, ran, extrapolation_pool);
    vprint("Extrapolated 'none' rows", verbose);
    
    // Take log of observed counts 
    count_log.resize(n_count_rows); 
    for (int r = 0; r < n_count_rows; r++) {
      count_log(r) = slog(count(r) + 1.0);
    }
    mean_count_log = vmean(count_log);
    vprint("Took log of observed counts", verbose);
    
    // Resize gamma dispersion matrix (... fill in next step)
    gamma_dispersion.resize(n_child, n_parent);
    gamma_dispersion.setZero();
    gd_child = IntegerVector(n_child);
    gd_parent = IntegerVector(n_parent);
    gd_child.names() = child_lvls;
    gd_parent.names() = parent_lvls;
    
    // Find mean observed ran effect per ran level
    for (int r = 1; r < n_ran; r++) {
      LogicalVector ran_mask = eq_left_broadcast(ran, ran_lvls[r]) & count_not_na_mask;
      IntegerVector ran_idx = Rwhich(ran_mask);
      for (int i = 0; i < ran_idx.size(); i++) {
        observed_mean_ran_eff[r - 1] += count_log(ran_idx(i)); 
      }
    }
    sdouble ran_count_log_mean = 0.0; 
    for (int r = 0; r < n_ran - 1; r++) {ran_count_log_mean += observed_mean_ran_eff[r];}
    ran_count_log_mean /= n_ran - 1.0; 
    for (int r = 0; r < n_ran - 1; r++) {
      observed_mean_ran_eff[r] /= ran_count_log_mean;
      observed_mean_ran_eff[r] -= 1.0;
    }
    vprint("Found mean observed ran effect per ran level", verbose);
    
    // Initialize matrix to hold degrees of each parent-child combination
    degMat = IntegerMatrix(n_child, n_parent);
    rownames(degMat) = child_lvls;
    colnames(degMat) = parent_lvls;
    
    // Compute running and filter window sizes for LRO change-point detection
    int ws = static_cast<int>(std::round(LROwindow_factor * (double)bin_num_i * buffer_factor.val()));
    int filter_ws = std::round(ws/2);
    int n_ran_trt = n_ran * treatment_num;
    
    // Estimate degree of each parent-child combination at baseline using LRO change-point detection 
    // ... store ref values in list
    List ref_values(n_parent);
    ref_values.names() = parent_lvls;
    // ... store rate effects in list
    List RtEffs(n_parent);
    RtEffs.names() = parent_lvls;
    // ... store tpoint effects in list
    List tpointEffs(n_parent);
    tpointEffs.names() = parent_lvls;
    // ... store change points in list
    List found_cp_list(n_parent);
    found_cp_list.names() = parent_lvls;
    
    // Loop through parent levels
    for (int p = 0; p < n_parent; p++) {
      
      LogicalVector parent_mask = eq_left_broadcast(parent, parent_lvls[p]);
      gd_parent[(String)parent_lvls[p]] = p;
      
      // Set up list for ref values
      ref_values[(String)parent_lvls[p]] = List(n_child);
      name_proxylist(ref_values[(String)parent_lvls[p]], child_lvls);
      // ... for rate effect values
      RtEffs[(String)parent_lvls[p]] = List(n_child);
      name_proxylist(RtEffs[(String)parent_lvls[p]], child_lvls);
      // ... for tpoint effect values 
      tpointEffs[(String)parent_lvls[p]] = List(n_child);
      name_proxylist(tpointEffs[(String)parent_lvls[p]], child_lvls);
      // ... for change points 
      found_cp_list[(String)parent_lvls[p]] = List(n_child);
      name_proxylist(found_cp_list[(String)parent_lvls[p]], child_lvls);
      
      // Loop through child levels
      for (int c = 0; c < n_child; c++) { 
        
        // Estimate dispersion of raw count (not log)
        LogicalVector child_mask = eq_left_broadcast(child, child_lvls[c]);
        LogicalVector pc_mask = parent_mask & child_mask & count_not_na_mask;
        sVec count_pc_masked = masked_vec(count, pc_mask); 
        sdouble count_pc_mean = vmean(count_pc_masked);
        sdouble count_pc_var = vvar(count_pc_masked);
        gd_child[(String)child_lvls[c]] = c;
        if (count_pc_var > count_pc_mean) {
          // Have: count_pc_var = count_pc_mean + count_pc_mean^2 * gdis
          // ... count_pc_var = count_pc_mean * (1 + count_pc_mean * gdis)
          // ... count_pc_var / count_pc_mean = 1 + count_pc_mean * gdis
          // ... (count_pc_var / count_pc_mean) - 1 = count_pc_mean * gdis
          gamma_dispersion(c, p) = ((count_pc_var / count_pc_mean) - 1) / count_pc_mean;
        }
        
        sMat count_masked_array(bin_num_i, n_ran_trt);
        count_masked_array.setZero();
        LogicalVector good_col(n_ran_trt);
        
        // Collect count values for each treatment-ran level interaction of this child-parent pair
        NumericMatrix count_avg_mat(bin_num_i, treatment_num);
        for (int t = 0; t < treatment_num; t++) {
          String trt = treatment_lvls[t];
         
          // Make mask for treatment rows of this parent-child pair
          LogicalVector trt_mask = eq_left_broadcast(treatment, trt);
          LogicalVector mask = pc_mask & trt_mask;
         
          // Grab average counts for this treatment level (used later to set initial values)
          dVec count_avg(bin_num_i); 
          LogicalVector mask0 = count_not_na_mask & parent_mask & child_mask & trt_mask;
          sVec count_trt_masked = masked_vec(count_log, mask0);
          sVec bin_trt_masked = masked_vec(bin, mask0);
          for (int b = 0; b < bin_num_i; b++) {
            LogicalVector mask_b = eq_left_broadcast(to_NumVec(bin_trt_masked), (double)b + 1.0);
            sVec count_b = masked_vec(count_trt_masked, mask_b);
            count_avg[b] = vmean(count_b).val();
          }
          count_avg_mat.column(t) = to_NumVec(count_avg);
          
          // Collect count values for each ran level and this treatment trt
          for (int r = 0; r < n_ran; r++) {
            // Make mask for ran level rows of this treatment (of this parent-child pair)
            LogicalVector ran_mask = eq_left_broadcast(ran, ran_lvls[r]);
            LogicalVector mask = mask0 & ran_mask;
            
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
          count_masked_array_good,    // series to scan
          ws,                         // running window size 
          filter_ws,                  // size of window for taking rolling mean
          LROcutoff                   // points more than this times sd considered outliers
        );
        
        // Estimate degree of this parent-child pair 
        int deg = found_cp_good.rows();
        int n_blocks = deg + 1;
        degMat(c, p) = deg;
        
        // Fill removed empty columns back into the found_cp matrix
        IntegerMatrix found_cp(deg, n_ran_trt);
        if (deg > 0) {
          for (int si = 0; si < good_col_idx.size(); si++) {
            int s = good_col_idx[si];
            found_cp.column(s) = found_cp_good.column(si);
          }
        }
        
        // Save
        assign_proxylist(found_cp_list[(String)parent_lvls[p]], (String)child_lvls[c], found_cp);
        
        // Extract treatment means 
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
       
        // Set initial parameters for fixed-effect treatments
        List placeholder_ref_values(mc_list.size());
        placeholder_ref_values.names() = mc_list;
        NumericMatrix placeholder_RtEffs(treatment_num, n_blocks);
        NumericMatrix placeholder_tpointEffs(treatment_num, deg);
        
        // Loop through model components (Rt, tslope, tpoint)
        for (String mc : mc_list) {
          
          // Mean rate per block
          if (mc == "Rt") {
            dVec Rt_est(n_blocks); 
            sMat RtVals(treatment_num, n_blocks);
           
            // Loop through treatments
            for (int t = 0; t < treatment_num; t++) {
             
              NumericVector count_avg = count_avg_mat.column(t);
              NumericVector found_cp_trt_t_num = found_cp_trt.column(t);
              IntegerVector found_cp_trt_t(deg);
              for (int d = 0; d < deg; d++) {
                found_cp_trt_t[d] = std::round(found_cp_trt_t_num[d]);
              }
              
              if (n_blocks == 1) { // case when deg = 0
                Rt_est[0] = vmean(count_avg);
                RtVals(t, 0) = (sdouble)Rt_est[0];
              } else {
                
                // Estimate rate values as mean of count values in each block
                for (int bk = 0; bk < n_blocks; bk++) {
                  int bin_start;
                  int bin_end;
                  if (bk == 0) {
                    bin_start = 0;
                    bin_end = found_cp_trt_t[bk];
                  } else if (bk == deg) {
                    bin_start = found_cp_trt_t[bk - 1] - 1;
                    bin_end = bin_num_i - 1;
                  } else {
                    bin_start = found_cp_trt_t[bk - 1] - 1;
                    bin_end = found_cp_trt_t[bk];
                  }
                  Rt_est[bk] = vmean_range(count_avg, bin_start, bin_end);
                  RtVals(t, bk) = (sdouble)Rt_est[bk];
                }
                
              }
              if (t == 0) {placeholder_ref_values[mc] = to_NumVec(Rt_est);}
             
            }
           
            // Estimate fixed effects from treatments for each block, rate
            for (int bk = 0; bk < n_blocks; bk++) {
              sVec beta_bk_Rt = weight_rows.fullPivLu().solve(RtVals.col(bk));
              placeholder_RtEffs.column(bk) = to_NumVec(beta_bk_Rt);
            }
            
          } else if (deg > 0) { 
            
            // Handle tpoint and tslope
            if (mc == "tpoint") {
              
              // Grab reference values 
              placeholder_ref_values[mc] = found_cp_trt.column(0);
              
              // Estimate fixed effects from treatments for each block, tpoint
              // ... put into Eigen for solving and make rows trts, cols degrees
              sMat found_cp_trt_transposed = to_sMat(found_cp_trt).transpose(); 
              for (int d = 0; d < deg; d++) {
                sVec beta_bk_tpoint = weight_rows.fullPivLu().solve(found_cp_trt_transposed.col(d));
                placeholder_tpointEffs.column(d) = to_NumVec(beta_bk_tpoint);
              }
              
            } else if (mc == "tslope") {
              
              // Grab reference values 
              NumericVector tslope_default(deg);
              for (int tp = 0; tp < deg; tp++) {tslope_default[tp] = tslope_initial;}
              placeholder_ref_values[mc] = tslope_default;
              
            }
          }
         
        }
        assign_proxylist(ref_values[(String)parent_lvls[p]], (String)child_lvls[c], placeholder_ref_values);
        assign_proxylist(RtEffs[(String)parent_lvls[p]], (String)child_lvls[c], placeholder_RtEffs);
        assign_proxylist(tpointEffs[(String)parent_lvls[p]], (String)child_lvls[c], placeholder_tpointEffs);
       
      }
    }
    // ... save found change points 
    change_points = found_cp_list;
    vprint("Estimated change points with LROcp and found initial parameter estimates for fixed-effect treatments", verbose); 
   
    // Build default fixed-effects matrices in shell
    List beta = build_beta_shell(mc_list, treatment_lvls, parent_lvls, child_lvls, ref_values, RtEffs, tpointEffs, degMat);
    vprint("Built initial beta (ref and fixed-effects) matrices", verbose); 
    
    // Initialize random effect warping factors 
    List wfactors = List(2);
    CharacterVector wfactors_names = {"point", "rate"};
    wfactors.names() = wfactors_names;
    double wf_initial_low = -1.0 * wf_initial;
    NumericVector wfs = dseq(wf_initial_low, wf_initial, n_ran - 1); 
    
    // Loop through warping factors (point and rate)
    for (String wf : wfactors_names) {
      
      // Initialize array to hold warping factors 
      NumericMatrix wf_array(n_ran, n_child);
      
      // Loop through child levels
      for (int c = 0; c < n_child; c++) {
        
        if (wf == "rate") {
          // Find mean count value for each random level
          dVec ran_lvl_means(n_ran - 1);
          LogicalVector child_mask = eq_left_broadcast(child, child_lvls[c]);
          for (int r = 1; r < n_ran; r++) { // skip "none" 
            LogicalVector r_mask = child_mask & eq_left_broadcast(ran, ran_lvls[r]);
            ran_lvl_means[r - 1] = vmean(masked_vec(count, r_mask)).val();
          }
          IntegerVector sorted_idx = Rorder(ran_lvl_means);
          NumericVector wfs_c = wfs[sorted_idx];
          wfs_c.push_front(0.0);            // add "none"
          wf_array.column(c) = wfs_c;
        } else { 
          NumericVector wfs_c(n_ran, 0.0);
          // ^ ... revise later to initiate in correct direction, as in rate??
          wf_array.column(c) = wfs_c;
        } 
        
      } 
      wfactors[wf] = wf_array;
    } 
    vprint("Initialized random effect warping factors", verbose);
    
    // Make and map parameter vector
    List params = make_parameter_vector(
      beta, wfactors,
      struc_values,
      parent_lvls, child_lvls, ran_lvls,
      mc_list, 
      treatment_lvls,
      struc_names,
      degMat
    );
    vprint("Made and mapped parameter vector", verbose);
    
    // Extract parameter vector information
    fitted_parameters = params["param_vec"];
    fitted_parameters_seed = fitted_parameters;
    param_names = params["param_names"];
    param_wfactor_point_idx = params["param_wfactor_point_idx"]; 
    param_wfactor_rate_idx = params["param_wfactor_rate_idx"];
    param_beta_Rt_idx = params["param_beta_Rt_idx"];
    param_beta_Rt_idx_no_ref = params["param_beta_Rt_idx_no_ref"];
    param_beta_tslope_idx = params["param_beta_tslope_idx"];
    param_beta_tslope_idx_no_ref = params["param_beta_tslope_idx_no_ref"];
    param_beta_tpoint_idx = params["param_beta_tpoint_idx"];
    param_beta_tpoint_idx_no_ref = params["param_beta_tpoint_idx_no_ref"];
    param_ref_values_tpoint_idx = params["param_ref_values_tpoint_idx"];
    param_struc_idx = params["param_struc_idx"]; 
    param_struc_idx.names() = struc_names;
    beta_idx = params["beta_idx"];
    wfactor_idx = params["wfactor_idx"];
    if (verbose) {vprint("Number of parameters: ", (int)fitted_parameters.size());}
   
    // Construct grouping variable ids as indexes for warping factor matrices, by count row
    gv_ranLr_int = IntegerVector(n_count_rows);
    gv_fixLr_int = IntegerVector(n_count_rows);
    for (int r = 0; r < n_count_rows; r++) {
      CharacterVector::iterator it_ran = std::find(ran_lvls.begin(), ran_lvls.end(), ran[r]);
      CharacterVector::iterator it_fix = std::find(child_lvls.begin(), child_lvls.end(), child[r]);
      gv_ranLr_int[r] = std::distance(ran_lvls.begin(), it_ran);
      gv_fixLr_int[r] = std::distance(child_lvls.begin(), it_fix);
    }
    vprint("Constructed grouping variable IDs", verbose);
    
    // Compute tpoint buffer
    tpoint_buffer = bin_num * buffer_factor;
    
    // Compute size of the parameter boundary vector 
    for (int r : idx_mc_unique) {
      // Grab degree for this row
      int c_num = Rwhich(eq_left_broadcast(child_lvls, child[r]))[0];
      int p_num = Rwhich(eq_left_broadcast(parent_lvls, parent[r]))[0];
      int deg = degMat(c_num, p_num);
      if (deg > 0){
        // Add slots for the tpoint boundary distance at each tpoint
        boundary_vec_size += deg + 1;
        // Add one slot for the R_sum boundary distance
        boundary_vec_size++;
      } else {
        // Add one slot for the R_sum boundary distance
        boundary_vec_size++;
      } 
    }
    boundary_vec_size += param_struc_idx.size() +  // Add slots for the structural parameter boundaries
      param_wfactor_point_idx.size()*2 +           // ... and for the warping factor boundaries
      param_wfactor_rate_idx.size()*2; 
    vprint("Computed size of boundary vector: " + std::to_string(boundary_vec_size), verbose);
    
    // Pre-compute bin ranges for bootstrap resampling 
    bs_bin_ranges = IntegerMatrix(n_count_rows, 3);
    for (int r : count_not_na_idx) {
      if (ran[r] != "none") {
        // ... find the treatment level and random level of this row
        int trt_num = Rwhich(eq_left_broadcast(treatment_lvls, treatment[r]))[0];
        int ran_num = Rwhich(eq_left_broadcast(ran_lvls, ran[r]))[0];
        // ... find the change points for the parent-child pair of this row
        List p_changepoints = change_points[(String)parent[r]];
        IntegerMatrix pc_changepoints = p_changepoints[(String)child[r]];
        NumericVector pc_changepoints_r = to_NumVec(pc_changepoints.column(trt_num * n_ran + ran_num));
        iVec block_range_idx = block_idx(pc_changepoints_r, bin(r).val());
        // ... find the range of bins in the same block as this row's bin
        double bin_low = 1;
        if (block_range_idx[0] >= 0) {
          bin_low = pc_changepoints_r[block_range_idx[0]];
        }
        double bin_high = bin_num.val();
        if (block_range_idx[1] >= 0) {
          bin_high = pc_changepoints_r[block_range_idx[1]];
        } 
        // ... narrow range by 1 in each direction and pack up range information
        int shift_back_max = static_cast<int>(std::round(bin(r).val() - bin_low)) - 1;
        int shift_up_max = static_cast<int>(std::round(bin_high - bin(r).val())) - 1;
        if (shift_back_max < 0) {shift_back_max = 0;}
        if (shift_up_max < 0) {shift_up_max = 0;}
        int shift_range = shift_up_max + shift_back_max + 1;
        IntegerVector back_up_range = {shift_back_max, shift_up_max, shift_range};
        bs_bin_ranges.row(r) = back_up_range;
      }
    }
    vprint("Pre-computed bin ranges for bootstrap resampling", verbose);
    
    // Initialize list to hold results from model fit
    optim_results = List::create(
      _["fitted_parameters"] = NumericVector(), 
      _["penalized_neg_loglik"] = NA_REAL,
      _["neg_loglik"] = NA_REAL, 
      _["success_code"] = NA_INTEGER,
      _["num_evals"] = NA_INTEGER,
      _["bs_times"] = NumericVector()
    );
    
    vprint("Finished initializing wspc object", verbose);
    
  }

// Destructor
wspc::~wspc() {}

// R copy 
wspc wspc::clone() const {
  wspc this_copy = wspc(*this);
  return this_copy;
}

/*
 * *************************************************************************
 * Computing predicted values from parameters
 */

// Compute model component values for rows of summed count data
// ... for Rt
sVec wspc::compute_mc_Rt_r(
    const int& r,
    const sVec& parameters,
    const sdouble& f_rw
  ) const {
    
    // Extract the parameter vector indexes for the current rate row, beta matrices
    List beta_idx_Rt = Rcpp::as<List>(beta_idx["Rt"]);
    List beta_idx_Rt_prt = Rcpp::as<List>(beta_idx_Rt[(String)parent[r]]);
    IntegerVector betas_Rt_idx = beta_idx_Rt_prt[Rcpp::as<std::string>(child[r])];
    
    // Grab degree
    int c_num = Rwhich(eq_left_broadcast(child_lvls, child[r]))[0];
    int p_num = Rwhich(eq_left_broadcast(parent_lvls, parent[r]))[0];
    int deg = degMat(c_num, p_num);
    
    // Re-roll beta matrices
    int idx_r = 0;
    sMat betas_Rt_r(treatment_num, deg + 1);
    for (int t = 0; t < deg + 1; t++) {
      for (int i = 0; i < treatment_num; i++) {
        betas_Rt_r(i, t) = parameters(betas_Rt_idx[idx_r]);
        idx_r++;
      }
    } 
    
    // Extract weight matrix row
    sVec weight_row = weights.row(r).transpose();
    
    // Compute unwarped Rt
    sVec Rt_vec = compute_mc(betas_Rt_r, weight_row);
    
    // Apply warping factor 
    for (int bt = 0; bt < Rt_vec.size(); bt++) {
      Rt_vec(bt) = rate_warp(Rt_vec(bt), f_rw);
    }
    
    // Send out
    return Rt_vec;
    
  }

// ... tslope
sVec wspc::compute_mc_tslope_r(
    const int& r,
    const sVec& parameters
  ) const {
    
    // Extract the parameter vector indexes for the current rate row, beta matrices
    List beta_idx_tslope = Rcpp::as<List>(beta_idx["tslope"]);
    List beta_idx_tslope_prt = Rcpp::as<List>(beta_idx_tslope[Rcpp::as<std::string>(parent[r])]);
    IntegerVector betas_tslope_idx = beta_idx_tslope_prt[Rcpp::as<std::string>(child[r])];
    
    // Grab degree
    int c_num = Rwhich(eq_left_broadcast(child_lvls, child[r]))[0];
    int p_num = Rwhich(eq_left_broadcast(parent_lvls, parent[r]))[0];
    int deg = degMat(c_num, p_num);
    
    // Re-roll beta matrices
    int idx_r = 0;
    sMat betas_tslope_r(treatment_num, deg);
    if (deg > 0) {
      for (int t = 0; t < deg; t++) {
        for (int i = 0; i < treatment_num; i++) {
          betas_tslope_r(i, t) = parameters(betas_tslope_idx[idx_r]);
          idx_r++;
        }
      } 
    }
    
    // Extract weight matrix row
    sVec weight_row = weights.row(r).transpose();
    
    // Compute tslope
    sVec tslope_vec = compute_mc(betas_tslope_r, weight_row);
    
    // Take exponential of slopes
    // ... slopes must be positive, but fit and reported as a normal parameter,
    //      so here the exp is taken.
    for (int i = 0; i < tslope_vec.size(); i++) {
      tslope_vec(i) = sexp(tslope_vec(i));
    }
    
    // Send out
    return tslope_vec;
    
  } 

// ... tpoint
sVec wspc::compute_mc_tpoint_r(
    const int& r,
    const sVec& parameters
  ) const {
    
    // Extract the parameter vector indexes for the current rate row, beta matrices
    List beta_idx_tpoint = Rcpp::as<List>(beta_idx["tpoint"]);
    List beta_idx_tpoint_prt = Rcpp::as<List>(beta_idx_tpoint[Rcpp::as<std::string>(parent[r])]);
    IntegerVector betas_tpoint_idx = beta_idx_tpoint_prt[Rcpp::as<std::string>(child[r])];
    
    // Grab degree
    int c_num = Rwhich(eq_left_broadcast(child_lvls, child[r]))[0];
    int p_num = Rwhich(eq_left_broadcast(parent_lvls, parent[r]))[0];
    int deg = degMat(c_num, p_num);
    
    // Re-roll beta matrices
    int idx_r = 0;
    sMat betas_tpoint_r(treatment_num, deg);
    if (deg > 0) {
      for (int t = 0; t < deg; t++) {
        for (int i = 0; i < treatment_num; i++) {
          betas_tpoint_r(i, t) = parameters(betas_tpoint_idx[idx_r]);
          idx_r++;
        }
      } 
    }
    
    // Extract weight matrix row
    sVec weight_row = weights.row(r).transpose();
    
    // Find tpoint and send out
    sVec tpoint_vec = compute_mc(betas_tpoint_r, weight_row);
    
    // Send out 
    return tpoint_vec;
    
  }  

// Predict log of rates
sVec wspc::predict_rates_log(
    const sVec& parameters,
    const bool& all_rows 
  ) const {
    
    /*
     * If "all_rows" is true, will compute all summed count rows, even with a count value of NA.
     * 
     * Function is called "log" because we assume the parameters are for counts
     *   that have been put into log space (plus 1). Hence, model parameters always 
     *   predict the log of the observed count rate.
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
    int f_pw_idx;
    int f_rw_idx;
    sdouble f_pw;
    sdouble f_rw;
    
    // Compute predicted rate for rows of the summed count data
    for (int r = 0; r < n_count_rows; r++) {
      
      // Skip rows with NA values in count, if requested
      if (count_not_na_mask[r] || all_rows) {
        
        // Grab warping factors
        if (gv_ranLr_int[r] == 0) {
          f_pw = 0.0;
          f_rw = 0.0; 
        } else {
          f_pw_idx = wfactor_idx_point(gv_ranLr_int[r], gv_fixLr_int[r]);
          f_rw_idx = wfactor_idx_rate(gv_ranLr_int[r], gv_fixLr_int[r]);
          f_pw = parameters(f_pw_idx);
          f_rw = parameters(f_rw_idx); 
        } 
        
        // Only update predicted model components if r begins a new batch of unique values 
        int cnt = std::count(mc_unique_rows.begin(), mc_unique_rows.end(), r);
        if (cnt > 0) { 
          
          // Compute model components for this row r
          Rt = compute_mc_Rt_r(r, parameters, f_rw);   // returned already warped
          tslope = compute_mc_tslope_r(r, parameters); // returned out of log space
          tpoint = compute_mc_tpoint_r(r, parameters);
          
        } 
        
        // Compute the predicted rate
        int c_num = Rwhich(eq_left_broadcast(child_lvls, child[r]))[0];
        int p_num = Rwhich(eq_left_broadcast(parent_lvls, parent[r]))[0];
        
        // Compute the actual poly-sigmoid!! 
        predicted_rates(r) = poly_sigmoid(
          point_warp(bin(r), bin_num, f_pw),
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

// Compute neg_loglik of the model under the given parameters
sdouble wspc::neg_loglik(
    const sVec& parameters
  ) const {
   
    // Initialize variable to hold (what will become) nll
    sdouble log_lik = 0.0;
    
    // Subset the wfactor and beta parameters
    sVec warping_factors_point = idx_vec(parameters, param_wfactor_point_idx);
    sVec warping_factors_rate = idx_vec(parameters, param_wfactor_rate_idx);
    sVec tslope_beta_values_no_ref = idx_vec(parameters, param_beta_tslope_idx_no_ref);
    sVec Rt_beta_values_no_ref = idx_vec(parameters, param_beta_Rt_idx_no_ref);
    sVec tpoint_beta_values_no_ref = idx_vec(parameters, param_beta_tpoint_idx_no_ref);
    
    // *****************************************************************************************************************
    // log-likelihood of warping factors (ensures warping factors fit a modelled beta distribution)
    
    // Compute the log-likelihood of the point warping factors, given the assumed beta distribution
    sdouble beta_shape_point = parameters[param_struc_idx["beta_shape_point"]];
    for (int i = 0; i < warping_factors_point.size(); i++) {
      sdouble wf = (warping_factors_point(i) + 1.0)/2.0; // rescale
      log_lik += slog(
        spower(wf, beta_shape_point - 1.0) * spower(1.0 - wf, beta_shape_point - 1.0) /
          (spower(stan::math::tgamma(beta_shape_point), 2.0) / stan::math::tgamma(2.0 * beta_shape_point))
      );
    }
    
    // Compute the log-likelihood of the rate warping factors, given the assumed beta distribution
    sdouble beta_shape_rate = parameters[param_struc_idx["beta_shape_rate"]];
    for (int i = 0; i < warping_factors_rate.size(); i++) {
      sdouble wf = (warping_factors_rate(i) + 1.0)/2.0; // rescale
      log_lik += slog(
        spower(wf, beta_shape_rate - 1.0) * spower(1.0 - wf, beta_shape_rate - 1.0) /
          (spower(stan::math::tgamma(beta_shape_rate), 2.0) / stan::math::tgamma(2.0 * beta_shape_rate))
      );
    }
    
    // *****************************************************************************************************************
    // log-likelihood of warping factors means (ensures warping factors tend to be centered around zero)
    
    // Compute the log-likelihood of the mean of point warping factors, given the expected normal distribution
    // ... the relevant probability distribution is a normal distribution of mean zero, sd = sqrt(1/((4*m)*(2*s + 1)))
    // ... as the actual warping factors are scaled to range from -1 to 1, the expected sd is twice the above formula
    sdouble pw_mean_p = vmean(warping_factors_point); 
    sdouble sd_pw_p = ssqrt(1.0 / ((4.0 * (sdouble)warping_factors_point.size()) * (2.0 * beta_shape_point + 1.0))) * 2.0;
    log_lik += log_dNorm(pw_mean_p, 0.0, sd_pw_p);
    
    // Compute the log-likelihood of the mean of rate warping factors, given the expected normal distribution
    // ... see above notes on the formula
    sdouble pw_mean_r = vmean(warping_factors_rate); 
    sdouble sd_pw_r = ssqrt(1.0 / ((4.0 * (sdouble)warping_factors_rate.size()) * (2.0 * beta_shape_rate + 1.0))) * 2.0;
    log_lik += log_dNorm(pw_mean_r, 0.0, sd_pw_r);
    
    // *****************************************************************************************************************
    // log-likelihood of the estimated overall rate warp of that warping factor,
    //  given the modelled child-specific rate warping factors for each random level, 
    // ... idea: warping factors for a random level can vary between child levels, but the mean of these warping factors 
    //      should still tend to line up with the observed mean random effect of that random level.
    
    // Compute the log-likelihood of the observed mean random rate effects, given these rate warping factors
    int n_child = child_lvls.size();
    int n_ran = ran_lvls.size(); 
    sdouble sqrt_n_ran = ssqrt((sdouble)n_ran - 1.0);
    // ... Grab warping indices and initiate variables to hold them
    NumericMatrix wfactor_idx_rate = wfactor_idx["rate"];
    sVec f_rw_row(n_child);
    for (int r = 1; r < n_ran; r++) {
      // ... Get rate-warp factors for this level (one per child)
      for (int c = 0; c < n_child; c++) {f_rw_row(c) = parameters[(int)wfactor_idx_rate(r, c)];}
      // ... find mean and scaled sd
      sdouble modeled_mean = vmean(f_rw_row); 
      sdouble modeled_sd = vsd(f_rw_row)/sqrt_n_ran; 
      // ... add log-likelihood
      log_lik += log_dNorm(observed_mean_ran_eff[r - 1], modeled_mean, modeled_sd);
    }
    
    // *****************************************************************************************************************
    // log-likelihood of beta values
    
    // Compute the log-likelihood of the Rt beta values, given the normal distribution implied by the beta-rate shape and fe_difference_ratio ratio
    sdouble expected_ran_effect_rate = 1.0 / ssqrt(4.0 * beta_shape_rate + 2.0);    // ... already includes the *= 2.0 rescaling to range from -1 to 1
    sdouble eff_mult_rate = fe_difference_ratio_Rt - 1.0;
    sdouble expected_Rt = mean_count_log;
    sdouble sd_Rt_effect = eff_mult_rate * expected_Rt * expected_ran_effect_rate;
    for (int i = 0; i < Rt_beta_values_no_ref.size(); i++) {
      log_lik += log_dNorm(Rt_beta_values_no_ref(i), 0.0, sd_Rt_effect);
    }
    
    // Compute the log-likelihood of the tslope beta values, given the assumed normal distribution
    // ... recall: slopes must be positive, but are fit and reported as a normal parameter
    int n_tslope = tslope_beta_values_no_ref.size();
    sdouble sd_tslope_effect = parameters[param_struc_idx["sd_tslope_effect"]];
    for (int i = 0; i < n_tslope; i++) {
      log_lik += log_dNorm(tslope_beta_values_no_ref(i), 0.0, sd_tslope_effect); 
    }
    
    // Compute the log-likelihood of the tpoint beta values, given the assumed normal distribution
    int n_tpoint = tpoint_beta_values_no_ref.size();
    sdouble expected_ran_effect_tpoint = 1.0 / ssqrt(4.0 * beta_shape_point + 2.0); // ... already includes the *= 2.0 rescaling to range from -1 to 1
    sdouble eff_mult_tpoint = fe_difference_ratio_tpoint + 1.0;                     // ... note the different sign
    sdouble expected_tpoint = bin_num / 2.0;
    sdouble sd_tpoint_effect = (eff_mult_tpoint * expected_tpoint * expected_ran_effect_tpoint) / (2.0 + eff_mult_tpoint * expected_ran_effect_tpoint);
    for (int i = 0; i < n_tpoint; i++) {
      log_lik += log_dNorm(tpoint_beta_values_no_ref(i), 0.0, sd_tpoint_effect);
    }
    
    // Check for zero likelihood
    if (std::isinf(log_lik) || std::isnan(log_lik)) {
      // A nan most likely means sd_tpoint_effect when negative
      return inf_;
    }
    
    // *****************************************************************************************************************
    // log-likelihood of observed counts, given predicted rates
    
    // Predict rates under these parameters
    sVec prate_log_var = predict_rates_log(
      parameters, // model parameters for generating prediction
      false       // compute all summed count rows, even with a count value of NA?
      );
    
    // Compute the log-likelihood of the count data, assuming a Poisson distribution with Gamma kernel for over-dispersion
    for (int r : count_not_na_idx) {
      
      if (std::isinf(prate_log_var(r)) || prate_log_var(r) < 0.0 || std::isnan(prate_log_var(r))) {
        return inf_;
      } else {
        
        // Find gamma variance for this row
        // ... grab parent and child index numbers
        int n_c = gd_child[(String)child[r]];
        int n_p = gd_parent[(String)parent[r]];
        // ... pull predicted rate from log space
        sdouble prate_var = sexp(prate_log_var(r)) - 1.0;
        // ... estimate variance of rate outside log space
        sdouble gamma_variance = prate_var + gamma_dispersion(n_c, n_p) * prate_var * prate_var;
        // ... estimate the corresponding variance of the rate back in log space 
        gamma_variance = delta_var_est(gamma_variance, prate_var);
        
        // Analytic solution to the log of the integral from 1 to positive infinity of the Poisson times Gamma densities
        log_lik += slog(poisson_gamma_integral(count_log(r), prate_log_var(r), gamma_variance));
        
      }
      
    }
    
    // Take negative and average, the latter so we're looking at values in the same range, regardless of data size
    sdouble negloglik = -log_lik / count_not_na_idx.size();
    
    if (std::isinf(negloglik) || negloglik > inf_) {
      negloglik = inf_;
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
    int f_pw_idx;
    int f_rw_idx;
    sdouble f_pw;
    sdouble f_rw;
    
    // Compute the boundary distance, for ...
    // ... transition points (enforces tpoint buffer), and 
    // ... Rsum (enforces positive predicted rates)
    for (int r : idx_mc_unique) {
      
      // Grab warping factors
      if (gv_ranLr_int[r] == 0) {
        f_pw = 0.0;
        f_rw = 0.0; 
      } else {
        f_pw_idx = wfactor_idx_point(gv_ranLr_int[r], gv_fixLr_int[r]);
        f_rw_idx = wfactor_idx_rate(gv_ranLr_int[r], gv_fixLr_int[r]);
        f_pw = parameters(f_pw_idx);
        f_rw = parameters(f_rw_idx); 
      } 
      
      // Compute Rt for this row r
      // ... note: Not computing rate for all rows, just ones with unique model components
      // ... this reduces number of rows to compute by a factor of bin_num, which is significant
      sVec Rt = compute_mc_Rt_r(r, parameters, f_rw);
      
      // Grab degree for this row
      int deg = Rt.size() - 1;
      
      // Compute t-point and R_sum boundary distances
      if (deg > 0){
        
        // Compute tslope and tpoint for this row r
        sVec tslope = compute_mc_tslope_r(r, parameters); 
        sVec tpoint = compute_mc_tpoint_r(r, parameters);
        
        // Find tpoint boundary distances
        // WARNING: this code is duplicated in test_tpoint
        sVec tpoint_ext(deg + 2);
        for (int bt = 0; bt < tpoint_ext.size(); bt++) {
          if (bt == 0) {
            tpoint_ext(bt) = 0.0;
          } else if (bt <= deg) {
            tpoint_ext(bt) = point_warp(tpoint(bt - 1), bin_num, f_pw);
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
        
        // Find R_sum boundary distance
        // ... Rates (Rt) must be positive, which happens iff Rt(0) - Rsum > 0.
        // ... Need to use this equivalent condition because we're not computing the 
        //      rate for every row, i.e., not computing for rows which differ only by 
        //      bin number, and this condition is independent of bin number. (Recall that 
        //      the variable "bin_num" is total number of bins, not the number of the bin 
        //      represented by a given row.) 
        sdouble Rsum = 0.0;
        for (int d = 0; d < deg; d++) {
          Rsum += (Rt(d) - Rt(d+1)) / (1.0 + sexp(-tslope(d)*(bin_num - tpoint(d))));
        } 
        sdouble R_sum_boundary_dist = Rt(0) - Rsum; 
        // ... need Rt(0) > Rsum, i.e., this difference should be positive
        boundary_dist_vec(ctr) = R_sum_boundary_dist;
        ctr++;
        
      } else {
        
        // ... trivial to check if rate (Rt) is positive
        boundary_dist_vec(ctr) = Rt(0);
        ctr++;
        
      }
      
    } 
    
    // Compute parameter boundary distance, w factors point
    for (int p : param_wfactor_point_idx) {
      // All warping factors must be between -1 and 1
      sdouble dist_low = parameters(p) + 1.0;
      boundary_dist_vec(ctr) = dist_low; 
      ctr++;
      sdouble dist_high = 1.0 - parameters(p);
      boundary_dist_vec(ctr) = dist_high;
      ctr++; 
    } 
    
    // Compute parameter boundary distance, w factors rate
    for (int p : param_wfactor_rate_idx) {
      // All warping factors must be between -1 and 1
      sdouble dist_low = parameters(p) + 1.0;
      boundary_dist_vec(ctr) = dist_low; 
      ctr++;
      sdouble dist_high = 1.0 - parameters(p);
      boundary_dist_vec(ctr) = dist_high;
      ctr++; 
    } 
    
    // Compute parameter boundary distance, structural parameters
    for (int p : param_struc_idx) {
      // All structural parameters must be > zero.
      sdouble dist_low = parameters(p);
      boundary_dist_vec(ctr) = dist_low;
      ctr++;
    }
    
    return boundary_dist_vec;
    
  }

// Test for tpoints below the buffer
bool wspc::test_tpoints(
    const sVec& parameters
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
        if (gv_ranLr_int[r] == 0) {
          f_pw = 0.0;
        } else {
          f_pw_idx = wfactor_idx_point(gv_ranLr_int[r], gv_fixLr_int[r]);
          f_pw = parameters(f_pw_idx);
        }
        
        // Find degree
        int c_num = Rwhich(eq_left_broadcast(child_lvls, child[r]))[0];
        int p_num = Rwhich(eq_left_broadcast(parent_lvls, parent[r]))[0];
        int deg = degMat(c_num, p_num);
        if (deg > 0) {
          
          // Compute tpoints for this row r
          tpoint = compute_mc_tpoint_r(r, parameters);
          
          // Find tpoint boundary distances
          // WARNING: this code is duplicated from boundary_dist
          sVec tpoint_ext(deg + 2);
          for (int bt = 0; bt < tpoint_ext.size(); bt++) {
            if (bt == 0) {
              tpoint_ext(bt) = 0.0;
            } else if (bt <= deg) {
              tpoint_ext(bt) = point_warp(tpoint(bt - 1), bin_num, f_pw);
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
    
    // Compute neg_loglik
    sdouble bnll = neg_loglik(parameters);
    
    // Compute boundary distance
    sVec bd = boundary_dist(parameters);
    
    // Add boundary penalty
    for (int i = 0; i < bd.size(); i++) {
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
 * Bootstrapping and model fitting, for statistical testing
 */

// Fit model using NLopt
void wspc::fit(
    const bool set_params, 
    const bool verbose
  ) { 
    
    // Set boundary-penalty coefficients 
    sVec initial_params_var = to_sVec(fitted_parameters);
    max_penalty_at_distance = neg_loglik(initial_params_var) * max_penalty_at_distance_factor;
    sdouble coefs_square = static_cast<double>(boundary_vec_size)/max_penalty_at_distance;
    bp_coefs = boundary_dist(initial_params_var);
    for (int i = 0; i < boundary_vec_size - 1; i++) {
      bp_coefs(i) = ssqrt(coefs_square)/bp_coefs(i);
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
      sVec ip = to_sVec(fitted_parameters);
      sdouble initial_nll = neg_loglik(ip);
      vprint("Initial neg_loglik: ", initial_nll);
      sdouble initial_obj = bounded_nll(ip);
      vprint("Initial neg_loglik with penalty: ", initial_obj);
      sdouble mb_dist = min_boundary_dist(ip);
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
    fitted_parameters = to_NumVec(x);
    optim_results["fitted_parameters"] = fitted_parameters; // for predicting log-linked count
    optim_results["penalized_neg_loglik"] = min_fx;
    optim_results["neg_loglik"] = final_nll.val();
    optim_results["success_code"] = success_code;
    optim_results["num_evals"] = opt.get_numevals();
    
  } 

// Fit model to bootstrap resample
dVec wspc::bs_fit(
    int bs_num,                  // A unique number for this resample
    bool clear_stan              // Recover stan memory at end?
  ) {
    
    // Set random number generator with unique seed based on resample ID
    unsigned int seed = rng_seed + bs_num;
    pcg32 rng(seed);
    
    // Resample (with replacement), re-sum token pool, and take logs
    for (int r : count_not_na_idx) {
      
      // ... but only for actual observations of random grouping variable
      if (ran[r] != "none") {
        
        // ... select a new row r2 from the block of r's bin
        IntegerVector back_up_range = bs_bin_ranges.row(r);
        int shift_back_max = back_up_range(0); 
        int shift_up_max = back_up_range(1); 
        int shift_range = back_up_range(2); 
        
        // ... randomly select integer between 0 and shift_range
        int shift_bins = 0;
        if (shift_range > 2) {
          int shift_idx = rng(shift_range); 
          IntegerVector shift = iseq(-shift_back_max, shift_up_max, shift_range);
          shift_bins = shift[shift_idx];
        }
        int r2 = r + shift_bins;
        
        // ... redraw randomly (with replacement) from the token pool of r2 and re-sum into r's count 
        IntegerVector token_pool_r = token_pool[r2];
        int resample_sz = token_pool_r.size();
        if (resample_sz < 1) {
          // ... ensure new point is viable
          token_pool_r = token_pool[r];
          resample_sz = token_pool_r.size();
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
    count = extrapolate_none(count, ran, extrapolation_pool);
    iVec r_rows = Rcpp::as<iVec>(count_row_nums[eq_left_broadcast(ran,"none")]);
    for (int r : r_rows) {
      count_log(r) = slog(count(r) + 1.0);
    }
    
    // Set fitted parameters to jitter of initial seed
    fitted_parameters = jitter_parameter_seed();
    
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
      stan::math::recover_memory();
    }
    
    return fitted_parameters;
    
  }

// Fork bootstraps (parallel processing)
Rcpp::NumericMatrix wspc::bs_batch(
    int bs_num_max,              // Number of bootstraps to perform
    int max_fork,                // Maximum number of forked processes per batch
    bool verbose
  ) {
    
    // Reset initial parameters to check their feasibility 
    set_parameters(fitted_parameters, true);
    
    // Fit full model from multiple jittered seed parameters
    vprint("Performing initial fit of full data", verbose);
    NumericVector fitted_parameters_best = Rcpp::clone(fitted_parameters);
    List optim_results_best = Rcpp::clone(optim_results);
    double least_pnll = inf_.val();
    for (int i = 0; i < initial_fits; i++) {
      
      // Report progress
      if (verbose) {
        Rcpp::Rcout << "Jittered initial position " << i + 1 << "/" << initial_fits << std::endl;
      }
      
      // Set fitted parameters to jitter of initial seed
      fitted_parameters = jitter_parameter_seed();
      
      // Check feasibility 
      List feasibility_results = check_parameter_feasibility(to_sVec(fitted_parameters), false); 
      fitted_parameters = Rcpp::as<NumericVector>(feasibility_results["parameters"]);
      
      // Save this seed
      NumericVector fitted_parameters_seed_this = Rcpp::clone(fitted_parameters);
      
      // Now fit
      fit(
        false, // don't bother setting parameters
        false  // don't print anything
      );
      
      // Grab penalized neg_loglik
      double pnll = optim_results["penalized_neg_loglik"]; 
      
      // Check if this is a better initial position 
      if (pnll < least_pnll) {
        least_pnll = pnll;
        // ... if so, save this seed
        fitted_parameters_seed = fitted_parameters_seed_this;
        // ... and save the fitted parameters
        fitted_parameters_best = fitted_parameters;
        // ... and optim results 
        optim_results_best = optim_results;
      }
      
    }
    // Set fitted parameters to best found 
    
    fitted_parameters = fitted_parameters_best;
    
    optim_results = optim_results_best;
    
    // Initiate variables to hold results
    const int c_num = fitted_parameters.size() + 4;
    const int r_num = bs_num_max + 1;
    NumericMatrix bs_results(r_num, c_num);
    
    // Save results from initial full fit
    dVec full_results = to_dVec(fitted_parameters);
    full_results.reserve(full_results.size() + 4);
    full_results.push_back(Rcpp::as<double>(optim_results["penalized_neg_loglik"]));
    full_results.push_back(Rcpp::as<double>(optim_results["neg_loglik"]));
    full_results.push_back(Rcpp::as<double>(optim_results["success_code"]));
    full_results.push_back(Rcpp::as<double>(optim_results["num_evals"]));
    
    // Save full fit results in last row of results matrix
    bs_results.row(bs_num_max) = to_NumVec(full_results);
    
    // Perform bootstrap fits in batches
    const int batch_num = std::round(bs_num_max / max_fork);
    for (int b = 0; b < batch_num; b++) {
      // Run in parallel with forking
      
      // Initiate timer and grab start time
      Timer batch_timer;
      batch_timer.step("start");  
      
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
      
      batch_timer.step("end");
      NumericVector batch_times(batch_timer);
      double batch_time = 1e-9 * (batch_times[1] - batch_times[0])/max_fork;
      if (verbose) {
        Rcpp::Rcout << "Batch " << b + 1 << "/" << batch_num << ", " << batch_time << " sec/bs" << std::endl;
      }
      
    }
    
    vprint("All complete!", verbose); 
    
    return bs_results;
    
  }

// Resample (demo, not used in package)
NumericMatrix wspc::resample(
    int n_resample               // total number of resamples to draw
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
      count_new = extrapolate_none(count_new, ran, extrapolation_pool);
      
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
     *  - Replaces structural parameters
     */
    
    // Ensure provided parameters are feasible (i.e., don't predict negative rates) 
    List feasibility_results = check_parameter_feasibility(to_sVec(parameters), verbose);
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
    
    // Update struc_values 
    for (int i = 0; i < param_struc_idx.size(); i++) {
      int p = param_struc_idx(i);
      struc_values[i] = fitted_parameters(p);
    }
    
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
    bool feasible = test_tpoints(parameters_var);
    if (verbose) {
      if (feasible) {
        vprint("... no tpoints below buffer");
      } else {
        vprint("... tpoints found below buffer");
      }
    }
    if (feasible) {
      
      // Predict rates 
      predicted_rates_log_var = predict_rates_log(
        parameters_var, // model parameters for generating prediction
        true            // compute all summed count rows, even with a count value of NA?
      );
      
      // Test if provided parameters produce any negative rates
      for (int i = 0; i < n_count_rows; i++) {
        if (predicted_rates_log_var(i) < 0 || std::isnan(predicted_rates_log_var(i))) {
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
    
    // If not feasible, attempt to find a feasible starting point
    if (!feasible) {
      
      // Find and report initial distance to boundary
      sdouble initial_dist = min_boundary_dist(parameters_var); 
      if (verbose) {
        vprint("Initial boundary distance (want to make >0): ", initial_dist);
      }
      
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
        predicted_rates_log_var = predict_rates_log(
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

// Jitter parameter seed, for use before fitting 
Rcpp::NumericVector wspc::jitter_parameter_seed() const {
    
    // Grab seed parameters
    NumericVector jittered_parameters = fitted_parameters_seed;
    
    // Jitter parameters
    double wf_cap = 0.75;        // Max magnitude for wf ... keep away from bounds
    double wf_sd = 0.2;          // Standard deviation for wf jitter
    double beta_Rt_sd = 0.1;     // Standard deviation for beta_Rt jitter ... be very conservative
    double beta_tpoint_sd = 0.5; // Standard deviation for beta_tpoint jitter ... be very conservative
    double beta_tslope_sd = 0.1; // Standard deviation for beta_tslope jitter ... be very conservative
    double struc_sd = 1.0;       // Standard deviation for structural parameter jitter 
    
    // Jitter warping factors, point
    for (int p : param_wfactor_point_idx) {
      double wfpj = safe_rnorm(jittered_parameters(p), wf_sd);
      if (wfpj < -wf_cap) {wfpj = -wf_cap;}
      if (wfpj > wf_cap) {wfpj = wf_cap;}
      jittered_parameters(p) = wfpj;
    }
    
    // Jitter warping factors, rate 
    for (int p : param_wfactor_rate_idx) {
      double wfrj = safe_rnorm(jittered_parameters(p), wf_sd);
      if (wfrj < -wf_cap) {wfrj = -wf_cap;}
      if (wfrj > wf_cap) {wfrj = wf_cap;}
      jittered_parameters(p) = wfrj;
    }
    
    // Jitter beta Rt parameters
    for (int p : param_beta_Rt_idx) {
      // might be sent into unfeasible position, but should be correctable
      double betaj = safe_rnorm(jittered_parameters(p), beta_Rt_sd);
      jittered_parameters(p) = betaj;
    }
    
    // Jitter beta tpoint parameters
    for (int p : param_beta_tpoint_idx) {
      // might be sent into unfeasible position, but should be correctable
      double betaj = safe_rnorm(jittered_parameters(p), beta_tpoint_sd);
      jittered_parameters(p) = betaj;
    }
    
    // Jitter beta tslope parameters
    for (int p : param_beta_tslope_idx) {
      // might be sent into unfeasible position, but should be correctable
      double betaj = safe_rnorm(jittered_parameters(p), beta_tslope_sd);
      jittered_parameters(p) = betaj;
    }
    
    // Jitter structural parameters
    for (int p : param_struc_idx) {
      double jit = safe_rnorm(1.0, struc_sd);
      double strucj = jittered_parameters(p) + jit * jit;
      jittered_parameters(p) = strucj;
    }
   
    return jittered_parameters;
   
  }

/*
 * *************************************************************************
 * Export/import data to/from R
 */

void wspc::import_fe_diff_ratio_Rt(
    const double& fe_diff_ratio,
    const bool& verbose
  ) {
    vprint("Imported fe_difference_ratio_Rt: " + std::to_string(fe_diff_ratio), verbose);
    fe_difference_ratio_Rt = (sdouble)fe_diff_ratio;
  }

void wspc::import_fe_diff_ratio_tpoint(
    const double& fe_diff_ratio,
    const bool& verbose
  ) {
    vprint("Imported fe_difference_ratio_tpoint: " + std::to_string(fe_diff_ratio), verbose);
    fe_difference_ratio_tpoint = (sdouble)fe_diff_ratio;
  }

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
    
    // Reformat parameter names for R
    int n_params = param_names.size();
    CharacterVector param_names_clean(n_params);
    for (int n = 0; n < n_params; n++) {
      CharacterVector name_comps = Rcpp::as<CharacterVector>(param_names[n]);
      int m = name_comps.size();
      if (m == 0) {Rcpp::stop("Empty parameter name!");}
      String name = name_comps[0]; 
      if (m > 1) {
        for (int i = 1; i < m; i++) {
          name += "_" + name_comps[i];
        }
      }
      param_names_clean[n] = name;
    }
    
    // Add parameter names to fitted parameter vector 
    fitted_parameters.names() = param_names_clean;
    
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
    
    // Add names to structure parameters
    struc_values.names() = struc_names;
    
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
    
    // Recompute sd_Rt_effect
    sdouble expected_ran_effect_rate = 1.0 / ssqrt(4.0 * struc_values["beta_shape_rate"] + 2.0);    // ... already includes the *= 2.0 rescaling to range from -1 to 1
    sdouble eff_mult_rate = fe_difference_ratio_Rt - 1.0;
    sdouble expected_Rt = mean_count_log;
    sdouble sd_Rt_effect = eff_mult_rate * expected_Rt * expected_ran_effect_rate;
    
    // Recompute sd_tpoint_effect
    sdouble expected_ran_effect_tpoint = 1.0 / ssqrt(4.0 * struc_values["beta_shape_point"] + 2.0); // ... already includes the *= 2.0 rescaling to range from -1 to 1
    sdouble eff_mult_tpoint = fe_difference_ratio_tpoint + 1.0;                                     // ... note the different sign
    sdouble expected_tpoint = bin_num / 2.0;
    sdouble sd_tpoint_effect = (eff_mult_tpoint * expected_tpoint * expected_ran_effect_tpoint) / (2.0 + eff_mult_tpoint * expected_ran_effect_tpoint);
    
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
      _["param.names"] = param_names_clean,
      _["fix"] = fixed_effects,
      _["treatment"] = treat,
      _["grouping.variables"] = grouping_variables,
      _["struc.params"] = struc_values,
      _["param.idx0"] = param_idx, // "0" to indicate this goes out w/ C++ zero-based indexing
      _["token.pool"] = token_pool_list,
      _["computed_sd_Rt_effect"] = sd_Rt_effect.val(),
      _["computed_sd_tpoint_effect"] = sd_tpoint_effect.val(),
      _["change.points"] = change_points,
      _["settings"] = model_settings
    );
    
    return results_list;
    
  }

/*
 * *************************************************************************
 * Testing and debugging in R
 */

// Wrap neg_loglik in form needed for R
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
    sdouble fx = bounded_nll(parameters_var);
    
    // Return the value of the neg_loglik
    return fx.val(); 
    
  } 

// Wrap neg_loglik in form needed for R
Rcpp::NumericVector wspc::grad_bounded_nll_debug(
    const dVec& x 
  ) { 
    
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
    
    // Return the value of the neg_loglik
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
    .method("bounded_nll_debug", &wspc::bounded_nll_debug)
    .method("grad_bounded_nll_debug", &wspc::grad_bounded_nll_debug)
    .method("fit", &wspc::fit)
    .method("bs_fit", &wspc::bs_fit)
    .method("bs_batch", &wspc::bs_batch)
    .method("resample", &wspc::resample)
    .method("import_fe_diff_ratio_Rt", &wspc::import_fe_diff_ratio_Rt)
    .method("import_fe_diff_ratio_tpoint", &wspc::import_fe_diff_ratio_tpoint)
    .method("results", &wspc::results);
  }
