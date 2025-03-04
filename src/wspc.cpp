
// wspc.cpp
#include "wspc.h"

// WSPmm class and methods *********************************************************************************************

/*
 * Object class to hold and fit Warped Sigmoidal Poisson-Process Mixed-Effect Model (WSPmm) model. 
 */

// Constructor
wspc::wspc(
    Rcpp::DataFrame count_data,
    Rcpp::List settings,
    bool verbose
  ) { 
    
    /*
     * Input data should be an R data frame with the columns specified below, and rows 
     *  as token count observations. 
     *  
     * Explanation of the difference between "summed" and "tokenized" count data:
     * 
     * Consider the following case. We have multiple parent levels, e.g. different cell types (glia, 
     *  excitatory neurons, inhibitory neurons, etc). In many cases, each bin will have multiple instances (tokens) 
     *  of that parent level. For example, a bin may have 3 glia cells, 5 excitatory neurons, and 2 inhibitory neurons.
     *  In this case, we have a decision to make: Either collapse those token together and put the total count (of whatever
     *  is being counted, e.g., gene transcripts) for a given bin in a single row (in the rate / count data frames), or 
     *  have individual rows for each token. Conceptually, one way to look at it is that in the rate data frame, 
     *  every row holds a unique combination of independent variables having an effect on expression rate (and so exhaust all
     *  we need to know for prediction), whereas in the count data frame, every row is a unique observation. Multiple rows 
     *  can have the same independent variables values. Practically speaking, only the rate data frame is relevant for 
     *  predicting rates, but the count data frame is the one relevant to computing measures of model fit, like likelihood. 
     */
    
    // Extract and save settings
    NumericVector struc_values_settings = settings["struc_values"];
    double buffer_factor_settings = settings["buffer_factor"];
    double ctol_settings = settings["ctol"];
    double max_penalty_at_distance_factor_settings = settings["max_penalty_at_distance_factor"];
    double LROcutoff_settings = settings["LROcutoff"];
    double tslope_initial_settings = settings["tslope_initial"];
    double wf_initial_settings = settings["wf_initial"];
    int max_evals_settings = settings["max_evals"];
    for (int i = 0; i < struc_values.size(); i++) {struc_values[i] = struc_values_settings[i];}
    buffer_factor = sdouble(buffer_factor_settings);
    ctol = ctol_settings;
    max_penalty_at_distance_factor = sdouble(max_penalty_at_distance_factor_settings);
    LROcutoff = LROcutoff_settings;
    tslope_initial = tslope_initial_settings;
    wf_initial = wf_initial_settings;
    max_evals = max_evals_settings;
    model_settings = Rcpp::clone(settings);
    
    // Check structure of input data
    CharacterVector col_names = count_data.names();
    CharacterVector required_cols = CharacterVector::create("count", "bin", "parent", "child", "ran");
    int n_cols = col_names.size();
    int r_cols = required_cols.size();
    for (int i = 0; i < r_cols; i++) {
      if (col_names[i] != required_cols[i]) {
        stop("Input data is missing required column (or out of order): " + required_cols[i]);
      }
    }
    if (tslope_initial < 1.0) {
      stop("tslope_initial must be greater than or equal to 1.0.");
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
    fix_names = CharacterVector(n_fix);  // names of fixed effect variables 
    fix_ref = CharacterVector(n_fix);    // reference level for each fixed effect
    fix_lvls.resize(n_fix);              // levels for each fixed effect
    fix_trt.resize(n_fix);               // treatments for each fixed effect
    for (int i = 0; i < n_fix; i++) {
      fix_names[i] = col_names[i + r_cols];
      CharacterVector lvls = Rcpp::sort_unique(Rcpp::as<CharacterVector>(count_data[i + r_cols]));
      if (lvls.size() < 2) {
        stop("Fixed effect " + fix_names[i] + " has less than 2 levels.");
      }
      fix_lvls[i] = lvls;
      fix_ref[i] = lvls[0];   // assume first level is reference
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
          // ref level must always have weight 1
        }
      }
    }
    vprint("Pre-computed weight matrix rows", verbose);
    
    // Extract grouping variables 
    parent_lvls = Rcpp::sort_unique(Rcpp::as<CharacterVector>(count_data["parent"]));
    child_lvls = Rcpp::sort_unique(Rcpp::as<CharacterVector>(count_data["child"]));
    ran_lvls = Rcpp::sort_unique(Rcpp::as<CharacterVector>(count_data["ran"]));
    ran_lvls.push_front("none"); // add "none" to represent no random effect (reference level)
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
    
    // Create summed count data, initializations
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
                token_pool[idx] = token_pool_idx; // save for bootstrap resampling
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
    
    // Estimate degree of each parent-child combination at baseline using LRO change-point detection 
    // ... store ref values in list
    List ref_values(n_parent);
    ref_values.names() = parent_lvls;
    // ... store rate effects in list
    List RtEffs(n_parent);
    RtEffs.names() = parent_lvls;
    // ... store change points in list
    List found_cp_list(n_parent);
    found_cp_list.names() = parent_lvls;
    for (int p = 0; p < n_parent; p++) {
      
      LogicalVector parent_mask = eq_left_broadcast(parent, parent_lvls[p]);
      
      // Set up list for ref values
      ref_values[(String)parent_lvls[p]] = List(n_child);
      name_proxylist(ref_values[(String)parent_lvls[p]], child_lvls);
      // ... for rate effect values
      RtEffs[(String)parent_lvls[p]] = List(n_child);
      name_proxylist(RtEffs[(String)parent_lvls[p]], child_lvls);
      // ... for change points 
      found_cp_list[(String)parent_lvls[p]] = List(n_child);
      name_proxylist(found_cp_list[(String)parent_lvls[p]], child_lvls);
      for (int c = 0; c < n_child; c++) {
        
        // Run LRO algorithm for each treatment level of this child-parent pair
        NumericMatrix count_avg_mat(bin_num_i, treatment_num);
        IntegerVector found_cp(0);
        for (int t = 0; t < treatment_num; t++) {
          String trt = treatment_lvls[t];
         
          // Make mask for treatment rows of this parent-child pair
          LogicalVector trt_mask = eq_left_broadcast(treatment, trt);
          LogicalVector mask = count_not_na_mask & parent_mask & trt_mask & eq_left_broadcast(child, child_lvls[c]);
         
          // Make masked copies of count_log and bin
          sVec count_masked = masked_vec(count_log, mask);  
          sVec bin_masked = masked_vec(bin, mask);
          dVec count_avg(bin_num_i); 
         
          // Compute average count for each bin
          for (int b = 0; b < bin_num_i; b++) {
            LogicalVector mask_b = eq_left_broadcast(to_NumVec(bin_masked), (double)b + 1.0);
            sVec count_b = masked_vec(count_masked, mask_b);
            count_avg[b] = vmean(count_b).val();
          }
         
          // Estimate change points from averaged count series
          int ws = static_cast<int>(std::round(2.0 * (double)bin_num_i * buffer_factor.val()));
          int filter_ws = std::round(ws/2);
          IntegerVector found_cp_trt = LROcp(
            count_avg,  // series to scan
            ws,         // running window size 
            filter_ws,  // size of window for taking rolling mean
            LROcutoff   // points more than this times sd considered outliers
          );
          
          // Merge with previously found change points 
          if (found_cp.size() == 0) {
            // ... if no cp found yet, replace with current treatment cp
            found_cp = found_cp_trt;
          } else if (found_cp_trt.size() > 0) {
            // ... if cp have been found and current treatment has cp, merge
            buffered_merge(found_cp, found_cp_trt, ws/2);
          }
         
          // Save averaged count values
          count_avg_mat.column(t) = to_NumVec(count_avg);
         
        }
        
        // Extract deg and blocks
        int deg = found_cp.size();
        degMat(c, p) = deg;
        int n_blocks = deg + 1;
        assign_proxylist(found_cp_list[(String)parent_lvls[p]], (String)child_lvls[c], to_NumVec(found_cp));
       
        // Set initial parameters for fixed-effect treatments
        List placeholder_ref_values(mc_list.size());
        placeholder_ref_values.names() = mc_list;
        NumericMatrix placeholder_RtEffs(treatment_num, n_blocks);
        for (String mc : mc_list) {
          
          // Mean rate per block
          if (mc == "Rt") {
            dVec Rt_est(n_blocks); 
            sMat RtVals(treatment_num, n_blocks);
           
            for (int t = 0; t < treatment_num; t++) {
             
              NumericVector count_avg = count_avg_mat.column(t);
              
              if (n_blocks == 1) { // case when deg = 0
                Rt_est[0] = vmean(count_avg);
                RtVals(t, 0) = (sdouble)Rt_est[0];
              } else {
                
                for (int bk = 0; bk < n_blocks; bk++) {
                  int bin_start;
                  int bin_end;
                  if (bk == 0) {
                    bin_start = 0;
                    bin_end = found_cp[bk];
                  } else if (bk == deg) {
                    bin_start = found_cp[bk - 1] - 1;
                    bin_end = bin_num_i - 1;
                  } else {
                    bin_start = found_cp[bk - 1] - 1;
                    bin_end = found_cp[bk];
                  }
                  Rt_est[bk] = vmean_range(count_avg, bin_start, bin_end);
                  RtVals(t, bk) = (sdouble)Rt_est[bk];
                }
                
              }
              if (t == 0) {placeholder_ref_values[mc] = to_NumVec(Rt_est);}
             
            }
           
            // Estimate fixed effects from treatments for each block
            for (int bk = 0; bk < n_blocks; bk++) {
              sVec beta_bk = weight_rows.fullPivLu().solve(RtVals.col(bk));
              placeholder_RtEffs.column(bk) = to_NumVec(beta_bk);
            }
            
          } else if (deg > 0) { // Handle tpoint and tslope
            if (mc == "tpoint") {
              placeholder_ref_values[mc] = Rcpp::as<NumericVector>(found_cp);
            } else if (mc == "tslope") {
              NumericVector tslope_default(deg);
              // ... but tslope into log space
              for (int tp = 0; tp < deg; tp++) {tslope_default[tp] = tslope_initial;}
              placeholder_ref_values[mc] = tslope_default;
            }
          }
         
        }
        assign_proxylist(ref_values[(String)parent_lvls[p]], (String)child_lvls[c], placeholder_ref_values);
        assign_proxylist(RtEffs[(String)parent_lvls[p]], (String)child_lvls[c], placeholder_RtEffs);
       
      }
    }
    // ... save found change points 
    change_points = found_cp_list;
    vprint("Estimated change points with LROcp and found initial parameter estimates for fixed-effect treatments", verbose); 
   
    // Build default fixed-effects matrices in shell
    List beta = build_beta_shell(mc_list, treatment_lvls, parent_lvls, child_lvls, ref_values, RtEffs, degMat);
    vprint("Built initial beta (ref and fixed-effects) matrices", verbose); 
    
    // Initialize random effect warping factors 
    List wfactors = List(2);
    CharacterVector wfactors_names = {"point", "rate"};
    wfactors.names() = wfactors_names;
    double wf_initial_low = -1.0 * wf_initial;
    NumericVector wfs = dseq(wf_initial_low, wf_initial, n_ran - 1); 
    for (String wf : wfactors_names) {
      NumericMatrix wf_array(n_ran, n_child);
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
          wfs_c = wfs_c * Rcpp::runif(n_ran - 1, 0.5, 1.0); // add noise
          wfs_c.push_front(0.0); // add "none"
          wf_array.column(c) = wfs_c;
        } else { 
          NumericVector wfs_c = Rcpp::runif(n_ran - 1, 0.0, 0.01); 
          wfs_c.push_front(0.0); // add "none"
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
   
    // Construct grouping variable ids as indexes for warping factor matrices, by count row
    gv_ranLr_int = IntegerVector(n_count_rows);
    gv_fixLr_int = IntegerVector(n_count_rows);
    for (int r = 0; r < n_count_rows; r++) {
      CharacterVector::iterator it_ran = std::find(ran_lvls.begin(), ran_lvls.end(), ran[r]);
      CharacterVector::iterator it_fix = std::find(child_lvls.begin(), child_lvls.end(), child[r]);
      gv_ranLr_int[r] = std::distance(ran_lvls.begin(), it_ran);
      gv_fixLr_int[r] = std::distance(child_lvls.begin(), it_fix);
    }
    vprint("Constructed grouping variable ids", verbose);
    
    // Compute tpoint buffer
    tpoint_buffer = bin_num * buffer_factor;
    
    // Compute size of the parameter boundary vector 
    for (int r : idx_mc_unique) {
      // Grab degree for this row
      int c_num = Rwhich(eq_left_broadcast(child_lvls, child[r]))[0];
      int p_num = Rwhich(eq_left_broadcast(parent_lvls, parent[r]))[0];
      int deg = degMat(c_num, p_num);
      // Compute t-point and R_sum boundary distances
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
    vprint("Computed size of rRsum boundary vector: " + std::to_string(boundary_vec_size), verbose);
    
    // Pre-compute bin ranges for bootstrap resampling 
    bs_bin_ranges = IntegerMatrix(n_count_rows, 3);
    for (int r : count_not_na_idx) {
      if (ran[r] != "none") {
        // ... find the change points for the parent-child pair of this row
        List p_changepoints = change_points[(String)parent[r]];
        NumericVector pc_changepoints = p_changepoints[(String)child[r]];
        iVec block_range_idx = block_idx(pc_changepoints, bin(r).val());
        // ... find the range of bins in the same block as this row's bin
        double bin_low = 1;
        if (block_range_idx[0] >= 0) {
          bin_low = pc_changepoints[block_range_idx[0]];
        }
        double bin_high = bin_num.val();
        if (block_range_idx[1] >= 0) {
          bin_high = pc_changepoints[block_range_idx[1]];
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
    const bool& all_rows // compute all summed count rows, even with a count value of NA?
  ) const {
    
    /*
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
   
    // Initialize variable to hold nll
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
    log_lik += log_dnorm(pw_mean_p, 0.0, sd_pw_p);
    
    // Compute the log-likelihood of the mean of rate warping factors, given the expected normal distribution
    // ... see above notes on the formula
    sdouble pw_mean_r = vmean(warping_factors_rate); 
    sdouble sd_pw_r = ssqrt(1.0 / ((4.0 * (sdouble)warping_factors_rate.size()) * (2.0 * beta_shape_rate + 1.0))) * 2.0;
    log_lik += log_dnorm(pw_mean_r, 0.0, sd_pw_r);
    
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
      log_lik += log_dnorm(observed_mean_ran_eff[r - 1], modeled_mean, modeled_sd);
    }
    
    // *****************************************************************************************************************
    // log-likelihood of beta values
    
    // Compute the log-likelihood of the Rt beta values, given the normal distribution implied by the beta-rate shape and fe_difference_ratio ratio
    sdouble expected_ran_effect = ssqrt(1.0 / (4.0 * (2.0 * beta_shape_rate + 1.0)));
    expected_ran_effect *= 2.0; // scale to range from -1 to 1
    sdouble eff_mult = fe_difference_ratio - 1.0;
    sdouble sd_Rt_effect = eff_mult * mean_count_log * expected_ran_effect;
    for (int i = 0; i < Rt_beta_values_no_ref.size(); i++) {
      log_lik += log_dnorm(Rt_beta_values_no_ref(i), 0.0, sd_Rt_effect);
    }
    
    // Compute the log-likelihood of the tslope beta values, given the assumed normal distribution
    int n_tslope = tslope_beta_values_no_ref.size();
    sdouble sd_tslope_effect = parameters[param_struc_idx["sd_tslope_effect"]];
    if (n_tslope != 0) {
      for (int i = 0; i < n_tslope; i++) {
        log_lik += log_dnorm(tslope_beta_values_no_ref(i), 0.0, sexp(sd_tslope_effect));
      }
    }
    
    // Compute the log-likelihood of the tpoint beta values, given the assumed normal distribution
    // TO DO: REWRITE SO THAT IT'S LIKE SD_RT_EFFECT, DERIVED FROM EXPECTED_RAN_EFFECT ON TPOINT
    int n_tpoint = tpoint_beta_values_no_ref.size();
    sdouble sd_tpoint_effect = parameters[param_struc_idx["sd_tpoint_effect"]];
    if (n_tpoint != 0) {
      for (int i = 0; i < n_tpoint; i++) {
        log_lik += log_dnorm(tpoint_beta_values_no_ref(i), 0.0, sd_tpoint_effect);
      }
    }
    
    // Check for zero likelihood
    if (std::isinf(log_lik) || std::isnan(log_lik)) {
      // A nan most likely means sd_tpoint_effect when negative
      return inf_;
    }
    
    // *****************************************************************************************************************
    // log-likelihood of observed counts, given predicted rates
    
    // Predict rates under these parameters
    sVec predicted_rates_log = predict_rates_log(
      parameters, // model parameters for generating prediction
      false       // compute all summed count rows, even with a count value of NA?
      );
    
    // Compute the log-likelihood of the count data
    for (int r : count_not_na_idx) {
      
      if (std::isinf(predicted_rates_log(r)) || predicted_rates_log(r) < 0.0 || std::isnan(predicted_rates_log(r))) {
        return inf_;
      } else {
        log_lik += count_log(r) * slog(predicted_rates_log(r)) - predicted_rates_log(r) - stan::math::lgamma(count_log(r) + 1);
        // ^ Hand-written function for log Poisson density
        //  ... could use stan implementation: log_lik += stan::math::poisson_lpmf(count_log(r), pred_rate_log);
        //  ... but seems identical in results and speed?
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
sVec wspc::neg_boundary_dist(
    const sVec& parameters
  ) const {
    
    // Initialize vector to hold negative boundary distances
    sVec neg_boundary_dist_vec(boundary_vec_size);
    int ctr = 0;
    
    // Grab warping indices and initiate variables to hold them
    NumericMatrix wfactor_idx_point = wfactor_idx["point"];
    NumericMatrix wfactor_idx_rate = wfactor_idx["rate"];
    int f_pw_idx;
    int f_rw_idx;
    sdouble f_pw;
    sdouble f_rw;
    
    // Compute the boundary distance
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
          neg_boundary_dist_vec(ctr) = -buffer_dist;
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
          // ... optimization wants to inflate Rt(deg), as that lowers Rsum; the same effect does not 
          //      happen for d <= deg, as that would also increase Rsum. 
        } 
        sdouble R_sum_boundary_dist = Rt(0) - Rsum; 
        // ... need Rt(0) > Rsum, i.e., this difference should be positive
        neg_boundary_dist_vec(ctr) = -R_sum_boundary_dist;
        ctr++;
        
      } else {
        
        // ... trivial to check if rate (Rt) is positive
        neg_boundary_dist_vec(ctr) = -Rt(0);
        ctr++;
        
      }
      
    } 
    
    // Compute parameter boundary distance, w factors point
    for (int p : param_wfactor_point_idx) {
      // All warping factors must be between -1 and 1
      sdouble dist_low = parameters(p) + 1.0;
      neg_boundary_dist_vec(ctr) = -dist_low; 
      ctr++;
      sdouble dist_high = 1.0 - parameters(p);
      neg_boundary_dist_vec(ctr) = -dist_high;
      ctr++; 
    } 
    
    // Compute parameter boundary distance, w factors rate
    for (int p : param_wfactor_rate_idx) {
      // All warping factors must be between -1 and 1
      sdouble dist_low = parameters(p) + 1.0;
      neg_boundary_dist_vec(ctr) = -dist_low; 
      ctr++;
      sdouble dist_high = 1.0 - parameters(p);
      neg_boundary_dist_vec(ctr) = -dist_high;
      ctr++; 
    } 
    
    // Compute parameter boundary distance, structural parameters
    for (int p : param_struc_idx) {
      // All structural parameters must be > zero.
      sdouble dist_low = parameters(p);
      neg_boundary_dist_vec(ctr) = -dist_low;
      ctr++;
    }
    
    return neg_boundary_dist_vec;
    
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
      
      // Grab warping factors
      if (gv_ranLr_int[r] == 0) {
        f_pw = 0.0;
      } else {
        f_pw_idx = wfactor_idx_point(gv_ranLr_int[r], gv_fixLr_int[r]);
        f_pw = parameters(f_pw_idx);
      }
      
      // Only update predicted model components if r begins a new batch of unique values 
      int cnt = std::count(mc_unique_rows.begin(), mc_unique_rows.end(), r);
      if (cnt > 0) { 
        // Find degree
        List beta_idx_Rt = Rcpp::as<List>(beta_idx["Rt"]);
        List beta_idx_Rt_prt = Rcpp::as<List>(beta_idx_Rt[Rcpp::as<std::string>(parent[r])]);
        IntegerVector betas_Rt_idx = beta_idx_Rt_prt[Rcpp::as<std::string>(child[r])];
        int c_num = Rwhich(eq_left_broadcast(child_lvls, child[r]))[0];
        int p_num = Rwhich(eq_left_broadcast(parent_lvls, parent[r]))[0];
        int deg = degMat(c_num, p_num);
        if (deg > 0) {
          // Compute tpoints for this row r
          tpoint = compute_mc_tpoint_r(r, parameters);
          // Check for negative tpoints 
          for (int p = 0; p < tpoint.size(); p++) {
            if (point_warp(tpoint(p), bin_num, f_pw) < tpoint_buffer) {
              return false;
            }
          }
        }
      }  
      
    }  
    
    return true;
    
  } 

// Compute min boundary penalty
sdouble wspc::max_neg_boundary_dist(
    const sVec& parameters
  ) const {
   
    // Compute neg_tRsum_boundary_dist and take max
    sVec bd = neg_boundary_dist(parameters);
    sdouble bd_max = smax(bd);
    
    return -bd_max;
    
  }

// Wrap neg_min_boundary_dist in form needed for NLopt constraint function
double wspc::max_neg_boundary_dist_NLopt(
    const dVec& x,
    dVec& grad,
    void* data
  ) { 
    
    // Grab model
    wspc* model = static_cast<wspc*>(data);
    
    // Convert dVec to Eigen with stan
    sVec parameters_var = to_sVec(x);
    
    // Compute max_neg_boundary_dist
    sdouble fx = model->max_neg_boundary_dist(parameters_var);
    
    // Compute gradient if needed
    if (!grad.empty()) {
      Eigen::VectorXd grad_eigen = model->grad_max_neg_boundary_dist(parameters_var);
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
    sVec nbd = neg_boundary_dist(parameters);
    
    // Add boundary penalty
    for (int i = 0; i < nbd.size(); i++) {
      bnll += boundary_penalty_transform(nbd(i), bp_coefs(i));
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

// Compute the gradient of the max_neg_boundary_dist function
// ... this is the gradient function used in the search for feasible parameters
Eigen::VectorXd wspc::grad_max_neg_boundary_dist(
    const sVec& p_
  ) const { 
    
    // Create nested autodiff context
    stan::math::nested_rev_autodiff nested;
    
    // Make copy
    sVec p = p_;
    
    // Initialize max_neg_boundary_dist variable
    sdouble mnbd = max_neg_boundary_dist(p);
    
    // Initialize variable to hold gradient
    Eigen::VectorXd gr_mnbd(p.size());
    
    // Compute max_neg_boundary_dist and its gradient
    stan::math::grad(mnbd, p, gr_mnbd);
    
    // Return max_neg_boundary_dist gradient
    return gr_mnbd;
    
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
    bp_coefs = neg_boundary_dist(initial_params_var);
    for (int i = 0; i < boundary_vec_size - 1; i++) {
      bp_coefs(i) = -1.0 * ssqrt(coefs_square)/bp_coefs(i);
    } 
    // ^ ... * -1.0 because starts as negative of boundary distance, and we want the actual distance
    
    // Grab initial parameters
    dVec x = to_dVec(fitted_parameters); 
    size_t n = x.size();
    
    // Set up NLopt optimizer
    nlopt::opt opt(nlopt::LD_LBFGS, n);
    opt.set_min_objective(wspc::bounded_nll_NLopt, this);
    opt.set_ftol_rel(ctol);       // stop when iteration changes objective fn value by less than this fraction 
    opt.set_maxeval(max_evals);   // Maximum number of evaluations to try
    
    // Find and print initial neg_loglik and total objective values
    if (verbose) {
      sVec ip = to_sVec(fitted_parameters);
      sdouble initial_nll = neg_loglik(ip);
      Rcpp::Rcout << "Initial neg_loglik: " << initial_nll << std::endl;
      sdouble initial_obj = bounded_nll(ip);
      Rcpp::Rcout << "Initial neg_loglik with penalty: " << initial_obj << std::endl;
      sdouble ff_dist = -1 * max_neg_boundary_dist(ip);
      Rcpp::Rcout << "Initial min boundary distance: " << ff_dist << std::endl;
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
      Rcpp::Rcout << "Final neg_loglik: " << final_nll << std::endl;
      Rcpp::Rcout << "Final neg_loglik with penalty: " << min_fx << std::endl;
      sdouble ff_dist = -1 * max_neg_boundary_dist(parameters_final);
      Rcpp::Rcout << "Final min boundary distance: " << ff_dist.val() << std::endl;
      int num_evals = opt.get_numevals();
      Rcpp::Rcout << "Number of evaluations: " << num_evals << std::endl;
      if (success_code == 0) {
        Rcpp::Rcout << "Warning: optimization did not converge." << std::endl;
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
    int bs_num,                 // A unique number for this resample
    bool clear_stan             // Recover stan memory at end?
  ) {
    
    // Set random number generator with unique seed based on resample ID
    unsigned int seed = 42u + bs_num;
    pcg32 rng(seed);
    
    // Resample (with replacement), re-sum token pool, and take logs
    for (int r : count_not_na_idx) {
      if (ran[r] != "none") {
        // ... select a new row r2 from the block of r's bin
        IntegerVector back_up_range = bs_bin_ranges.row(r);
        int shift_back_max = back_up_range(0); 
        int shift_up_max = back_up_range(1); 
        int shift_range = back_up_range(2); 
        // ... randomly select integer between 0 and shift_range
        int shift_idx = rng(shift_range); 
        IntegerVector shift = iseq(-shift_back_max, shift_up_max, shift_range);
        int r2 = r + shift[shift_idx];
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
    
    // Extrapolate none's and take logs
    count = extrapolate_none(count, ran, extrapolation_pool);
    iVec r_rows = Rcpp::as<iVec>(count_row_nums[eq_left_broadcast(ran,"none")]);
    for (int r : r_rows) {
      count_log(r) = slog(count(r) + 1.0);
    }
    
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
    
    // Fit full model
    fit(
      false, // don't bother setting parameters
      false  // don't print anything
    );
    
    // Initiate variables to hold results
    const int c_num = fitted_parameters.size() + 4;
    const int r_num = bs_num_max + 1;
    NumericMatrix bs_results(r_num, c_num);
    
    dVec full_results = to_dVec(fitted_parameters);
    full_results.reserve(full_results.size() + 4);
    full_results.push_back(Rcpp::as<double>(optim_results["penalized_neg_loglik"]));
    full_results.push_back(Rcpp::as<double>(optim_results["neg_loglik"]));
    full_results.push_back(Rcpp::as<double>(optim_results["success_code"]));
    full_results.push_back(Rcpp::as<double>(optim_results["num_evals"]));
    
    // Save full fit results in last row of results matrix
    bs_results.row(bs_num_max) = to_NumVec(full_results);
    
    // Perform bootstrap fits in batches
    const int batch_num = bs_num_max / max_fork;
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
        Rcpp::Rcout << "Batch " << b << "/" << batch_num << ", " << batch_time << " sec/bs" << std::endl;
      }
      
    }
    
    if (verbose) {
      Rcpp::Rcout << "All complete!" << std::endl;
    }
    
    return bs_results;
    
  }

// Resample (demo)
NumericMatrix wspc::resample(
    int n_resample                 // total number of resamples to draw
  ) {
    
    NumericMatrix resamples(n_count_rows, n_resample);
    IntegerVector NA_idx = Rwhich(!count_not_na_mask);
    
    for (int j = 0; j < n_resample; j++) {
      
      // Set random number generator with unique seed based on resample ID
      unsigned int seed = 42u + j;
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
    
    // Convert parameters to math stan 
    sVec parameters_var = to_sVec(parameters);
    
    // Ensure provided parameters are feasible (i.e., don't predict negative rates) 
    List feasibility_results = check_parameter_feasibility(parameters_var, verbose);
    parameters_var = to_sVec(Rcpp::as<NumericVector>(feasibility_results["parameters_var"]));
    sVec predicted_rates_log_var = to_sVec(Rcpp::as<NumericVector>(feasibility_results["predicted_rates_log_var"]));
    bool feasible = feasibility_results["feasible"];
    
    // Stop if not feasible 
    if (!feasible) {
      Rcpp::stop("Negative or nan rates predicted!");
    }
    
    // Convert back to doubles, remove from log space, and save
    predicted_rates_log = NumericVector(n_count_rows);
    predicted_rates = NumericVector(n_count_rows); 
    for (int i = 0; i < n_count_rows; i++) {
      predicted_rates_log[i] = predicted_rates_log_var(i).val();
      predicted_rates[i] = exp(predicted_rates_log_var(i).val()) - 1.0;
    }
    
    // Update saved parameter vector 
    fitted_parameters = parameters;
    optim_results["fitted_parameters"] = parameters;
    
    // Update struc_values 
    for (int i = 0; i < param_struc_idx.size(); i++) {
      int p = param_struc_idx(i);
      struc_values[i] = parameters(p);
    }
    
  } 

// Check and correct parameter feasibility
Rcpp::List wspc::check_parameter_feasibility(
    const sVec& parameters_var,
    const bool verbose
  ) { 
    
    if (verbose) {
      Rcpp::Rcout << "Checking feasibility of provided parameters" << std::endl;
    }
    
    // Initialize vectors to return 
    sVec feasible_parameters_var = parameters_var; 
    sVec predicted_rates_log_var;
    
    // Test if any tpoints are below the buffer 
    bool feasible = test_tpoints(parameters_var);
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
      
    }
    
    if (verbose) {
      if (feasible) {
        Rcpp::Rcout << "Provided parameters are feasible" << std::endl;
      } else {
        Rcpp::Rcout << "Provided parameters not feasible, searching nearby" << std::endl;
      }
    }
    
    // If not feasible, attempt to find a feasible starting point
    if (!feasible) {
      
      // Find and report initial distance to boundary
      sdouble initial_dist = -1 * max_neg_boundary_dist(parameters_var); 
      Rcpp::Rcout << "Initial boundary distance: " << initial_dist << std::endl;
      
      // Variables for optimization
      dVec x = to_dVec(to_NumVec(parameters_var));
      size_t n = x.size();
      
      // Set up NLopt optimizer
      nlopt::opt opt(nlopt::LD_LBFGS, n);
      opt.set_min_objective(wspc::max_neg_boundary_dist_NLopt, this);
      opt.set_ftol_rel(ctol);               // stop when iteration changes objective fn value by less than this fraction 
      opt.set_maxeval(max_evals);           // Maximum number of evaluations to try 
      opt.set_stopval(-0.001);              // ensure negative boundary distance is at least just under zero
      double min_fx;
      int success_code = 0;
      
      // Optimize
      try {
        nlopt::result sc = opt.optimize(x, min_fx);
        success_code = static_cast<int>(sc);
      } catch (std::exception& e) {
        if (verbose) {
          Rcpp::Rcout << "Optimization failed: " << e.what() << std::endl;
        }
        success_code = 0;
      }
      
      // Find and report final distance to boundary
      if (verbose) {
        Rcpp::Rcout << "Numer of evals: " << opt.get_numevals() << std::endl;
        Rcpp::Rcout << "Success code: " << success_code << std::endl;
        Rcpp::Rcout << "Final boundary distance: " << -1 * min_fx << std::endl;
      }
      
      // Final check of negative boundary distance
      if (min_fx >= 0) {
        success_code = 0;
      }
      
      // Check for success
      if (success_code == 0) {
        if (verbose) {
          Rcpp::Rcout << "Could not find a nearby feasible parameters, returning provided ones" << std::endl;
        }
      } else { // found a feasible starting point, save
        if (verbose) {
          Rcpp::Rcout << "Nearby feasible parameters found!" << std::endl;
        }
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
      _["predicted_rates_log_var"] = to_NumVec(predicted_rates_log_var),
      _["parameters_var"] = to_NumVec(feasible_parameters_var)
    );
    
  } 

/*
 * *************************************************************************
 * Export data to R
 */

void wspc::import_fe_diff_ratio(
    const double& fe_diff_ratio,
    const bool& verbose
  ) {
    vprint("Imported fe_difference_ratio: " + std::to_string(fe_diff_ratio), verbose);
    fe_difference_ratio = (sdouble)fe_diff_ratio;
  }

Rcpp::List wspc::results() {
  
  NumericVector predicted_rates_out(n_count_rows);
  NumericVector predicted_rates_log_out(n_count_rows);
  if (predicted_rates.size() == n_count_rows) {
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
  sdouble expected_ran_effect = ssqrt(1.0 / (4.0 * (2.0 * struc_values["beta_shape_rate"] + 1.0)));
  expected_ran_effect *= 2.0; // scale to range from -1 to 1
  sdouble eff_mult = fe_difference_ratio - 1.0;
  sdouble sd_Rt_effect = eff_mult * mean_count_log * expected_ran_effect;
  
  // Make final list to return 
  List results_list = List::create(
    _["model.component.list"] = mc_list,
    _["count.data.summed"] = count_data_summed,
    _["fitted.parameters"] = fitted_parameters,
    _["param.names"] = param_names_clean,
    _["fix"] = fixed_effects,
    _["treatment"] = treat,
    _["grouping.variables"] = grouping_variables,
    _["struc.params"] = struc_values,
    _["param.idx0"] = param_idx, // "0" to indicate this goes out w/ C++ zero-based indexing
    _["token.pool"] = token_pool_list,
    _["computed_sd_Rt_effect"] = sd_Rt_effect.val(),
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
    .method("import_fe_diff_ratio", &wspc::import_fe_diff_ratio)
    .method("results", &wspc::results);
  }
