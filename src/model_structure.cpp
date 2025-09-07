
// model_structure.cpp
#include "wspc.h"

// Treatment combinations **********************************************************************************************

// Generate all combinations of j indices from {0, ..., n-1}
void combinations(
    int n, int j, int start, 
    std::vector<int>& current, 
    std::vector<std::vector<int>>& result
  ) {
    if (current.size() == j) {
      result.push_back(current);
      return;
    }
    for (int i = start; i < n; ++i) {
      current.push_back(i);
      combinations(n, j, i + 1, current, result);
      current.pop_back();
    }
  }

// Generate Cartesian product of chosen CharacterVectors
void cartesian_product_CharVec(
    const std::vector<CharacterVector>& selected_vectors, 
    std::vector<std::vector<String>>& results, 
    std::vector<String>& current, 
    int depth
  ) {
    if (depth == selected_vectors.size()) { 
      results.push_back(current);
      return;
    }
    for (String val : selected_vectors[depth]) {
      current.push_back(val);  
      cartesian_product_CharVec(selected_vectors, results, current, depth + 1);
      current.pop_back();
    }
  }

// Make treatment interactions
std::vector<CharacterVector> make_treatments(
    std::vector<CharacterVector> fix_trt
  ) {
    int n = fix_trt.size();
    
    std::vector<CharacterVector> treatments;
    
    for (int j = 1; j <= n; j++) {
      
      // Generate all index combinations of size j from {0, ..., n-1}
      std::vector<iVec> index_combinations;
      iVec temp; 
      combinations(n, j, 0, temp, index_combinations);
      
      // Generate all possible selections for each combination
      for (iVec indices : index_combinations) {
        std::vector<CharacterVector> selected_vectors;
        
        // Collect the selected fix_trt vectors
        for (int idx : indices) {
          selected_vectors.push_back(fix_trt[idx]);
        }
        
        // Compute Cartesian product of selected vectors
        std::vector<std::vector<String>> draws;
        std::vector<String> current;
        cartesian_product_CharVec(selected_vectors, draws, current, 0);
        
        // Convert each draw (vector<String>) into a CharacterVector and store in treatments
        for (std::vector<String> draw : draws) {
          CharacterVector Rdraw(draw.begin(), draw.end());
          treatments.push_back(Rdraw);  // Convert std::vector<String> to CharacterVector
        }
        
      }
      
    }
    
    return treatments;
    
  }

// Model parameter management ******************************************************************************************

// Build default fixed-effects matrices in shell
List build_beta_shell(
    const CharacterVector& mc_names,
    const CharacterVector& treatment_lvls,
    const CharacterVector& context_lvls,
    const CharacterVector& species_lvls,
    const List& ref_values,
    const List& RtEffs,
    const List& tpointEffs,
    const List& tslopeEffs,
    const IntegerMatrix& degs
  ) {
    
    List beta(mc_names.size());
    beta.names() = mc_names;
    
    int treat_num = treatment_lvls.size(); 
    int context_num = context_lvls.size();
    int species_num = species_lvls.size();
    
    // Loop through model components 
    for (String mc : mc_names) {
      List beta_mc(context_num);
      beta_mc.names() = context_lvls;
     
      // Loop through context values
      for (int c = 0; c < context_num; c++) {
        List beta_mc_cxt(species_num);
        beta_mc_cxt.names() = species_lvls;
       
        // Loop through species values
        for (int s = 0; s < species_num; s++) {
          
          int deg = degs(s, c);
          int col_num = deg; 
          if (mc == "Rt") {col_num++;}
          NumericMatrix bta(treat_num, col_num);
          rownames(bta) = treatment_lvls;
          
          if (deg > 0 || mc == "Rt") {
         
            List RtEffs_cxt = RtEffs[c];
            NumericMatrix RtEffs_Mat = RtEffs_cxt[s];
            
            List tpointEffs_cxt = tpointEffs[c];
            NumericMatrix tpointEffs_Mat = tpointEffs_cxt[s];
            
            List tslopeEffs_cxt = tslopeEffs[c];
            NumericMatrix tslopeEffs_Mat = tslopeEffs_cxt[s];
           
            for (int t = 0; t < treat_num; t++) {
              for (int i = 0; i < col_num; i++) {
                if (t == 0) {
                  List ref_values_cxt = ref_values[c];
                  List ref_values_sps = ref_values_cxt[s];
                  List ref_values_mc = ref_values_sps[mc];
                  bta(t,i) = ref_values_mc(i);  // set ref level
                } else {
                  if (mc == "Rt") {
                    // use estimated rate effect for non-ref treatment levels
                    bta(t, i) = RtEffs_Mat(t, i);
                  } else if (mc =="tpoint") {
                    // use estimated tpoint effect for non-ref treatment levels 
                    bta(t, i) = tpointEffs_Mat(t, i);
                  } else if (mc == "tslope") {
                    // Seed with zero effect for tslope
                    bta(t, i) = tslopeEffs_Mat(t, i); 
                  }
                }
              }
            }
          }
          beta_mc_cxt[s] = bta;
        } 
        beta_mc[c] = beta_mc_cxt;
      } 
      beta[mc] = beta_mc;
    } 
    
    return beta;
    
  }

// Function for converting structured encoding of parameters into a vector
List make_parameter_vector(
    const List& beta, 
    const List& wfactor, 
    const CharacterVector& cxt_lvls, 
    const CharacterVector& sps_lvls, 
    const CharacterVector& rn_lvls, 
    const CharacterVector& mc_list, 
    const CharacterVector& treatment_lvls,
    const IntegerMatrix& degs
  ) {
    
    /*
     * Initializations
     */
    
    // Initialize vector to hold all parameter values (faster to use std and wrap at end)
    dVec param_vector(0);
    
    // Initialize List to hold names
    List param_names;
    
    // Initialize vectors to hold indexes (faster to use std and wrap at end)
    int idx = 0;
    iVec param_wfactor_point_idx; 
    iVec param_wfactor_rate_idx;
    iVec param_wfactor_slope_idx;
    iVec param_beta_Rt_idx;        // ... Excluding ref level (same below)
    iVec param_beta_tslope_idx;
    iVec param_beta_tpoint_idx;
    iVec param_baseline_idx;
    List beta_idx = clone(beta);
    List wfactor_idx = clone(wfactor); 
    
    int n_context = cxt_lvls.size();
    int n_species = sps_lvls.size();
    int n_ran = rn_lvls.size();
   
    /*
     * Beta parameters (ref level and fixed effects)
     */
    
    // Loop through model components 
    for (String mc : mc_list) {
      
      // For make
      List beta_mc = beta[mc];
      // For map
      List beta_idx_mc = clone(beta_mc);
      
      // Loop through context values
      for (int c = 0; c < n_context; c++) {
        String cxt = cxt_lvls[c];
        
        // For make
        List beta_mc_cxt = beta_mc[cxt];
        // For map
        List beta_idx_mc_cxt = clone(beta_mc_cxt);
        
        // Loop through species values
        for (int s = 0; s < n_species; s++) {
          String sps = sps_lvls[s];
          
          // For make
          int deg = degs(s, c);
          // For map
          iVec beta_idx_mc_cxt_sps;
          
          if (deg > 0 || mc == "Rt") { 
            NumericMatrix beta_mc_cxt_sps = beta_mc_cxt[sps];
            
            // Unroll the effects matrices
            for (int t = 0; t < beta_mc_cxt_sps.ncol(); t++) {
              
              // For map
              String t_name = "Tns/Blk" + std::to_string(t + 1);
              
              for (int i = 0; i < treatment_lvls.length(); i++) {
                
                // Map
                CharacterVector param_name;
                if (i == 0) {
                  param_name = CharacterVector::create("baseline", cxt, mc, sps, t_name);
                } else {
                  param_name = CharacterVector::create("beta", mc, cxt, sps, treatment_lvls[i], "X", t_name);
                }
                
                // Add name
                param_names.push_back(param_name);
               
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
                }
                beta_idx_mc_cxt_sps.push_back(idx);
                idx++;  
                
                // Make
                param_vector.push_back(beta_mc_cxt_sps(i,t));
                
              } 
            }
          } 
          beta_idx_mc_cxt[sps] = beta_idx_mc_cxt_sps;
        } 
        beta_idx_mc[cxt] = beta_idx_mc_cxt;
      }
      beta_idx[mc] = beta_idx_mc;
    }
    
    /*
     * Random-effect parameters
     */
    
    // Unroll the warping factors
    
    // For make
    NumericMatrix wfactor_point = wfactor["point"];
    NumericMatrix wfactor_rate = wfactor["rate"];
    NumericMatrix wfactor_slope = wfactor["slope"];
    // For map
    NumericMatrix wfactor_idx_point = clone(wfactor_point); 
    NumericMatrix wfactor_idx_rate = clone(wfactor_rate); 
    NumericMatrix wfactor_idx_slope = clone(wfactor_slope);
    
    for (int s = 0; s < n_species; s++) {
      // Intentionally skip the first level, which is "none"
      for (int r = 1; r < n_ran; r++) {
        
        // Map
        String c_name = sps_lvls[s];
        String r_name = rn_lvls[r];
        // ... Point warp
        CharacterVector param_name_point = CharacterVector::create("wfactor", "point", r_name, "X", c_name);
        param_names.push_back(param_name_point);
        param_wfactor_point_idx.push_back(idx);
        wfactor_idx_point(r, s) = idx;
        idx++;
        // ... Rate warp
        CharacterVector param_name_rate = CharacterVector::create("wfactor", "rate", r_name, "X", c_name);
        param_names.push_back(param_name_rate);
        param_wfactor_rate_idx.push_back(idx);
        wfactor_idx_rate(r, s) = idx; 
        idx++;
        // ... Slope warp
        CharacterVector param_name_slope = CharacterVector::create("wfactor", "slope", r_name, "X", c_name);
        param_names.push_back(param_name_slope);
        param_wfactor_slope_idx.push_back(idx);
        wfactor_idx_slope(r, s) = idx;
        idx++;
       
        // Make
        // ... Point warp
        param_vector.push_back(wfactor_point(r, s));
        // ... Rate warp
        param_vector.push_back(wfactor_rate(r, s));
        // ... Slope warp
        param_vector.push_back(wfactor_slope(r, s));
        
      }
    } 
    wfactor_idx["point"] = wfactor_idx_point;
    wfactor_idx["rate"] = wfactor_idx_rate; 
    wfactor_idx["slope"] = wfactor_idx_slope;
    
    // Reformat parameter names as single strings
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
    
    // Pack up parameter vector and mappings
    List params = List::create(
      _["param_vec"] = Rcpp::wrap(param_vector),
      _["param_names"] = param_names_clean,
      _["param_wfactor_point_idx"] = wrap(param_wfactor_point_idx),
      _["param_wfactor_rate_idx"] = wrap(param_wfactor_rate_idx),
      _["param_wfactor_slope_idx"] = wrap(param_wfactor_slope_idx),
      _["param_beta_Rt_idx"] = wrap(param_beta_Rt_idx),
      _["param_beta_tslope_idx"] = wrap(param_beta_tslope_idx),
      _["param_beta_tpoint_idx"] = wrap(param_beta_tpoint_idx),
      _["param_baseline_idx"] = wrap(param_baseline_idx),
      _["beta_idx"] = beta_idx,
      _["wfactor_idx"] = wfactor_idx
    );
    
    return params;
    
  }

// Extrapolating random-effect free counts ("none's") ******************************************************************
std::vector<IntegerVector> make_extrapolation_pool(
    const sVec& bin,
    const sVec& count, 
    const CharacterVector& context,
    const CharacterVector& species,
    const CharacterVector& ran, 
    const CharacterVector& treatment,
    bool verbose
  ) {
    
    // Construct nan and none masks
    LogicalVector nan_mask = !Rcpp::is_na(to_NumVec(count));
    LogicalVector none_mask = !eq_left_broadcast(ran, "none");
    
    // Find indexes of none rows
    IntegerVector r_idx = Rwhich(!none_mask);
    
    // Initialze vector to hold extrapolation pools
    int n_count = count.size();
    std::vector<IntegerVector> extrapolation_pool(n_count);
    
    // Convert bin vector to numeric and get max bin number
    NumericVector bin_NumVec = to_NumVec(bin);
    int max_bin = Rcpp::max(bin_NumVec);
    
    // Pre-make masks for bin, context, species, and treatment levels
    CharacterVector context_lvls = Rcpp::unique(context);
    CharacterVector species_lvls = Rcpp::unique(species);
    CharacterVector treatment_lvls = Rcpp::unique(treatment);
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
      
      // Attempt to find interpolation pool
      bool found_pool = false;
      int bin_range = 1;
      
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
      if (any_true(eq_left_broadcast(tracker, ri))) {
        vprint("row: " + std::to_string(ri + 1) + "/" + std::to_string(n_rows), verbose);
      } 
      
    }
    
    return extrapolation_pool;
    
  }

sVec extrapolate_none(
    const sVec& count,
    const CharacterVector& ran, 
    const std::vector<IntegerVector>& extrapolation_pool,
    const bool& log_transform
  ) {
    
    sVec count_out = count;
    IntegerVector r_idx = Rwhich(eq_left_broadcast(ran, "none"));
    
    for (int r : r_idx) {
      
      
      IntegerVector extrapolation_pool_r = extrapolation_pool[r];
      int extrapolation_sz = extrapolation_pool_r.size();
      if (extrapolation_sz == 0 || extrapolation_pool_r[0] < 1) {
        // If no extrapolation pool, set to zero
        count_out[r] = 0.0;
      } else {
        // If extrapolation pool found, take mean of pool
        sdouble running_sum = 0.0;
        
        // Take running sum
        if (log_transform) {
          for (int i = 0; i < extrapolation_sz; i++) {
            running_sum += slog(count[extrapolation_pool_r[i]] + 1.0);
          }
        } else {
          for (int i = 0; i < extrapolation_sz; i++) {
            running_sum += count[extrapolation_pool_r[i]];
          }
        }
        
        // Divide by size of pool and round
        sdouble roundedMean;
        if (log_transform) {
          roundedMean = running_sum / (sdouble)extrapolation_sz;
          roundedMean = sexp(roundedMean) - 1.0; // Convert back from log space
          roundedMean = stan::math::round(roundedMean);
        } else {
          roundedMean = stan::math::round(running_sum / (sdouble)extrapolation_sz);
        }
        count_out[r] = roundedMean;
        
      }
      
    }
    
    return count_out;
    
  }
