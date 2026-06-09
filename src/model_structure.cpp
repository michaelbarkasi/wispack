
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
   
    // Initialize vector to hold all treatment combinations
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

// Extrapolating random-effect free counts ("none's") ******************************************************************

sVec extrapolate_none(
    const sVec& count,
    const CharacterVector& ran, 
    const std::vector<IntegerVector>& extrapolation_pool,
    const bool& log_transform,
    const bool& round_count
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
        
        // Divide by size of pool
        sdouble count_out_r;
        if (log_transform) {
          count_out_r = running_sum / (sdouble)extrapolation_sz;
          count_out_r = sexp(count_out_r) - 1.0; // Convert back from log space
        } else {
          count_out_r = running_sum / (sdouble)extrapolation_sz;
        }
        
        // Round, if specified
        if (round_count) {count_out_r = stan::math::round(count_out_r);}
        
        // Set
        count_out[r] = count_out_r;
        
      }
      
    }
    
    return count_out;
    
  }
