
// model_structure.cpp
#include "wspc.h"

// Treatment combinations **********************************************************************************************

// Generate all combinations of j indices from {0, ..., n-1}
void combinations(
    int n, int j, int start, 
    iVec& current, 
    std::vector<iVec>& result
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



