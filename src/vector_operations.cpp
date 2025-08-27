
// vector_operations.cpp
#include "wspc.h"

// Vectorized means ****************************************************************************************************

// Mean of vector elements
double vmean(
    const dVec& x
  ) {
    double sum = 0.0;
    int ctr = 0;
    for (double xi : x) {
      if (!std::isnan(xi)) {
        sum += xi;
        ctr++;
      }
    }
    return sum / (double)ctr; 
  }

// ... overload 
sdouble vmean(
    const sVec& x
  ) {
    sdouble sum = 0.0;
    int ctr = 0;
    for (sdouble xi : x) {
      if (!std::isnan(xi)) {
        sum += xi;
        ctr++;
      }
    }
    return sum / (sdouble)ctr; 
  }

// ... overload 
double vmean(
    const NumericVector& x
  ) {
    double sum = 0.0;
    int ctr = 0;
    for (int i = 0; i < x.size(); i++) {
      if (!std::isnan(x[i])) {
        sum += x[i];
        ctr++;
      }
    }
    return sum / (double)ctr; 
  }

// ... overload 
int vmean(
    const iVec& x
  ) {
    int sum = 0;
    int ctr = 0;
    for (int xi : x) {
      if (!std::isnan(xi)) {
        sum += xi;
        ctr++;
      }
    }
    return std::round(sum / ctr); 
  }

// Mean of vector elements within a range
double vmean_range(
    const dVec& x,
    const int& start,
    const int& end
  ) {
    if (start == end) {return x[start];}
    double sum = 0.0;
    int ctr = 0;
    for (int i = start; i < end; i++) {
      if (!std::isnan(x[i])) {
        sum += x[i];
        ctr++;
      }
    }
    return sum / (double)ctr; 
  }

// ... overload 
sdouble vmean_range(
    const sVec& x,
    const int& start,
    const int& end
  ) {
    if (start == end) {return x[start];}
    sdouble sum = 0.0;
    int ctr = 0;
    for (int i = start; i < end; i++) {
      if (!std::isnan(x[i].val())) {
        sum += x[i];
        ctr++;
      }
    }
    return sum / (sdouble)ctr; 
  }

// ... overload
double vmean_range(
    const NumericVector& x,
    const int& start,
    const int& end
  ) {
    if (start == end) {return x[start];}
    double sum = 0.0;
    int ctr = 0;
    for (int i = start; i < end; i++) {
      if (!std::isnan(x[i])) {
        sum += x[i];
        ctr++;
      }
    }
    return sum / (double)ctr; 
  }

// Rolling mean
dVec roll_mean(
    const dVec& series,    // 1D vector of points to take rolling mean
    int filter_ws          // Size of window for taking rolling mean
  ) {
    int n = series.size();
    dVec series_out(n);
    for (int i = 0; i < n; i++) {
      int ws = std::min(i + 1, filter_ws);
      // Ensure iterator bounds stay within range
      auto start = series.begin() + (i - ws + 1);
      auto end = series.begin() + (i + 1);
      series_out[i] = std::accumulate(start, end, 0.0) / static_cast<double>(ws);
    } 
    return series_out;
  }

// Variance of vector elements
sdouble vvar(
    const sVec& x
  ) {
    sdouble mean = vmean(x);
    sdouble sum = 0.0;
    int ctr = 0;
    for (sdouble xi : x) {
      if (!std::isnan(xi)) {
        sum += (xi - mean) * (xi - mean);
        ctr++;
      }
    }
    return sum / (sdouble)ctr;
  }

// Standard deviation of vector elements 
sdouble vsd(
    const sVec& x
  ) {
    return ssqrt(vvar(x));
  }

// ... overload
double vsd(
    const dVec& x
  ) {
    double mean = vmean(x);
    double sum = 0.0;
    int ctr = 0;
    for (double xi : x) {
      if (!std::isnan(xi)) {
        sum += (xi - mean) * (xi - mean);
        ctr++;
      }
    }
    return std::sqrt(sum / (double)ctr);
  }

// ... overload
double vsd(
    const NumericVector& x
  ) {
    double mean = vmean(x);
    double sum = 0.0;
    int ctr = 0;
    for (double xi : x) {
      if (!std::isnan(xi)) {
        sum += (xi - mean) * (xi - mean);
        ctr++;
      }
    }
    return std::sqrt(sum / (double)ctr);
  }

// Component-wise operations *******************************************************************************************

// Component-wise multiplication
dVec vmult(
    const dVec& x,
    const dVec& y
  ) {
    int n = x.size();
    dVec out(n);
    for (int i = 0; i < n; i++) {
      out[i] = x[i] * y[i];
    }
    return out;
  }

// ... overload for NumericVector
NumericVector vmult(
    const NumericVector& x,
    const NumericVector& y
  ) {
    int n = x.size();
    NumericVector out(n);
    for (int i = 0; i < n; i++) {
      out[i] = x[i] * y[i];
    }
    return out;
  }

// Component-wise division
dVec vdivide(
    const dVec& x,
    const dVec& y
  ) {
    int n = x.size();
    dVec out(n);
    for (int i = 0; i < n; i++) {
      out[i] = x[i] / y[i];
    }
    return out;
  }

// ... overload for NumericVector
NumericVector vdivide(
    const NumericVector& x,
    const NumericVector& y
  ) {
    int n = x.size();
    NumericVector out(n);
    for (int i = 0; i < n; i++) {
      out[i] = x[i] / y[i];
    }
    return out;
  }

// Component-wise addition
dVec vadd(
    const dVec& x,
    const dVec& y
  ) {
    int n = x.size();
    dVec out(n);
    for (int i = 0; i < n; i++) {
      out[i] = x[i] + y[i];
    }
    return out;
  }

// ... overload for NumericVector
NumericVector vadd(
    const NumericVector& x,
    const NumericVector& y
  ) {
    int n = x.size();
    NumericVector out(n);
    for (int i = 0; i < n; i++) {
      out[i] = x[i] + y[i];
    }
    return out;
  }

// Component-wise subtraction 
dVec vsubtract(
    const dVec& x,
    const dVec& y
  ) {
    int n = x.size();
    dVec out(n);
    for (int i = 0; i < n; i++) {
      out[i] = x[i] - y[i];
    }
    return out;
  }

// ... overload for NumericVector
NumericVector vsubtract(
    const NumericVector& x,
    const NumericVector& y
  ) {
    int n = x.size();
    NumericVector out(n);
    for (int i = 0; i < n; i++) {
      out[i] = x[i] - y[i];
    }
    return out;
  }

// Vectorized sqrt
NumericVector vsqrt(
    const NumericVector& x
  ) {
    int n = x.size();
    NumericVector out(n);
    for (int i = 0; i < n; i++) {
      out[i] = std::sqrt(x[i]);
    }
    return out;
  }

// Masks and indexes ***************************************************************************************************

// Convert boolean masks to integer indexes
IntegerVector Rwhich(
    const LogicalVector& x
  ) {
    iVec indices;  // Use std::vector for efficient dynamic resizing
    for (int i = 0; i < x.size(); ++i) {
      if (x[i]) {
        indices.push_back(i);
      }
    }
    return wrap(indices);  // Convert std::vector to IntegerVector
  }

// Return input's indexes in the ascending order of its elements
IntegerVector Rorder(
    const NumericVector& x
  ) {
    IntegerVector indices = Rcpp::seq(0, x.size() - 1);
    std::sort(indices.begin(), indices.end(), [&x](const int i1, const int i2) { return x[i1] < x[i2]; });
    return indices;
  }

// ... overload
IntegerVector Rorder(
    const dVec& x
  ) {
    IntegerVector indices = Rcpp::seq(0, x.size() - 1);
    std::sort(indices.begin(), indices.end(), [&x](const int i1, const int i2) { return x[i1] < x[i2]; });
    return indices;
  }

// For subsetting vectors with Rcpp mask 
sVec masked_vec(
    sVec vec, 
    Rcpp::LogicalVector mask
  ) {
    int m = mask.size();
    if (vec.size() != m) {Rcpp::stop("Vector and mask must have the same length.");}
    int n = 0;
    for (int i = 0; i < m; i++) {if (mask[i]) {n++;}}
    sVec mvec = sVec(n);
    int j = 0;
    for (int i = 0; i < m; i++) {if (mask[i]) {mvec[j++] = vec[i];}}
    return mvec;
  } 

// ... overload 
dVec masked_vec(
    dVec vec, 
    Rcpp::LogicalVector mask
  ) {
    int m = mask.size();
    if (vec.size() != m) {Rcpp::stop("Vector and mask must have the same length.");}
    int n = 0;
    for (int i = 0; i < m; i++) {if (mask[i]) {n++;}}
    dVec mvec = dVec(n);
    int j = 0;
    for (int i = 0; i < m; i++) {if (mask[i]) {mvec[j++] = vec[i];}}
    return mvec;
  }

// ... overload
iVec masked_vec(
    iVec vec, 
    Rcpp::LogicalVector mask
  ) {
    int m = mask.size();
    if (vec.size() != m) {Rcpp::stop("Vector and mask must have the same length.");}
    int n = 0;
    for (int i = 0; i < m; i++) {if (mask[i]) {n++;}}
    iVec mvec = iVec(n);
    int j = 0;
    for (int i = 0; i < m; i++) {if (mask[i]) {mvec[j++] = vec[i];}}
    return mvec;
  }

// ... overload
IntegerVector masked_vec(
    IntegerVector vec, 
    Rcpp::LogicalVector mask
  ) {
    int m = mask.size();
    if (vec.size() != m) {Rcpp::stop("Vector and mask must have the same length.");}
    int n = 0;
    for (int i = 0; i < m; i++) {if (mask[i]) {n++;}}
    IntegerVector mvec = IntegerVector(n);
    int j = 0;
    for (int i = 0; i < m; i++) {if (mask[i]) {mvec[j++] = vec[i];}}
    return mvec;
  }

// For subsetting vectors with Rcpp IntegerVector
sVec idx_vec(
    sVec vec,
    Rcpp::IntegerVector idx
  ) {
    int m1 = Rcpp::max(idx);
    if (vec.size() <= m1) {Rcpp::stop("Index out of bounds.");}
    int m2 = Rcpp::min(idx); 
    if (m2 < 0) {Rcpp::stop("Index out of bounds.");}
    int n = idx.size();
    sVec mvec = sVec(n);
    for (int i = 0; i < n; i++) {mvec[i] = vec[idx[i]];}
    return mvec;
  }

// ... overload 
dVec idx_vec(
    dVec vec,
    Rcpp::IntegerVector idx
  ) {
    int m1 = Rcpp::max(idx);
    if (vec.size() <= m1) {Rcpp::stop("Index out of bounds.");}
    int m2 = Rcpp::min(idx); 
    if (m2 < 0) {Rcpp::stop("Index out of bounds.");}
    int n = idx.size();
    dVec mvec = dVec(n);
    for (int i = 0; i < n; i++) {mvec[i] = vec[idx[i]];}
    return mvec;
  }

// ... overload 
iVec idx_vec(
    iVec vec,
    Rcpp::IntegerVector idx
  ) {
    int m1 = Rcpp::max(idx);
    if (vec.size() <= m1) {Rcpp::stop("Index out of bounds.");}
    int m2 = Rcpp::min(idx); 
    if (m2 < 0) {Rcpp::stop("Index out of bounds.");}
    int n = idx.size();
    iVec mvec = iVec(n);
    for (int i = 0; i < n; i++) {mvec[i] = vec[idx[i]];}
    return mvec;
  }

// ... overload
CharacterVector idx_vec(
    CharacterVector vec,
    Rcpp::IntegerVector idx
  ) {
    int m1 = Rcpp::max(idx);
    if (vec.size() <= m1) {Rcpp::stop("Index out of bounds.");}
    int m2 = Rcpp::min(idx); 
    if (m2 < 0) {Rcpp::stop("Index out of bounds.");}
    int n = idx.size();
    CharacterVector mvec = CharacterVector(n);
    for (int i = 0; i < n; i++) {mvec[i] = vec[idx[i]];}
    return mvec;
  }

// Find matches in character vector
IntegerVector grep_cpp(
    CharacterVector V, 
    std::string pattern
  ) {
    IntegerVector idx;
    for (int i = 0; i < V.size(); i++) {
      std::string str = Rcpp::as<std::string>(V[i]);
      if (str.find(pattern) != std::string::npos) {
        idx.push_back(i); 
      }
    }
    return idx;
  }

// Find matches in string
bool pattern_match(
    std::string pattern,
    std::string test
  ) {
    bool match = false;
    if (test.find(pattern) != std::string::npos) {
      match = true; 
    }
    return match;
  }

// Sequence generation *************************************************************************************************

// Sequence of doubles
NumericVector dseq(
    const double& start, 
    const double& end, 
    const int& lengthout
  ) {
    if (lengthout < 2) {Rcpp::stop("lengthout must be at least 2.");}
    double by = (end - start) / (lengthout - 1.0);
    NumericVector seq(lengthout);
    for (int i = 0; i < lengthout; i++) {
      seq[i] = start + i * by;
    }
    return seq;
  }

// Sequence of integers
IntegerVector iseq(
    const int& start,
    const int& end, 
    const int& lengthout
  ) {
    if (lengthout < 2) {Rcpp::stop("lengthout must be at least 2.");}
    int by = (end - start) / (lengthout - 1);
    IntegerVector seq(lengthout);
    for (int i = 0; i < lengthout; i++) {
      seq[i] = start + i * by;
    }
    return seq;
  }

// List assignments ****************************************************************************************************

// Name proxy list 
void name_proxylist(
    List list,
    const CharacterVector& new_names
  ) {
    List out = List(list);
    out.names() = new_names;
    list = out;
  }

// Assign to proxy list 
void assign_proxylist(
    List list,
    String element, 
    List assignment
  ) {
    List out = List(list);
    out[element] = assignment;
    list = out;
  }

// ... overload
void assign_proxylist(
    List list,
    String element, 
    NumericVector assignment
  ) {
    List out = List(list);
    out[element] = assignment;
    list = out;
  }

// ... overload 
void assign_proxylist(
    List list,
    String element, 
    IntegerMatrix assignment
  ) {
    List out = List(list);
    out[element] = assignment;
    list = out;
  }

// Misc ****************************************************************************************************************

// Merge two integer vectors 
IntegerVector buffered_merge(
    const IntegerVector& a, 
    const IntegerVector& b,
    const int& buffer
  ) {
    IntegerVector out;
    IntegerVector a_range = Rcpp::range(a);
    IntegerVector b_range = Rcpp::range(b);
    int max_val = a_range[1];
    if (b_range[1] > a_range[1]) {max_val = b_range[1];}
    IntegerVector test_range = Rcpp::seq(1, max_val);
    for (int i : test_range) {
      if (any_true(eq_left_broadcast(a, i)) || any_true(eq_left_broadcast(b, i))) {
        int os = out.size();
        if (os > 0) {
          if (i - out[os - 1] > buffer) {out.push_back(i);}
        } else {
          out.push_back(i);
        }
      }
    }
    return out;
  }

// Return the indices of the elements of vec of which x is between
iVec block_idx(
    const NumericVector& vec,
    const double& x
  ) {
    iVec out = {-1, -1};
    int vec_size = vec.size();
    if (vec_size > 0) {
      for (int i = 0; i < vec_size; i++) {
        if (x < vec[i]) {
          out[1] = i;
          if (i > 0) {
            if (x > vec[i - 1]) {
              out[0] = i - 1;
            }
          }
          break;
        } else if (i == vec_size - 1) {
          out[0] = i;
        }
      }
    }
    return out;
  }

// Return differences between elements 
IntegerVector vdiff(
    const IntegerVector& x
  ) {
    int n = x.size();
    if (n < 2) {Rcpp::stop("Vector must have at least two elements for diff.");}
    IntegerVector out(n - 1);
    for (int i = 0; i < n - 1; i++) {
      out[i] = x[i + 1] - x[i];
    }
    return out;
  }

#include <Rcpp.h>
using namespace Rcpp;

// Remove row from IntegerMatrix
IntegerMatrix remove_row(
    IntegerMatrix M, 
    int i
  ) {
    int nrow = M.nrow();
    int ncol = M.ncol();
    
    if (i < 0 || i >= nrow) stop("Row index out of range");
   
    IntegerMatrix out(nrow - 1, ncol);
   
    // Copy rows before i
    for (int r = 0; r < i; r++) {
      for (int c = 0; c < ncol; c++) {
        out(r, c) = M(r, c);
      }
    } 
    // Copy rows after i
    for (int r = i + 1; r < nrow; r++) {
      for (int c = 0; c < ncol; c++) {
        out(r - 1, c) = M(r, c);
      }
    } 
    
    return out;
  }
 