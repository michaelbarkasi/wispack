
// wspmm.h
#ifndef WSPMM_H
#define WSPMM_H

// Rcpp
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(StanHeaders)]]
#include <stan/math.hpp>                      // pulls in everything from rev/ and prim/
#include <stan/math/memory/stack_alloc.hpp>   // for better/faster stack memory allocation
#include <Rcpp.h>
#include <RcppEigen.h>
#include <Rcpp/Benchmark/Timer.h>
#include <nlopt.hpp>
#include "pcg/pcg_random.hpp"
#include <unistd.h> // for forking on unix
using namespace Rcpp;

// [[Rcpp::plugins(cpp14)]]

// Define aliases for convenience
using sdouble = stan::math::var;
using sMat = Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic>;
using sVec = Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>;
using dVec = std::vector<double>; 
using iVec = std::vector<int>;
constexpr sdouble (*spower)(const sdouble&, const sdouble&) = stan::math::pow;
constexpr sdouble (*smax)(const sVec&) = stan::math::max;
constexpr sdouble (*smin)(const sVec&) = stan::math::min;
constexpr sdouble (*slog)(const sdouble&) = stan::math::log;
constexpr sdouble (*sexp)(const sdouble&) = stan::math::exp;
constexpr sdouble (*ssqrt)(const sdouble&) = stan::math::sqrt;

// Constants ***********************************************************************************************************

const sdouble inf_ = 1e100;   // pseudo-infinity value for optimization

// Helper functions, printing ******************************************************************************************

// Function for printing character messages
void vprint(const std::string& message, bool verbose);
void vprintV(const CharacterVector& message, bool verbose);

// Function for printing vectors
void print_Vec(const NumericVector& vec);
// ... overload for IntegerVector
void print_Vec(const IntegerVector& vec);
// ... overload for std::vector
void print_Vec(const dVec& vec);
// ... overload for sdouble
void print_Vec(const sVec& vec);

// Vector-type conversions *********************************************************************************************

// Convert to Eigen::Matrix with sdouble elements
// ... from NumericVector
sVec to_sVec(const NumericVector& vec);
// ... overload, from std::vector with doubles 
sVec to_sVec(const dVec& vec);
// ... overload, from std::vector with int
sVec to_sVec(const iVec& vec);
// ... overload, from IntegerVector
sVec to_sVec(const IntegerVector& vec);

// Convert to NumericVector 
// ... from Eigen::Matrix with sdouble elements
NumericVector to_NumVec(const sVec& vec);
// ... overload, from std::vector with doubles
NumericVector to_NumVec(const dVec& vec);

// Convert to IntegerVector
// ... from std::vector with int
IntegerVector to_IntVec(const iVec& vec);

// Convert to std::vector with doubles 
// ... from NumericVector
dVec to_dVec(const NumericVector& vec);

// Convert to std::vector with int
// ... from IntegerVector
iVec to_iVec(const IntegerVector& vec);

// Vector operations ***************************************************************************************************

// Means
// -----

// Mean of vector elements
double vmean(const dVec& x);
// ... overload 
sdouble vmean(const sVec& x);
// ... overload 
double vmean(const NumericVector& x);

// Mean of vector elements within a range
double vmean_range(const dVec& x, const int& start, const int& end);
// ... overload 
sdouble vmean_range(const sVec& x, const int& start, const int& end);
// ... overload
double vmean_range(const NumericVector& x, const int& start, const int& end);

// Rolling mean
dVec roll_mean(const dVec& series, int filter_ws);

// Component-wise operations
// -------------------------

// Subtraction
dVec vsubtract(const dVec& x, const dVec& y);

// Masks and indexes
// -----------------

// Convert boolean masks to integer indexes
IntegerVector Rwhich(const LogicalVector& x);

// Return input's indexes in the ascending order of its elements
IntegerVector Rorder(const NumericVector& x);
// ... overload
IntegerVector Rorder(const dVec& x);

// For subsetting vectors with Rcpp mask 
sVec masked_vec(sVec vec, Rcpp::LogicalVector mask);
// ... overload 
dVec masked_vec(dVec vec, Rcpp::LogicalVector mask);
// ... overload
iVec masked_vec(iVec vec, Rcpp::LogicalVector mask);

// For subsetting vectors with Rcpp IntegerVector
sVec idx_vec(sVec vec, Rcpp::IntegerVector idx);
// ... overload 
dVec idx_vec(dVec vec, Rcpp::IntegerVector idx);
// ... overload 
iVec idx_vec(iVec vec, Rcpp::IntegerVector idx);

// Sequence generation
// -------------------

// Sequence of doubles
NumericVector dseq(const double& start, const double& end, const int& lengthout);

// Sequence of integers
IntegerVector iseq(const int& start, const int& end, const int& lengthout);

// List assignments 
// ----------------

// Name proxy list 
void name_proxylist(List list, const CharacterVector& new_names);

// Assign to proxy list 
void assign_proxylist(List list, String element, List assignment);
// ... overload
void assign_proxylist(List list, String element, NumericVector assignment);

// Misc
// ____

// merge two integer vectors 
IntegerVector buffered_merge(
    const IntegerVector& a, 
    const IntegerVector& b,
    const int& buffer
  );

// Vectorized logic ****************************************************************************************************

// Return logical vector giving elements of left which match right
LogicalVector eq_left_broadcast(const CharacterVector& left, const String& right);
// ... overload 
LogicalVector eq_left_broadcast(const IntegerVector& left, const int& right);
// ... overload
LogicalVector eq_left_broadcast(const NumericVector& left, const double& right);

// Quantifiers
bool any_true(const LogicalVector& x);
bool all_true(const LogicalVector& x);

// Model structure *****************************************************************************************************

// Treatment combinations 
// ----------------------

// Generate all combinations of j indices from {0, ..., n-1}
void combinations(
  int n, int j, int start, 
  std::vector<int>& current, 
  std::vector<std::vector<int>>& result
  );

// Generate Cartesian product of chosen CharacterVectors
void cartesian_product_CharVec(
  const std::vector<CharacterVector>& selected_vectors, 
  std::vector<std::vector<String>>& results, 
  std::vector<String>& current, 
  int depth
  );

// Make treatment interactions
std::vector<CharacterVector> make_treatments(
  std::vector<CharacterVector> fix_trt
  );

// Model parameter management
// --------------------------

// Build default fixed-effects matrices in shell
List build_beta_shell(
  const CharacterVector& mc_names,
  const CharacterVector& treatment_lvls,
  const CharacterVector& parent_lvls,
  const CharacterVector& child_lvls,
  const List& ref_values,
  const List& RtEffs,
  const IntegerMatrix& degs
  );

// Function for converting structured encoding of parameters into a vector
List make_parameter_vector(
  const List& beta, 
  const List& wfactor, 
  const NumericVector& struc,
  const CharacterVector& prt_lvls, 
  const CharacterVector& cld_lvls, 
  const CharacterVector& rn_lvls, 
  const CharacterVector& mc_list, 
  const CharacterVector& treatment_lvls,
  const CharacterVector& struc_names,
  const IntegerMatrix& degs
  );

// Extrapolating random-effect free counts
// ---------------------------------------

// Find row numbers ("pool") of observed counts to use for extrapolation of "none's"
 std::vector<IntegerVector> make_extrapolation_pool(
  const sVec& bin,
  const sVec& count, 
  const CharacterVector& parent,
  const CharacterVector& child,
  const CharacterVector& ran, 
  const CharacterVector& treatment,
  bool verbose
  );

// Use pools to extrapolate counts
sVec extrapolate_none(
  const sVec& count,
  const CharacterVector& ran, 
  const std::vector<IntegerVector>& extrapolation_pool
  );

// Model fitting *******************************************************************************************************

// Numerically stable implementation of sigmoid function
sdouble sigmoid_stable(
  const sdouble& x
  );

// Core poly-sigmoid function of the model
sdouble poly_sigmoid(
  const sdouble& b,        // input variable
  const int& deg,          // degree of the poly-sigmoid, i.e., number of transitions between blocks
  const sVec& Rt,          // vector containing the height ("rate") of each block
  const sVec& tslope,      // vector containing the slope of each transition between blocks
  const sVec& tpoint       // vector containing the point of each transition in the bin space
  );

// Warping function for random effects on transition points
sdouble point_warp(
  const sdouble& b,        // Bin coordinate
  const sdouble& bin_num,  // Max bin coordinate
  const sdouble& f         // warping factor
  );

// Warping function for random effects on count rate
sdouble rate_warp(
  const sdouble& rate,     // rate output by the poly-sigmoid
  const sdouble& f         // warping factor
  );

// Function for computing model components
sVec compute_mc(
  const sMat& betas,       // effects matrix, rows as interactions, columns as blocks/transitions
  const sVec& weight_row   // row of weight matrix, elements as fixed effects
  );

// Inverse quadratic ramp function for boundary penalty
sdouble boundary_penalty_transform(
  const sdouble& x,
  const sdouble& a
  );

// Calculate rolling-window negloglik of a series being generated by a given rate, with or without change-point
dVec series_nll(
  const dVec& series0,          // 1D vector of points for which to take negative log-likelihood of a change-point
  const int& ws,                // Running window size
  const int& filter_ws,         // Size of window for taking rolling mean
  const bool& null              // If true, compute likelihood of data assuming no transitions; otherwise, assuming transition
  );

// Likelihood ratio outlier change-point detection
IntegerVector LROcp(
  const dVec& series,          // 1D vector of points to test for change points
  const int& ws,                // Running window size
  const int& filter_ws,         // Size of window for taking rolling mean
  const double& out_mult        // Outlier multiplier
  );

#endif