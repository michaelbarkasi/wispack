
// wspc.h
#ifndef WSPC_H
#define WSPC_H

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

// Main class **********************************************************************************************************

/*
 * Object class to hold and fit Warped Sigmoidal Poisson-Process Mixed-Effect Model (WSPmm) model. 
 */

class wspc {
  
  public:
    
    // Fields **********************************************************************************************************
    
    // Data sizes
    int n_count_rows;       // number of rows in summed count data frame
    int treatment_num;      // number of treatment combinations (including "none", i.e., the reference)
    sdouble bin_num;        // max bin (i.e., number of bins)
    IntegerVector count_row_nums;  // sequence of integers from 0 to n_count_rows - 1
    
    // Data variables, stan
    sVec bin;               // bin column of summed data
    sVec count;             // count column of summed data
    sVec count_log;         // log of observed counts
    sVec count_tokenized;   // tokenized count data
    
    // Data variables, Rcpp ("factors")
    CharacterVector parent;                 // parent column of summed data
    CharacterVector child;                  // child column of summed data
    CharacterVector ran;                    // random effect column of summed data
    CharacterVector treatment;              // treatment column of summed data
    
    // Model predictions, Rcpp 
    NumericVector predicted_rates_log;
    NumericVector predicted_rates; 
    
    // Fixed effect variables
    CharacterVector fix_names;              // names of fixed effect variables
    std::vector<CharacterVector> fix_lvls;  // levels for each fixed effect
    std::vector<CharacterVector> fix_trt;   // treatment levels for each fixed effect
    std::vector<CharacterVector> treatment_components;   // all possible treatment combinations, level components
    CharacterVector treatment_lvls;         // all possible treatment combinations, levels as single-string name
    CharacterVector fix_ref;                // reference level for each fixed effect
    
    // Grouping variables
    CharacterVector parent_lvls;            // levels of parent grouping variable (fixed-effects)
    CharacterVector child_lvls;             // levels of child grouping variable (fixed-effects)
    CharacterVector ran_lvls;               // levels of random effect grouping variable
    
    // Variables to help with data manipulation
    sMat weights;                           // weight matrix, rows as rows of summed count data, columns as treatments (first column is reference)
    IntegerVector idx_mc_unique;            // count data rows at which model component values will change
    std::vector<IntegerVector> token_pool;  // list of token count indexes associated with each summed count row
    std::vector<IntegerVector> extrapolation_pool;  // list of summed-count indexes giving summed count rows from which to extrapolate
    IntegerVector count_not_na_idx;         // indexes of non-NA rows in summed count data
    LogicalVector count_not_na_mask;        // mask of non-NA rows in summed count data
    
    // Variables related to model parameters
    CharacterVector mc_list = {"Rt", "tslope", "tpoint"}; // list of model components, rates, transition slopes, transition points
    IntegerMatrix degMat;                   // matrix of degrees for each parent (column) - child (rows) pair
    NumericVector fitted_parameters;        // vector holding the model parameters
    List param_names;                       // list holding the names of the model parameters as they appear in fitted_parameters
    IntegerVector param_wfactor_point_idx;  // ... indexes in parameter vector for different kinds of model parameters
    IntegerVector param_wfactor_rate_idx;
    IntegerVector param_beta_Rt_idx;
    IntegerVector param_beta_Rt_idx_no_ref;
    IntegerVector param_beta_tslope_idx;
    IntegerVector param_beta_tslope_idx_no_ref;
    IntegerVector param_beta_tpoint_idx;
    IntegerVector param_beta_tpoint_idx_no_ref;
    IntegerVector param_ref_values_tpoint_idx;
    IntegerVector param_struc_idx;
    List beta_idx;                          // lists giving the structured array indices for named parameters
    List wfactor_idx; 
    IntegerVector gv_ranLr_int;             // indices (row and column) for random effect arrays 
    IntegerVector gv_fixLr_int;  
    NumericVector struc_values = {1.0, 1.0, 1.0, 1.0};
    CharacterVector struc_names = {
      "beta_shape_point", 
      "beta_shape_rate",
      "sd_tpoint_effect",
      "sd_tslope_effect"
    };
    sdouble buffer_factor = 0.05;           // scaling factor for buffer value, the minimum distance between transition points 
    sdouble tpoint_buffer;                  // min number of bins between transition points (immutable structural parameter)
    double LROcutoff = 2.0;                 // cutoff (x sd) for likelihood ratio outlier detection
    double tslope_initial = 1.0;            // initial value for transition slope
    double wf_initial = 0.5;                // initial value for warping factor ... any sensible magnitude > 0.1 and < 0.75 should do? 
    
    // Optimization settings
    int max_evals = 200;  // max number of evaluations
    double ctol = 5e-4;   // convergence tolerance
    
    // Other settings 
    List model_settings;
    
    /*
     * Basic idea: there will be a number of relevant distances to a boundary. Want to transform
     *  these distances so that they go to infinity as they approach zero, yet while they are 
     *  sufficiently far from zero ("at distance"), they are not too large.
     */
    
    // Boundary penalty variables
    int boundary_vec_size = 0;                     // number of boundary components
    sVec bp_coefs;                                 // Coefficients used to scale boundary penality so it's negligible when far from boundary, and infinity at boundary
    sdouble max_penalty_at_distance_factor = 0.01; // as fraction of initial neg_loglik, the max penalty when far from the boundary
    sdouble max_penalty_at_distance = 1;           // variable to hold the max penalty value, once computed
    
    // Variables for holding results 
    List optim_results;                            // results from optimization
   
    // Methods *********************************************************************************************************
   
    // Constructor
    wspc(Rcpp::DataFrame count_data, Rcpp::List settings, bool verbose);
    // Destructor
    ~wspc();
    // R copy 
    wspc clone() const;
    
    // ... computing predicted values from parameters 
    
    // Compute model component values for rows of summed count data
    // ... for Rt
    sVec compute_mc_Rt_r(
        const int& r,
        const sVec& parameters,
        const sdouble& f_rw
    ) const;
    
    // ... tslope
    sVec compute_mc_tslope_r(
        const int& r,
        const sVec& parameters
    ) const;
    
    // ... tpoint
    sVec compute_mc_tpoint_r(
        const int& r,
        const sVec& parameters
    ) const;
    
    // Predict log of rates
    sVec predict_rates_log(
        const sVec& parameters,
        const bool& all_rows 
    ) const;
    
    // ... computing objective function (i.e., fitting model and parameter boundary distances)
    
    // Compute neg_loglik of the model under the given parameters
    sdouble neg_loglik(
        const sVec& parameters
    ) const;
    
    // Compute boundary distances
    sVec neg_boundary_dist(
        const sVec& parameters
    ) const;
    
    // Test for tpoints below the buffer
    bool test_tpoints(
        const sVec& parameters
    ) const;
    
    // Compute min boundary penalty
    sdouble max_neg_boundary_dist(
        const sVec& parameters
    ) const;
    
    // Wrap neg_min_boundary_dist in form needed for NLopt constraint function
    static double max_neg_boundary_dist_NLopt(
        const dVec& x,
        dVec& grad,
        void* data
    );
    
    // Compute neg_loglik plus boundary penalty (main objective function) 
    sdouble bounded_nll(
        const sVec& parameters
    ) const;
    
    // Wrap bounded_nll in form needed for NLopt objective function
    static double bounded_nll_NLopt(
        const dVec& x, 
        dVec& grad, 
        void* data
    );
    
    // ... Computing gradients with stan reverse-mode autodiff
    
    // Compute the gradient of the bounded_nll function
    // ... this is the gradient function used in model optimization
    Eigen::VectorXd grad_bounded_nll(
        const sVec& p_
    ) const;
    
    // Compute the gradient of the max_neg_boundary_dist function
    // ... this is the gradient function used in the search for feasible parameters
    Eigen::VectorXd grad_max_neg_boundary_dist(
        const sVec& p_
    ) const;
    
    // ... Bootstrapping and model fitting, for statistical testing
    
    // Fit model using NLopt
    void fit(
        const bool set_params, 
        const bool verbose
    );
    
    // Fit model to bootstrap resample
    dVec bs_fit(
        int bs_num,                 // A unique number for this resample
        bool clear_stan             // Recover stan memory at end?
    );
    
    // Fork bootstraps (parallel processing)
    Rcpp::NumericMatrix bs_batch(
        int bs_num_max,              // Number of bootstraps to perform
        int max_fork,                // Maximum number of forked processes per batch
        bool verbose
    );
    
    // Resample (demo)
    Rcpp::NumericMatrix resample(
        int n_resample                 // total number of resamples to draw
    );
    
    // ... Setting parameters
    
    // Set the model with the given parameters
    void set_parameters(
        const NumericVector& parameters,
        const bool verbose
    );
    
    // Check and correct parameter feasibility
    Rcpp::List check_parameter_feasibility(
        const sVec& parameters_var,
        const bool verbose
    );
    
    // ... export data to R
    Rcpp::List results(); 
    
    // ... misc and debugging 
    // Wrap neg_loglik in form needed for R
    double bounded_nll_debug(
        const dVec& x
    );
    
    // Wrap neg_loglik in form needed for R
    Rcpp::NumericVector grad_bounded_nll_debug(
        const dVec& x 
    );
  
};

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

// Log of density of normal distribution centered on zero
sdouble log_dnorm_centered(
    const sdouble& x,        // value to evaluate
    const sdouble& sd        // standard deviation
  );

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

#endif // WSPC_H