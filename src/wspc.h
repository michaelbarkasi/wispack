
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
#include <unistd.h>                           // for forking on unix
#include <random>                             // for thread-safe random number generation
#include <sys/wait.h>                         // for Ubuntu C++ compiler to recognize waitpid
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
typedef int IndexType;

// Constants ***********************************************************************************************************

static const double inf_ = 1e100;          // pseudo-infinity value for optimization
static const double rt_lower_bound = -0.01; // Technically zero, but need wiggle room for "rebound"

/*
 * "Rebound"
 * 
 *  Run the following in R to see why we need a little wiggle room:
 *  
 *  x <- c(1:100)
 *  Rt <- c(0, 0.122661, 0, 0)
 *  tslope <- c(0.25157, 1.3417, 4.02511)
 *  tpoint <- c(36.8081, 44.7962, 50.7942) 
 *  f_pw <- -0.090846
 *  f_rw <- 0.069678
 *  f_sw <- 0.079237
 * 
 *  Rt_w <- wisp.warp(Rt, 1e4, f_rw)
 *  tslope_w <- wisp.warp(tslope, 1e4, f_sw)
 *  tpoint_w <- wisp.warp(tpoint, 100, f_pw)
 *  y_w <- wisp.sigmoid(x, Rt_w, tslope_w, tpoint_w)
 *  y <- wisp.sigmoid(x, Rt, tslope, tpoint)
 *  plot(x, y_w, type='l', lwd=2, col='blue', xlab='x', ylab='y', main='Warped Sigmoid Function')
 *  lines(x, y, col='black', lwd=2)
 *  abline(v=tpoint_w, col='red', lty=2)
 *  abline(h=0, col='gray', lty=2)
 *  min(y)
 *  
 */

// Main class **********************************************************************************************************

/*
 * Object class to hold and fit Warped Sigmoidal Poisson-Process Mixed-Effect Model (WSPmm, aka "WiSP") model. 
 */

class wspc {
  
  public:
    
    // Fields **********************************************************************************************************
    
    // Data sizes
    int n_count_rows;                       // number of rows in summed count data frame
    int treatment_num;                      // number of treatment combinations (including the reference, i.e. no treatment)
    sdouble bin_num;                        // max bin (i.e., number of bins)
    IntegerVector count_row_nums;           // sequence of integers from 0 to n_count_rows - 1
    
    // Data variables, stan
    sVec bin;                               // bin column of summed data
    sVec count;                             // count column of summed data
    sVec count_log;                         // log of observed counts
    sVec count_tokenized;                   // tokenized count data
    
    // Data variables, Rcpp ("factors")
    CharacterVector parent;                 // parent column of summed data
    CharacterVector child;                  // child column of summed data
    CharacterVector ran;                    // random effect column of summed data
    CharacterVector treatment;              // treatment column of summed data
    
    // Model predictions, Rcpp 
    NumericVector predicted_rates_log;      // log of values predicted by model
    NumericVector predicted_rates;          // values predicted by model
    
    // Fixed and random effect variables
    CharacterVector fix_names;              // names of fixed effect variables
    CharacterVector wfactors_names = {
      "point", 
      "rate", 
      "slope"
      };
    std::vector<CharacterVector> fix_lvls;  // levels for each fixed effect
    std::vector<CharacterVector> fix_trt;                // treatment levels for each fixed effect
    std::vector<CharacterVector> treatment_components;   // all possible treatment combinations, level components
    CharacterVector treatment_lvls;                      // all possible treatment combinations, levels as single-string name
    CharacterVector fix_ref;                // reference level for each fixed effect
    
    // Grouping variables
    CharacterVector parent_lvls;            // levels of parent grouping variable (fixed-effects)
    CharacterVector child_lvls;             // levels of child grouping variable (fixed-effects)
    CharacterVector ran_lvls;               // levels of random effect grouping variable
    
    // Variables to help with data manipulation
    sMat weights;                           // weight matrix, rows as rows of summed count data, columns as treatments (first column is reference)
    sMat weight_rows;                       // ... for making weights matrix and initial effects estimates
    IntegerVector idx_mc_unique;            // count data rows at which model component values will change
    std::vector<IntegerVector> token_pool;               // list of token count indexes associated with each summed count row
    std::vector<IntegerVector> extrapolation_pool;       // list of summed-count indexes giving summed count rows from which to extrapolate
    IntegerVector count_not_na_idx;                      // indexes of non-NA rows in summed count data
    LogicalVector count_not_na_mask;        // mask of non-NA rows in summed count data
    
    // Variables related to model parameters
    CharacterVector mc_list = {             // list of model components
      "Rt",                                 // ... rates
      "tslope",                             // ... transition slopes
      "tpoint"                              // ... transition points
      }; 
    IntegerMatrix degMat;                   // matrix of degrees for each parent (column) -- child (rows) pair
    List found_cp_list;                     // list of found change points (IntegerMatrix) for each parent-child combination
    List found_cp_trt_list;                 // list of found change points (NumericMatrix) for each parent-child combination, averaged across treatments
    List count_log_avg_mat_list;            // list of average log counts (NumericMatrix) for each parent-child combination
    NumericVector fitted_parameters;        // vector holding the model parameters
    CharacterVector param_names;            // list holding the names of the model parameters as they appear in fitted_parameters
    sdouble buffer_factor = 0.05;           // scaling factor for buffer value, the minimum distance between transition points 
    sdouble tpoint_buffer;                  // min number of bins between transition points (immutable structural parameter)
    sdouble warp_precision;                 // precision surviving in calculations of warping
    sdouble inf_warp;                       // pseudo-infinity value for warping (representing no warp boundary)
    sVec warp_bounds;                       // warping bounds for each model component
    IntegerVector warp_bounds_idx = IntegerVector::create(0, 1, 2);
    
    // Indices for managing parameters vector
    IntegerVector param_wfactor_point_idx;  // ... indexes of parameter vector for quick access of different kinds of model parameters
    IntegerVector param_wfactor_rate_idx;
    IntegerVector param_wfactor_slope_idx;
    IntegerVector param_beta_Rt_idx;
    IntegerVector param_beta_tslope_idx;
    IntegerVector param_beta_tpoint_idx;
    IntegerVector param_baseline_idx;
    List beta_idx;                          // ... lists giving the structured array indices for named parameters
    List wfactor_idx; 
    IntegerVector gv_ran_idx;               // ... indices (row and column) for random effect arrays 
    IntegerVector gv_fix_idx;  
    
    // Variables for initial degree estimation
    double LROcutoff = 2.0;                 // cutoff (x sd) for likelihood ratio outlier detection
    double LROwindow_factor = 2.0;          // factor for window size in likelihood ratio outlier detection (bigger is bigger window)
    double rise_threshold_factor = 0.8;     // amount of detected rise as fraction of total required to end run
    double min_initialization_slope = 0.25; // minimum slope for initialization of transition slopes
    int ws;                                 // running window size for likelihood ratio outlier detection (high-pass filter)
    
    // Model settings and results 
    List model_settings;
    List optim_results; 
    
    // Variables for optimization
    int max_evals = 1000;                       // max number of evaluations
    double ctol = 5e-6;                         // convergence tolerance
    unsigned int rng_seed = 42u;                // seed for random number generator
    sMat gamma_dispersion;                      // dispersion terms for "kernel" of gamma-Poisson model
    IntegerVector gd_child_idx;                 // indexes of child levels in gamma_dispersion
    IntegerVector gd_parent_idx;                // indexes of parent levels in gamma_dispersion
    
    // Boundary penalty variables
    int boundary_vec_size = 0;                           // number of boundary components
    sVec bp_coefs;                                       // coefficients used to scale boundary penality so it's negligible when far from boundary, and infinity at boundary
    sdouble max_penalty_at_distance_factor = 0.01;       // the max penalty when far from the boundary, as fraction of initial neg_loglik
    
    /*
     * Basic idea of boundary penalty: there will be a number of relevant distances to a boundary. Want to transform
     *  these distances so that they go to infinity as they approach zero, yet while they are 
     *  sufficiently far from zero ("at distance"), they are not too large.
     */
    
    // Methods *********************************************************************************************************
   
    // Constructor
    wspc(Rcpp::DataFrame count_data, Rcpp::List settings, bool verbose);
    // Destructor
    ~wspc();
    // R copy 
    wspc clone() const;
    
    // Clear Stan
    void clear_stan_mem();
    
    // ***** computing predicted values from parameters 
    
    // Compute model component values for rows of summed count data
    sVec compute_warped_mc(
        const String& mc,          // Model component for which to compute values
        const int& r,              // Row of summed count data for which to compute model component 
        const sVec& parameters,    // Parameters to use in computation
        const sdouble& wf          // Warping factor to apply 
    ) const;
    
    // Predict log of rates
    sVec predict_rates(
        const sVec& parameters,
        const bool& all_rows 
    ) const;
    
    // ***** computing objective function (i.e., fitting model and parameter boundary distances)
    
    // Compute weighted total neg_loglik of observations under the given parameters
    sdouble neg_loglik(
        const sVec& parameters
    ) const;
    
    // Compute boundary distances
    sVec boundary_dist(
        const sVec& parameters
    ) const;
    
    // Test for tpoints below the buffer
    bool test_tpoints(
        const sVec& parameters,
        const bool& verbose
    ) const;
    
    // Compute min boundary penalty
    sdouble min_boundary_dist(
        const sVec& parameters
    ) const;
    
    // Wrap neg_min_boundary_dist in form needed for NLopt constraint function
    static double min_boundary_dist_NLopt(
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
    
    // ***** Computing gradients with stan reverse-mode autodiff
    
    // Compute the gradient of the bounded_nll function
    // ... this is the gradient function used in model optimization
    Eigen::VectorXd grad_bounded_nll(
        const sVec& p_
    ) const;
    
    // Compute the gradient of the min_boundary_dist function
    // ... this is the gradient function used in the search for feasible parameters
    Eigen::VectorXd grad_min_boundary_dist(
        const sVec& p_
    ) const;
    
    // ***** Bootstrapping and model fitting, for statistical testing
    
    // Fit model using NLopt
    void fit(
        const bool set_params, 
        const bool verbose
    );
    
    // Fit model to bootstrap resample
    dVec bs_fit(
        int bs_num,                         // A unique number for this resample
        bool clear_stan                     // Recover stan memory at end?
    );
    
    // Fork bootstraps (parallel processing)
    Rcpp::NumericMatrix bs_batch(
        int bs_num_max,                     // Number of bootstraps to perform
        int max_fork,                       // Maximum number of forked processes per batch
        bool verbose
    );
    
    // Markov-chain Monte Carlo (MCMC) simulation
    Rcpp::NumericMatrix MCMC(
        int n_steps,                        // Number of steps to take in random walk
        int neighbor_filter,                // Keep only ever neighbor_filter step
        double step_size,                   // Step size for random walk
        double prior_sd,                    // Standard deviation of prior distribution
        bool verbose
    );
    
    // Resample (demo)
    Rcpp::NumericMatrix resample(
        int n_resample                      // total number of resamples to draw
    );
    
    // ***** Setting parameters
    
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
    
    // ***** export data to R
    Rcpp::List results(); 
    
    // ***** misc and debugging 
    
    // Wrap neg_loglik in form needed for R
    double neg_loglik_debug(
        const dVec& x
    );
    
    // Wrap bounded_nll in form needed for R
    double bounded_nll_debug(
        const dVec& x
    );
    
    // Wrap grad_bounded_nll_debug in form needed for R
    Rcpp::NumericVector grad_bounded_nll_debug(
        const dVec& x 
    ) const;
  
};

// Helper functions, printing ******************************************************************************************

// Function for printing character messages
void vprint_header(const std::string& message, bool verbose);
void vprint_header(const std::string& message);
void vprint(const std::string& message, bool verbose);
void vprint(const std::string& message);
void vprint(const std::string& message, sdouble val);
void vprint(const std::string& message, double val);
void vprint(const std::string& message, int val);
void vprintV(const CharacterVector& message, bool verbose);

// Function for printing vectors
void vprintV(const NumericVector& vec, bool verbose);
// ... overload for IntegerVector
void vprintV(const IntegerVector& vec, bool verbose);
// ... overload for std::vector
void vprintV(const dVec& vec, bool verbose);
// ... overload for sdouble
void vprintV(const sVec& vec, bool verbose);

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
// ... overload, from IntegerVector
NumericVector to_NumVec(const IntegerVector& vec);
// ... convert to NumericMatrix
NumericMatrix to_NumMat(const sMat& mat);

// Convert to IntegerVector
// ... from std::vector with int
IntegerVector to_IntVec(const iVec& vec);

// Convert to std::vector with doubles 
// ... from NumericVector
dVec to_dVec(const NumericVector& vec);
// ... from Eigen::Matrix with sdouble elements
dVec to_dVec(const sVec& vec);

// Convert to std::vector with int
// ... from IntegerVector
iVec to_iVec(const IntegerVector& vec);

// Misc 
// ... convert to Eigen::Matrix with sdouble elements
sMat to_sMat(const IntegerMatrix& mat);
// ... overload
sMat to_sMat(const NumericMatrix& mat);

// Vector operations ***************************************************************************************************

// Means
// -----

// Mean of vector elements
double vmean(const dVec& x);
// ... overload 
sdouble vmean(const sVec& x);
// ... overload 
double vmean(const NumericVector& x);
// ... overload 
int vmean(const iVec& x);

// Mean of vector elements within a range
double vmean_range(const dVec& x, const int& start, const int& end);
// ... overload 
sdouble vmean_range(const sVec& x, const int& start, const int& end);
// ... overload
double vmean_range(const NumericVector& x, const int& start, const int& end);

// Rolling mean
dVec roll_mean(const dVec& series, int filter_ws);

// Variance of vector elements 
sdouble vvar(const sVec& x);

// Standard deviations of vector elements 
double vsd(const dVec& x); 
// ... overload 
sdouble vsd(const sVec& x);
// ... overload 
double vsd(const NumericVector& x); 

// Component-wise operations
// -------------------------

// Multiplication 
dVec vmult(const dVec& x, const dVec& y);
// ... overload
NumericVector vmult(const NumericVector& x, const NumericVector& y);

// Division
dVec vdivide(const dVec& x, const dVec& y);
// ... overload
NumericVector vdivide(const NumericVector& x, const NumericVector& y);

// Addition 
dVec vadd(const dVec& x, const dVec& y);
// ... overload
NumericVector vadd(const NumericVector& x, const NumericVector& y);

// Subtraction
dVec vsubtract(const dVec& x, const dVec& y);
// ... overload 
NumericVector vsubtract(const NumericVector& x, const NumericVector& y);

// Vectorized sqrt
NumericVector vsqrt(const NumericVector& x);

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
// ... overload 
IntegerVector masked_vec(IntegerVector vec, Rcpp::LogicalVector mask);

// For subsetting vectors with Rcpp IntegerVector
sVec idx_vec(sVec vec, Rcpp::IntegerVector idx);
// ... overload 
dVec idx_vec(dVec vec, Rcpp::IntegerVector idx);
// ... overload 
iVec idx_vec(iVec vec, Rcpp::IntegerVector idx);
// ... overload
CharacterVector idx_vec(CharacterVector vec, Rcpp::IntegerVector idx);

// Find matches in character vector
IntegerVector grep_cpp(CharacterVector V, std::string pattern);

// Find matches in string
bool pattern_match(std::string pattern, std::string test);

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
// ... overload
void assign_proxylist(List list, String element, IntegerMatrix assignment);

// Misc
// ____

// Merge two integer vectors 
IntegerVector buffered_merge(
    const IntegerVector& a, 
    const IntegerVector& b,
    const int& buffer
  );

// Return the indices of the elements of vec of which x is between
iVec block_idx(
    const NumericVector& vec,
    const double& x
  );

// Return differences between elements 
IntegerVector vdiff(
    const IntegerVector& x
  );

// Remove row from IntegerMatrix
IntegerMatrix remove_row(
    IntegerMatrix M, 
    int i
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
  const List& tpointEffs,
  const List& tslopeEffs,
  const IntegerMatrix& degs
  );

// Function for converting structured encoding of parameters into a vector
List make_parameter_vector(
  const List& beta, 
  const List& wfactor, 
  const CharacterVector& prt_lvls, 
  const CharacterVector& cld_lvls, 
  const CharacterVector& rn_lvls, 
  const CharacterVector& mc_list, 
  const CharacterVector& treatment_lvls,
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

// Use extrapolation pool to extrapolate counts
sVec extrapolate_none(
  const sVec& count,
  const CharacterVector& ran, 
  const std::vector<IntegerVector>& extrapolation_pool,
  const bool& log_transform
  );

// Model fitting *******************************************************************************************************

// Thread-safe normal distribution function
double safe_rnorm(
    double mean, 
    double sd
  );

// Better normal distribution function, with PCG and Box-Muller
double pcg_rnorm(
    double mean, 
    double sd,
    pcg32& rng
  );

// Log of density of normal distribution centered on zero
sdouble log_dNorm(
    const sdouble& x,              // value to evaluate
    const sdouble& mu,             // mean
    const sdouble& sd              // standard deviation
  );

// ... overload
double log_dNorm(
    const double& x,               // value to evaluate
    const double& mu,              // mean
    const double& sd               // standard deviation
  );

// Log of density of normal distribution, normalized so highest value is 0
double log_dNorm0(
    const double& x,               // value to evaluate
    const double& mu,              // mean
    const double& sd               // standard deviation
  );

// Log of density of beta distribution, assuming equal shape parameters 
sdouble log_dbeta1(
    const sdouble& x,        // value to evaluate
    const sdouble& shape     // shape parameter (same for both)
  );

// Log of density of gamma distribution
sdouble log_dGamma(
    const sdouble& x,              // value to evaluate
    const sdouble& expected_value, // expected value   
    const sdouble& variance        // variance
  );

// Density of gamma distribution
sdouble dGamma(
    const sdouble& x,              // value to evaluate
    const sdouble& expected_value, // expected value   
    const sdouble& variance        // variance
  );

// Log of density of Poisson distribution
sdouble log_dPois(
    const sdouble& x,              // value to evaluate
    const sdouble& lambda          // rate parameter
  );

// Density of Poisson distribution
sdouble dPois(
    const sdouble& x,              // value to evaluate
    const sdouble& lambda          // rate parameter
  );

// Integral of Poisson-Gamma distribution, from 1 to positive infinity
// ... analytic solution
sdouble poisson_gamma_integral(
    sdouble y, 
    sdouble r, 
    sdouble v
  );

// Estimate variation after x -> log(x + 1) transform 
// ... critical (!!) for Gaussian kernel of Poisson distribution
sdouble delta_var_est(
    const sdouble& var, 
    const sdouble& mu
  );

// Formula to calculate gamme dispersion factor from mean and variance of counts
sdouble gamma_dispersion_formula(
    const sdouble& count_pc_mean, // mean of counts for parent-child combination
    const sdouble& count_pc_var   // variance of counts for parent-child combination
  );

// Recomputes gamma_dispersion matrix, but not gd_child_idx and gd_parent_idx
List compute_gamma_dispersion(
    const sVec& count,                             // count data vector (raw, not log)
    const LogicalVector& count_not_na_mask,        // mask for non-NA counts
    const CharacterVector& parent,                 // parent column of summed data
    const CharacterVector& child,                  // child column of summed data
    const CharacterVector& parent_lvls,            // levels of parent grouping variable (fixed-effects)
    const CharacterVector& child_lvls              // levels of child grouping variable (fixed-effects)
  );

// Function to set warping ratios 
sdouble warp_ratio(
    const sdouble& basis,    // parameterizing coordinate to set the warp
    const sdouble& b,        // bound on this value 
    const sdouble& w         // warping parameter
  );

// Warping function for model components 
sdouble warp_mc(
    const sdouble& x,        // value to warp
    const sdouble& b,        // bound on this value 
    const sdouble& w         // warping parameter
  );

// Numerically stable implementation of sigmoid function
sdouble sigmoid_stable(
  const sdouble& x
  );

// Core poly-sigmoid function of the model
sdouble poly_sigmoid(
  const sdouble& b,                // input variable
  const int& deg,                  // degree of the poly-sigmoid, i.e., number of transitions between blocks
  const sVec& Rt,                  // vector containing the height ("rate") of each block
  const sVec& tslope,              // vector containing the slope of each transition between blocks
  const sVec& tpoint               // vector containing the point of each transition in the bin space
  );

// Inverse quadratic ramp function for boundary penalty
sdouble boundary_penalty_transform(
  const sdouble& x,
  const sdouble& a
  );

// Calculate rolling-window negloglik of a series being generated by a given rate, with or without change-point
dVec series_loglik(
  const dVec& series0,             // 1D vector of points for which to take negative log-likelihood of a change-point
  const int& ws,                   // Running window size
  const bool& null                 // If true, compute likelihood of data assuming no transitions; otherwise, assuming transition
  );

// Likelihood ratio outlier change-point detection
IntegerVector LROcp_find(
    const dVec& nll_ratio,         // 1D vector of points to test for change points
    const int& ws,                 // Running window size
    const double& out_mult         // Outlier multiplier
  );

// ... overload
IntegerVector LROcp_find(
    const NumericMatrix& loglik_ratio_mat,     // NumericMatrix of vectors (columns) to test for change points
    const int& ws,                             // Running window size
    const double& out_mult                     // Outlier multiplier
  );

// Compute likelihood ratios of change points for a series
dVec LROcp_logRatio(
    const dVec& series,           // 1D vector of points to test for change points
    const int& ws                 // Running window size
  );

// Likelihood ratio outlier change-point detection, array input and output
IntegerMatrix LROcp_array(
    const sMat& series_array,     // 2D matrix of points to test for change points
    const int& ws,                // Running window size
    const double& out_mult,       // Outlier multiplier
    const double& cp_buffer       // Minimum distance between two change points
  );

// Function to estimate block rate and transition slopes from count series and change points
std::vector<dVec> est_bkRates_tRuns(
    const int& n_blocks,                // number of blocks
    const NumericVector& count_series,  // count series
    const IntegerVector& cp_series,     // found change points
    const double& rise_threshold_factor // amount of detected rise as fraction of total required to end run
  );

// Function to estimate change points
List estimate_change_points(
    const sVec& bin,                               // bin data vector 
    const sVec& count_log,                         // log of count data vector
    const LogicalVector& count_not_na_mask,        // mask for non-NA counts
    const double& cp_buffer,                       // minimum distance between two change points
    const int& ws,                                 // running window size
    const int& bin_num_i,                          // number of bins in the count data
    const double& LROcutoff,                       // cutoff for outlier detection in change-point detection
    const CharacterVector& parent,                 // parent column of summed data
    const CharacterVector& child,                  // child column of summed data
    const CharacterVector& ran,                    // random effect column of summed data
    const CharacterVector& treatment,              // treatment column of summed data
    const CharacterVector& parent_lvls,            // levels of parent grouping variable (fixed-effects)
    const CharacterVector& child_lvls,             // levels of child grouping variable (fixed-effects)
    const CharacterVector& ran_lvls,               // levels of random effect grouping variable (random-effects)
    const CharacterVector& treatment_lvls,                      // all possible treatment combinations, levels as single-string name
    const std::vector<CharacterVector>& treatment_components    // all possible treatment combinations, level components
  );

// Function to take means of count_log 
List find_count_log_means(
    const sVec& bin,                               // bin data vector 
    const sVec& count_log,                         // log of count data vector
    const LogicalVector& count_not_na_mask,        // mask for non-NA counts
    const int& bin_num_i,                          // number of bins in the count data
    const CharacterVector& parent,                 // parent column of summed data
    const CharacterVector& child,                  // child column of summed data
    const CharacterVector& treatment,              // treatment column of summed data
    const CharacterVector& parent_lvls,            // levels of parent grouping variable (fixed-effects)
    const CharacterVector& child_lvls,             // levels of child grouping variable (fixed-effects)
    const CharacterVector& treatment_lvls,                      // all possible treatment combinations, levels as single-string name
    const std::vector<CharacterVector>& treatment_components    // all possible treatment combinations, level components
  );

// Function to estimate change points
List estimate_initial_parameters(
    const sVec& count,                             // count data vector (raw, not log)
    const LogicalVector& count_not_na_mask,        // mask for non-NA counts
    const int& treatment_num,                      // number of treatments
    const double& rise_threshold_factor,           // amount of detected rise as fraction of total required to end run
    const double& min_initialization_slope,        // minimum slope for initial parameterization
    const sMat& weight_rows,                       // ... for making weights matrix and initial effects estimates
    const CharacterVector& mc_list,                // list of model component types to estimate initial parameters for
    const CharacterVector& parent,                 // parent column of summed data
    const CharacterVector& child,                  // child column of summed data
    const CharacterVector& parent_lvls,            // levels of parent grouping variable (fixed-effects)
    const CharacterVector& child_lvls,             // levels of child grouping variable (fixed-effects)
    const IntegerMatrix& degMat,                   // matrix of degrees of each parent-child combination
    const List& found_cp_list,                     // list of found change points (IntegerMatrix) for each parent-child combination
    const List& found_cp_trt_list,                 // list of found change points (NumericMatrix) for each parent-child combination, averaged across treatments
    const List& count_log_avg_mat_list             // list of average log counts (NumericMatrix) for each parent-child combination
  );

// Function to initialize random effect warping factors 
List make_initial_random_effects(
    const CharacterVector& wfactors_names,  // names of warping factors to initialize
    const int& n_ran,                       // number of random effects
    const int& n_child                      // number of child levels
  );

#endif // WSPC_H