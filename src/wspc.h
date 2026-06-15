
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
using eMat = Eigen::MatrixXd;
using eVec = Eigen::VectorXd;
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

// Indexes of parameter vector for quick access of different kinds of model parameters
struct ParamIdx {
  iVec wfactor_point;  
  iVec wfactor_rate;
  iVec wfactor_slope;
  iVec beta_Rt;
  iVec beta_tslope;
  iVec beta_tpoint;
  iVec baseline_Rt;
  iVec baseline_tslope;
  iVec baseline_tpoint;
  iVec baseline;
};

class wspc {
  
  public:
    
    // Fields **********************************************************************************************************
    
    // Data sizes
    int n_count_rows;                       // number of rows in summed count data frame
    int n_treatments;                       // number of treatment combinations (including the reference, i.e. no treatment)
    int n_bin;                              // max bin (i.e., number of bins)
    iVec factor_sizes;                      // lengths of the vectors holding factor levels, in order of summed count construction: <n_ran, n_context, n_species, n_treatments, n_bin>
    
    // Data variables, numeric
    dVec bin;                               // bin column of summed data
    dVec count;                             // count column of summed data
    dVec count_log;                         // log of summed counts
    dVec count_tokenized;                   // token (non-summed) count data
    
    // Data variables, factors
    iVec ran_num;                           // numeric encoding of ran factor, with 0 as reference level
    iVec context_num;                       // numeric encoding of context factor
    iVec species_num;                       // numeric encoding of species factor
    iVec treatment_num;                     // numeric encoding of treatment factor
    
    // Model predictions, Rcpp 
    NumericVector predicted_rates_log;      // log of values predicted by model
    NumericVector predicted_rates;          // values predicted by model
    
    // Fixed and random effect variables
    CharacterVector fix_names;              // names of fixed effect variables
    CharacterVector wfactors_names = {"point", "rate", "slope"};
    std::vector<CharacterVector> fix_lvls;  // levels for each fixed effect
    std::vector<CharacterVector> fix_trt;                // treatment levels for each fixed effect
    std::vector<CharacterVector> treatment_components;   // all possible treatment combinations, level components
    CharacterVector treatment_lvls;                      // all possible treatment combinations, levels as single-string name
    CharacterVector fix_ref;                // reference level for each fixed effect
    CharacterVector trtKO;                  // Names of treatments to "knock out" (remove from model), if any. 
    IntegerVector timeseries_rank;          // integer sequence the elements of which are ranks, element names as elements from the special "timeseries" fixed effect variable, giving order of time series
    
    // Grouping variables
    CharacterVector context_lvls;           // levels of context grouping variable (fixed-effects)
    CharacterVector species_lvls;           // levels of species grouping variable (fixed-effects)
    CharacterVector ran_lvls;               // levels of random effect grouping variable
    
    // Variables to help with data manipulation
    eMat weights;                           // weight matrix, rows as rows of summed count data, columns as treatments (first column is reference)
    eMat weight_rows;                       // ... for making weights matrix and initial effects estimates
    IntegerVector idx_mc_unique;            // count data rows at which model component values will change
    std::vector<iVec> token_pool;           // list of token count indexes associated with each summed count row
    std::vector<iVec> extrapolation_pool;   // list of summed-count indexes giving summed count rows from which to extrapolate
    IntegerVector r_idx;                    // indexes of rows with "none" in count
    IntegerVector count_not_na_idx;         // indexes of non-NA rows in summed count data
    LogicalVector count_not_na_mask;        // mask of non-NA rows in summed count data
    LogicalVector not_none_mask;            // mask of rows in summed count data with "none" in ran
    bool round_none;                        // whether to round extrapolated "none" counts to nearest integer
    
    // Variables related to model parameters
    CharacterVector mc_list = {"Rt", "tslope", "tpoint"}; 
    IntegerMatrix degMat;                   // matrix of degrees for each context (column) -- species (rows) pair
    std::vector<std::vector<IntegerMatrix>> found_cp;       // list of found change points (IntegerMatrix) for each context-species combination
    std::vector<std::vector<NumericMatrix>> found_cp_trt;   // list of found change points (NumericMatrix) for each context-species combination, averaged across treatments
    std::vector<std::vector<NumericMatrix>> count_log_avg_mat;   // list of average log counts (NumericMatrix) for each context-species combination
    NumericVector fitted_parameters;        // vector holding the model parameters
    CharacterVector param_names;            // list holding the names of the model parameters as they appear in fitted_parameters
    double buffer_factor;                   // scaling factor for buffer value, the minimum distance between transition points 
    double tpoint_buffer;                   // min number of bins between transition points (immutable structural parameter)
    double warp_precision;                  // precision surviving in calculations of warping
    double inf_warp;                        // pseudo-infinity value for warping (representing no warp boundary)
    eVec warp_bounds;                       // warping bounds for each model component
    
    // Indices for managing parameters vector
    ParamIdx param_idx; // ... indexes of parameter vector for quick access of different kinds of model parameters
    std::vector<std::vector<std::vector<iVec>>> beta_idx;
    std::vector<iVec> wfactor_idx; 
    
    // Variables for initial degree estimation
    String LRO_cost;                        // Cost measure for LRO_grid_search: AIC or BIC
    NumericMatrix LRO_grid_search_results;  // matrix to hold results of LRO_grid_search, with columns for settings and cost measure
    double LROcutoff;                       // cutoff (x sd) for likelihood ratio outlier detection
    double LROwindow_factor;                // factor for window size in likelihood ratio outlier detection (bigger is bigger window)
    double rise_threshold_factor;           // amount of detected rise as fraction of total required to end run
    double min_initialization_slope = 0.25; // minimum slope for initialization of transition slopes
    int ws;                                 // running window size for likelihood ratio outlier detection (high-pass filter)
    
    // Model settings and results 
    List model_settings;
    
    // Variables for optimization
    int max_evals;                              // max number of evaluations
    double ctol;                                // convergence tolerance
    unsigned int rng_seed;                      // seed for random number generator
    std::vector<dVec> gamma_dispersion;         // dispersion terms for "kernel" of gamma-Poisson model
    double bnll;
    double nll;
    int success_code;
    int num_evals;
    NumericVector bs_times;
    
    // Boundary penalty variables
    int boundary_vec_size;                               // number of boundary components
    eVec bp_coefs;                                       // coefficients used to scale boundary penality so it's negligible when far from boundary, and infinity at boundary
    double max_penalty_at_distance_factor;               // the max penalty when far from the boundary, as fraction of initial nll
    
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
    
    // ***** Initial setup
    
    // Get row indices of count data from the given factor levels
    iVec fetch_count_idx(iVec I);
    
    // Computes gamma_dispersion matrix and related vectors
    void compute_gamma_dispersion();
    
    // Find row numbers ("pool") of observed counts to use for extrapolation of "none's"
    void make_extrapolation_pool(bool verbose);
    
    // Use extrapolation pool to extrapolate counts
    void extrapolate_none();
    
    // Function to take means of count_log 
    void find_count_log_means();
    
    // Function to estimate change points
    void estimate_change_points();
    
    // Estimate change points and initial parameters for model fitting
    void LRO_initial_param_ests(
      bool verbose = false,
      double LROwf = 0.0,
      double LROco = 0.0
    );
    
    // Search for best LRO change-point detection settings
    void LRO_grid_search(bool verbose);
    
    // ***** computing predicted values from parameters 
    
    // Compute model component values for rows of summed count data
    sVec compute_warped_mc(
        int mc,                    // Model component for which to compute values
        int r,                     // Row of summed count data for which to compute model component 
        const sVec& parameters,    // Parameters to use in computation
        const sdouble& wf          // Warping factor to apply 
    ) const;
    
    // Predict log of rates
    sVec predict_rates(
        const sVec& parameters,
        const bool& all_rows 
    ) const;
    
    // Predict log of rates, R wrapper 
    NumericVector predict_rates_R(
        const NumericVector& parameters_R,
        const bool& all_rows 
    );
    
    // ***** computing objective function (i.e., fitting model and parameter boundary distances)
    
    // Compute weighted total nll of observations under the given parameters
    sdouble compute_nll(
        const sVec& parameters
    ) const;
    
    // Compute boundary distances
    sVec boundary_dist(
        const sVec& parameters
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
    
    // Compute nll plus boundary penalty (main objective function) 
    sdouble compute_bnll(
        const sVec& parameters
    ) const;
    
    // Wrap compute_bnll in form needed for NLopt objective function
    static double compute_bnll_NLopt(
        const dVec& x, 
        dVec& grad, 
        void* data
    );
    
    // ***** Computing gradients with stan reverse-mode autodiff
    
    // Compute the gradient of the compute_bnll function
    // ... this is the gradient function used in model optimization
    Eigen::VectorXd grad_compute_bnll(
        const sVec& p_
    ) const;
    
    // Compute the gradient of the min_boundary_dist function
    // ... this is the gradient function used in the search for feasible parameters
    Eigen::VectorXd grad_min_boundary_dist(
        const sVec& p_
    ) const;
    
    // ***** Bootstrapping and model fitting, for statistical testing
    
    // Fit model using NLopt
    void fit(const bool verbose);
    
    // Resample counts for bootstrapping
    void resample(unsigned int seed);
    
    // Fit model to bootstrap resample
    dVec bs_fit(
        int bs_num,                         // A unique number for this resample
        bool clear_stan                     // Recover stan memory at end?
    );
    
    // Fork bootstraps (parallel processing)
    Rcpp::NumericMatrix bs_batch(
        int bs_num_max,           // Number of bootstraps to perform
        int max_fork,             // Maximum number of forked processes per batch
        bool use_median,
        bool verbose
    );
    
    // Markov-chain Monte Carlo (MCMC) simulation
    Rcpp::NumericMatrix MCMC(
        int n_burnin,             // Number of initial steps to take (to find optimal parameters)
        int n_steps,              // Number of steps to take in random walk (post-burnin)
        int neighbor_filter,      // Keep only ever neighbor_filter step
        double step_size,         // Step size for random walk
        double prior_sd,          // standard deviation to use in prior
        bool start_from_fit,      // Start from parameters found with gradient descent? 
        bool use_median,
        bool verbose
    );
    
    // ***** Setting parameters
    
    // Estimate initial parameters for model fitting
    void estimate_initial_parameters();
    
    // Check and correct parameter feasibility
    void check_parameter_feasibility(const sVec& parameters_var);
    
    // ***** export data to R
    Rcpp::List results(); 
    
    // ***** misc and debugging 
    
    // Wrap compute_nll in form needed for R
    double compute_nll_debug(
        const dVec& x
    );
    
    // Wrap compute_bnll in form needed for R
    double compute_bnll_debug(
        const dVec& x
    );
    
    // Wrap grad_compute_bnll_debug in form needed for R
    Rcpp::NumericVector grad_compute_bnll_debug(
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

// Convert to Eigen::Matrix with double elements
// ... from NumericVector
eVec to_eVec(const NumericVector& vec);
// ... overload, from std::vector with doubles 
eVec to_eVec(const dVec& vec);
// ... overload, from std::vector with int
eVec to_eVec(const iVec& vec);
// ... overload, from IntegerVector
eVec to_eVec(const IntegerVector& vec);
// ... overload, from sVec
eVec to_eVec(const sVec& vec);

// Convert to NumericVector 
// ... from Eigen::Matrix with sdouble elements
NumericVector to_NumVec(const sVec& vec);
// ... overload, from std::vector with doubles
NumericVector to_NumVec(const dVec& vec);
// ... overload, from IntegerVector
NumericVector to_NumVec(const IntegerVector& vec);
// ... convert to NumericMatrix
NumericMatrix to_NumMat(const sMat& mat);
// ... overload, from vector<dVec>
NumericMatrix to_NumMat(const std::vector<dVec>& mat);

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
// ... convert to Eigen::Matrix with double elements
eMat to_eMat(const IntegerMatrix& mat);
// ... overload
eMat to_eMat(const NumericMatrix& mat);

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
// ... overload 
double vmean(const eVec& x);

// Mean of vector elements within a range
double vmean_range(const dVec& x, const int& start, const int& end);
// ... overload 
sdouble vmean_range(const sVec& x, const int& start, const int& end);
// ... overload
double vmean_range(const NumericVector& x, const int& start, const int& end);

// Rolling mean
dVec roll_mean(const dVec& series, int filter_ws);

// Variance of vector elements 
double vvar(const dVec& x);
// ... overload 
sdouble vvar(const sVec& x);
// ... overload 
double vvar(const eVec& x);

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
LogicalVector eq_left_broadcast(const iVec& left, const int& right);
// ... overload
LogicalVector eq_left_broadcast(const NumericVector& left, const double& right);

// Quantifiers
bool any_true(const LogicalVector& x);
bool all_true(const LogicalVector& x);

// Treatment conditions ************************************************************************************************

// Generate all combinations of j indices from {0, ..., n-1}
void combinations(
  int n, int j, int start, 
  iVec& current, 
  std::vector<iVec>& result
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

// Model fitting *******************************************************************************************************

// Random draw from uniform distribution, integers
int rUnifi(
    int min,
    int max,
    std::mt19937& rng
  );

// Random draw from uniform distribution, reals
double rUnifr(
    double min,
    double max,
    std::mt19937& rng
  );

// Random draw from normal distribution
double rNorm(
    double mean, 
    double sd,
    std::mt19937& rng
  );

// Random draw from normal distribution, thread-safe
double safe_rNorm(
    double mean, 
    double sd
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
double gamma_dispersion_formula(
    const double& count_cs_mean, // mean of counts for context-species combination
    const double& count_cs_var   // variance of counts for context-species combination
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
    const double& b,         // bound on this value 
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
  const double& a
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
    const eMat& series_array,     // 2D matrix of points to test for change points
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

#endif // WSPC_H
