
// model_fit.cpp
#include "wspc.h"

// Empirical cdf 
sdouble emp_cdf(
    const sdouble& val,    // value to evaluate
    const sVec& vec        // vector of values
  ) {
    // Count the number of elements less than or equal to vec
    int count = 0;
    for (int i = 0; i < vec.size(); i++) {
      if (vec(i) <= val) {
        count++;
      }
    }
    // Return the empirical cdf value
    return sdouble(count) / sdouble(vec.size());
  }

// Empirical pdf from normalized sigmoid fit to empirical cdf
sdouble emp_pdf(
    const sdouble& val,     // value to evaluate
    const sVec& ecdf_params // empirical cdf parameters
  ) {
    /*
     * Basic idea: First fit to an empirical cdf a normalized sigmoid of the form: 
     *    F_X(u) = N(G(u)), for 
     *    G(u) = h + 1/(1+exp(-s*(u+p)))
     *    N(v) = (v - min(v))/(max(v) - min(v))
     * The pdf of this cdf is dF_X/du, which is: 
     *    dF_X/du = (1/(max(v) - min(v))) * (s * exp(-s*(u+p)))/(1 + exp(-s*(u+p)))^2 
     */
    sdouble min_val = ecdf_params(3);
    sdouble max_val = ecdf_params(4);
    sdouble p_val = ecdf_params(0); 
    sdouble s_val = ecdf_params(2); 
    // ... h is the second element of ecdf_params, but not needed 
    sdouble min_max = 1/(max_val - min_val); 
    sdouble exp_val = sexp(-s_val * (val + p_val));
    sdouble dGdu = (s_val * exp_val) / spower(1 + exp_val, 2.0); 
    sdouble pdf = min_max * dGdu;
    return pdf;
  }

// Thread-safe normal distribution function
double safe_rnorm(
    double mean, 
    double sd
  ) {
    // Create a thread-local random engine
    static thread_local std::mt19937 generator(std::random_device{}());
    
    // Create a normal distribution with the given mean and standard deviation
    std::normal_distribution<double> distribution(mean, sd);
    
    // Generate and return a random number
    return distribution(generator);
  }

// Density of normal distribution 
double dNorm(
    const double& x,        // value to evaluate
    const double& mu,       // mean (expected value)
    const double& sd        // standard deviation
  ) {
    return std::exp(-((x - mu) * (x - mu)) / (2.0 * sd * sd)) / (std::sqrt(2.0 * M_PI) * sd);
  }

// Log of density of normal distribution
sdouble log_dNorm(
    const sdouble& x,        // value to evaluate
    const sdouble& mu,       // mean (expected value)
    const sdouble& sd        // standard deviation
  ) {
    /*
     * Have: log_dNorm = log(exp(-((x - mu)^2) / (2 * sd^2)) / (sqrt(2 * pi) * sd))
     * ... since log(x/y) = log(x) - log(y), and log(exp(x)) = x, then:
     * log_dNorm = -((x - mu)^2) / (2 * sd^2) - log(sqrt(2 * pi) * sd)
     */
    return (-spower(x - mu, 2.0) / (2.0 * spower(sd, 2.0))) - slog(2.0 * M_PI * spower(sd, 2.0))/2.0;
  }

// ... overload
double log_dNorm(
    const double& x,        // value to evaluate
    const double& mu,       // mean (expected value)
    const double& sd        // standard deviation
  ) {
    return (-((x - mu) * (x - mu)) / (2.0 * sd * sd)) - std::log(2.0 * M_PI * sd * sd)/2.0;
  }

// Log of density of normal distribution, normalized so highest value is 0
double log_dNorm0(
    const double& x,        // value to evaluate
    const double& mu,       // mean (expected value)
    const double& sd        // standard deviation
  ) {
    return -((x - mu) * (x - mu)) / (2.0 * sd * sd);
  }

// Log of density of beta distribution, assuming equal shape parameters 
sdouble log_dbeta1(
    const sdouble& x,        // value to evaluate
    const sdouble& shape     // shape parameter (same for both)
  ) {
    return slog((spower(x, shape - 1.0) * spower(1.0 - x, shape - 1.0)) / (spower(stan::math::tgamma(shape), 2.0) / stan::math::tgamma(2.0 * shape)));
  }

// Log of density of gamma distribution
sdouble log_dGamma(
    const sdouble& x,              // value to evaluate
    const sdouble& expected_value, // expected value   
    const sdouble& variance        // variance
  ) {
    // Have: 
    // rate = shape / expected_value;
    // shape = variance / (rate * rate);
    sdouble shape = (expected_value * expected_value) / variance;
    sdouble rate = shape / expected_value;
    return slog(
      (spower(rate, shape) * spower(x, shape - 1.0) * sexp(-rate * x)) / stan::math::tgamma(shape)
    );
  }

// Density of gamma distribution
sdouble dGamma(
    const sdouble& x,              // value to evaluate
    const sdouble& expected_value, // expected value   
    const sdouble& variance        // variance
  ) {
    // Have: 
    // rate = shape / expected_value;
    // shape = variance / (rate * rate);
    sdouble shape = (expected_value * expected_value) / variance;
    sdouble rate = shape / expected_value;
    return (
      (spower(rate, shape) * spower(x, shape - 1.0) * sexp(-rate * x)) / stan::math::tgamma(shape)
    );
  }

// Log of density of Poisson distribution
sdouble log_dPois(
    const sdouble& x,        // value to evaluate
    const sdouble& lambda    // rate parameter
  ) {
    return x * slog(lambda) - lambda - stan::math::lgamma(x + 1.0);
    // ^ Hand-written function for log Poisson density
    //  ... could use stan implementation: stan::math::poisson_lpmf(count_log(r), pred_rate_log);
    //  ... but seems identical in results and speed?
  }

// Density of Poisson distribution
sdouble dPois(
    const sdouble& x,        // value to evaluate
    const sdouble& lambda    // rate parameter
  ) {
    return sexp(log_dPois(x, lambda));
  }

// Integral of Poisson-Gamma distribution, from 1 to positive infinity
sdouble poisson_gamma_integral(
    sdouble y, // observed count value
    sdouble r, // expected process rate drawn from the gamma distribution
    sdouble v  // variance of the gamma distribution
  ) {
   
    // convert more intuitive rate and variance parameters into the standard shape-rate parameters 
    // ... for the gamma distribution
    sdouble s = r * r / v; // Gamma distribution "shape" parameter
    sdouble R = s / r;     // Gamma distribution "rate" parameter
    
    // Idea: This is an analytic solution to the integral of dPois(y, lambda) * dGamma(lambda, r, v) from 1 to positive infinity. 
    //        ... The solution takes the form of a ratio num/denom which is subject to overflow/underflow, and so instead of 
    //            computing this ratio directly, we compute the log of the numerator and denominator separately, and then exponentiate the difference.
    sdouble log_num = s * slog(R) + stan::math::lgamma(y + s) + stan::math::log1m(stan::math::gamma_p(y + s, R + 1.0));
    sdouble log_denom = stan::math::lgamma(y + 1.0) + stan::math::lgamma(s) + (y + s) * slog(R + 1.0);
    sdouble integral = sexp(log_num - log_denom);
   
    // note: this return value is *not* the log of the density! It's the density itself. 
    return integral;
    
  }

// Estimate variation after x -> log(x + 1) transform 
// ... critical (!!) for Gaussian kernel of Poisson distribution
sdouble delta_var_est(
    const sdouble& var,
    const sdouble& mu
  ) {
    // For a random variable X with variance var and mean mu, 
    // ... estimate the variance of log(X + 1) using the delta method 
    //      applied to the first-order Taylor series of g(x) = log(x + 1) around mu.
    return var / ((mu + 1) * (mu + 1));
  }

// Generic warping function with bound 
sdouble warp_ratio(
    const sdouble& basis,    // parameterizing coordinate to set the warp
    const sdouble& b,        // bound on this value 
    const sdouble& w         // warping parameter
  ) {
    return(1.0 - spower(sexp(w*w), -basis/b));
  }

// Warping function for model components 
sdouble warp_mc(
    const sdouble& s,        // value to warp
    const sdouble& b,        // bound on this value 
    const sdouble& w         // warping parameter
  ) {
    sdouble basis;           // parameterizing coordinate to set the warp
    sdouble magnitude;       // value by which to warp
    sdouble direction;       // direction of the warp
    if (w > 0.0) {
      basis = s;
      magnitude = b - s;
      direction = 1.0;
    } else {
      basis = b - s;
      magnitude = s;
      direction = -1.0;
    }
    return(s + direction * warp_ratio(basis, b, w) * magnitude);
  }

// Wrapper for R
// [[Rcpp::export]]
double warp_mc_R(
    const double& s,        // value to warp
    const double& b,        // bound on this value 
    const double& w         // warping parameter
  ) {
    return(warp_mc(sdouble(s), sdouble(b), sdouble(w)).val());
  }

// Numerically stable implementation of sigmoid function
sdouble sigmoid_stable(
    const sdouble& x
  ) {
    // Implementation of sigmoid function that's numerically stable
    if (x.val() > 0) {
      return 1.0 / (1.0 + sexp(-x));
    } else {
      sdouble exp_x = sexp(x);
      return exp_x / (1.0 + exp_x);
    }
  }

// Wrapper for R
// [[Rcpp::export]]
double sigmoid_stable_R(
    const double& x
  ) {
    return(sigmoid_stable(sdouble(x)).val());
  }

// Core poly-sigmoid function of the model
sdouble poly_sigmoid(
    const sdouble& b,        // input variable
    const int& deg,          // degree of the poly-sigmoid, i.e., number of transitions between blocks
    const sVec& Rt,          // vector containing the height ("rate") of each block
    const sVec& tslope,      // vector containing the slope of each transition between blocks
    const sVec& tpoint       // vector containing the point of each transition in the bin space
  ) {
    
    /*
    * Core function of the WSPmm framework. 
    * ... This function parameterizes the 1D distribution of rates which 
    *      stay constant in "blocks" separated by a fixed number of transition
    *      points. The number of transition points is determined by the "degree"
    *      of the function, which the degree is the number of sigmoid functions
    *      summed together. 
    *      
    *  Model components: 
    *   - Rt: the rate at each block
    *   - tslope: the slope of the transition between blocks (not literal rise over run, but a scaler of sorts)
    *   - tpoint: the point in the bin space where the transition occurs
    */
    
    if ( deg == 0 ) {
      return Rt(0);
    } else {
      sdouble ps = Rt(0);
      for (int i = 0; i < deg; i++) {
        if (std::isinf(Rt(i+1))) {
          return Rt(i+1);
        } else {
          sdouble sig_term = tslope(i)*(b - tpoint(i));
          ps += (Rt(i+1) - Rt(i)) * sigmoid_stable(sig_term);
        }
      } 
      return ps;
    } 
    
  }

// Wrapper for R
// [[Rcpp::export]]
double poly_sigmoid_R(
    const double& b,            // input variable
    const int& deg,             // degree of the poly-sigmoid, i.e., number of transitions between blocks
    const NumericVector& Rt,           // vector containing the height ("rate") of each block
    const NumericVector& tslope,       // vector containing the slope of each transition between blocks
    const NumericVector& tpoint        // vector containing the point of each transition in the bin space
  ) {
    
    return poly_sigmoid(
      sdouble(b), 
      deg, 
      to_sVec(Rt), 
      to_sVec(tslope), 
      to_sVec(tpoint)
    ).val();
    
  }

// Inverse quadratic ramp function for boundary penalty
sdouble boundary_penalty_transform(
    const sdouble& x,
    const sdouble& a
  ) {
    sdouble x_ = x;
    if (x_ <= 0) {
      x_ = 1e-30;
    } 
    return 1.0 / spower(a*x_, 2.0);
  }

// Calculate rolling-window log-likelihood of a series being generated by a given rate, with or without change-point
dVec series_loglik(
    const dVec& series0,          // 1D vector of points for which to take likelihood of a change-point
    const int& ws,                // Running window size
    const bool& null              // If true, compute likelihood of data assuming no transitions; otherwise, assuming transition
  ) {
    // Assumes points are normally distributed about their process rate
    
    int n = series0.size();
    dVec series = series0; 
    dVec loglik(n); // length out should equal length in
    
    // Compute joint log-likelihood of observations within each window
    for (int i = 0; i <= (n - ws); i++) {
     
      loglik[i] = 0.0;
      // ... series values in window
      dVec series_i(series.begin() + i, series.begin() + (i + ws));
      // ... series values in first half of window
      dVec series_i1(series_i.begin(), series_i.begin() + (ws / 2));
      // ... series values in second half of window
      dVec series_i2(series_i.begin() + (ws / 2), series_i.end());
      
      if (null) {
        double mean_i = vmean(series_i);
        double sd_i = vsd(series_i);
        for (int j = 0; j < ws; j++) {
          loglik[i] += log_dNorm(series_i[j], mean_i, sd_i);
        }
      } else {
        double mean_i1 = vmean(series_i1);
        double sd_i1 = vsd(series_i1);
        for (int j = 0; j < ws/2; j++) {
          loglik[i] += log_dNorm(series_i[j], mean_i1, sd_i1);
        }
        double mean_i2 = vmean(series_i2);
        double sd_i2 = vsd(series_i2);
        for (int j = ws/2; j < ws; j++) {
          loglik[i] += log_dNorm(series_i[j], mean_i2, sd_i2);
        }
      }
      
      // Check for nan's and replace with 1
      if (std::isnan(loglik[i])) {loglik[i] = 1.0;}
     
    }
   
    // Fill in the last part of the vector with the mean
    double loglik_mean = vmean(loglik);
    for (int i = n - ws + 1; i < n; i++) {
      loglik[i] = loglik_mean;
    }
    
    return loglik;
    
  }

// Likelihood ratio outlier change-point detection
IntegerVector LROcp_find(
    const dVec& loglik_ratio,     // 1D vector of points to test for change points
    const int& ws,                // Running window size
    const double& out_mult        // Outlier multiplier
  ) {
    
    //Find change points
    iVec fcp;
    int last_cp = 0;
    double loglik_ratio_mean = vmean(loglik_ratio);
    double loglik_ratio_sd = vsd(loglik_ratio);
    double loglik_max = loglik_ratio_mean + out_mult * loglik_ratio_sd;
    int n_nll = loglik_ratio.size();
    for (int i = 0; i < n_nll; i++) {
      if (
          loglik_ratio[i] > loglik_max && 
            i < (loglik_ratio.size() - ws/2) &&   // can't be too close to end
            i > ws/2                              // can't be too close to beginning
      ) {
        if (i - last_cp > int(ws/2)) {            // can't be too close to last change point
          fcp.push_back(i);
          last_cp = i;
        } else if (loglik_ratio[i] > loglik_ratio[last_cp]) {
          // If too close to last change point, check if this ratio is larger
          fcp[fcp.size() - 1] = i;
        }
      }
    }
    
    // Make found_cp vector
    IntegerVector found_cp;
    if (fcp.size() > 0) {
      // Grab the indexes of the found change points and add 1 (first bin is 1)
      found_cp = to_IntVec(fcp) + 1;
    } 
    
    return found_cp;
    
  } 

// Compute likelihood ratios of change points for a series
dVec LROcp_logRatio(
    const dVec& series,           // 1D vector of points to test for change points
    const int& ws                 // Running window size
  ) {
    
    // Compute the null likelihood at each bin (no change point)
    dVec loglik_null = series_loglik(series, ws, true);
    
    // Compute the change-point likelihood at each bin
    dVec loglik_cp = series_loglik(series, ws, false);
    
    // Find the loglik ratio at each bin 
    // ... Expect lik_null to be smaller than lik_cp (i.e., less likely, a worse fit), 
    // ... so expect loglik_null to be smaller (e.g., a negative number of greater magnitude) 
    // ... than loglik_cp. So, loglik_cp - loglik_null should be positive, with a minimum of 
    // ... zero, i.e., the log of the ratio 1/1, which is 0. Importantly, the *larger* this 
    // ... number, the *better* the cp hypothesis fits the data. So, a larger ratio is still 
    // ... better, even though we're in log space. 
    dVec loglik_ratio = vsubtract(loglik_cp, loglik_null);
    
    return loglik_ratio;
    
  }

// Likelihood ratio outlier change-point detection, arrays
IntegerMatrix LROcp_array(
    const sMat& series_array,     // 2D matrix of points to test for change points
    const int& ws,                // Running window size
    const double& out_mult        // Outlier multiplier
  ) {
    // The test here will treat each column of the matrix as a separate series. 
    //  it further assumes that the series are dependent on each other, and share 
    //  change points, plus some variance between them on the exact location. 
    
    // Control random number generator
    Rcpp::RNGScope scope;
    
    int n_series = series_array.cols();
    int n_samples = series_array.rows();
    NumericMatrix loglik_ratio_array(n_samples, n_series);
    
    // Convert series_array of raw (or log) count values into an array of likelihood ratios
    for (int s = 0; s < n_series; s++) {
      
      // Get the series for this column
      dVec series = to_dVec(series_array.col(s));
      
      // Find the likelihood ratios of change points for this series
      dVec loglik_ratio = LROcp_logRatio(series, ws);
      
      // Save in array 
      loglik_ratio_array.column(s) = to_NumVec(loglik_ratio);
      
    }
    
    // Find centroid
    NumericVector centroid(n_samples);
    for (int i = 0; i < n_samples; i++) {
      centroid[i] = vmean(loglik_ratio_array.row(i));
    }
    
    // Use LROcp on this series to estimate change points 
    IntegerVector found_cp = LROcp_find(to_dVec(centroid), ws, out_mult);
    
    if (found_cp.size() > 0) {
      
      // Project found_cp back to original series 
      // ... these will be one-based indices that need to be subtracted back to zero-based
      Function project_cp("project_cp");
      IntegerMatrix found_cp_array = project_cp(found_cp, centroid, loglik_ratio_array);
      // ... ^ in this matrix, rows are change points (by deg), columns are trt x ran interactions
     
      for (int i = 0; i < found_cp_array.ncol(); i++) {
        found_cp_array.column(i) = found_cp_array.column(i) - 1;
        for (int k = 0; k < found_cp_array.nrow(); k++) {
          if (found_cp_array(k, i) < ws) {found_cp_array(k, i) = ws;}
          if (found_cp_array(k, i) > n_samples - ws) {found_cp_array(k, i) = n_samples - ws;}
        }
      }
      return found_cp_array;
      
    } else {
      // If no change points found, return empty matrix
      return IntegerMatrix(0, 0);
    }
    
  } 

// Function to estimate block rate and transition slopes from count series and change points
std::vector<dVec> est_bkRates_tRuns(
    const int& n_blocks,                // number of blocks
    const NumericVector& count_series,  // count series
    const IntegerVector& cp_series,     // found change points
    const double& rise_threshold_factor // amount of detected rise as fraction of total required to end run
  ) {
    // Grab bin number and series length
    int deg = n_blocks - 1;
    int bin_num_i = count_series.size();
    // Initialize vectors for block rate and run estimates
    dVec Rt_est(n_blocks); 
    dVec run_estimates(n_blocks - 1);
    NumericVector block_sizes(n_blocks);
    for (int bk = 0; bk < n_blocks; bk++) {
      // ... find bins which bound this block
      int bin_start;
      int bin_end;
      if (bk == 0) {
        bin_start = 0;
        bin_end = cp_series[bk];
      } else if (bk == deg) {
        bin_start = cp_series[bk - 1] - 1;
        bin_end = bin_num_i - 1;
      } else {
        bin_start = cp_series[bk - 1] - 1;
        bin_end = cp_series[bk];
      }
      block_sizes[bk] = bin_end - bin_start;
      // ... take mean of count values in this range
      Rt_est[bk] = vmean_range(count_series, bin_start, bin_end);
      // ... estimate run, number of bins needed to complete rate transition
      if (bk > 0) {
        int bb = 0; // count of bins stepped back
        int bf = 0; // count of bins stepped forward
        double rise_found = 0.0;
        double rise_threshold = (Rt_est[bk] - Rt_est[bk - 1]) * rise_threshold_factor;
        double sign = 1.0;
        if (rise_threshold < 0) {sign = -1.0;}
        bool go_forward = true;
        while (
            rise_found * sign < rise_threshold * sign && 
            (bb < block_sizes[bk - 1] || bf < block_sizes[bk]) 
        ) {
          double Rt_back = vmean_range(count_series, bin_start - bb, bin_start);
          double Rt_forward = vmean_range(count_series, bin_start, bin_start + bf);
          rise_found = Rt_forward - Rt_back;
          if (go_forward) {
            if (bf < block_sizes[bk]) {bf++;}
          } else {
            if (bb < block_sizes[bk - 1]) {bb++;}
          }
          go_forward = !go_forward;
        }
        int run_length = bb + bf;
        if (run_length == 0) {run_length++;} // ... ensure at least one value in run
        run_estimates[bk - 1] = (double)run_length;
      }
    }
    std::vector<dVec> out(2); 
    out[0] = Rt_est;
    out[1] = run_estimates;
    return(out);
    
  }

// Function to estimate the ratio (FixEff + RanEff) / RanEff 
sdouble effect_ratio_est(
    const sdouble& bl, // Expected baseline spatial parameter value
    const sdouble& fe, // Expected fixed effect parameter value
    const sdouble& re, // Expected random effect parameter value
    const sdouble& bd  // Warp boundary
  ) {
    sdouble z = bl + fe; 
    sdouble B = bd / (bd - z);
    sdouble Z = 2.0 * z * re;
    return (B + re*re + Z) / Z;
  }