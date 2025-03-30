
// model_fit.cpp
#include "wspc.h"

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

// Calculate rolling-window likelihood of a series being generated by a given rate, with or without change-point
dVec series_lik(
    const dVec& series0,          // 1D vector of points for which to take likelihood of a change-point
    const int& ws,                // Running window size
    const int& filter_ws,         // Size of window for taking rolling mean
    const bool& null              // If true, compute likelihood of data assuming no transitions; otherwise, assuming transition
  ) {
    // Assumes points are normally distributed about their process rate
    
    int n = series0.size();
    dVec series = roll_mean(series0, filter_ws);
    dVec lik(n); // length out should equal length in
    
    // Compute mean likelihood of an observation within each window
    for (int i = 0; i <= (n - ws); i++) {
     
      lik[i] = 0.0;
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
          lik[i] += dNorm(series_i[j], mean_i, sd_i);
        }
      } else {
        double mean_i1 = vmean(series_i1);
        double sd_i1 = vsd(series_i1);
        for (int j = 0; j < ws/2; j++) {
          lik[i] += dNorm(series_i[j], mean_i1, sd_i1);
        }
        double mean_i2 = vmean(series_i2);
        double sd_i2 = vsd(series_i2);
        for (int j = ws/2; j < ws; j++) {
          lik[i] += dNorm(series_i[j], mean_i2, sd_i2);
        }
      }
      
      lik[i] = lik[i] / (double)ws;
      
      // Check for nan's and replace with 1
      if (std::isnan(lik[i])) {lik[i] = 1.0;}
     
    }
   
    // Fill in the last part of the vector with the mean
    double lik_mean = vmean(lik);
    for (int i = n - ws + 1; i < n; i++) {
      lik[i] = lik_mean;
    }
    
    return lik;
    
  }

// Likelihood ratio outlier change-point detection
IntegerVector LROcp_find(
    const dVec& lik_ratio,        // 1D vector of points to test for change points
    const int& ws,                // Running window size
    const double& out_mult        // Outlier multiplier
  ) {
    
    //Find change points
    iVec fcp;
    int last_cp = 0;
    double lik_ratio_mean = vmean(lik_ratio);
    double lik_ratio_sd = vsd(lik_ratio);
    double lik_max = lik_ratio_mean + out_mult * lik_ratio_sd;
    int n_nll = lik_ratio.size();
    for (int i = 0; i < n_nll; i++) {
      if (
          lik_ratio[i] > lik_max && 
            i < (lik_ratio.size() - ws/2) &&  // can't be too close to end
            i > ws/2                          // can't be too close to beginning
      ) {
        if (i - last_cp > ws/2) {             // can't be too close to last change point
          fcp.push_back(i);
          last_cp = i;
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
dVec LROcp_ratio(
    const dVec& series,           // 1D vector of points to test for change points
    const int& ws,                // Running window size
    const int& filter_ws          // Size of window for taking rolling mean
  ) {
    
    // Compute the null likelihood at each bin (no change point)
    dVec lik_null = series_lik(series, ws, filter_ws, true);
    
    // Compute the change-point likelihood at each bin
    dVec lik_cp = series_lik(series, ws, filter_ws, false);
    
    // Find the nll ratio at each bin 
    // ... expected value is 1. Expect lik_null to be 
    // ... smaller than lik_cp (i.e., less likely, a worse fit), 
    // ... so expect this value to be 1 or larger. 
    dVec lik_ratio = vdivide(lik_cp, lik_null);
    
    return lik_ratio;
    
  }

// Likelihood ratio outlier change-point detection, single-series
IntegerVector LROcp(
    const dVec& series,           // 1D vector of points to test for change points
    const int& ws,                // Running window size
    const int& filter_ws,         // Size of window for taking rolling mean
    const double& out_mult        // Outlier multiplier
  ) {
    
    // Find the likelihood ratios of change points for the series
    dVec lik_ratio = LROcp_ratio(series, ws, filter_ws);
    
    // Use LROcp on this series to estimate change points 
    IntegerVector found_cp = LROcp_find(series, ws, out_mult);
    
    return found_cp;
    
  }

// Likelihood ratio outlier change-point detection, arrays
IntegerMatrix LROcp_array(
    const sMat& series_array,     // 2D matrix of points to test for change points
    const int& ws,                // Running window size
    const int& filter_ws,         // Size of window for taking rolling mean
    const double& out_mult        // Outlier multiplier
  ) {
    // The test here will treat each column of the matrix as a separate series. 
    //  it further assumes that the series are dependent on each other, and share 
    //  change points, plus some variance between them on the exact location. 
    
    // Control random number generator
    Rcpp::RNGScope scope;
    
    int n_series = series_array.cols();
    int n_samples = series_array.rows();
    NumericMatrix lik_ratio_array(n_samples, n_series);
    
    // Convert series_array of raw (or log) count values into an array of likelihood ratios
    for (int s = 0; s < n_series; s++) {
      
      // Get the series for this column
      dVec series = to_dVec(series_array.col(s));
      
      // Find the likelihood ratios of change points for this series
      dVec lik_ratio = LROcp_ratio(series, ws, filter_ws);
      
      // Square the likelihood ratios so DTW-DBA algorithm more heavily weights peaks when aligning series
      //NumericVector lik_ratio2 = to_NumVec(vmult(lik_ratio, lik_ratio));
      
      // Save in array 
      lik_ratio_array.column(s) = to_NumVec(lik_ratio);
      
    }
    
    // Find the centroid of lik_ratio_array
    // ... calling function from dtwclust package in R
    //int dtw_window = ws / 2;
    //Function find_centroid("find_centroid");   
    //NumericVector centroid = find_centroid(lik_ratio_array, dtw_window);
    
    NumericVector centroid(n_samples);
    for (int i = 0; i < n_samples; i++) {
      centroid[i] = vmean(lik_ratio_array.row(i));
    }
    
    // Take centroid out of squared space
    //centroid = vsqrt(centroid);
    
    // Use LROcp on this series to estimate change points 
    IntegerVector found_cp = LROcp_find(to_dVec(centroid), ws, out_mult);
    
    if (found_cp.size() > 0) {
      
      // Project found_cp back to original series 
      // ... these will be one-based indices that need to be subtracted back to zero-based
      Function project_cp("project_cp");
      IntegerMatrix found_cp_array = project_cp(found_cp, centroid, lik_ratio_array);
      // ... ^ in this matrix, rows are change points (by deg), columns are trt x ran interactions
     
      for (int i = 0; i < found_cp_array.ncol(); i++) {
        found_cp_array.column(i) = found_cp_array.column(i) - 1;
        for (int k = 0; k < found_cp_array.nrow(); k++) {
          if (found_cp_array(k, i) < filter_ws) {found_cp_array(k, i) = filter_ws;}
          if (found_cp_array(k, i) > n_samples - filter_ws) {found_cp_array(k, i) = n_samples - filter_ws;}
        }
      }
      return found_cp_array;
      
    } else {
      // If no change points found, return empty matrix
      return IntegerMatrix(0, 0);
    }
    
  } 

// Function to derive sd of fixed effect from variance of random effect
sdouble get_beta_sd(
    const sdouble& beta_shape, // shape parameter of the beta distribution
    const sdouble& diff_ratio, // estimated ratio of fixed effect to random effect
    const sdouble& expected_value 
  ) {
    sdouble expected_ran_effect = 1.0 / ssqrt(4.0 * beta_shape + 2.0);       // ... already includes the *= 2.0 rescaling to range from -1 to 1
    sdouble eff_mult = diff_ratio - 1.0;                                     // ... note the different sign
    sdouble sd_beta_effect = (eff_mult * expected_value * expected_ran_effect) / (2.0 + eff_mult * expected_ran_effect);
    return sd_beta_effect;
  }

// Function to estimate block rate and transition slopes from count series and change points
std::vector<dVec> est_bkRates_tRuns(
    const int& n_blocks,                // number of blocks
    const NumericVector& count_series,  // count series
    const IntegerVector& cp_series,     // found change points
    const double& rise_threshold_factor // amount of detected rise as fraction of total required to end run
  ) {
    int deg = n_blocks - 1;
    int bin_num_i = count_series.size();
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
      // ... estimate run, number of blocks needed to complete rate transition
      if (bk > 0) {
        int bb = 0;
        int bf = 0;
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

