
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

// Warping function for random effects on transition points
sdouble point_warp(
    const sdouble& b,        // Bin coordinate
    const sdouble& bin_num,  // Max bin coordinate
    const sdouble& f         // warping factor
  ) {
   
    /*
    * The main poly-sigmoid itself is "warped" via two functions: 
    *  - a point-warp function which warps bin coordinates, and 
    *  - a rate-warp function which warps the rate value output by the poly-sigmoid.
    *  
    *  In terms of modelling, the fixed-effects are built into the poly-sigmoid function, 
    *   while the random effects are modelled via "warping factors" which control these two warping functions. 
    */
   
    return b - f * b * (1.0 - (b / bin_num));
    
  } 

// Warping function for random effects on count rate
sdouble rate_warp(
    const sdouble& rate,     // rate output by the poly-sigmoid
    const sdouble& f         // warping factor
  ) {
    // See description in point-warp function
    sdouble warped_rate;
    if (std::isinf(rate)) {
      warped_rate = rate;
    } else { 
      warped_rate = rate + (f*rate);
    }
    return warped_rate;
  }

// Function for computing model components
sVec compute_mc(
    const sMat& betas,       // effects matrix, rows as interactions, columns as blocks/transitions
    const sVec& weight_row   // row of weight matrix, elements as fixed effects
  ) {
    
    /*
    * This function computes the model components, starting with the baseline value 
    *  and then adding the fixed effects based on the inputted weights. 
    */
    
    int i_max = betas.rows();
    int num_blocks = betas.cols();
    
    // Initialize mc_value as zero
    sVec mc_value(num_blocks);
    mc_value.setZero();
    
    // Calculate value from betas (ref and fixed-effects)
    for (int i = 0; i < i_max; i++) {
      sdouble wi = weight_row(i);
      sVec betai = betas.row(i); // Effects column
      for (int bt = 0; bt < num_blocks; bt++) {
        mc_value(bt) += wi*betai(bt);
      }
    }
    
    // Return updated model-component values 
    return mc_value;
    
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

// Calculate rolling-window negloglik of a series being generated by a given rate, with or without change-point
dVec series_nll(
    const dVec& series0,          // 1D vector of points for which to take negative log-likelihood of a change-point
    const int& ws,                // Running window size
    const int& filter_ws,         // Size of window for taking rolling mean
    const bool& null              // If true, compute likelihood of data assuming no transitions; otherwise, assuming transition
  ) {
    
    int n = series0.size();
    dVec series = roll_mean(series0, filter_ws);
    dVec nlls(n); // length out should equal length in
    
    // Compute nll for each window
    for (int i = 0; i <= (n - ws); i++) {
     
      nlls[i] = 0.0;
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
          nlls[i] += dNorm(series_i[j], mean_i, sd_i);
        }
      } else {
        double mean_i1 = vmean(series_i1);
        double sd_i1 = vsd(series_i1);
        for (int j = 0; j < ws/2; j++) {
          nlls[i] += dNorm(series_i[j], mean_i1, sd_i1);
        }
        double mean_i2 = vmean(series_i2);
        double sd_i2 = vsd(series_i2);
        for (int j = ws/2; j < ws; j++) {
          nlls[i] += dNorm(series_i[j], mean_i2, sd_i2);
        }
      }
      
      // Make negative
      //nlls[i] *= -1;
      
      // Check for nan's and replace with 1
      if (std::isnan(nlls[i])) {nlls[i] = 1.0;}
     
    }
   
    // Fill in the last part of the vector with the mean
    double nlls_mean = vmean(nlls);
    for (int i = n - ws + 1; i < n; i++) {
      nlls[i] = nlls_mean;
    }
    
    return nlls;
    
  }

// Likelihood ratio outlier change-point detection
IntegerVector LROcp_find(
    const dVec& nll_ratio,        // 1D vector of points to test for change points
    const int& ws,                // Running window size
    const double& out_mult        // Outlier multiplier
  ) {
    
    //Find change points
    iVec fcp;
    int last_cp = 0;
    double nll_ratio_mean = vmean(nll_ratio);
    double nll_ratio_sd = vsd(nll_ratio);
    double nll_max = nll_ratio_mean + out_mult * nll_ratio_sd;
    int n_nll = nll_ratio.size();
    for (int i = 0; i < n_nll; i++) {
      if (
          nll_ratio[i] > nll_max && 
            i < (nll_ratio.size() - ws/2) && // can't be too close to end
            i > ws/2                         // can't be too close to beginning
      ) {
        if (i - last_cp > ws/2) {            // can't be too close to last change point
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

// Likelihood ratio outlier change-point detection, single-series
IntegerVector LROcp(
    const dVec& series,           // 1D vector of points to test for change points
    const int& ws,                // Running window size
    const int& filter_ws,         // Size of window for taking rolling mean
    const double& out_mult        // Outlier multiplier
  ) {
    
    // Compute the null negative log-likelihood at each bin (no change point)
    dVec nll_null = series_nll(series, ws, filter_ws, true);
    
    // Compute the change-point negative log-likelihood at each bin
    dVec nll_cp = series_nll(series, ws, filter_ws, false);
    
    // Find the nll ratio at each bin
    dVec nll_ratio = vsubtract(nll_cp, nll_null);
    
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
    NumericMatrix nll_ratio_array(n_samples, n_series);
    
    // Convert series_array of raw (or log) count values into an array of likelihood ratios
    for (int s = 0; s < n_series; s++) {
      
      // Get the series for this column
      dVec series = to_dVec(series_array.col(s));
      
      // Compute the null negative log-likelihood at each bin (no change point)
      dVec nll_null = series_nll(series, ws, filter_ws, true);
      
      // Compute the change-point negative log-likelihood at each bin
      dVec nll_cp = series_nll(series, ws, filter_ws, false);
      
      // Find the nll ratio at each bin 
      // ... expected value is 1. Expect nll_null to be 
      // ... larger than nll_cp (i.e., a worse fit), so expect
      // ... this value to be 1 or larger. 
      dVec nll_ratio = vdivide(nll_cp, nll_null);
      
      // Save in array 
      nll_ratio_array.column(s) = to_NumVec(nll_ratio);
      
    }
    
    // Find the centroid of nll_ratio_array
    // ... calling function from dtwclust package in R
    Function find_centroid("find_centroid");   
    NumericVector centroid = find_centroid(nll_ratio_array);
    
    // Use LROcp on this series to estimate change points 
    IntegerVector found_cp = LROcp_find(to_dVec(centroid), ws, out_mult);
    
    if (found_cp.size() > 0) {
      
      // Project found_cp back to original series 
      // ... these will be one-based indices that need to be subtracted back to zero-based
      Function project_cp("project_cp");
      IntegerMatrix found_cp_array = project_cp(found_cp, centroid, nll_ratio_array);
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

