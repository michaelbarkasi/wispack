
// model_fit.cpp
#include "wspc.h"

// Log of density of normal distribution
sdouble log_dNorm(
    const sdouble& x,        // value to evaluate
    const sdouble& mu,       // mean (expected value)
    const sdouble& sd        // standard deviation
  ) {
    return slog(
      sexp(-spower(x - mu, 2.0) / (2.0 * spower(sd, 2.0))) / (ssqrt(2.0 * M_PI) * sd)
    );
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
    sdouble y, 
    sdouble r, 
    sdouble v
  ) {
    sdouble s = r * r / v;
    sdouble R = s / r;
    
    //sdouble num = spower(R, s) * stan::math::tgamma(y + s) * (1.0 - stan::math::gamma_q(y + s, R + 1.0));
    sdouble log_num = s * slog(R) + stan::math::lgamma(y + s) + stan::math::log1m(stan::math::gamma_p(y + s, R + 1.0));
    sdouble num = exp(log_num);
    sdouble denom = (sexp(stan::math::lgamma(y + 1.0)) * stan::math::tgamma(s) * spower(R + 1.0, y + s));
    return num / denom;
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
    sdouble x_ = -x;
    // x is negative boundary distance, so want to flip back to boundary distance 
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
    
    for (int i = 0; i <= (n - ws); i++) {
    
      nlls[i] = 0.0;
      dVec series_i(series.begin() + i, series.begin() + (i + ws));
      dVec series_i1(series_i.begin(), series_i.begin() + (ws / 2));
      dVec series_i2(series_i.begin() + (ws / 2), series_i.end());
      
      if (null) {
        double mean_i = vmean(series_i);
        double sd_i = Rcpp::sd(to_NumVec(series_i));
        for (int j = 0; j < ws; j++) {
          nlls[i] += R::dnorm(series_i[j], mean_i, sd_i, true);
        }
      } else {
        double mean_i1 = vmean(series_i1);
        double sd_i1 = Rcpp::sd(to_NumVec(series_i1));
        for (int j = 0; j < ws/2; j++) {
          nlls[i] += R::dnorm(series_i1[j], mean_i1, sd_i1, true);
        }
        double mean_i2 = vmean(series_i2);
        double sd_i2 = Rcpp::sd(to_NumVec(series_i2));
        for (int j = ws/2; j < ws; j++) {
          nlls[i] += R::dnorm(series_i2[j], mean_i2, sd_i2, true);
        }
      }
      
      nlls[i] *= -1;
    
    }
   
    // Fill in the last part of the vector with 1.0
    for (int i = n - ws + 1; i < n; i++) {
      nlls[i] = 1.0;
    }
    
    return nlls;
    
  }

// Likelihood ratio outlier change-point detection
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
    
    //Find change points
    iVec fcp;
    int last_cp = 0;
    double nll_ratio_mean = vmean(nll_ratio);
    double nll_ratio_sd = vsd(nll_ratio);
    double nll_max = nll_ratio_mean + out_mult*nll_ratio_sd;
    int n_nll = nll_ratio.size();
    for (int i = 0; i < n_nll; i++) {
      if (
          nll_ratio[i] > nll_max && 
          i < (series.size() - ws/2) && // can't be too close to end
          i > ws/2                      // can't be too close to beginning
      ) {
        if (i - last_cp > ws/2) {       // can't be too close to last change point
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

// ... overload for arrays
IntegerMatrix LROcp_array(
    const sMat& series_array,     // 2D matrix of points to test for change points
    const int& ws,                // Running window size
    const int& filter_ws,         // Size of window for taking rolling mean
    const double& out_mult        // Outlier multiplier
  ) {
    // The test here will treat each column of the matrix as a separate series. 
    //  it further assumes that the series are dependent on each other, and share 
    //  change points, plus some variance between them on the exact location. 
    
    int n_series = series_array.cols();
    int series_length = series_array.rows();
    int window_steps = series_length - ws; 
    sMat nll_ratio_array(series_length, n_series);
    sMat alignment_values(series_length, n_series);
    sMat divisor_values(series_length, n_series);
    alignment_values.setZero();
    divisor_values.setOnes();
    
    for (int s = 0; s < n_series; s++) {
      
      // Get the series for this column
      dVec series = to_dVec(series_array.col(s));
      
      // Compute the null negative log-likelihood at each bin (no change point)
      dVec nll_null = series_nll(series, ws, filter_ws, true);
      
      // Compute the change-point negative log-likelihood at each bin
      dVec nll_cp = series_nll(series, ws, filter_ws, false);
      
      // Find the nll ratio at each bin
      dVec nll_ratio = vsubtract(nll_cp, nll_null);
      
      // Save in array 
      nll_ratio_array.col(s) = to_sVec(nll_ratio);
      
    }
    
    // Align the individual series so that peak nll ratios are at the same time within the window
    for (int w = 0; w < window_steps; w++) {
      
      int idx0 = w;
      int idx1 = w + ws - 1;
      sMat nll_ratio_array_block = nll_ratio_array.block(idx0, 0, ws, n_series);
      iVec idx_max(n_series); 
      
      // Find max nll ratio in this window block, for each series
      for (int s = 0; s < n_series; s++) {nll_ratio_array_block.col(s).maxCoeff(&idx_max[s]);}
      
      // Find mean position of the max nll ratio 
      int idx_max_mean = vmean(idx_max);
      
      // Compute and save alignment values 
      for (int s = 0; s < n_series; s++) {alignment_values(w, s) += (sdouble)(idx_max_mean - idx_max[s]);}
      
      // Perform alignment
      for (int s = 0; s < n_series; s++) {
        for (int i = idx0; i < idx1 + 1; i++) {
          int new_idx = i + (int)alignment_values(w, s).val();
          if (new_idx < 0) {new_idx = 0;}
          if (new_idx >= series_length) {new_idx = series_length - 1;}
          divisor_values(new_idx, s) += 1.0;
          nll_ratio_array(new_idx, s) += nll_ratio_array(i, s);
          nll_ratio_array(new_idx, s) /= divisor_values(new_idx, s);
        }
      }
      
    }
    
    // Collapse along rows
    dVec nll_ratio = dVec(series_length);
    for (int i = 0; i < series_length; i++) {
      nll_ratio[i] = nll_ratio_array.row(i).mean().val();
    }
    
    //Find change points
    iVec fcp;
    int last_cp = 0;
    double nll_ratio_mean = vmean(nll_ratio);
    double nll_ratio_sd = vsd(nll_ratio);
    double nll_max = nll_ratio_mean + out_mult*nll_ratio_sd;
    for (int i = 0; i < series_length; i++) {
      if (
          nll_ratio[i] > nll_max && 
            i < (series_length - ws/2) && // can't be too close to end
            i > ws/2                      // can't be too close to beginning
      ) {
        if (i - last_cp > ws/2) {         // can't be too close to last change point
          fcp.push_back(i);
          last_cp = i;
        }
      }
    }
    
    // Make (1D) found_cp vector
    IntegerVector found_cp;
    if (fcp.size() > 0) {
      // Grab the indexes of the found change points and add 1 (first bin is 1)
      found_cp = to_IntVec(fcp) + 1;
    } 
    
    // Expand back out to array
    int n_cp = found_cp.size();
    IntegerMatrix found_cp_array(n_cp, n_series);
    for (int s = 0; s < n_series; s++) {
      for (int i = 0; i < n_cp; i++) {
        int align_value = (int)alignment_values(found_cp(i), s).val();
        found_cp_array(i, s) = found_cp(i) - align_value;
        if (found_cp_array(i, s) < 1) {found_cp_array(i, s) = 1;}
        if (found_cp_array(i, s) > series_length) {found_cp_array(i, s) = series_length;}
      }
    }
    
    return found_cp_array;
    
  } 