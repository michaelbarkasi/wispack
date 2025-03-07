
// model_fit.cpp
#include "wspc.h"

// Log of density of normal distribution
sdouble log_dnorm(
    const sdouble& x,        // value to evaluate
    const sdouble& mu,       // mean (expected value)
    const sdouble& sd        // standard deviation
  ) {
    return slog(
      sexp(-spower(x - mu, 2.0) / (2.0 * spower(sd, 2.0))) / (ssqrt(2.0 * M_PI) * sd)
    );
  }

// Log of density of gamma distribution
sdouble log_dgamma(
    const sdouble& x,        // value to evaluate
    const sdouble& shape,    // shape parameter
    const sdouble& rate      // rate parameter
  ) {
    return slog(
      (spower(rate, shape) * spower(x, shape - 1.0) * sexp(-rate * x)) / stan::math::tgamma(shape)
    );
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
    dVec nlls(n);
    
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
   
    for (int i = n - ws + 1; i < n; i++) {
      nlls[i] = 1.0;
    }
    
    return nlls;
    
  }

// Likelihood ratio outlier change-point detection
IntegerVector LROcp(
    const dVec& series,          // 1D vector of points to test for change points
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
    double nll_ratio_sd = Rcpp::sd(to_NumVec(nll_ratio));
    double nll_max = nll_ratio_mean + out_mult*nll_ratio_sd;
    int n_nll = nll_ratio.size();
    LogicalVector outlier_mask(n_nll);
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
    
    IntegerVector found_cp;
    if (fcp.size() > 0) {
      // Grab the indexes of the found change points and add 1 (first bin is 1)
      found_cp = to_IntVec(fcp) + 1;
    } 
    
    return found_cp;
    
  }