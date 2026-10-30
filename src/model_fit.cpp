
// model_fit.cpp
#include "wspc.h"
#include <cmath>
#include <limits>
#include <algorithm>
#include <utility>

// Random draw from uniform distribution, integers
int rUnifi(
    int min,
    int max,
    std::mt19937& rng
  ) {
    std::uniform_int_distribution<int> dist(min, max);
    return dist(rng);
  }

// Random draw from uniform distribution, reals
double rUnifr(
    double min,
    double max,
    std::mt19937& rng
  ) {
    std::uniform_real_distribution<double> dist(min, max);
    return dist(rng);
  }

// Random draw from normal distribution
double rNorm(
    double mean, 
    double sd,
    std::mt19937& rng
  ) {
    std::normal_distribution<double> norm(mean, sd);
    return norm(rng);
  }

// Random draw from normal distribution, thread-safe
double safe_rNorm(
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

// Integral of Poisson-Gamma distribution, from 0 to positive infinity
sdouble poisson_gamma_integral(
    sdouble y, // observed count value
    sdouble r, // expected process rate drawn from the gamma distribution
    sdouble v  // variance of the gamma distribution
  ) {
   
    // convert more intuitive rate and variance parameters into the standard shape-rate parameters 
    // ... for the gamma distribution
    sdouble s = r * r / v; // Gamma distribution "shape" parameter
    sdouble R = s / r;     // Gamma distribution "rate" parameter
    
    // Idea: This is an analytic solution to the integral of dPois(y, lambda) * dGamma(lambda, r, v) from 0 to positive infinity. 
    //        ... The solution takes the form of a ratio num/denom which is subject to overflow/underflow, and so instead of 
    //            computing this ratio directly, we compute the log of the numerator and denominator separately, and then exponentiate the difference.
    sdouble log_num = s * slog(R) + stan::math::lgamma(y + s); 
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

// Formula to calculate gamme dispersion factor from mean and variance of counts
double gamma_dispersion_formula(
    const double& count_cs_mean, // mean of counts for context-species combination
    const double& count_cs_var   // variance of counts for context-species combination
  ) {
    
    if (count_cs_var > count_cs_mean) {
      // Have: count_cs_var = count_cs_mean + count_cs_mean^2 * gdis
      // ... count_cs_var = count_cs_mean * (1 + count_cs_mean * gdis)
      // ... count_cs_var / count_cs_mean = 1 + count_cs_mean * gdis
      // ... (count_cs_var / count_cs_mean) - 1 = count_cs_mean * gdis
      return ((count_cs_var / count_cs_mean) - 1.0) / count_cs_mean;
    } else {
      return 0.0; // no dispersion if variance is less than or equal to mean
    }
  }

// Generic warping function with bound 
sdouble warp_ratio(
    const sdouble& basis,    // parameterizing coordinate to set the warp
    const sdouble& b,        // bound on this value 
    const sdouble& w         // warping parameter
  ) {
    return 1.0 - spower(sexp(w*w), -basis/b);
  }

// Warping function for model components 
sdouble warp_mc(
    const sdouble& s,        // value to warp
    const double& b,         // bound on this value 
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
    return s + direction * warp_ratio(basis, b, w) * magnitude;
  }

// Wrapper for R
// [[Rcpp::export]]
double warp_mc_R(
    const double& s,        // value to warp
    const double& b,        // bound on this value 
    const double& w         // warping parameter
  ) {
    return warp_mc(sdouble(s), b, sdouble(w)).val();
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
    return sigmoid_stable(sdouble(x)).val();
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
      if (Rt.size() < deg + 1) {Rcpp::stop("Error in poly_sigmoid: Rt.size() < deg + 1");}
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
    const double& b,                   // input variable
    const int& deg,                    // degree of the poly-sigmoid, i.e., number of transitions between blocks
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
    const double& a
  ) {
    sdouble z = a * x;
    if (z <= 1e-30) {z = 1e-30;} 
    return 1.0 / (z * z);
  }

// Calculate rolling-window log-likelihood of a series being generated by a given rate, with or without change-point
dVec series_loglik(
    const dVec& series0,          // 1D vector of points for which to take likelihood of a change-point
    const int& ws,                // Running window size
    const bool& null              // If true, compute likelihood of data assuming no transitions; otherwise, assuming transition
  ) {
    // Assumes points are normally distributed about their process rate
    
    int n = series0.size();
    
    // Guard: ws must be >= 2 (the algorithm splits the window into halves)
    // and must not exceed n; return zeros if the window is degenerate.
    if (ws < 2 || ws > n) {
      return dVec(n, 0.0);
    }
    
    dVec series = series0; 
    dVec loglik(n); // length out should equal length in
    
    // Compute joint log-likelihood of observations within each window
    for (int i = 0; i <= (n - ws); i++) {
      
      int i_mid = i + int(ws/2);
     
      loglik[i_mid] = 0.0;
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
          loglik[i_mid] += log_dNorm(series_i[j], mean_i, sd_i);
        }
      } else {
        double mean_i1 = vmean(series_i1);
        double sd_i1 = vsd(series_i1);
        for (int j = 0; j < ws/2; j++) {
          loglik[i_mid] += log_dNorm(series_i[j], mean_i1, sd_i1);
        }
        double mean_i2 = vmean(series_i2);
        double sd_i2 = vsd(series_i2);
        for (int j = ws/2; j < ws; j++) {
          loglik[i_mid] += log_dNorm(series_i[j], mean_i2, sd_i2);
        }
      }
      
      // Check for nan's and replace with 1
      if (std::isnan(loglik[i_mid])) {loglik[i_mid] = 1.0;}
     
    }
   
    // Fill in the last part of the vector with the mean
    double loglik_mean = vmean(loglik);
    for (int i = 0; i < int(ws/2); i++) {
      loglik[i] = loglik_mean;
    }
    for (int i = n - ws + int(ws/2) + 1; i < n; i++) {
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

// ... overload
IntegerVector LROcp_find(
    const NumericMatrix& loglik_ratio_mat,     // NumericMatrix of vectors (columns) to test for change points
    const int& ws,                             // Running window size
    const double& out_mult                     // Outlier multiplier
  ) {
    
    //Find change points
    iVec fcp;
    int last_cp = 0;
    double last_nll = 0.0;
    int n_series = loglik_ratio_mat.cols();
    dVec loglik_max(n_series);
    for (int s = 0; s < n_series; s++) {
      dVec loglik_ratio = to_dVec(loglik_ratio_mat.column(s));
      loglik_max[s] = vmean(loglik_ratio) + out_mult * vsd(loglik_ratio);
    }
    
    int n_nll = loglik_ratio_mat.rows();
    for (int i = 0; i < n_nll; i++) {
      bool search = true;
      for (int s = 0; s < n_series && search; s++) {
        if (
            loglik_ratio_mat(i,s) > loglik_max[s] && 
              i < (n_nll - ws/2) &&                 // can't be too close to end
              i > ws/2                              // can't be too close to beginning
        ) {
          if (i - last_cp > int(ws/2)) {            // can't be too close to last change point
            fcp.push_back(i);
            last_cp = i;
            last_nll = loglik_ratio_mat(i,s);
            search = false; 
          } else if (loglik_ratio_mat(i,s) > last_nll) {
            // If too close to last change point, check if this ratio is larger
            fcp[fcp.size() - 1] = i;
            last_cp = i;
            last_nll = loglik_ratio_mat(i,s);
            search = false;
          }
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

// Dynamic-time-warp alignment under the Sakoe-Chiba symmetric step pattern of slope 0.5 (symmetricP05;
// Sakoe & Chiba 1978), used by project_cp() below to project a change point found on a cross-series
// centroid back onto one individual series. Ported from R (this package's former R-level project_cp(),
// which called dtw::dtw(..., step.pattern = dtw::symmetricP05)) so the compiled model-fitting code no
// longer needs to call back into R for this. The recursion and backtrace mirror dtw:::globalCostMatrix()'s
// interpreted fallback and dtw:::backtrack() line for line (confirmed by reading the dtw package's R
// source directly), validated against dtw's own output on this package's test fixtures before this
// replaced the R callback (see tests/testthat/test-dtw-port.R):
//
//   g[i,j] = min(
//     g[i-1,j-3] + 2*d[i,  j-2] +   d[i,  j-1] +   d[i,j],
//     g[i-1,j-2] + 2*d[i,  j-1] +   d[i,j],
//     g[i-1,j-1] + 2*d[i,j],
//     g[i-2,j-1] + 2*d[i-1,j  ] +   d[i,j],
//     g[i-3,j-1] + 2*d[i-2,j  ] +   d[i-1,j] + d[i,j]
//   )
//
// with local cost d[i,j] = |query[i] - ref[j]| (dtw()'s default dist.method = "Euclidean" reduces to
// absolute difference for 1-D series), global alignment (no open begin/end -- dtw()'s defaults, left
// unchanged by the original call), and no window restriction (dtw()'s default window.type = "none").
void dtw_symmetricP05_align(
    const dVec& query,     // series to align, length n
    const dVec& ref,       // reference/centroid series, length m
    iVec&       index1,    // out: aligned positions in query (1-indexed, dtw()'s index1 convention)
    iVec&       index2     // out: aligned positions in ref (1-indexed, dtw()'s index2 convention)
  ) {

    int n = query.size();
    int m = ref.size();

    auto d = [&](int i, int j) { return std::fabs(query[i - 1] - ref[j - 1]); }; // 1-indexed access

    // Global cost matrix (1-indexed; row/col 0 unused) and chosen-pattern matrix (0 = unreachable).
    // Only cm[1][1] is seeded -- every other cell is computed only if reachable via one of the 5
    // patterns below, exactly mirroring dtw:::globalCostMatrix()'s cm[1,1] <- lm[1,1] seed with every
    // other cell left NA until (and unless) a valid predecessor is found.
    const double NA_ = std::numeric_limits<double>::quiet_NaN();
    std::vector<std::vector<double>> cm(n + 1, std::vector<double>(m + 1, NA_));
    std::vector<std::vector<int>>    sm(n + 1, std::vector<int>(m + 1, 0));
    cm[1][1] = d(1, 1);

    for (int j = 1; j <= m; j++) {
      for (int i = 1; i <= n; i++) {
        if (i == 1 && j == 1) continue;
        double best = NA_;
        int bestp = 0;
        if (i - 1 >= 1 && j - 3 >= 1 && !std::isnan(cm[i-1][j-3])) {
          double c = cm[i-1][j-3] + 2*d(i,j-2) + d(i,j-1) + d(i,j);
          if (std::isnan(best) || c < best) { best = c; bestp = 1; }
        }
        if (i - 1 >= 1 && j - 2 >= 1 && !std::isnan(cm[i-1][j-2])) {
          double c = cm[i-1][j-2] + 2*d(i,j-1) + d(i,j);
          if (std::isnan(best) || c < best) { best = c; bestp = 2; }
        }
        if (i - 1 >= 1 && j - 1 >= 1 && !std::isnan(cm[i-1][j-1])) {
          double c = cm[i-1][j-1] + 2*d(i,j);
          if (std::isnan(best) || c < best) { best = c; bestp = 3; }
        }
        if (i - 2 >= 1 && j - 1 >= 1 && !std::isnan(cm[i-2][j-1])) {
          double c = cm[i-2][j-1] + 2*d(i-1,j) + d(i,j);
          if (std::isnan(best) || c < best) { best = c; bestp = 4; }
        }
        if (i - 3 >= 1 && j - 1 >= 1 && !std::isnan(cm[i-3][j-1])) {
          double c = cm[i-3][j-1] + 2*d(i-2,j) + d(i-1,j) + d(i,j);
          if (std::isnan(best) || c < best) { best = c; bestp = 5; }
        }
        if (bestp > 0) { cm[i][j] = best; sm[i][j] = bestp; }
      }
    }

    if (std::isnan(cm[n][m])) {
      Rcpp::stop("dtw_symmetricP05_align: no warping path exists between the given series");
    }

    // Backtrack from (n, m) to (1, 1): replay each winning pattern's intermediate points
    // (nearest-to-farthest -- dtw:::.extractpattern()'s row order after its reversal), then jump to
    // its anchor cell and repeat, exactly mirroring dtw:::backtrack().
    index1.clear(); index2.clear();
    index1.push_back(n); index2.push_back(m);
    int i = n, j = m;
    while (!(i == 1 && j == 1)) {
      int p = sm[i][j];
      if (p == 0) break; // mirrors dtw:::backtrack()'s is.na(s) guard; shouldn't occur for a valid path
      std::vector<std::pair<int,int>> steps; // (di, dj) offsets from (i,j); last entry is the anchor
      switch (p) {
        case 1: steps = {{0,1},{0,2},{1,3}}; break;
        case 2: steps = {{0,1},{1,2}};       break;
        case 3: steps = {{1,1}};             break;
        case 4: steps = {{1,0},{2,1}};       break;
        case 5: steps = {{1,0},{2,0},{3,1}}; break;
      }
      for (auto& st : steps) {
        index1.push_back(i - st.first);
        index2.push_back(j - st.second);
      }
      i -= steps.back().first;
      j -= steps.back().second;
    }
    std::reverse(index1.begin(), index1.end());
    std::reverse(index2.begin(), index2.end());

  }

// Round to nearest integer, ties to even (matches R's round(), IEC 60559). Exact for the small
// integer-valued averages computed in project_cp() below (series index positions), where
// floating-point representation of a x.5 halfway value is exact at these magnitudes.
int round_half_even(
    double x
  ) {
    double fl = std::floor(x);
    double frac = x - fl;
    int fi = (int) fl;
    if (frac < 0.5) return fi;
    if (frac > 0.5) return fi + 1;
    return (fi % 2 == 0) ? fi : fi + 1;
  }

// Project change points found on the across-series centroid back onto each individual series, via DTW
// alignment of each series to the centroid (dtw_symmetricP05_align() above). Ports what was formerly
// an R-level function (project_cp(), which did the same thing by calling dtw::dtw(...,
// step.pattern = dtw::symmetricP05)) to C++, so the compiled model-fitting code no longer needs to call
// back into R for this -- and project_cp need not be an exported R function any more.
IntegerMatrix project_cp(
    const IntegerVector& found_cp,     // change-point indices on centroid (1-indexed)
    const NumericVector& centroid,     // across-series centroid, length n_samples
    const NumericMatrix& X_lik         // per-series data the centroid was computed from (n_samples x n_series)
  ) {

    int n_cp      = found_cp.size();
    int n_series  = X_lik.ncol();
    int n_samples = X_lik.nrow();

    if (n_cp == 0 || centroid.size() == 0 || n_series == 0) {
      IntegerMatrix err(1, 1);
      err(0, 0) = -1;
      return err;
    }

    dVec centroid_v = to_dVec(centroid);
    IntegerMatrix out(n_cp, n_series);

    for (int s = 0; s < n_series; s++) {

      // Series for this column, with any NA replaced by 1 (matches R's Xi[is.na(Xi)] <- 1)
      dVec series(n_samples);
      for (int k = 0; k < n_samples; k++) {
        double v = X_lik(k, s);
        series[k] = std::isnan(v) ? 1.0 : v;
      }

      iVec index1, index2;
      dtw_symmetricP05_align(series, centroid_v, index1, index2);

      for (int c = 0; c < n_cp; c++) {
        int cp = found_cp[c];
        double sum = 0.0;
        int count = 0;
        for (size_t k = 0; k < index2.size(); k++) {
          if (index2[k] == cp) { sum += index1[k]; count++; }
        }
        out(c, s) = round_half_even(sum / count);
      }

    }

    return out;

  }

// Internal-only test hook for project_cp()'s DTW port (tests/testthat/test-dtw-port.R), accessed via
// wispack:::.project_cp_test_hook() -- not in NAMESPACE's export() list, so not part of the public API.
// [[Rcpp::export(name = ".project_cp_test_hook")]]
IntegerMatrix project_cp_test_hook(
    const IntegerVector& found_cp,
    const NumericVector& centroid,
    const NumericMatrix& X_lik
  ) {
    return project_cp(found_cp, centroid, X_lik);
  }

// Likelihood ratio outlier change-point detection, arrays
IntegerMatrix LROcp_array(
    const eMat& series_array,     // 2D matrix of points to test for change points
    const int& ws,                // Running window size
    const double& out_mult,       // Outlier multiplier
    const double& cp_buffer       // Minimum distance between two change points
  ) {
    // Test treats each column of the matrix as a separate series. 
    //  It further assumes that the series are dependent on each other, and share 
    //  change points, plus some variance between them on the exact location. 
    
    // Control random number generator
    Rcpp::RNGScope scope;
    
    int n_series = series_array.cols();
    int n_samples = series_array.rows();
    NumericMatrix loglik_ratio_array_highpass(n_samples, n_series);
    NumericMatrix loglik_ratio_array_midpass(n_samples, n_series);
    NumericMatrix loglik_ratio_array_lowpass(n_samples, n_series);
    
    // Convert series_array of (raw or log) count values into an array of likelihood ratios
    for (int s = 0; s < n_series; s++) {
      
      // Get the series for this column
      dVec series = to_dVec(series_array.col(s));
      
      // Find the likelihood ratios of change points for this series
      dVec loglik_ratio_highpass = LROcp_logRatio(series, ws);
      dVec loglik_ratio_midpass  = LROcp_logRatio(series, ws*2);
      dVec loglik_ratio_lowpass  = LROcp_logRatio(series, ws*4);
      
      // Save in array 
      loglik_ratio_array_highpass.column(s) = to_NumVec(loglik_ratio_highpass);
      loglik_ratio_array_midpass.column(s) = to_NumVec(loglik_ratio_midpass);
      loglik_ratio_array_lowpass.column(s) = to_NumVec(loglik_ratio_lowpass);
      
    }
    
    // Find centroid
    NumericMatrix centroid(n_samples, 3);
    for (int i = 0; i < n_samples; i++) {
      centroid(i,0) = vmean(loglik_ratio_array_highpass.row(i));
      centroid(i,1) = vmean(loglik_ratio_array_midpass.row(i));
      centroid(i,2) = vmean(loglik_ratio_array_lowpass.row(i));
    }
    
    // Use LROcp on this series to estimate change points 
    IntegerVector found_cp = LROcp_find(centroid, ws, out_mult);
    
    if (found_cp.size() > 0) {
      
      // Project found_cp back to original series
      IntegerMatrix found_cp_array = project_cp(found_cp, NumericVector(centroid.column(1)), loglik_ratio_array_midpass);
      // ... ^ in this matrix, rows are change points (by deg), columns are trt x ran interactions
      if (found_cp_array(0,0) == -1) {Rcpp::stop("Problem running project_cp");}
     
      // Ensure change points are separated by at least cp_buffer (and correct indexes to zero-indexing)
      IntegerVector collision_rows;
      for (int i = 0; i < found_cp_array.ncol(); i++) {
        // Convert one-based indices from project_cp to zero-based indices
        found_cp_array.column(i) = found_cp_array.column(i) - 1;
        for (int k = 0; k < found_cp_array.nrow(); k++) {
          if (k == 0) {
            if (found_cp_array(k, i) <= cp_buffer) {
              found_cp_array(k, i) = cp_buffer + 1.0;
            }
          } else {
            int cp_gap = found_cp_array(k, i) - found_cp_array(k - 1, i);
            if (cp_gap <= 0 || found_cp_array(k, i) >= n_samples - cp_buffer - 1.0) {
              // Change point collision! 
              if (k == found_cp_array.nrow() - 1 && (n_samples - cp_buffer - 1.0 - found_cp_array(k - 1, i) > cp_buffer)) {
                found_cp_array(k, i) = n_samples - cp_buffer - 1;
              } else {
                // Mark row for removal
                collision_rows.push_back(k);
                // Advance counter k
                k++;
              }
            } else if (cp_gap <= cp_buffer) {
              // Change points too close together, so move this one forward
              found_cp_array(k, i) = found_cp_array(k - 1, i) + cp_buffer + 1.0;
            }
          }
        }
      }
      
      // Remove collision rows
      if (collision_rows.size() > 0) {
        collision_rows = Rcpp::sort_unique(collision_rows); 
        int n_collisions = collision_rows.size();
        for (int i = 0; i < n_collisions; i++) {
          int r = collision_rows[i];
          found_cp_array = remove_row(found_cp_array, r);
          collision_rows = collision_rows - 1; // after removing a row, all subsequent rows shift up by one
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
    return out;
    
  }
