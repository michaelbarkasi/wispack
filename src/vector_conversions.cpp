
// vector_conversions.cpp
#include "wspc.h"

// Convert to Eigen::Matrix with sdouble elements **********************************************************************

// ... from NumericVector
sVec to_sVec(
    const NumericVector& vec
  ) {
    sVec stan_vec(vec.size());
    for (int i = 0; i < vec.size(); i++) {
      stan_vec(i) = sdouble(static_cast<double>(vec[i]));
    }
    return stan_vec;
  }

// ... overload, from std::vector with doubles 
sVec to_sVec(
    const dVec& vec
  ) {
    sVec stan_vec(vec.size());
    for (int i = 0; i < vec.size(); i++) {
      stan_vec(i) = sdouble(vec[i]);
    }
    return stan_vec;
  }

// ... overload, from std::vector with int
sVec to_sVec(
    const iVec& vec
  ) {
    sVec stan_vec(vec.size());
    for (int i = 0; i < vec.size(); i++) {
        stan_vec(i) = sdouble((double)vec[i]);
    }
    return stan_vec;
  }

// ... overload, from IntegerVector
sVec to_sVec(
    const IntegerVector& vec
  ) {
    NumericVector vec_ = as<NumericVector>(vec);
    sVec stan_vec(vec_.size());
    for (int i = 0; i < vec_.size(); i++) {
      stan_vec(i) = sdouble(static_cast<double>(vec_[i]));
    }
    return stan_vec;
  }

// Convert to NumericVector ********************************************************************************************

// ... from Eigen::Matrix with sdouble elements
NumericVector to_NumVec(
    const sVec& vec
  ) {
    NumericVector num_vec(vec.size());
    for (int i = 0; i < vec.size(); i++) {
      if (std::isnan(vec(i))) {
        num_vec(i) = NA_REAL;
      } else {
        num_vec(i) = vec(i).val();
      }
    }
    return num_vec;
  }

// ... overload, from std::vector with doubles
NumericVector to_NumVec(
    const dVec& vec
  ) {
    return Rcpp::wrap(vec);
  }

// ... overload, from IntegerVector
NumericVector to_NumVec(
    const IntegerVector& vec
  ) {
    NumericVector num_vec(vec.size());
    for (int i = 0; i < vec.size(); i++) {
      num_vec(i) = vec[i];
    }
    return num_vec;
  }

// Convert to IntegerVector ********************************************************************************************

// ... from std::vector with int
IntegerVector to_IntVec(
    const iVec& vec
  ) {
    return Rcpp::wrap(vec);
  }

// Convert to std::vector with doubles *********************************************************************************

// ... from NumericVector
dVec to_dVec(
    const NumericVector& vec
  ) {
    return Rcpp::as<dVec>(vec);
  }

// Convert to std::vector with int *************************************************************************************

// ... from IntegerVector
iVec to_iVec(
    const IntegerVector& vec
  ) {
    return Rcpp::as<iVec>(vec);
  }

