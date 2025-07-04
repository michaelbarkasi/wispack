
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

// ... to NumericMatrix
NumericMatrix to_NumMat(
    const sMat& mat
  ) {
    NumericMatrix num_mat(mat.rows(), mat.cols());
    for (int i = 0; i < mat.rows(); i++) {
      for (int j = 0; j < mat.cols(); j++) {
        if (std::isnan(mat(i, j))) {
          num_mat(i, j) = NA_REAL;
        } else {
          num_mat(i, j) = mat(i, j).val();
        }
      }
    }
    return num_mat;
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

// ... from Eigen::Matrix with sdouble elements
dVec to_dVec(
    const sVec& vec
  ) {
    dVec dvec(vec.size());
    for (int i = 0; i < vec.size(); i++) {
      dvec[i] = vec(i).val();
    }
    return dvec;
  }

// Convert to std::vector with int *************************************************************************************

// ... from IntegerVector
iVec to_iVec(
    const IntegerVector& vec
  ) {
    return Rcpp::as<iVec>(vec);
  }

// Misc ****************************************************************************************************************

// Convert to Eigen::Matrix with sdouble elements
sMat to_sMat(
    const IntegerMatrix& mat
  ) {
    sMat stan_mat(mat.nrow(), mat.ncol());
    for (int i = 0; i < mat.nrow(); i++) {
      for (int j = 0; j < mat.ncol(); j++) {
        stan_mat(i, j) = sdouble(static_cast<double>(mat(i, j)));
      }
    }
    return stan_mat;
  }

// ... overload
sMat to_sMat(
    const NumericMatrix& mat
  ) {
    sMat stan_mat(mat.nrow(), mat.ncol());
    for (int i = 0; i < mat.nrow(); i++) {
      for (int j = 0; j < mat.ncol(); j++) {
        stan_mat(i, j) = sdouble(mat(i, j));
      }
    }
    return stan_mat;
  }

