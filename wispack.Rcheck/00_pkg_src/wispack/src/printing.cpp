
// printing.cpp
#include "wspc.h"

// Function for printing character messages
void vprint(
    const std::string& message,
    bool verbose
  ) {
    if (verbose) {
      Rcpp::Rcout << message << std::endl;
    }
  }

// ... overload for no verbose flag
void vprint(
    const std::string& message
  ) {
    Rcpp::Rcout << message << std::endl;
  }

// ... overload for sdouble variant
void vprint(
    const std::string& message, 
    sdouble val
  ) {
    Rcpp::Rcout << message << val << std::endl;
  }

// ... overload for double variant
void vprint(
    const std::string& message, 
    double val
  ) {
    Rcpp::Rcout << message << val << std::endl;
  }

// ... overload for int variant
void vprint(
    const std::string& message, 
    int val
  ) {
    Rcpp::Rcout << message << val << std::endl;
  }

// Function for printing vectors
void vprintV(
    const CharacterVector& message,
    bool verbose
  ) {
    if (verbose) {
      Rcpp::Rcout << message << std::endl;
    }
  }

// ... overload for NumericVector
void vprintV(
    const NumericVector& vec,
    bool verbose
  ) {
    if (vec.size() > 20) {
      Rcpp::Rcout << "vec (first 20): ";
      for (int i = 0; i < 20; ++i) {
        Rcpp::Rcout << vec(i) << " ";
      }
      Rcpp::Rcout << std::endl;
    } else {
      Rcpp::Rcout << "vec: ";
      for (int i = 0; i < vec.size(); ++i) {
        Rcpp::Rcout << vec(i) << " ";
      }
      Rcpp::Rcout << std::endl;
    }
  }

// ... overload for IntegerVector
void vprintV(
    const IntegerVector& vec,
    bool verbose
  ) {
    if (vec.size() > 20) {
      Rcpp::Rcout << "vec (first 20): ";
      for (int i = 0; i < 20; ++i) {
        Rcpp::Rcout << vec(i) << " ";
      }
      Rcpp::Rcout << std::endl;
    } else { 
      Rcpp::Rcout << "vec: ";
      for (int i = 0; i < vec.size(); ++i) {
        Rcpp::Rcout << vec(i) << " ";
      } 
      Rcpp::Rcout << std::endl;
    }
  } 

// ... overload for std::vector
void vprintV(
    const dVec& vec,
    bool verbose
  ) {
    if (vec.size() > 20) {
      Rcpp::Rcout << "vec (first 20): ";
      for (int i = 0; i < 20; ++i) {
        Rcpp::Rcout << vec[i] << " ";
      }
      Rcpp::Rcout << std::endl;
    } else { 
      Rcpp::Rcout << "vec: ";
      for (int i = 0; i < vec.size(); ++i) {
        Rcpp::Rcout << vec[i] << " ";
      } 
      Rcpp::Rcout << std::endl;
    } 
  }
  
// ... overload for sdouble
void vprintV(
    const sVec& vec,
    bool verbose
  ) {
    if (vec.size() > 20) {
      Rcpp::Rcout << "vec (first 20): ";
      for (int i = 0; i < 20; ++i) {
        Rcpp::Rcout << vec(i).val() << " ";
      }
      Rcpp::Rcout << std::endl;
    } else { 
      Rcpp::Rcout << "vec: ";
      for (int i = 0; i < vec.size(); ++i) {
        Rcpp::Rcout << vec(i).val() << " ";
      }
      Rcpp::Rcout << std::endl;
    } 
  }
