
// printing.cpp
#include "wspmm.h"

// Function for printing character messages
void vprint(
    const std::string& message,
    bool verbose
  ) {
    if (verbose) {
      Rcpp::Rcout << message << std::endl;
    }
  }

void vprintV(
    const CharacterVector& message,
    bool verbose
  ) {
    if (verbose) {
      Rcpp::Rcout << message << std::endl;
    }
  }

// Function for printing vectors
void print_Vec(
    const NumericVector& vec
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
void print_Vec(
    const IntegerVector& vec
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
void print_Vec(
    const dVec& vec
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
void print_Vec(
    const sVec& vec
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