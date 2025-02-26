
// logic.cpp
#include "wspc.h"

// Return logical vector giving elements of left which match right
LogicalVector eq_left_broadcast(
    const CharacterVector& left,
    const String& right
  ) {
    int n = left.size();
    LogicalVector out(n);
    for (int i = 0; i < n; i++) {
      out[i] = left[i] == right;
    }
    return out;
  }

// ... overload 
LogicalVector eq_left_broadcast(
    const IntegerVector& left,
    const int& right
  ) {
    int n = left.size();
    LogicalVector out(n);
    for (int i = 0; i < n; i++) {
      out[i] = left[i] == right;
    }
    return out;
  }

// ... overload
LogicalVector eq_left_broadcast(
    const NumericVector& left,
    const double& right
  ) {
    int n = left.size();
    LogicalVector out(n);
    for (int i = 0; i < n; i++) {
      out[i] = left[i] == right;
    }
    return out;
  }

// Quantifiers

bool any_true(
    const LogicalVector& x
  ) {
    for (int i = 0; i < x.size(); i++) {
      if (x[i]) {return true;}
    }
    return false;
  }

bool all_true(
    const LogicalVector& x
) {
  for (int i = 0; i < x.size(); i++) {
    if (!x[i]) {return false;}
  }
  return true;
}