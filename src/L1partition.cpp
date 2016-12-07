#include <Rcpp.h>
using namespace Rcpp;
#include "cutil.h"

// [[Rcpp::export]]
NumericVector l1partition(NumericVector x, double epsilon, double ratio, long seed) {
  int n = x.length();
  NumericVector Rhist(n);

  int *X = new int[n];
  int *hist = new int[n];
  for (int i = 0; i < n; i++) {
    X[i] = (int)x[i];
  }

  L1partition(hist, n, X, n, epsilon, ratio, seed);
  for (int i = 0; i < n; i++) {
    Rhist[i] = hist[i];
  }
  delete [] X;
  delete [] hist;

  return(Rhist);
}

// [[Rcpp::export]]
NumericVector l1partition_approx(NumericVector x, double epsilon, double ratio, long seed) {
  int n = x.length();
  NumericVector Rhist(n);

  int *X = new int[n];
  int *hist = new int[n];
  for (int i = 0; i < n; i++) {
    X[i] = (int)x[i];
  }

  L1partition_approx(hist, n, X, n, epsilon, ratio, seed);
  for (int i = 0; i < n; i++) {
    Rhist[i] = hist[i];
  }
  delete [] X;
  delete [] hist;

  return(Rhist);
}
