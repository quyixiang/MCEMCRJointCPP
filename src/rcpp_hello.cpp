#include "rtruncatenormal.h"

#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

// This is a simple function using Rcpp that creates an R list
// containing a character vector and a numeric vector.
//
// Learn more about how to use Rcpp at:
//
//   https://www.rcpp.org/
//   https://adv-r.hadley.nz/rcpp.html
//
// and browse examples of code using Rcpp at:
//
//   https://gallery.rcpp.org/
//

// [[Rcpp::export]]
List rcpp_hello() {
  CharacterVector x = CharacterVector::create("foo", "bar");
  NumericVector y   = NumericVector::create(0.0, 1.0);
  List z            = List::create(x, y);
  return z;
}

// // [[Rcpp::export]]
// arma::vec useTruncNorm() {
//   int n = 100;
//   arma::vec a(n), b(n), mean(n), sd(n);
//
//   // Initialize vectors
//   std::fill(a.begin(), a.end(), 0);    // Lower bounds initialized to 0
//   std::fill(b.begin(), b.end(), 5);    // Upper bounds initialized to 5
//   std::fill(mean.begin(), mean.end(), 2);  // Means initialized to 2
//   std::fill(sd.begin(), sd.end(), 1);  // Standard deviations initialized to 1
//
//   // Assuming rtruncnorm is correctly implemented and available
//   return rtruncnorm(n, a, b, mean, sd);
// }
