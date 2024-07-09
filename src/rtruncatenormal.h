#ifndef RTRUNCATENORMAL_H
#define RTRUNCATENORMAL_H

#include <RcppArmadillo.h>

// Declare the function
arma::vec rtruncnorm_cpp(int n, arma::vec a, arma::vec b, arma::vec mean, arma::vec sd);
#endif // RTRUNCATENORMAL_H
