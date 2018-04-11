//[[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' @title mydgamma
//' @description The function \code{mydgamma} is meant for computing the matrix
//' as the input for the function \code{mysample} where the family of densities 
//' are gamma. See also \code{mysample}.
//' @param x a vector of values at which the gamma densities are evaluated.  
//' @param s a vector of shape parameters.
//' @param r a vector of rate parameters.
//' @return A matrix with (i,j) component equals the j-th gamma density evaluated
//' at x_i.
//' @examples
//' N <- 500
//' K <- 3
//' alpha <- c(2.928, 2.979, 2.754)
//' beta <- c(0.124, 0.238, 0.045)
//' theta <-  rgamma(N, shape=2.6, rate=0.3)
//' d.ordinates <- mydgamma(theta, alpha, beta)
//' @export 
//[[Rcpp::export]]
arma::mat mydgamma( NumericVector x, NumericVector s, NumericVector r){
    // Environment st("package:stats");
    // Function dg = st["dgamma"];
    arma::mat res(x.size(), s.size()); 
    int n = s.size() ;
    for( int i=0; i<n; i++) {
        NumericVector temp = Rcpp::dgamma(x, s[i], 1/r[i], 0);
        res.col(i) = as<arma::vec>(temp);
    }
    return res ;
}

