//The function mydgamma defined here are meant for computing the quantities which
//the (conditional) posterior probabilities of z_1, ..., z_n is proportional to

//[[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

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

