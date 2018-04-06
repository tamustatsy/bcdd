//The function mysample is used to update/sample the indicators z_1, ..., z_n

//[[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
arma::vec mysample(NumericVector K, const arma::mat& ord, const arma::colvec& pvec,
                       bool replace){
    //RNGScope scope;             
    const int size = 1;
    arma::vec z(ord.n_rows);
    NumericVector prob; 
    NumericVector temp;
    
    for(int ii = 0; ii < ord.n_rows; ii++){
        prob = pvec%trans(ord.row(ii));
        temp = Rcpp::RcppArmadillo::sample(K, size, replace, prob);
        z(ii) = temp(0);
    }
    
    return z;
}