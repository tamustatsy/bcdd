//[[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;
using namespace arma;

//' @title sample_pos_z
//' @description The function \code{sample_pos_z} is used to update/sample the 
//' indicators z_1, ..., z_n. The core function is the basic sample function. The 
//' component probability equals a density evaluated at the latent variable coexisting 
//' with z and a probability that is associated with the (latent) true prior of 
//' these component probability.  
//' @param K a value that gives the total category that z varies from.
//' @param ord a matrix with (i,j) component equals the prespecified density evaluated
//' at the latent variable coexisting with z_i.
//' @param pvec a vector of probabilities that is associated with the (latent) 
//' true prior of these component probability.
//' @return A vector containing the posterior samples z_1,...,z_n.
//' @examples
//' set.seed(123)
//' N <- 500
//' K <- 3
//' alpha <- c(2.928, 2.979, 2.754)
//' beta <- c(0.124, 0.238, 0.045)
//' p <- c(0.069, 0.758, 0.173)
//' theta <-  rgamma(N, shape=2.6, rate=0.3)
//' d.ordinates <- matrix(0, nrow=N, ncol=K)
//' for(jj in 1:K){
//'     d.ordinates[,jj] <- dgamma(theta, shape=alpha[jj], rate=beta[jj])
//' }
//' z <- sample_pos_z(1:K, d.ordinates, p) 
//' @export
//[[Rcpp::export]]
arma::vec sample_pos_z(NumericVector K, const arma::mat& ord, const arma::colvec& pvec){
    //RNGScope scope;             
    const int size = 1;
    arma::vec z(ord.n_rows);
    NumericVector prob; 
    NumericVector temp;
    
    for(int ii = 0; ii < ord.n_rows; ii++){
        prob = pvec%trans(ord.row(ii));
        temp = Rcpp::RcppArmadillo::sample(K, size, prob);
        z(ii) = temp(0);
    }
    
    return z;
}