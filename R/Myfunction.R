#' Calculate the integrand for the marginal distribution of X over theta.
#' 
#' The (posterior) marginal distribution of X equals the integral of the conditional 
#' distribution of X given theta, uniform over [-theta, theta], weighted by the 
#' (posterior) distribution of theta, a mixture of Gammas. The function fun_integ
#' explicitly computes the integrand as a function of theta.
#' @param x a number (vector) specifies the value of the domain argument. 
#' @param pi a vector containing the component probabilities of the Gammas 
#' @param shape, rate numeric vectors corresponding to shape and rate parameters
#' of Gammas
#' @return The (posterior) marginal distribution of X (over theta) at the input 
#' parameter x. 
#' @examples 
#' pi <- c(0.2,0.5,0.3)
#' shape <- c(2.3,3.3,4.3) 
#' rate <- c(1.5,2.5,3.5)
#' theta <- 0.1 
#' fun_integ(theta, pi, shape, rate)
#' theta <- c(0.1,0.5,3)
#' fun_integ(theta, pi, shape, rate)
#' @export 
fun_integ <- function(x, pi, shape, rate){
    temp <- numeric(length(x))
    temp <- 0
    for(ii in 1:length(pi))
        temp <- temp + (2*x)^{-1}*pi[ii]*dgamma(x, shape=shape[ii], rate=rate[ii])
    return(temp)
}






