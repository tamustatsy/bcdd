###########################################################
#The function fun_integ is to calculate the integrand for #
#the marginal distribution of X over theta                #
###########################################################

fun_integ <- function(x, pi, shape, rate){
    temp <- numeric(length(x))
    temp <- 0
    for(ii in 1:length(pi))
        temp <- temp + (2*x)^{-1}*pi[ii]*dgamma(x, shape=shape[ii], rate=rate[ii])
    return(temp)
}






