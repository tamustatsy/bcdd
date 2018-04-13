#' The MCMC algorithm for Constrained Bayes density deconvolution estimator.
#' 
#' The main function \code{bcdd_mcmc} in this package can be used to deconvolve 
#' a density with a unimodal and symmetric shape. The method is coined the name 
#' Constrained Bayes and is essentially a nonparametric approach using Bayesian
#' hierarchical models. The speed is bumped by utilizing functions written in Rcpp. 
#'
#' @param w A vector containing the proxy variables.
#' @param sd_u A vector containing the standard deviation of error, if the error
#' is homoscedastic, then it should be a scalar.
#' @param n.burnin A number (default 1000) specifies the number of burning steps. 
#' @param n.MCMC A number (default 5000) specifies the number of total MCMC steps.
#' @param hpar A list containing the hyperparameters in the algorithm, the 
#' default values are chosen if left empty. m is the concentration parameter and 
#' K is the number of categories of Dirichlet distribution prior for p. Xi_1, Xi_2 
#' are the shape and rate parameters of the Gamma prior for beta. lambda is the 
#' rate parameter of the exponential prior for alpha. shape.lambda is one parameter
#' in the proposal distribution for the MH step of alpha. tt corresponds to the 
#' lower bound that is prespecified for alpha.
#' @return A list containing the expected posterior density of the latent true 
#' variable (y) evaluated on a grid (x) chosen automatically in the function.
#' @examples 
#' df.t <- 5
#' n <- 500
#' x <- rt(n, df = df.t)
#' sd_u <- 1.29
#' u <- rnorm(n, mean = 0, sd = sd_u)               
#' w <- x + u
#' bcdd_mcmc(w, sd_u, n.burnin=100, n.MCMC=500)


bcdd_mcmc <- function(w, sd_u, n.burnin=1000, n.MCMC=5000, hpar=list(m=20, K=8, 
                        Xi_1=1.0, Xi_2=4.0, lambda=2, shape.lambda=2,tt=2.5)){
m <- hpar$m
K <- hpar$K
Xi_1 <- hpar$Xi_1
Xi_2 <- hpar$Xi_2
lambda <- hpar$lambda
shape.lambda <- hpar$shape.lambda
tt <- hpar$tt
n <- length(w)

################
# MCMC storage #
################
    
d.ordinates <- matrix(0, nrow = n, ncol = K)
n.kk.z <- numeric(K)
s.kk.theta <- numeric(K)
s.kk.ltheta <- numeric(K)
mh.ratio <- numeric(K)
log.mh.ratio <- numeric(K)
alpha.proposed <- numeric(K) 

##################
# Initialization #
##################

beta <- rgamma(K, shape = Xi_1, rate = Xi_2)
alpha <- truncdist::rtrunc(K, "exp", a = tt, b = Inf, rate = lambda)
p <- MCMCpack::rdirichlet(1, rep(m/K, K))
z <- sample(1:K, size = n, prob = p, replace = TRUE)
theta <- numeric(n)
for(ii in 1:n){
    theta[ii] <- rgamma(1, shape = alpha[z[ii]], rate = beta[z[ii]])
}
x <- runif(n, min = -1, max = 1)
x <- x*theta
x.init <- x

##############
# Start MCMC #
##############

density.x <- rep(0, lfg)
for(iii in 1:n.MCMC){
    if(iii%%100==0)
        print(iii)
    #-------update x_1, ..., x_n from a truncated normal distribution--------#
    ############################################################################
    #X_i = W_i - U_i                                                           #
    #1. generate U_i in a vectorized way:                                      #
    #the cutoff for U_i is [W_i - \theta_i, W_i + \theta_i]                    #
    #generate uniform variable between N(W_i - \theta_i) and N(W_i + \theta_i) #
    #get the quantiles corresponding to these uniform variables                #
    ############################################################################
    uu_cut_low <- w - theta
    uu_cut_up <- w + theta
    puu_cut_low <- pnorm(uu_cut_low, mean = 0, sd = sd_u)
    puu_cut_up <- pnorm(uu_cut_up, mean = 0, sd = sd_u)
    
    uu <- runif(n)
    uu_ext <- (puu_cut_up - puu_cut_low)*uu + puu_cut_low
    uu_ext[uu_ext == 1] <- 1 - 1E-5
    uu_ext[uu_ext == 0] <- 1E-5
    x <- w - qnorm(uu_ext, mean = 0, sd = sd_u)
    
    #----------update z_1, ..., z_n----------#
    d.ordinates = cand_dgamma(theta, s=alpha, r=beta)          #see Mydist.cpp
    z <- sample_pos_z(1:K, d.ordinates, p)                     #see Mysample.cpp
    
    #------update theta_1,...,theta_n from a truncated gamma distribution------#
    sg <- alpha[z] - 1
    rg <- beta[z]
    puu_cut_low <- pgamma(abs(x), shape = sg, rate = rg)
    uu <- runif(n)
    uu_ext <- (1 - puu_cut_low)*uu + puu_cut_low
    uu_ext[uu_ext == 1] <- 1 - 1E-5
    uu_ext[uu_ext == 0] <- 1E-5
    theta <- qgamma(uu_ext, shape = sg, rate = rg)
    
    ltheta <- log(theta)
    for(kk in 1:K){
        index <- (z == kk)
        n.kk.z[kk] <- sum(index)              #represents group sum of z
        s.kk.theta[kk] <- sum(theta[index])
        s.kk.ltheta[kk] <- sum(ltheta[index])
    }
    
    #----------update p_1, ..., p_K----------#
    p = MCMCpack::rdirichlet(1,m/K + n.kk.z)
    
    #----------update beta_1, ..., beta_K----------#
    beta <- rgamma(rep(1, K), shape = alpha*n.kk.z + Xi_1,
                   rate = s.kk.theta + Xi_2)
    
    #----------update alpha_1, ..., alpha_K----------#
    
    for(kk in 1:K){
        alpha.proposed[kk] <- truncdist::rtrunc(1,"gamma",a=tt,b=Inf,shape = shape.lambda,
                                     rate = shape.lambda/alpha[kk])
        tp_proposed <- 1 - pgamma(tt,shape.lambda,rate = shape.lambda/alpha[kk])
        tp_current <- 1 - pgamma(tt,shape.lambda,rate = shape.lambda/alpha.proposed[kk])
        term_1 <- shape.lambda*(alpha[kk]/alpha.proposed[kk]-alpha.proposed[kk]/alpha[kk])
        term_2 <- (alpha.proposed[kk] - alpha[kk])*(lambda - n.kk.z[kk]*
                    log(beta[kk]) - s.kk.ltheta[kk])
        log.mh.ratio[kk] <- (2*shape.lambda-1)*(log(alpha[kk])-log(alpha.proposed[kk]))
                        + n.kk.z[kk]*(lgamma(alpha[kk])-lgamma(alpha.proposed[kk]))
                        -term_1-term_2 + log(tp_proposed)-log(tp_current)
        mh.ratio[kk] = exp(log.mh.ratio[kk])
    }
    
    mh.ratio[which(is.nan(mh.ratio)==T)] = 0
    acc.prob = runif(K)
    inds.to.replace = (1:K)[acc.prob<mh.ratio]
    alpha[inds.to.replace] = alpha.proposed[inds.to.replace]
    
    lfg <- 200                                        #indicates length for grid
    integ_grid <- numeric(lfg)
    x.grid <- seq(0, max(abs(w)), length = lfg)
    if(iii > n.burnin){
        integ_grid <- 0
        p.vec <- as.vector(p)
        for(gg in 1:(lfg - 1)){
            integ_grid[gg] <- integrate(fun_integ, lower=x.grid[gg], upper=x.grid[gg+1],
                         p=p.vec, shape=alpha, rate=beta, subdivisions=2000)$value
        }
        integ_grid[lfg] <- integrate(fun_integ, lower = x.grid[lfg], upper = Inf,
                        p=p.vec, shape=alpha, rate=beta, subdivisions=2000)$value
        density.x <- density.x + rev(cumsum(rev(integ_grid)))
    }
}#MCMC iteration

density.x <- density.x/(n.MCMC - n.burnin)
return(list(x=x.grid, y=density.x))

}#bcdd_mcmc
