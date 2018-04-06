#####################################################################
#ddsc_mcmc is the main function to conduct the hybrid Gibbs sampler #
#####################################################################

ddsc_mcmc <- function(w, n.burnin, n.MCMC, n, K, sd_u, z, theta, beta, alpha, p,
                      m, lambda, shape.lambda, tt){
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
    d.ordinates = mydgamma(theta, s=alpha, r=beta)             #see Mydist.cpp
    z <- mysample(1:K, d.ordinates, p, TRUE)                   #see Mysample.cpp
    
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
    p = rdirichlet(1,m/K + n.kk.z)
    
    #----------update beta_1, ..., beta_K----------#
    beta <- rgamma(rep(1, K), shape = alpha*n.kk.z + Xi_1,
                   rate = s.kk.theta + Xi_2)
    
    #----------update alpha_1, ..., alpha_K----------#
    
    for(kk in 1:K){
        alpha.proposed[kk] <- rtrunc(1,"gamma",a=tt,b=Inf,shape = shape.lambda,
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
return(density.x)

}#ddsc_mcmc
