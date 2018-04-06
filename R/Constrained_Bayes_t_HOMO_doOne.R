###########################################################################################################################
#The current file is built for batch runs. It can be run on its own by specifying lines 40-42 and comment out lines 37-38 # 
###########################################################################################################################


##################################
# Add source files and libraries #
##################################

library(coda)                                                              #required by MCMCpack
library(MCMCpack)
library(evd)                                                               #required by truncdist
library(truncdist)
library(Rcpp)
library(RcppArmadillo)
library(VGAM)                                                              #required by Mysample.cpp
library(batch)                                                             #required if conduct batch code for parallel

#--files containing functions for Constrained Bayes--#
source("../Myfunction.R")                                                  #this is used in computing the marginal density of X 
source("../ddsc_mcmc_prior1_cpp.R")                                        #this contains the main function to conduct the hybrid Gibbs sampler
sourceCpp("../Mysample.cpp")                                               #this is used in updating the group indicator variables z
sourceCpp("../Mydist.cpp")                                                 #also used in updating the group indicator variables z
#--files containing functions for Kernel deconvolution estimator--#
source("../PI_deconvUknownth4.r")
source("../fdecUknown.r")
source("../phiK2.r")
source("../rlap.R")
source("../outerop.r")
source("../kernel_decon_known_u_homo.R")


#####################
#Global parameters  #
#####################

parseCommandArgs()           #seed, n
set.seed(seed)

# seed <- 1100
# set.seed(seed)
# n <- 5000

#####################
# Tuning parameters #
#####################

m <- 20.0                    #prior on mixing probabilities
K <- 8                       #number of mixture components
Xi_1 <- 1.0                  #Xi_1, Xi_2 are for prior on rate of theta
Xi_2 <- 4.0
lambda <- 2                  #prior on shape of theta
shape.lambda <- 2            #proposal distribution, in the paper, it is fixed at 2 and not treated as hyperparameter
tt <- 2.5                    #truncated value for alpha_k and its proposal, the larger it is, the smoother
                             #(or flatter) the peak is   

########################################
# Parameters that are used in the MCMC #
########################################

n.burnin <- 1000
n.MCMC <- 5000

#--length of grid to plot density--#
lfg <- 200                                 #indicates length for grid
integ_grid <- numeric(lfg)


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

#####################################################
# Generate data x using a sym and unimodal density #
#####################################################

df.t <- 5
x <- rt(n, df = df.t)
sd_u <- 1.29
u <- rnorm(n, mean = 0, sd = sd_u)               
w <- x + u

x.grid <- seq(0, 11, length = lfg)

##################
# Initialization #
##################

beta <- rgamma(K, shape = Xi_1, rate = Xi_2)
alpha <- rtrunc(K, "exp", a = tt, b = Inf, rate = lambda)
p <- rdirichlet(1, rep(m/K, K))
z <- sample(1:K, size = n, prob = p, replace = TRUE)
theta <- numeric(n)
for(ii in 1:n){
  theta[ii] <- rgamma(1, shape = alpha[z[ii]], rate = beta[z[ii]])
}
x <- runif(n, min = -1, max = 1)                                                      #first generate uniform[-1,1] then scale by theta
x <- x*theta
x.init <- x

density.CB <- ddsc_mcmc(w, n.burnin, n.MCMC, n, K, sd_u, z, theta, beta, alpha, p,
                        m, lambda, shape.lambda, tt)

density.kernel <- kernel_homo(n = n, W = w, sigU = sd_u, xs = x.grid, lfg = 2*lfg)

#--compute the mean integrated absolute error--#
imae_CB <- sum(abs(density.CB - dt(x.grid, df = df.t))*(x.grid[2]-x.grid[1]))*2        #factor 2 for the other half
imse_CB <- sqrt(sum((density.CB - dt(x.grid, df = df.t))^2*(x.grid[2]-x.grid[1])))*2

density.Ker <- density.kernel$y2
xx.Ker <- density.kernel$xx

imse_Kernel <- sqrt(sum((density.Ker - dt(xx.Ker, df = df.t))^2*(xx.Ker[2]-xx.Ker[1])))
imae_Kernel <- sum(abs(density.Ker - dt(xx.Ker, df = df.t))*(xx.Ker[2]-xx.Ker[1]))

#--write the density estimates and other parameters in csv file--#
write.table(data.frame(seed, imae_CB, imse_CB, imae_Kernel, imse_Kernel, x.grid, density.CB, xx.Ker, density.Ker), 
            file=paste("./", n, "/seed",seed,"_t_HOMO.csv", sep = ""), append = TRUE, row.names = FALSE, col.names = FALSE, sep = ",")








