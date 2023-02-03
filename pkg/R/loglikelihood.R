#######################################################
### Authors: Giulia Marcon and Simone Padoan        ###
### Emails: giulia.marcon@phd.unibocconi.it,        ###
### simone.padoan@unibocconi.it                     ###
### Institution: Department of Decision Sciences,   ###
### University Bocconi of Milan                     ###
### File name: loglikelihood.r                      ###
### Description:                                    ###
### This file provides the log-likelihood function  ###
### corresponding to the bivariate max-stable       ###
### distribution, using the angular distribution    ###
### or the Pickands dependence function,            ###
### in Bernstein form.                              ###
### See eq. 3.16 in Marcon et al. (2016)            ###
### Last change: 15/08/2016                         ###
#######################################################

###################### START LOG-LIKELIHOOD ####################################

################################################################################
# INPUT:                                                                     ###
# coef is a vector of the eta coefficients                                   ###
# k is the polynomial order                                                  ###
# data is a list of:                                                         ###
#      z = 1/x + 1/y (Frechet scale)                                         ###
#      w = angular of data                                                   ###
#      r = radius of the data                                                ###
#      w2 = w^2                                                              ###
#      r2 = r^2                                                              ###
#      den = x^2, y^2                                                        ###
# pm is a vector of point masses at zero and one                             ###
# bpb, bpb1, bpb2 are matrices of the Bernstein polynomial basis             ###
# nsim is the number of the iterations of the chain                          ###
# if approx = TRUE, only coef, k and bp1 are needed.                         ###
################################################################################ 

llik <- function(coef, k, data = NULL, pm = NULL, bpb = NULL, bpb1, bpb2 = NULL, approx = FALSE){
  # z: 1/x + 1/y (Frechet scale)
  # w angular of data
  # r: radius of the data
  # w2 <- w^2
  # r2 <- r^2
  # den <- x^2, y^2
  
  if(approx){
    # compute the angular density:
    h <- k*c(bpb1 %*% diff(coef))
    llik <- log(h)
  }
  
  else{
    p0 <- pm$p0
    p1 <- pm$p1
    # derive the betas coefficients
    beta <- net(coef, from='H')$beta
  
    # compute the difference of the coefficients:
    beta1 <- diff(beta); beta2 <- diff(beta1)
  
    # compute the Pickands and its derivatives:
    A <- c(bpb %*% beta)
    A1 <- c((k+1) * (bpb1 %*% beta1))
    A2 <- c(k * (k+1) * (bpb2 %*% beta2))
  
    # GEV distribution
    lG <- -data$z*A
  
    C <- data$w*A1
    A1.s <- (A - C[,1])*(A + C[,2]) / data$den[,2]
    A2.s <- A2*data$w2[,1]/data$r
    llik <- lG + log(A1.s+A2.s) - log(data$den[,1])
  }
    return(sum(llik))
}

###################### END LOG-LIKELIHOOD ####################################


###################### START CONSTRAINED MAXIMUM LOG-LIKELIHOOD FITTING ####################################
###################### FOR MAXIMA AND THRESHOLD EXCEEDANCES             ####################################
constmle <- function(data, k, w, start=NULL, type="maxima", r0=NULL, q=NULL){
  ################################
  # START auxilary functions
  
  llik_pickands_min <- function(coef){
    # compute the difference of the coefficients:
    beta1 <- diff(coef); beta2 <- diff(beta1)
    
    # compute the Pickands and its derivatives:
    A <- c(bpb %*% coef)
    A1 <- c(k * (bpb1 %*% beta1))
    A2 <- c(k * (k-1) * (bpb2 %*% beta2))
    
    # GEV distribution
    lG <- -pseudo$z*A
    
    C <- pseudo$w*A1
    A1.s <- (A - C[,1])*(A + C[,2]) / pseudo$den[,2]
    A2.s <- A2*pseudo$w2[,1]/pseudo$r
    if(any((A1.s+A2.s)<0)){
      llik <- 1e10
    }else{
      llik <- lG + log(A1.s+A2.s) - log(pseudo$den[,1])
      llik <- -sum(llik)  
    }
    
    if(is.na(llik)) llik <- 1e10
    return(llik)
  }
  
  llik_h_beta <- function(coef){
    # define the eta coefficients
    #eta <- net(beta=coef,from='A')$eta
    eta <- 1/2 + k*diff(coef)/2
    # compute the angular density:
    h <- (k-1)*c(bpb2 %*% diff(eta))
    # check if it is the same to dh in CODE_Simulation...
    # is it faster?
    if(any(h<0)){
      llik <- 1e100
    }else{
      llik <- log(h)
      llik <- -sum(llik)      
    }
  
    if(is.na(llik)) llik <- 1e100
    return(llik)
  }
  
  confun <- function(x) return(constraint$r - constraint$R%*%x)
  
  # END auxilary functions
  ################################
  kp <- k+1
  km <- k-1
  # define the linear constraints
  constraint <- constraints(k)
  # define the setup for the estimation
  pseudo <- list(z=rowSums(1/data), r=rowSums(data), w=data/rowSums(data), 
                 r2=(rowSums(data))^2, w2=(data/rowSums(data))^2, den=data^2)
  # define options
  #opts <- list(algorithm="NLOPT_LN_NELDERMEAD", xtol_rel=1e-25, maxeval=100000)
  opts <- list(algorithm="NLOPT_LN_COBYLA", xtol_rel=1e-25, maxeval=100000)
  #opts <- list(algorithm="NLOPT_LN_NEWUOA_BOUND", xtol_rel=1e-25, maxeval=100000)
  #opts <- list(algorithm="NLOPT_LN_SBPLX", xtol_rel=1e-08)
  # Compute starting values
  #if(is.null(start)){
  #  eta0 <- prior_eta_sampler(km)$eta
  #  beta0 <- net(eta=eta0, from='H')$beta}
  #else beta0 <- start
  if(is.null(start)) beta0 <- rcoef(km)$beta
  else beta0 <- start

  # maximize the constraint log-likelihood function
  if(type=="maxima"){
    bpb <- bpb1 <- bpb2 <- NULL
    bpb <- bp2d(pseudo$w, k)
    bpb1 <- bp2d(pseudo$w, km)
    bpb2 <- bp2d(pseudo$w, km-1)
    fit <- nloptr(beta0, eval_f=llik_pickands_min, eval_g_ineq=confun, 
                  lb=rep(0,kp), ub=rep(1,kp), opts=opts)
  }
  if(type=="rawdata"){
    bpb2 <- NULL 
    bpb2 <- bp2d(pseudo$w[pseudo$r>r0, ], km-1)
    fit <- list()
    fvalue <- c(Inf, numeric(14))
    fit[[1]] <- nloptr(beta0, eval_f=llik_h_beta, eval_g_ineq=confun, 
                       lb=rep(0,kp), ub=rep(1,kp), opts=opts)
    for(i in 2:15){
      fit[[i]]<- nloptr(fit[[i-1]]$solution[kp:1], eval_f=llik_h_beta, eval_g_ineq=confun, 
                        lb=rep(0,kp), ub=rep(1,kp), opts=opts)
      fvalue[i] <- fit[[i]]$objective}
    
    ind <- c(which(fvalue==min(fvalue)))[1]
    fit <- fit[[ind]]
  }
  beta <- fit$solution[kp:1]
  A <- beed(data, cbind(w, 1-w), 2, 'md', 'emp', k, beta=beta, plot=FALSE)
  
  return(list(beta=beta, A=A$A, status=fit$status, start=beta0, iterations=fit$iterations, message=fit$message))
}

###################### FOR MAXIMA AND THRESHOLD EXCEEDANCES             ####################################
###################### END CONSTRAINED MAXIMUM LOG-LIKELIHOOD FITTING ####################################
