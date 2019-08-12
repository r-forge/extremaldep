#######################################################
### Authors: Giulia Marcon and Simone Padoan        ###
### Emails: giulia.marcon@phd.unibocconi.it,        ###
### simone.padoan@unibocconi.it                     ###
### Institution: Department of Decision Sciences,   ###
### University Bocconi of Milan                     ###
### File name: prior.r                              ###
### Description:                                    ###
### This file provides the random generation        ###
### function for the Bernstein representation       ###
### coefficients and the transformation function to ###
### move from eta to beta, and vice versa.          ###
### See Corollary 3.4 and Proposition 3.2 in        ### 
### Marcon et al. (2016)                            ###
### Last change: 15/08/2016                         ###
#######################################################

###################### START RANDOM GENERATION COEFFICIENTS ####################

################################################################################
# INPUT:                                                                     ###
# k is the polynomial order                                                  ###
# pm is a vector of point masses at zero and one                             ###
################################################################################ 

rcoef <- function(k, pm){
  if(missing(pm)){
    pm <- NULL
    pm$p0 <- pm$p1 <- 0
    warning('point masses p0 and p1 set equal to zero.')
  }
    p0 <- pm$p0
    p1 <- pm$p1

  kp <- k+1
  km <- k-1
  sample <- inf <- sup <- numeric(kp)
  
  sample[1] <- p0
  # i = 1,...,k-1
  for(i in 1:km){
    ip <- i+1
    im <- i-1
    ss <- sum(sample[1:i])
    inf[ip] <- max( sample[i], kp/2 - ss +(k-i)*(p1-1))
    sup[ip] <- min(1-p1, 1/(k-i) * ( kp/2 -1 - ss + p1) )
    if( (-1e-10 < inf[ip]-sup[ip]) & (1e-10 > inf[ip]-sup[ip]) ) sample[ip] <- sup[ip]
    else sample[ip] <- runif(1,inf[ip],sup[ip])
  }
  sample[kp] <- 1-p1
  
  beta <- net(sample, from = 'H')$beta
    
  return(list(eta=sample, beta=beta))
}
  
###################### END RANDOM GENERATION COEFFICIENTS ####################

###################### START COEFFICIENTS TRANSFORMATION ####################

################################################################################
# INPUT:                                                                     ###
# coef is a vector of the beta or eta coefficients, depending which extremal ###
# dependence function is specified, A or H, respectively.                    ###
################################################################################ 

net <- function (coef, from=c('A','H')) 
{
  
  if(from=='A'){
    k <- length(coef) - 2
    kp <- k+1
    
    # From A(w) to H(w)
    # beta: j=0,...,k+1
    # eta: j=0,...,k
    eta <- numeric(kp)
    eta <- 1/2 + kp/2 * diff(coef)
    
    return(list(eta=eta))
  }
  
  if(from=='H'){
    k <- length(coef) - 1
    kp <- k+1
    
    # From H(w) to A(w)
    # eta: j=0,...,k
    # beta: j=0,...,k+1
    beta <- numeric(kp+1)
    beta[1] <- 1
    beta[2:(kp+1)] <- 1/kp * ( 2 * cumsum(coef[1:kp]) +kp-0:k-1 )
    
    return(list(beta=beta))
  }
}

###################### END COEFFICIENTS TRANSFORMATION ####################
