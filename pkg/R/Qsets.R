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

ellipse <- function(center=c(0,0), alpha=c(0,0), sigma=diag(2), df=1, prob=0.01, npoints=250, pos=FALSE)
{
  es <- eigen(sigma)
  e1 <- es$vec %*% diag(sqrt(es$val))
  
  if(!all(alpha==0)){
    h <- 2*log(1+exp(-1.544/sqrt(alpha%*%sigma%*%t(alpha))))
    r1 <- sqrt(qchisq(prob, 2))-h
  }else{
    r1 <- sqrt(2*qf(prob, 2, df))
  }
  
  theta <- seq(0, 2*pi, len=npoints)
  v1 <- cbind(r1 * cos(theta), r1 * sin(theta))
  pts <- t(center - (e1 %*% t(v1)))
  if(pos) pts <- pts[pts[,1]>=0 & pts[,2]>=0, ]
  return(pts)
  
}  

# Routine for hill-based estimator of the tail index

"Moment" <- function(data, 
                     k = 5:(sum(data>0)-1), 
                     CI.type = "Wald", 
                     CI.p = NULL, plot = TRUE, 
                     test = "s", alpha = 0.5,
                     ...) {
  
  X <- sort(data[data>0], decreasing=TRUE)
  n <- length(X)
  if (!all(k < n)) {
    k <- k[k < n]
    warning("Only those k for which X_{n-k:n} is positive are retained.", call. = FALSE)
  }
  std.err <- numeric(length(k))
  
  # --- Moment estimates
  
  l <- log(X[1:(max(k)+1)])
  s1 <- cumsum(l[1:max(k)])[k]
  s2 <-	cumsum((l[1:max(k)])^2)[k]
  M1 <- s1 / k - l[k+1]
  M2 <- s2 / k - 2 * l[k+1] * s1 / k + (l[k+1])^2
  Moment <- M1 + 1 - 0.5 / (1 - M1^2 / M2)
  
  # --- standard errors
  
  if (any(Moment >= 0)) {
    I <- Moment >= 0
    g <- Moment[I]
    std.err[I] <- 1 + g^2
  }
  if (any(Moment < 0)) {
    I <- Moment < 0
    g <- Moment[I]
    std.err[I] <- (1-g)^2 * (1-2*g) * (6*g^2 - g + 1) / ((1-3*g) * (1-4*g))
  }
  std.err <- sqrt(std.err/k)
  
  # --- Confidence intervals
  
  if (is.null(CI.p)) {
    CI <- NULL
  }
  else if (!is.numeric(CI.p)) {
    CI <- NULL
    warning("CI.p should be NULL or a number")
  }
  else if ((CI.p-0.5) * (1-CI.p) <= 0) {
    CI <- NULL
    warning("CI.p should be between 0.5 and 1")
  }
  else {
    z <- qnorm((CI.p+1)/2)
    CI <- array(0, dim=c(length(k),2))
    CI[,1] <- Moment - std.err * z
    CI[,2] <- Moment + std.err * z
  }
  
  
  # --- output list
  
  out <- list(n = n, k = k, threshold = X[k+1], estimate = Moment, 
              CI = CI, CI.type = CI.type, CI.p = CI.p, 
              std.err = std.err,
              data = deparse(substitute(data)),
              quantity = "gamma",
              method = "Moment")
  class <- "EVI"

  out <- structure(out, class = class)
  
  # --- plot
  if (plot) {
    if (is.na(charmatch("plot.EVI", ls()))) source("plot.EVI.R")
    plot(out, ...)
  }
  
  invisible(out)
  
}

