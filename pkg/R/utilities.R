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

### random generation function for the Bernstein representation coefficients
### See Corollary 3.4 in Marcon et al. (2016)  

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

### transformation function to move from eta to beta, and vice versa.
### See Proposition 3.2 in Marcon et al. (2016)        

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

