fGEV <- function(data, par.start, method="Frequentist", u, cov, 
                 optim.method="BFGS", optim.trace=0,  sig0, nsim){
  
  methods <- c("Frequentist", "Bayesian")
  if(!any(method ==  methods)){ stop("estimation method wrongly specified")}
  
  if(!is.vector(data)){stop("Data must be a vector")}
  if(missing(cov)){cov <- as.matrix(rep(1,length(data)))}
  if(!is.matrix(cov)){stop("'cov' must be a matrix")}
  
  # Handling NAs
  index.na <- is.na(data)
  if(any(index.na)){
    data <- data[!index.na]
    cov <- as.matrix(cov[!index.na,])
  }
  
  for(j in 1:ncol(cov)){
    index.na <- is.na(cov[,j])
    if(any(index.na)){
      data <- data[!index.na]
      cov <- as.matrix(cov[!index.na,])
    }
  }
  
  if(missing(par.start) && ncol(cov)==1){
    momest <- MomEst(data, length(data))
    location <- momest$location
    scale <- momest$scale
    shape <- momest$shape
    if(shape >= 0) shape <- .1
    else shape <- -.1
    par.start <- c(location, scale, shape)
  }else if(missing(par.start) && ncol(cov)>1){
    stop("A vector of starting values needs to be specified")
  }
   
  if(missing(u)) u <- NULL
  
  if(method=="Frequentist"){
    
    if(is.null(optim.method)){stop("Need to specify an optimisation method in frequentist setting")}
    
    optimfun <- function(para){
        return(llik.cens(data = data, cov=cov, param = para, u = u))
    }
    
    fit <- optim(par.start, optimfun, method=optim.method,
                    control=list(fnscale=-1, reltol=1e-14, maxit=1e8, trace=optim.trace), hessian=TRUE)
    
    varcov <- try(-solve(fit$hessian), silent = TRUE)
    if(!is.matrix(varcov)){
      varcov <- 'none'
      stderr <- 'none'
    }else{
      stderr <- diag(varcov)
      if(any(stderr < 0))
        stderr <- 'none'
      else
        stderr <- sqrt(stderr)
    }
    
    return(list(est=fit$par, varcov=varcov, stderr=stderr))
  
  }else if(method=="Bayesian"){
    
    mcmc <- RWMH.gev(data=data, cov=cov, param0=par.start, U=u, p=0.234, sig0=sig0, nsim=nsim)
    return(mcmc)
  }
  
}

##########################################
##########################################
##########################################
### Internal functions ###################
##########################################
##########################################
##########################################

#
# Method="Frequentist", u=NULL
#

GevLogLik <- function(data, numdata, param)
{
  ### Compute the log-likelihood:
  result <- .C('GevLogLik', as.double(data), as.integer(numdata), as.double(param),
               res=double(1), DUP=TRUE, NAOK=TRUE)$res
  return(result)
}

MomEst <- function(data, n)
{
  # Scale estimate:
  scale <- sqrt(6 * var(data)) / pi
  # Location estimate:
  location <- mean(data) - 0.58 * scale
  
  k <- round(.5 * n)
  mind <- min(data)
  if(mind < 0) data <- data - mind
  
  data <- sort(data)
  delta <- log(data[n-0:(k-1)] / data[n-k])
  M1 <- sum(delta) / k
  M2 <- sum(delta^2) / k
  # Shape estimate:
  shape <- M1 + 1 - .5 * (M2 / (M2 - M1^2))
  
  return(list(location=location, scale=scale, shape=shape))
}

#
# Method="Frequentist", u!=NULL
#

llik.cens <- function(data, cov, param, u){
  
  if(!is.vector(data)){stop("Data must be a vector")}
  if(length(param) != (ncol(cov)+2)){stop("Wrong length of parameter vector")}
  mu <- cov %*% param[1:(ncol(cov))]
  sigma <- param[ncol(cov)+1]; gamma <- param[ncol(cov)+2];
  
  if(is.null(u)){ #Block maxima approach
    
    if(sigma <=0){return(-1e300)}
    #if( any(data < (mu - sigma/gamma)) ){ return(-1e300)}
    
    #return( sum(dgev(x=data, loc=mu, scale=sigma, shape=gamma, log=TRUE)) )
    if(ncol(cov) == 1){
      return( GevLogLik(data=data, numdata=length(data), param=c(mu[1],sigma,gamma)) )
    }else{
      llik <- 0
      for(i in 1:nrow(cov)){
        llik <- llik + GevLogLik(data=data, numdata=length(data), param=c(mu[i],sigma,gamma))
      }
      return(llik)
    }

    #GevLogLik(data=data, numdata=length(data), param=c(mu[i],sigma,gamma))
    
  }else{ #threshold exceedances
    
    if(any(c(mu, sigma, gamma)<=0)){return(-1e300)} # Focus on mean > 0  for now
    # Not considering Weibull yet
    
    if( any(u < (mu - sigma/gamma)) ){ return(-1e300)}
    kn <- sum(data>u)/length(data)
    
    part1 <- sum( log(kn) - log(sigma) + kn * log(pGEV(q=data[data>u], loc=mu[data>u], scale=sigma, shape=gamma)) - (1/gamma + 1) * log(1 + gamma * (data[data>u] - mu[data>u]) /sigma) )
    part2 <- kn * sum(log(pGEV(q=u, loc=mu[data<=u], scale=sigma, shape=gamma)))
    
    return(part1 + part2)
  }
  
}

#
# Method="Bayesian"
#

RWMH.gev <- function(data, cov, param0, U, p, sig0, nsim){
  
  alpha  <- -qnorm(p/2)
  d <- length(param0) # Dimension of the vector of parameters
  if(d != (ncol(cov)+2)){stop("Wrong length of parameter vector")}
  sig <- sig0 # Initial value of sigma
  sig.vec <- sig
  sigMat <-diag(d)  # Initial proposal covarince matrix
  acc.vec <- rep(NA, nsim) # Vector of acceptances
  accepted <- rep(0, nsim)
  straight.reject <- rep(0, nsim) # Monitor the number of proposed parameters that don't respect the constraints
  sig.start<- sig
  sig.restart<- sig
  
  n0 = round(5/(p*(1-p)))
  iMax=100 # iMax is the max number of iterations before the last restart
  Numbig=0
  Numsmall=0
  param <- param0
  
  # set the chains
  sparam <- param
  # create a progress bar
  pb = txtProgressBar(min = 0, max = nsim, initial = 0, style=3) 
  
  for(i in 1:nsim){
    
    # simulation from the proposal (MVN):
    param_new <- as.vector(rmvnorm(1, mean = param, sigma = sig^2 * sigMat))
    
    mu <- cov %*% param_new[1:(ncol(cov))]
    
    if(is.null(U)){ # Block Maxima approach
      if( (any(data < (mu-param_new[ncol(cov)+1]/param_new[ncol(cov)+2]) ) && param_new[ncol(cov)+2]>0)
          || (any(data > (mu-param_new[ncol(cov)+1]/param_new[ncol(cov)+2]) ) && param_new[ncol(cov)+2]<0)
          || any(param_new[ncol(cov)+ 1]<0)
      ){
        #cat(paste("Staright reject with params:", param_new[1], ", ",  param_new[2], ", and ", param_new[3],"\n"))
        straight.reject[i] <- 1
        acc.vec[i] <- 0
      }
    }else{ # Threshold Exceedances
      if( any(U < (mu-param_new[ncol(cov)+1]/param_new[ncol(cov)+2]) ) || any(param_new[ncol(cov)+ 1:2]<0) ){
        straight.reject[i] <- 1
        acc.vec[i] <- 0
      }
    }
    
    if(straight.reject[i] != 1){
      # compute the acceptance probability
      llik.new <- llik.cens(data=data, cov=cov, param = param_new, u = U)
      llik.cur <- llik.cens(data = data, cov=cov, param = param, u = U)
      #cat(paste("Proposal:", param_new, "and likelihood: ", llik.new,"\n"))
      ratio <- min( exp(llik.new - llik.cur + log(param[ncol(cov)+ 1]) - log(param_new[ncol(cov)+ 1]) ), 1)
      acc.vec[i] <- ratio
      
      u <- runif(1)
      if(u<ratio){
        param <- param_new
        accepted[i] <- 1
      }
    }
    
    sparam <- rbind( sparam, param)
    
    # update covariance matrix with adaptive MCMC
    if (i > 100) {
      if (i==101) {
        sigMat=cov(sparam)
        thetaM=apply(sparam, 2, mean)
      } else
      {
        tmp=update.cov(sigMat = sigMat, i = i, thetaM = thetaM, theta = param, d = d)
        sigMat=tmp$sigMat
        thetaM=tmp$thetaM
      }
    }
    
    # Update sigma
    if (i>n0) {
      sig <- update.sig(sig = sig, acc = ratio, d = d, p = p, alpha = alpha, i = i)
      sig.vec <- c(sig.vec, sig)
      if ((i <= (iMax+n0)) && (Numbig<5 || Numsmall<5)) {
        Toobig<- (sig > (3*sig.start))
        Toosmall<-(sig < (sig.start/3))
        
        if (Toobig || Toosmall) {
          # restart the algorithm
          # message("\n restart the program at", i, "th iteration", "\n")
          sig.restart <- c(sig.restart, sig)
          Numbig <- Numbig + Toobig
          Numsmall <- Numsmall + Toosmall
          i <- n0
          sig.start <- sig
        }
      } #end iMax
    }
    
    print.i <- seq(0, nsim, by=100)
    if(i %in% print.i) setTxtProgressBar(pb, i)
    
  }
  
  close(pb)
  
  return(list(param_post=sparam, accepted=accepted, straight.reject=straight.reject , nsim=nsim, sig.vec=sig.vec, sig.restart=sig.restart, acc.vec=acc.vec ))
}

# Internal function to adjust the value of sigma (=log(theta)) in each iteration
update.sig <- function(sig, acc, d = d, p = p, alpha = alpha, i) {
  c = ((1-1/d) * sqrt(2*pi) * exp(alpha^2 / 2) / (2 * alpha) + 1 / (d*p*(1-p)) ) # Eq(7) of Garthwaite, Fan & Sisson (2016)
  Theta = log(sig)
  # Theta = Theta + c * (acc-p) / max(200, i/d)
  Theta = Theta + c * (acc-p) / i
  return(exp(Theta))
}

# Internal function to update the covariance matrix
update.cov<-function(sigMat, i, thetaM, theta, d){
  epsilon=1/i
  thetaM2=((thetaM*i)+theta)/(i+1)
  sigMat=(i-1)/i*sigMat + thetaM%*%t(thetaM)-(i+1)/i*thetaM2%*%t(thetaM2)+1/i*theta%*%t(theta) + epsilon*diag(d)
  return(list(sigMat=sigMat, thetaM=thetaM2))
}





