#######################################################
### Authors: Boris Beraner and Simone Padoan        ###
### Emails: borisberanger@gmail.com,                ###
### simone.padoan@unibocconi.it                     ###
### Institutions:                                   ###
### (BB) School of Mathematics and Statistics       ###
###      UNSW Sydney, Australia                     ###
### (SP) Department of Decision Sciences            ###
###      University Bocconi of Milan                ###
### File name: UnivExtremeQ.R                       ###
### Description:                                    ###
### This file enables to compute univairate extreme ###
### quantiles as proposed in Beranger et al. (2019) ###
### Last change: 19/07/2019                         ###
#######################################################

UniExtQ <- function(data, P=NULL, method="bayesian", U=NULL, cov=as.matrix(rep(1,length(data))), QatCov=NULL, param0=NULL, sig0=NULL, nsim=NULL, p=0.234, optim.meth="BFGS", control=NULL, ...){
  # Computes univariate extreme quantiles
  #
  #
  # If Bayesian need the parameters: data, cov, param0, U, p, sig0, nsim
  #
  # If frequentist need the parameters: data, cov, param0, U, optim.meth,... (extra parameters of optim function)
  
  if(missing(data)){stop("Missing data input")}
  if(method != "bayesian" && method != "frequentist" ){stop("Method should be `bayesian` or `frequentist` ")}
  if(!is.vector(data)){stop("Data must be a vector")}
  if(is.null(P)){stop("Must specify the probability associated with the quantile(s)")}
  if(nrow(cov) != length(data)){stop("Wrong dimensions of the covariates")}
  if(is.null(QatCov) && ncol(cov)>1){stop("QatCov argument required")}
  if(!is.null(QatCov) && (ncol(cov)>1) && (ncol(QatCov) != (ncol(cov)-1))){stop("Wrong dimensions of QatCov")}
  if(ncol(cov)==1){
    Cov <- as.matrix(1)
  }else{
    Cov <- cbind(rep(1, nrow(QatCov)),QatCov)
  }
  if(is.null(U)){ 
    U <- quantile(data, probs=0.9, type=3)
    cat("U set to 90% quantile by default \n")
  }
  if(is.null(param0)){stop("Need to specify param0 parameter")} # param0 is the starting value for both methods
  if(length(param0) != (ncol(cov)+2)){stop("Wrong length of initial parameter vector")}
  
  kn <- sum(data>U)/length(data)
  
  if(method=="bayesian"){
    
    if(is.null(sig0) || is.null(nsim) || is.null(p)){stop("Missing initial values")}
    if(is.null(nsim)){stop("Missing number of replicates in algorithm")}
    
    mcmc <- RWMH.gev(data = data, cov = cov, U=U, param0 = param0, p=p, sig0 = sig0, nsim = nsim)
    meanacc<-rep(NA, length(mcmc$acc.vec))
    for (j in c(1:length(mcmc$acc.vec))) {
      meanacc[j] = mean(mcmc$acc.vec[round(j/2) :j])
    }

    index <- c(1:2000,seq(2001,length(mcmc$sig.vec), by=100))

    par(mar=c(4.4,4.8,0.5,0.5))
    plot( cbind(index, mcmc$sig.vec[index]^2) , type="l", col=3, ylim=c(0, max(mcmc$sig.vec^2)), ylab=expression(tau^2), xlab="Iterations", cex.lab=1.8, cex.axis=1.8, lwd=2)
    abline(h=0, lwd=2)
    plot(cbind(index, meanacc[index]), type="l", col=2, ylim=c(0,1), ylab="Acceptance Probability", xlab="Iterations", cex.lab=1.8, cex.axis=1.8, lwd=2)
    abline(h=p, lwd=2)

    burn <- readline("What is the burn-in period? \n")
    burn <- as.numeric(unlist(strsplit(burn, ",")))
    if(burn<0 || burn>nsim){
      burn <- readline(paste("The burn-in period needs to be between 0 and", nsim, ", new value? \n",sep=""))
      burn <- as.numeric(unlist(strsplit(burn, ",")))
    }
  
    post_sample <- mcmc$param_post[-c(1:burn),] 
    mu.ps <-  post_sample[,1:(ncol(cov))] %*% t(Cov)
    sig.ps <- post_sample[, ncol(cov)+1]
    gam.ps <- post_sample[, ncol(cov)+2]
    
    if(ncol(cov)==1){
      Q.est <- matrix(ncol=length(P), nrow=nrow(post_sample))
      for(j in 1:length(P)){
        Q.est[,j] <- mu.ps + sig.ps * ((kn/P[j])^gam.ps-1) / gam.ps
      }  
    }else{
      Q.est <- list(length=nrow(QatCov))
      for(j in 1:length(P)){
        Q.est[[j]] <- mu.ps + sig.ps * ((kn/P[j])^gam.ps-1) / gam.ps  
      }  
    }
    
    return(list(Q.est=Q.est, post_sample=post_sample, burn=burn, straight.reject=mcmc$straight.reject[-c(1:(burn-1))], sig.vec=mcmc$sig.vec))
    
  }else if(method=="frequentist"){
    
    if(is.null(optim.meth)){stop("Need to specify an optimisation method in frequentist setting")}
    
    optimfun <- function(para){
      return(llik.cens(data = data, cov=cov, param = para, u = U))
    }
    
    if(!("control" %in% names(sys.call())) ){
      control <- list(fnscale=-1)
    }else{
      if(!("fnscale" %in% names(control))){
        control[["fnscale"]] <- -1  
      }
    }
    
    est.para <- optim(param0, optimfun, method=optim.meth, control=control,...)
    mu <- est.para$par[1:(ncol(cov))] %*% t(Cov)
    sigma <- est.para$par[ncol(cov)+1] 
    gamma <- est.para$par[ncol(cov)+2]
    
    if(ncol(cov)==1){
      Q.est <- as.vector(mu) + sigma * ((kn/P)^gamma-1) / gamma  
    }else{
      Q.est <- matrix(nrow=length(P), ncol=nrow(Cov))
      for(j in 1:length(P)){
        Q.est[j,] <- mu + sigma * ((kn/P[j])^gamma-1) / gamma
      }  
    }
    
    if("hessian" %in% names(est.para)){
      return(list(Q.est=Q.est, kn=kn, est = est.para$par, VarCov = solve(-est.para$hessian) ))}
    else{
      return(list(Q.est=Q.est, kn=kn, est = est.para$par ))
    }
    
  }
  
  
}

### Hidden functions

# Random Walk Metropolis-Hastings

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
    
    if( any(U < (mu-param_new[ncol(cov)+1]/param_new[ncol(cov)+2]) ) || any(param_new[ncol(cov)+ 1:2]<0) ){
      straight.reject[i] <- 1
      acc.vec[i] <- 0
    }else{
      # compute the acceptance probability
      ratio <- min( exp(llik.cens(data=data, cov=cov, param = param_new, u = U) - llik.cens(data = data, cov=cov, param = param, u = U) + log(param[ncol(cov)+ 1]) - log(param_new[ncol(cov)+ 1]) ), 1)
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
          #restart the algorithm
          cat("restart the program at", i, "th iteration", "\n")
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
  
  return(list(param_post=sparam, accepted=accepted, straight.reject=straight.reject , nsim=nsim, sig.vec=sig.vec, sig.restart=sig.restart, acc.vec=acc.vec ))
}

# Censored likelihood
llik.cens <- function(data, cov, param, u){
  
  if(!is.vector(data)){stop("Data must be a vector")}
  if(length(param) != (ncol(cov)+2)){stop("Wrong length of parameter vector")}
  mu <- cov %*% param[1:(ncol(cov))]
  sigma <- param[ncol(cov)+1]; gamma <- param[ncol(cov)+2];
  
  if(any(c(mu, sigma, gamma)<=0)){return(-1e300)} # Focus on mean > 0  for now
  # Not considering Weibull yet
  if( (gamma > 0) && (u < (mu - sigma/gamma) ) ){ return(-1e300)}
  
  kn <- sum(data>u)/length(data)
  
  part1 <- sum( log(kn) - log(sigma) + kn * log(pgev(q=data[data>u], loc=mu[data>u], scale=sigma, shape=gamma)) - (1/gamma + 1) * log(1 + gamma * (data[data>u] - mu[data>u]) /sigma) )
  part2 <- kn * sum(log(pgev(q=u, loc=mu[data<=u], scale=sigma, shape=gamma)))
  
  return(part1 + part2)
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



