diagnostics <- function(mcmc){
  
  dim <- 1
  if("eta" %in% names(mcmc)) dim <- 2
   
  if(dim == 1){
    return(diagnostics1d(mcmc))
  }else if(dim==2){
    return(diagnostics2d(mcmc))
  }
    
}


###############################################################################
###############################################################################
## Hidden functions
###############################################################################
###############################################################################

# Function to display the diagnostics plots from the univariate adaptive MCMC scheme.

diagnostics1d <- function(mcmc){
  
  # mcmc should be the ouput fom the RWMH.gev function
  # or at least a list containing acc.vec, sig.vec and nsim
  
  index <- c(1:2000,seq(2001,mcmc$nsim, by=100))  
  
  meanacc<-rep(NA, mcmc$nsim)
  for (j in c(1:mcmc$nsim)) {
    meanacc[j] = mean(mcmc$acc.vec[round(j/2) :j])
  }
  
  par(mfrow=c(1,2), mar=c(4.4,4.8,0.5,0.5))
  plot(index, mcmc$sig.vec[index]^2 , type="l", col=3, ylim=c(0, max(mcmc$sig.vec^2)), ylab=expression(tau^2), xlab="Iterations", cex.lab=1.8, cex.axis=1.8, lwd=2)
  abline(h=0, lwd=2)
  plot(index, meanacc[index], type="l", col=2, ylim=c(0,1), ylab="Acceptance Probability", xlab="Iterations", cex.lab=1.8, cex.axis=1.8, lwd=2)
  abline(h=0.234, lwd=2)
  par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1)) # reset graphical parameters
}

# Function to display the diagnostics plots from the bivariate adaptive MCMC scheme.

diagnostics2d <- function(mcmc){
  
  mar.fit <- FALSE
  if("mar1" %in% names(mcmc)) mar.fit <- TRUE

  index <- c(1:2000,seq(2001,mcmc$nsim, by=100))  
  
  if(mar.fit){
    
    meanacc <- meanacc.mar1 <- meanacc.mar2 <- rep(NA, mcmc$nsim)
    for (j in c(1:mcmc$nsim)) {
      meanacc.mar1[j] = mean(mcmc$acc.vec.mar1[round(j/2) :j])
      meanacc.mar2[j] = mean(mcmc$acc.vec.mar2[round(j/2) :j])
      meanacc[j] = mean(mcmc$acc.vec[round(j/2) :j])
    }
    
    oldpar1 <- par(mfrow=c(2,3), mar=c(4.4,4.8,0.5,0.5))
    on.exit(par(oldpar1))
    
    plot( cbind(index, mcmc$sig1.vec[index]^2) , type="l", col=3, ylim=c(0, max(mcmc$sig1.vec^2)), ylab=expression(tau[1]), xlab="Iterations", cex.lab=1.8, cex.axis=1.8, lwd=2)
    abline(h=0, lwd=2)
    plot( cbind(index, mcmc$sig2.vec[index]^2) , type="l", col=3, ylim=c(0, max(mcmc$sig2.vec^2)), ylab=expression(tau[2]), xlab="Iterations", cex.lab=1.8, cex.axis=1.8, lwd=2)
    abline(h=0, lwd=2)
    plot( cbind(index, mcmc$k[index]), type="l", col=3, ylim=c(0, max(mcmc$k)+0.5), ylab=expression(kappa), xlab="Iterations", cex.lab=1.8, cex.axis=1.8, lwd=2 )
    
    plot(cbind(index, meanacc.mar1[index]), type="l", col=2, ylim=c(0,1), ylab="Acceptance Probability", xlab="Iterations", cex.lab=1.8, cex.axis=1.8, lwd=2)
    abline(h=0.234, lwd=2)
    plot(cbind(index, meanacc.mar2[index]), type="l", col=2, ylim=c(0,1), ylab="Acceptance Probability", xlab="Iterations", cex.lab=1.8, cex.axis=1.8, lwd=2)
    abline(h=0.234, lwd=2)
    plot(cbind(index, meanacc[index]), type="l", col=2, ylim=c(0,1), ylab="Acceptance Probability", xlab="Iterations", cex.lab=1.8, cex.axis=1.8, lwd=2)
    
  }else{
    
    meanacc <- rep(NA, mcmc$nsim)
    for (j in c(1:mcmc$nsim)) {
      meanacc[j] = mean(mcmc$acc.vec[round(j/2) :j])
    }
    
    oldpar1 <- par(mfrow=c(1,2), mar=c(4.4,4.8,0.5,0.5))
    on.exit(par(oldpar1))
    
    plot( cbind(index, mcmc$k[index]), type="l", col=3, ylim=c(0, max(mcmc$k)+0.5), ylab=expression(kappa), xlab="Iterations", cex.lab=1.8, cex.axis=1.8, lwd=2 )
    plot(cbind(index, meanacc[index]), type="l", col=2, ylim=c(0,1), ylab="Acceptance Probability", xlab="Iterations", cex.lab=1.8, cex.axis=1.8, lwd=2)

  }
  
}  

