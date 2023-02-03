####################### START RETURN VALUES ####################################

################################################################################
# P(Y_1 > y_1,Y_2 > y_2) = int_0^1 min(w/y_1, 1-w/y_2) h(w) dw               ###
# INPUT:                                                                     ###
# out is the output of fExtDep.np function                                   ###
# summary.mcmc is the output of summary.bbeed function                       ###
# y = c(y_1, y_2) vector/matrix of thresholds                                ###
################################################################################

returns <- function(out, summary.mcmc, y, plot=FALSE, labels=NULL, data=NULL){
  # P(X>x,Y>y) = int_0^1 min(w/x,1-w/y) h(w) dw
  # y = c(y_1, y_2) vector of thresholds
  
  # Check if margins were fitted and transform to unit Frechet if needed
  if("mar1" %in% names(out)){
    z <- cbind( sapply(y[,1], function(x) trans.UFrech(c(x,summary.mcmc$mar1.mean)) ),
                sapply(y[,2], function(x) trans.UFrech(c(x,summary.mcmc$mar2.mean)) ) )
  }else{
    z <- y
  }
  
  # Subroutine
  returns.int <- function(y, eta){
    k <- length(eta) - 1
    p0 <- eta[1]
    p1 <- 1-eta[k+1]
    
    v <- y[1]/sum(y)
    j <- 1:k
    P<- diff(eta) * (j * pbeta(v, j+1, k-j+1) / y[1] + 
                       (k-j+1) * pbeta(1-v, k-j+2, j) / y[2])
    
    return( 2*sum(P)/(k+1))
  }
  
  id1 <- which(out$k==summary.mcmc$k.median)
  idb1 <- id1[which(id1>=summary.mcmc$burn)]
  
  etalow <- apply(out$eta[idb1,1:(summary.mcmc$k.median+1)],2,quantile,.05)
  etaup <- apply(out$eta[idb1,1:(summary.mcmc$k.median+1)],2,quantile,.95)
  etamean <- colMeans(out$eta[idb1,1:(summary.mcmc$k.median+1)])
  etamedian <- apply(out$eta[idb1,1:(summary.mcmc$k.median+1)], 2, quantile,.5)
  
  uno <- due <- tre <- quattro <- nrow(y)
  for(i in 1:nrow(y)){
    uno[i] <- returns.int(y=z[i,], etalow)
    due[i] <- returns.int(y=z[i,], etaup)
    tre[i] <- returns.int(y=z[i,], etamean)
    quattro[i] <- returns.int(y=z[i,], etamedian)
  }

  rrmcmc <- rowMeans(cbind(uno,due,tre,quattro))
  
  if(plot) plot_ExtDep.np(type = "returns", out=out, summary.mcmc=summary.mcmc, 
                          y=y, probs=rrmcmc, labels=labels, data=data)
    
  return(rrmcmc)
}

#################### END RETURN VALUES #########################################


