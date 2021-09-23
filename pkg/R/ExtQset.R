#######################################################
### Authors: Boris Beraner and Simone Padoan        ###
### Emails: borisberanger@gmail.com,                ###
### simone.padoan@unibocconi.it                     ###
### Institutions:                                   ###
### (BB) School of Mathematics and Statistics       ###
###      UNSW Sydney, Australia                     ###
### (SP) Department of Decision Sciences            ###
###      University Bocconi of Milan                ###
### File name: ExtQset.R                            ###
### Description:                                    ###
### This file enables to compute bivairate extreme  ###
### quantile regions as proposed in Beranger et     ###
### al. (2019)                                      ###
###                                                 ###
### Last change: 22/07/2019                         ###
#######################################################

ExtQset <- function(data, P=NULL, method="bayesian", U=NULL, 
                    cov1=as.matrix(rep(1,nrow(data))), cov2=as.matrix(rep(1,nrow(data))),
                    QatCov1=NULL, QatCov2=NULL, mar=TRUE, par10=c(1,2,1), par20=c(1,2,1), sig10=1, sig20=1,
                    param0=NULL, k0=NULL, pm0=NULL, prior.k="nbinom", prior.pm="unif", 
                    hyperparam = list(mu.nbinom = 3.2, var.nbinom = 4.48), nsim=NULL,
                    lo=NULL, up=NULL, d=5){

  if(missing(data) || !is.matrix(data) || ncol(data) != 2){stop("Input data missing or not defined as a 2-column matrix")}
  if(is.null(P)){stop("Must specify the probability associated with the quantile(s)")}
  # if(method != "bayesian" && method != "frequentist" ){stop("Method should be `bayesian` or `frequentist` ")}
  if(method=="bayesian"){
    if(is.null(U) || !is.vector(U) ){ 
      U <- apply(data,2,function(x) quantile(x, probs=0.9, type=3))
      message("U set to 90% quantile by default on both margins \n")
    }
    if(is.vector(U) && length(U)!=2){
      U <- rep(U[1],2)
      message(paste("Warning, U incorrectly specified and set to U=(", U[1], ",", U[1],") \n",sep=""))
    }
    if(nrow(cov1) != nrow(data) || nrow(cov2) != nrow(data)){stop("Wrong dimensions of the covariates")}
    if(is.null(QatCov1) && ncol(cov1)>1){stop("QatCov1 argument required")}
    if(is.null(QatCov2) && ncol(cov2)>1){stop("QatCov2 argument required")}
    if(!is.null(QatCov1) && (ncol(cov1)>1) && (ncol(QatCov1) != (ncol(cov1)-1))){stop("Wrong dimensions of QatCov1")}
    if(!is.null(QatCov2) && (ncol(cov2)>1) && (ncol(QatCov2) != (ncol(cov2)-1))){stop("Wrong dimensions of QatCov2")}
    if(ncol(cov1)==1){
      Cov1 <- as.matrix(1)
    }else{
      Cov1 <- cbind(rep(1, nrow(QatCov1)),QatCov1)
    }
    if(ncol(cov2)==1){
      Cov2 <- as.matrix(1)
    }else{
      Cov2 <- cbind(rep(1, nrow(QatCov2)),QatCov2)
    }
    if(is.null(par10)){stop("Need to specify par10 parameter")} 
    if(length(par10) != (ncol(cov1)+2)){stop("Wrong length of initial first marginal parameter vector")}
    if(is.null(par20)){stop("Need to specify par20 parameter")} 
    if(length(par20) != (ncol(cov2)+2)){stop("Wrong length of initial second marginal parameter vector")}
    
    kn <- c(sum(data[,1]>U[1]), sum(data[,2]>U[2])) / nrow(data)
    
  }else if(method=="EDhK"){
    if(is.null(lo) || is.null(up)){stop("Arguments 'lo' and 'up' are required (for hill estimator of shape parameters) ")}
    if(lo>=up){stop("Argument 'lo' should be strictly less than 'up' ")}
    if(up>nrow(data)){stop("Argument 'up' cannot be greater than the number of observations ")}
    if(lo<1){stop("Argument 'lo' cannot be less than 1 ")}
  }else if(method=="frequentist"){
    if(is.null(lo) || is.null(up)){stop("Arguments 'lo' and 'up' are required (for hill estimator of shape parameters) ")}
    if(lo>=up){stop("Argument 'lo' should be strictly less than 'up' ")}
    if(up>nrow(data)){stop("Argument 'up' cannot be greater than the number of observations ")}
    if(lo<1){stop("Argument 'lo' cannot be less than 1 ")}
    if(missing(d)){stop("Argument d should be passed for frequentist estimation.")}
  }else{
    stop("Method should be `bayesian`, `EDhK` or `frequentist` ")
  }
  
  nw=100
  
  # Define the density function for the exponent measure:
  den <- function(w, beta, gamma){
    return(2 * dh(w, beta) * (2 * dh(w, beta) * (w)^(1-gamma[1]) * (1-w)^(1-gamma[2]) / gamma[1] / gamma[2])^(-1/(1+gamma[1]+gamma[2])))
  }
  
  # Bound function
  ghat_fun <- function(w, hhat, gamma1, gamma2){
    return((2*hhat * (1-w)^(1-gamma1) * (w)^(1-gamma2) / gamma1 / gamma2)^(1/(1+gamma1+gamma2)))
  }
    
  # 
  # Bayesian estimation
  #
  
  if(method=="bayesian"){
    
    if(is.null(sig10) || is.null(sig20)){stop("Missing initial values")}
    if(is.null(nsim)){stop("Missing number of replimessagees in algorithm")}
    
    if(mar){
      message("\n Preliminary on margin 1 \n")
      rwmh.mar1 <- RWMH.gev(data = data[,1], cov = cov1, param0 = par10, U=U[1], p=0.234, sig0 = sig10, nsim = 5e+4)
      message("\n Preliminary on margin 2 \n")
      rwmh.mar2 <- RWMH.gev(data = data[,2], cov = cov2, param0 = par20, U=U[2], p=0.234, sig0 = sig20, nsim = 5e+4)
      
      par10 <- apply(rwmh.mar1$param_post[-c(1:30000),], 2, mean)
      par20 <- apply(rwmh.mar2$param_post[-c(1:30000),], 2, mean)
      sig10 <- tail(rwmh.mar1$sig.vec,1)
      sig20 <- tail(rwmh.mar2$sig.vec,1)
    }
    
    ###############################################################
    # START Estimation of the extremal dependence (angular density)
    ###############################################################
    
    # Set the grid point for angles:
    w <- seq(0.00001, .99999, length=nw)
    
    m <- length(P)
    Qhat <- Qhat_up <- Qhat_low <- array(0, c(nw,2,m))
    
    # pm0 <- list(p0=0, p1=0); 
    # If specifying pm0 then need to also give the log densities for p0 and p1
    # p0 is unif(0,1/2) but p1 is uniform depending on k and p0
    
    message("\n Estimation of the extremal dependence and margins \n")
    mcmc <- cens.bbeed(data=data, cov1=cov1, cov2=cov2, pm0=pm0, param0=param0, k0=k0, par10=par10, par20=par20, sig10=sig10, sig20=sig20, kn=kn,
                       prior.k = prior.k, prior.pm=prior.pm, hyperparam=hyperparam, nsim=nsim, U=U, p.star=0.234)
    
    # Plots to decide on the burn-in period
    index <- c(1:2000,seq(2001,length(mcmc$sig1.vec), by=100))
    meanacc <- meanacc.mar1 <- meanacc.mar2 <- rep(NA, length(mcmc$acc.vec))
    for (j in c(1:length(mcmc$acc.vec))) {
      meanacc.mar1[j] = mean(mcmc$acc.vec.mar1[round(j/2) :j])
      meanacc.mar2[j] = mean(mcmc$acc.vec.mar2[round(j/2) :j])
      meanacc[j] = mean(mcmc$acc.vec[round(j/2) :j])
    }
    
    oldpar1 <- par(mfrow=c(2,3), mar=c(4.4,4.8,0.5,0.5))
    on.exit(par(oldpar1))

    plot( cbind(index, mcmc$sig1.vec[index]^2) , type="l", col=3, ylim=c(0, max(mcmc$sig1.vec^2)), ylab=expression(sigma[1]^2), xlab="Iterations", cex.lab=1.8, cex.axis=1.8, lwd=2)
    abline(h=0, lwd=2)
    plot( cbind(index, mcmc$sig2.vec[index]^2) , type="l", col=3, ylim=c(0, max(mcmc$sig2.vec^2)), ylab=expression(sigma[2]^2), xlab="Iterations", cex.lab=1.8, cex.axis=1.8, lwd=2)
    abline(h=0, lwd=2)
    plot( cbind(index, mcmc$k[index]), type="l", col=3, ylim=c(0, max(mcmc$k)+0.5), ylab=expression(k), xlab="Iterations", cex.lab=1.8, cex.axis=1.8, lwd=2 )
    
    plot(cbind(index, meanacc.mar1[index]), type="l", col=2, ylim=c(0,1), ylab="Acceptance Probability", xlab="Iterations", cex.lab=1.8, cex.axis=1.8, lwd=2)
    abline(h=0.234, lwd=2)
    plot(cbind(index, meanacc.mar2[index]), type="l", col=2, ylim=c(0,1), ylab="Acceptance Probability", xlab="Iterations", cex.lab=1.8, cex.axis=1.8, lwd=2)
    abline(h=0.234, lwd=2)
    plot(cbind(index, meanacc[index]), type="l", col=2, ylim=c(0,1), ylab="Acceptance Probability", xlab="Iterations", cex.lab=1.8, cex.axis=1.8, lwd=2)
    abline(h=0.234, lwd=2)
    
    # Set burn-in period
    
    burn <- readline("What is the burn-in period? \n")
    burn <- as.numeric(unlist(strsplit(burn, ",")))
    if(burn<0 || burn>nsim){
      burn <- readline(paste("The burn-in period needs to be between 0 and", nsim, ", new value? \n",sep=""))
      burn <- as.numeric(unlist(strsplit(burn, ",")))
    }
    
    stat_mcmc <- summary.bbeed(object=w , mcmc, burn) # Should look at diagnostic plots to decide on the value of burn
    
    hhat <- rbind(stat_mcmc$h.low,stat_mcmc$h.mean,stat_mcmc$h.up)
    
    muhat1_post <- stat_mcmc$mar1_post[,1:(ncol(Cov1))] %*% t(Cov1) # A npost by nrow(Cov1) matrix
    muhat2_post <- stat_mcmc$mar2_post[,1:(ncol(Cov2))] %*% t(Cov2) # A npost by nrow(Cov2) matrix
    sighat1_post <- stat_mcmc$mar1_post[,ncol(Cov1)+1]
    sighat2_post <- stat_mcmc$mar2_post[,ncol(Cov2)+1]  
    gamhat1_post <- stat_mcmc$mar1_post[,ncol(Cov1)+2]
    gamhat2_post <- stat_mcmc$mar2_post[,ncol(Cov2)+2]
    npost <- nsim-burn+1
    
    ###############################################################
    # START Estimation of the basic-set S and measure nu(S) :
    ###############################################################
    
    # Estimation of the bound
    ghat_post <- matrix(nrow = npost, ncol = nw)
    for(i in 1:nw) {
      # ghat_post[,i] <- (2*stat_mcmc$h_post[,i] * (w[i])^(1-gamhat1_post) * (1-w[i])^(1-gamhat2_post) / gamhat1_post / gamhat2_post)^(1/(1+gamhat1_post+gamhat2_post))
      ghat_post[,i] <- ghat_fun(w[i], hhat=stat_mcmc$h_post[,i], gamma1 = gamhat1_post, gamma2 = gamhat2_post)
    }                                                                                                                        
    ghat <- apply(ghat_post,2, func)
    
    # Estimation of the basic-set S:
    Shat_post <- array(0, c(nw,2,npost))
    for(j in 1:npost) Shat_post[,,j] <- cbind(ghat_post[j,]*w, ghat_post[j,]*(1-w)) 
    
    Shat <- array(0, c(nw,2,3))
    for(j in 1:3) Shat[,,j] <- cbind(ghat[j,]*w, ghat[j,]*(1-w)) 
    
    #Estimation of the measure nu(S):
    nuShat_post <- numeric(npost)
    nuShat <- NULL
    
    for(i in 1:npost) nuShat_post[i] <- integrate(Vectorize(function(t){den(t,beta=stat_mcmc$beta_post[[i]], gamma=c(gamhat1_post[i],gamhat2_post[i]))}), 0,1)$value
    nuShat <- func(nuShat_post)
    
    for(i in 1:length(P)){
      for(z in 1:nrow(Cov1)){
        
        tmp <- array(0, c(nw,2,npost))
        tmp2 <- array(0, c(nw,2,3)) # 3 because we compute the mean and 90% cred. regions
        
        for(j in 1:nw){
          
          tmp <- rbind( muhat1_post[,z] + sighat1_post / gamhat1_post  * (( kn[1] * nuShat_post*Shat_post[j,1,]/P[i])^(gamhat1_post) -1) ,
                        muhat2_post[,z] + sighat2_post / gamhat2_post  * (( kn[2] * nuShat_post*Shat_post[j,2,]/P[i])^(gamhat2_post) -1) )
          
          tmp2[j,,] <- cbind(c(tmp[1,which(tmp[1,]==func(tmp[1,])[1])[1]],tmp[2,which(tmp[2,]==func(tmp[2,])[1])[1]]), 
                             apply(tmp,1,mean),
                             c(tmp[1,which(tmp[1,]==func(tmp[1,])[3])[1]],tmp[2,which(tmp[2,]==func(tmp[2,])[3])[1]])) 
        }
        
        assign(paste("Qset_P",i,"_CovNum_",z,"_post",sep=""), tmp)
        assign(paste("Qset_P",i,"_CovNum_",z,sep=""), tmp2) # 3 because we compute the mean and 90% cred. regions
        
      }
      
    }
    
    return( c(mget(ls(pattern = "Qset_P")), list(ghat=ghat, Shat=Shat, nuShat=nuShat, burn=burn)))
    
  } # end Bayesian estimation  
  
  # 
  # EDhK estimator
  #
  
  if(method=="EDhK"){
    
    gamma1 = mean(Moment(data[,1], k=lo:up, plot=FALSE)$estimate)
    gamma2 = mean(Moment(data[,2], k=lo:up, plot=FALSE)$estimate)
    
    n <- nrow(data)
    k = c(1:n) 
    U = apply(data, 2, sort)
    f1 = function(k){(U[n-k+1,1])*(k/n)^gamma1}
    f2 = function(k){(U[n-k+1,2])*(k/n)^gamma2}
    c1 = mean(f1(k[lo:up])) 
    c2 = mean(f2(k[lo:up]))
    
    ## estimate the spectral measure and its density
    fi = AngularMeasure(data[,1], data[,2], k=nw, method="c", plot=FALSE)
    w = fi$weights
    loc = fi$angles
    lok = pi/2 - loc
    m = length(loc)
    h = 0.45
    
    ## Kernels
    kernK = function(x){if(x >= -1 && x <= 1){ret = (15/16)*((1-x^2)^2)}
      if(x < -1 || x > 1){ret = 0}
      ret}
    
    kernK1 = function(x){if(x >= -1 && x <= 1){ret = x*(15/16)*((1-x^2)^2)}
      if(x < -1 || x > 1){ret = 0}
      ret}
    
    kernK2 = function(x){if(x >= -1 && x <= 1){ret = (x^2)*(15/16)*((1-x^2)^2)}
      if(x < -1 || x > 1){ret = 0}
      ret}
    
    kernK = Vectorize(kernK)
    kernK1 = Vectorize(kernK1)
    kernK2 = Vectorize(kernK2)
    
    a0 = function(y){(15/16)*(y^5/5 - 2*y^3/3 + y + 8/15)}
    a1 = function(y){(15/16)*(y^6/6 -   y^4/2 + y^2/2 - 1/6)}
    a2 = function(y){(15/16)*(y^7/7 - 2*y^5/5 + y^3/3 + 8/105)}
    a0 = Vectorize(a0)
    a1 = Vectorize(a1)
    a2 = Vectorize(a2)
    
    kernBL = function(y){
      ly = length(y)    
      ret = c(1:ly)
      p = (y+loc/h)[1]
      i = 1
      while(i <= ly){       
        ret[i] = (  (  a2(p)-a1(p)*y[i]   )*kernK(y[i])  )/( a0(p)*a2(p)-(a1(p))^2 )
        i = i+1  } # end of while 
      ret} # end of function
    
    kernBR = function(y){
      ly = length(y)
      ret = c(1:ly)
      p = ((y*h+lok)/h)[1]
      i = 1
      while(i<=ly){
        ret[i] = (  (  a2(p)-(a1(p))*y[i]  )*kernK(y[i])  )/( a0(p)*a2(p)-(a1(p))^2 )
        i = i+1  } # end of while 
      ret} # end of function
    
    fi_hat = function(x){
      if(x >= h && x <= pi/2 - h){ret = sum((w/h)*kernK((x-loc)/h))}
      if(x >= 0 && x < h){ret = sum((w/h)*kernBL((x-loc)/h))}
      if(x >= pi/2 - h && x <= pi/2){ret = sum((w/h)*kernBR((pi/2-x-lok)/h))}
      ret}
    fi_hat = Vectorize(fi_hat)
    
    
    fi_hat0 = function(x){
      ret = (w/h)*kernK((x-loc)/h)
      ret = sum(ret)
      ret}
    fi_hat0 = Vectorize(fi_hat0)
    
    g = function(theta){
      (  ((2/pi  )*((cos(min(theta,pi/2-h)))^(1-gamma1)  )*(  (sin(max(theta,h)))^(1-gamma2)  ))/(gamma1*gamma2))^(-1/(gamma1+gamma2+1))
    }
    
    ## \nu(S)
    nuS_integrand = function(x){g(x)*fi_hat(x)}
    nuS_integrand = Vectorize(nuS_integrand)
    nuS_hat = integrate(nuS_integrand, 0, pi/2)$value
    
    ## S
    S_r = function(theta){ 1/g(theta)}
    S_r = Vectorize(S_r)
    
    theta = seq(0, pi/2, pi/200)
    Srr = S_r(theta)
    
    xS = Srr*cos(theta)
    yS = Srr*sin(theta)
    
    for(i in 1:length(P)){
      assign(paste("xn_hat",i,sep=""), c1*(((nuS_hat*xS/P[i]))^gamma1) )
      assign(paste("yn_hat",i,sep=""), c2*(((nuS_hat*yS/P[i]))^gamma2) )
    }
  
    return(mget(ls(pattern = "n_hat")))
    
  } # End EDhK estimator
  
  # 
  # Frequentist estimation
  #
  
  if(method=="frequentist"){
    
    gamhat1 = mean(Moment(data[,1], k=lo:up, plot=FALSE)$estimate)
    gamhat2 = mean(Moment(data[,2], k=lo:up, plot=FALSE)$estimate)
    
    n <- nrow(data)
    k = c(1:n) 
    U = apply(data, 2, sort)
    f1 = function(k){(U[n-k+1,1])*(k/n)^gamhat1}
    f2 = function(k){(U[n-k+1,2])*(k/n)^gamhat2}
    c1 = mean(f1(k[lo:up])) 
    c2 = mean(f2(k[lo:up]))
    
    ## estimate the spectral measure and its density
    
    if(is.null(U) || !is.vector(U) ){ 
      q=0.9
    }else{
      q <- sum(data[,1]<=U[1])/nrow(data)
    }
    
    # Compute empirical distributions: 
    y1 <- ecdf(data[,1])
    y2 <- ecdf(data[,2])
    
    # Transform marginal distributions into unit-frechet:
    fdata <- cbind(1/(-log(y1(data[,1]))), 1/(-log(y2(data[,2]))))
    rdata <- rowSums(fdata) # Radial components
    wdata <- fdata[,1]/rdata # Angular components
    r0 <- quantile(rdata, probs=q) # Set high threhold
    extdata <- fdata[rdata>=r0,] # Extract data from the tail
    
    # Set the grid point for angles:
    w <- seq(0.00001, .99999, length=nw)
    
    m <- length(P)
    Qhat <- Qhat_up <- Qhat_low <- array(0, c(nw,2,m))
    
    # Compute starting value:
    Aest_proj <- beed(extdata, cbind(w, 1-w), 2, "md", "emp", d) 
    
    # Optimize the angular density:
    Aest_mle <- constmle(fdata, d, w, start=Aest_proj$beta, type="rawdata", r0=r0)
    
    # Compute the angular density:
    hhat <- sapply(1:nw, function(i, w, beta) dh(w[i], beta), w, Aest_mle$beta)
    
    # START Estimation of the basic-set S and measure nu(S) :
    
    # Estimation of the bound
    ghat <- ghat_fun(w=w, hhat=hhat, gamma1=gamhat1, gamma2=gamhat2)
    
    # Estimation of the basic-set S:
    xS <- ghat*w
    yS <- ghat*(1-w)
    Shat <- cbind(xS,yS)
    
    #Estimation of the measure nu(S):
    nuShat <- adaptIntegrate(den, 0 , 1, beta=Aest_mle$beta, gamma=c(gamhat1,gamhat2))$integral
    
    # START Estimation of quantile regions:
    for(i in 1:m)
      Qhat[,,i] <- cbind(c1*(((nuShat*xS/P[i]))^gamhat1),
                         c2*(((nuShat*yS/P[i]))^gamhat2))  
    
    return(list(hhat=hhat, ghat=ghat, Shat=Shat, nuShat=nuShat, 
                Qhat=Qhat, gamhat=c(gamhat1,gamhat2), uhat=c(c1,c2)))
    
  } # End frequentist estimation
  
}

###
### Hidden functions for Bayesian estimation
###

cens.bbeed <- function(data, cov1, cov2, pm0, param0, k0, par10, par20, sig10, sig20, kn, prior.k = c('pois','nbinom'), prior.pm = c('unif', 'beta'), hyperparam, nk=70, nsim, U=NULL, p.star){
  
  if(nrow(data) != nrow(cov1) || nrow(data) != nrow(cov2)){stop("Different number of observations/covariates")}
  
  if(length(par10) != (ncol(cov1)+2)){stop("Wrong parameter length for margin 1")}
  if(length(par20) != (ncol(cov2)+2)){stop("Wrong parameter length for margin 2")} 
  if(ncol(cov1)==1 && all(cov1==1)){
    Par10 <- t(as.matrix(par10))
  }else{
    Par10 <- cbind( cov1 %*% par10[1:(ncol(cov1))], rep(par10[ncol(cov1)+1], nrow(cov1) ), rep(par10[ncol(cov1)+2], nrow(cov1) )  )
  }
  if(ncol(cov2)==1 && all(cov2==1)){
    Par20 <- t(as.matrix(par20))
  }else{
    Par20 <- cbind( cov2 %*% par20[1:(ncol(cov2))], rep(par20[ncol(cov2)+1], nrow(cov2) ), rep(par20[ncol(cov2)+2], nrow(cov2) )  )
  }  
  
  if(is.null(U)){
    U <- apply(data,2,function(x) quantile(x, probs=0.9, type=3))
    message("U set to 90% quantile by default on both margins \n")
  }
  
  data.cens <- data
  data.cens[data.cens[,1]<=U[1],1] <- U[1]
  data.cens[data.cens[,2]<=U[2],2] <- U[2]
  
  if(ncol(cov1)==1 && all(cov1==1) && ncol(cov2)==1 && all(cov2==1)){
    data.censored <- apply(data.cens,1, function(x) ((x[1]==U[1]) && (x[2]==U[2])) )
    censored <- which(data.censored==1)[1]
    n.censored <- sum(data.censored)
    uncensored <- which(data.censored==0)
    data.cens <- cbind(data.cens[c(censored,uncensored),], c(n.censored, rep(1,length(uncensored))))
    xy <- as.matrix(data[c(censored,uncensored),])
  }else{
    xy <- as.matrix(data)
  }  
  
  if(ncol(cov1)==1 && all(cov1==1)){
    data.cens.uv <- t( apply(data.cens, 1, function(val) c(marg(val[1], kn=kn[1], par = Par10, der=FALSE), marg(val[2], kn=kn[2], par = Par20, der=FALSE)) ) )
  }else{
    data.cens.uv <- t( apply(cbind(data.cens, Par10, Par20), 1, function(val) c(marg(val[1], par = val[3:5], kn = kn[1], der=FALSE), marg(val[2], par = val[6:8], kn = kn[2], der=FALSE)) ) )
  }
  r.cens.uv <- rowSums(data.cens.uv)
  w.cens.uv <- data.cens.uv/r.cens.uv
  data0 <- list(r.cens.uv=r.cens.uv, w.cens.uv=w.cens.uv, xy = as.matrix(xy), data.cens = as.matrix(data.cens), data.cens.uv = as.matrix(data.cens.uv) )
  
  # derive the polynomial bases:
  bpb <- bp(cbind(w.cens.uv[,2], w.cens.uv[,1]), k0 + 1)
  bpb1 <- bp(cbind(w.cens.uv[,2], w.cens.uv[,1]), k0)
  bpb2 <- bp(cbind(w.cens.uv[,2], w.cens.uv[,1]), k0 - 1)
  
  # Checks:
  Check <- check.bayes(prior.k = prior.k, prior.pm = prior.pm, hyperparam = hyperparam, pm0 = pm0, param0 = param0, k0 = k0)
  prior.k <- Check$prior.k
  prior.pm <- Check$prior.pm
  hyperparam <- Check$hyperparam
  pm0 <- Check$pm0
  param0 <- Check$param0
  k0 <- Check$k0
  a <- Check$a
  b <- Check$b
  mu.pois <- Check$mu.pois
  pnb <- Check$pnb
  rnb <- Check$rnb
  
  n0 = round(5/(p.star*(1-p.star)))
  iMax=100 # iMax is the max number of iterations before the last restart
  Numbig1 <- Numbig2 <- 0
  Numsmall1 <- Numsmall2 <- 0
  
  n <- nrow(data)

  acc.vec.mar1 <- acc.vec.mar2 <- acc.vec <- rep(NA, nsim) # vector of acceptances
  accepted.mar1 <- accepted.mar2 <- accepted <- rep(0, nsim)
  straight.reject1 <- straight.reject2 <- rep(0,nsim) # to monitor the number of straight rejections of the marginal parameters (need to fulfil conditions)
  
  smar1 <- par1 <- par10
  smar2 <- par2 <- par20
  Par1 <- Par10
  Par2 <- Par20
  pm <- pm0
  param <- param0
  spm <- c(pm$p0, pm$p1) 
  sk <- k <- k_new <- k0
  seta <- array(0, dim=c(nsim+1,nk+2))
  seta[1,1:(k+1)] <- param$eta
  
  if(k==3) q <- 0.5 else q <- 1
  
  alpha  <- -qnorm(p.star/2);
  d <- length(par1); # dimension of the vector of marginal parameters
  sig1 <- sig10; sig2 <- sig20; # initial value of sigma
  sig1.vec <- sig1; sig2.vec <- sig2;
  sig1Mat <- sig2Mat <- diag(d)  #initial proposal covariance matrix
  sig1.start<- sig1; sig2.start<- sig2;
  sig1.restart<- sig1; sig2.restart<- sig2; 
  
  # to store the polynomial bases:
  bpb_mat <- NULL 
  
  pb = txtProgressBar(min = 0, max = nsim, initial = 0, style=3) #creates a progress bar
  print.i <- seq(0, nsim, by=100)
  
  for(i in 1:nsim){
    
    ###
    ### Margin 1
    ###
    
    par1_new <- as.vector(rmvnorm(1, mean = par1, sigma = sig1^2 * sig1Mat))
    if(ncol(cov1)==1 && all(cov1==1)){
      Par1_new <- t(as.matrix(par1_new))
    }else{
      Par1_new <- cbind( cov1 %*% par1_new[1:(ncol(cov1))], rep(par1_new[ncol(cov1)+1],nrow(cov1)), rep(par1_new[ncol(cov1)+2],nrow(cov1))  )
    }
        
    if( any(U[1] < (Par1_new[,1]-Par1_new[,2]/Par1_new[,3]) ) || any(Par1_new[,3]<0) ){
      straight.reject1[i] <- 1
      acc.vec.mar1[i] <- 0
    }else{
      # Marginalised data
      if(ncol(cov1)==1 && all(cov1==1)){
        data.cens.uv_new <- cbind( sapply(data.cens[,1], function(val) marg(val, par = Par1_new, kn = kn[1] , der=FALSE)), data0$data.cens.uv[,2] )
      }else{
        data.cens.uv_new <- cbind( apply(cbind(data.cens[,1], Par1_new), 1, function(val) marg(val[1], par = val[2:4], kn = kn[1] , der=FALSE)), data0$data.cens.uv[,2] )
      }
      r.cens.uv_new <- rowSums(data.cens.uv_new)
      w.cens.uv_new <- data.cens.uv_new/r.cens.uv_new
      
      data_new <- list(xy = xy, data.cens = as.matrix(data.cens), data.cens.uv = data.cens.uv_new, w.cens.uv = w.cens.uv_new, r.cens.uv = r.cens.uv_new)
      
      # derive the polynomial bases:
      bpb_new <- bp(cbind(w.cens.uv_new[,2], w.cens.uv_new[,1]), k + 1)
      bpb1_new <- bp(cbind(w.cens.uv_new[,2], w.cens.uv_new[,1]), k)
      bpb2_new <- bp(cbind(w.cens.uv_new[,2], w.cens.uv_new[,1]), k - 1)
      
      ratio.mar1 <- min( exp(cens.llik(data = data_new, coef = param$eta, bpb = bpb_new, bpb1 = bpb1_new, bpb2 = bpb2_new, thresh = U, par1 = Par1_new, par2 = Par2, kn = kn)
                             - cens.llik(data = data0, coef = param$eta, bpb = bpb, bpb1 = bpb1, bpb2 = bpb2, thresh = U, par1 = Par1, par2 = Par2, kn = kn) 
                             + log(Par1[2]) - log(Par1_new[2])), 1)
      acc.vec.mar1[i] <- ratio.mar1
      
      u.mar1 <- runif(1)
      
      if(u.mar1 < ratio.mar1){
        data0 <- data_new
        bpb <- bpb_new
        bpb1 <- bpb1_new
        bpb2 <- bpb2_new
        par1 <- par1_new
        Par1 <- Par1_new
        accepted.mar1[i] <- 1
      }
    }
    
    smar1 <- rbind(smar1, par1)
    
    ###
    ### Margin 2
    ###
    
    par2_new <- as.vector(rmvnorm(1, mean = par2, sigma = sig2^2 * sig2Mat))
    if(ncol(cov2)==1 && all(cov2==1)){
      Par2_new <- t(as.matrix(par2_new))
    }else{
      Par2_new <- cbind( cov2 %*% par2_new[1:(ncol(cov2))], rep(par2_new[ncol(cov2)+1],nrow(cov2)), rep(par2_new[ncol(cov2)+2],nrow(cov2))  )
    }
    
    if( any(U[2] < (Par2_new[,1]-Par2_new[,2]/Par2_new[,3]) ) || any(Par2_new[,3]<0) ){
      straight.reject2[i] <- 1
      acc.vec.mar2[i] <- 0
    }else{
      # Marginalised data
      if(ncol(cov2)==1 && all(cov2==1)){
        data.cens.uv_new <- cbind( data0$data.cens.uv[,1], sapply(data.cens[,2], function(val) marg(val[1], par = Par2_new, kn = kn[2], der=FALSE)) )
      }else{
        data.cens.uv_new <- cbind( data0$data.cens.uv[,1], apply(cbind(data.cens[,2], Par2_new), 1, function(val) marg(val[1], par = val[2:4], kn = kn[2], der=FALSE)) )
      }
      
      r.cens.uv_new <- rowSums(data.cens.uv_new)
      w.cens.uv_new <- data.cens.uv_new/r.cens.uv_new
      data_new <- list(xy = xy, data.cens = as.matrix(data.cens), data.cens.uv = data.cens.uv_new, w.cens.uv = w.cens.uv_new, r.cens.uv = r.cens.uv_new)
      
      # derive the polynomial bases:
      bpb_new <- bp(cbind(w.cens.uv_new[,2], w.cens.uv_new[,1]), k + 1)
      bpb1_new <- bp(cbind(w.cens.uv_new[,2], w.cens.uv_new[,1]), k)
      bpb2_new <- bp(cbind(w.cens.uv_new[,2], w.cens.uv_new[,1]), k - 1)
      
      # compute the acceptance probability of 
      ratio.mar2 <- min( exp(cens.llik(data = data_new, coef = param$eta, bpb = bpb_new, bpb1 = bpb1_new, bpb2 = bpb2_new, thresh = U, par1 = Par1, par2 = Par2_new, kn= kn) 
                             - cens.llik(data = data0, coef = param$eta, bpb = bpb, bpb1 = bpb1, bpb2 = bpb2, thresh = U, par1 = Par1, par2 = Par2, kn = kn) 
                             + log(Par2[2]) - log(Par2_new[2])), 1)
      acc.vec.mar2[i] <- ratio.mar2
      
      u.mar2 <- runif(1)
      
      if(u.mar2 < ratio.mar2){
        data0 <- data_new
        bpb <- bpb_new
        bpb1 <- bpb1_new
        bpb2 <- bpb2_new
        par2 <- par2_new
        Par2 <- Par2_new
        accepted.mar2[i] <- 1
      }
      
    }  
    
    smar2 <- rbind(smar2, par2)
    
    ###
    ### Dependence structure
    ###
    
    # polynomial order
    k_new <- prior_k_sampler(k)
    
    # derive the polynomial bases:
    bpb_new <- bp(cbind(w.cens.uv_new[,2], w.cens.uv_new[,1]), k_new + 1)
    bpb1_new <- bp(cbind(w.cens.uv_new[,2], w.cens.uv_new[,1]), k_new)
    bpb2_new <- bp(cbind(w.cens.uv_new[,2], w.cens.uv_new[,1]), k_new - 1)
    
    if(prior.k=="pois"){p <- dpois(k_new-3,mu.pois)/dpois(k-3,mu.pois)}
    if(prior.k=="nbinom"){p <- dnbinom(k_new-3,prob=pnb,size=rnb)/dnbinom(k-3,prob=pnb,size=rnb)}
    
    if(k_new==nk){
      message('maximum value k reached')
      break
    }
    
    # point masses
    pm_new <- prior_p_sampler(a=a, b=b, prior.pm=prior.pm)
    while(check.p(pm_new,k_new)) pm_new <- prior_p_sampler(a=a, b=b, prior.pm=prior.pm)
    
    # polynomial coefficients
    param_new <- rcoef(k=k_new, pm=pm_new)
    
    # compute the acceptance probability of 
    ratio <- min( exp(cens.llik(data = data0, coef = param_new$eta, bpb = bpb_new, bpb1 = bpb1_new, bpb2 = bpb2_new, thresh = U, par1 = Par1, par2 = Par2, kn = kn) + log(q) -
                        cens.llik(data = data0, coef = param$eta, bpb = bpb, bpb1 = bpb1, bpb2 = bpb2, thresh = U, par1 = Par1, par2 = Par2, kn = kn) + log(p)), 1)
    acc.vec[i] <- ratio
    
    u <- runif(1)
    if(u < ratio){
      pm <- pm_new
      param <- param_new
      bpb <- bpb_new
      bpb1 <- bpb1_new
      bpb2 <- bpb2_new
      k <- k_new
      if(k==3) q <- 0.5 else q <- 1
      accepted[i] <- 1
    }
    
    sk <- c(sk,k)
    spm <- rbind(spm, c(pm$p0, pm$p1))
    seta[i+1,1:(k+1)] <- param$eta
    
    if(i %in% print.i) setTxtProgressBar(pb, i) 
    
    ###
    ### Margins: update covariance matrices with adaptive MCMC
    ###
    
    if (i > 100) {
      if (i==101) {
        # 1st margin
        sig1Mat <- cov(smar1)
        thetaM1 <- apply(smar1, 2, mean)
        # 2nd margin
        sig2Mat <- cov(smar2)
        thetaM2 <- apply(smar2, 2, mean)
      } else
      {
        # 1st margin
        tmp1 <- update.cov(sigMat = sig1Mat, i = i, thetaM = thetaM1, theta = par1, d = d)
        sig1Mat <- tmp1$sigMat
        thetaM1 <- tmp1$thetaM
        # 2nd margin
        tmp2 <- update.cov(sigMat = sig2Mat, i = i, thetaM = thetaM2, theta = par2, d = d)
        sig2Mat <- tmp2$sigMat
        thetaM2 <- tmp2$thetaM
      }
    }
    
    ###
    ### Margins: update sigmas (sig1 and sig2)
    ###
    
    if (i>n0) {
      
      # Margin 1
      sig1 <- update.sig(sig = sig1, acc = ratio.mar1, d = d, p = p.star, alpha = alpha, i = i)
      sig1.vec <- c(sig1.vec, sig1)
      
      if ((i <= (iMax+n0)) && (Numbig1<5 || Numsmall1<5) ) {
        
        Toobig1 <- (sig1 > (3*sig1.start))
        Toosmall1 <-(sig1 < (sig1.start/3))
        
        if (Toobig1 || Toosmall1) {
          #restart the algorithm
          message("\n restart the program at", i, "th iteration", "\n")
          sig1.restart <- c(sig1.restart, sig1)
          Numbig1 <- Numbig1 + Toobig1
          Numsmall1 <- Numsmall1 + Toosmall1
          i <- n0
          sig1.start <- sig1
        }
        
      } #end iMax mar 1
      
      # Margin 2
      sig2 <- update.sig(sig = sig2, acc = ratio.mar2, d = d, p = p.star, alpha = alpha, i = i)
      sig2.vec <- c(sig2.vec, sig2)
      
      if ((i <= (iMax+n0)) && (Numbig2<5 || Numsmall2<5) ) {
        
        Toobig2 <- (sig2 > (3*sig2.start))
        Toosmall2 <-(sig2 < (sig2.start/3))
        
        if (Toobig2 || Toosmall2) {
          #restart the algorithm
          message("\n restart the program at", i, "th iteration", "\n")
          sig2.restart <- c(sig2.restart, sig2)
          Numbig2 <- Numbig2 + Toobig2
          Numsmall2 <- Numsmall2 + Toosmall2
          i <- n0
          sig2.start <- sig2
        }
        
      } #end iMax mar 2
      
    }
    
  }
  
  return(list(pm=spm, eta=seta, k=sk, mar1=smar1, mar2=smar2, 
              accepted.mar1=accepted.mar1, accepted.mar2=accepted.mar2, accepted=accepted,
              straight.reject1=straight.reject1, straight.reject2=straight.reject2,
              acc.vec=acc.vec, acc.vec.mar1=acc.vec.mar1, acc.vec.mar2=acc.vec.mar2, nsim=nsim, 
              sig1.vec=sig1.vec, sig1.restart=sig1.restart,
              sig2.vec=sig2.vec, sig2.restart=sig2.restart, threshold=U ))
}

# Marginal transformation to Exponential and its derivative (der=TRUE)
# Corresponds to the u(x) function

marg <- function(x, kn, par, der=FALSE){
  # par = c(mu, sigma, gamma)
  if(length(x) != 1){stop(" 'x' should be of length 1")}
  if(length(par) != 3){stop("The vector of marginal parameters 'par' should be of length 3")}
  
  q <- 1 + par[3] *  (x - par[1]) / par[2]
  if(!der){
    return( kn * q^(-1/par[3]) )
  }else{
    return( - kn * q^(-1/par[3]-1) / par[2] )      
  }  
}


# Considering exponential margins, exp(-u(x))
cens.llik <- function(coef, data = NULL, bpb = NULL, bpb1 = NULL, 
                      bpb2 = NULL, thresh, par1, par2, kn){
  
  # coef:             A vector of the eta coefficients
  # data:             A list containing:
  #                      xy: the original data, 
  #                      data.cens: original censored data
  #                      data.cens.uv: the censored data, marginally transformed to unit Frechet, 
  #                      w.cens.uv: angular data of data.cens.uv
  #                      r.cens.uv: radius of data.cens.uv 
  # bpb, bpb1, bpb2:  Matrices of the Bernstein polynomial basis
  # thresh:           Some high threhsolds for each variable
  # par1, par2:       Vectors of length 3 representing the marginal parameters as (mu, sigma, gamma)
  
  if(ncol(par1) != 3 || ncol(par2) != 3){stop("Wrong length of marginal parameters")}
  if( length(thresh) != 2){stop("Wrong length of threshold parameters")}
  
  if( any(par1[,1]<=0) || any(par2[1]<=0)){return(-1e300)}
  
  if( any(par1[,3] > 0) && any( thresh[1] < (par1[,1] - par1[,2]/par1[,3]) ) ){return(-1e300)}
  if( any(par1[,3] < 0) && any( thresh[1] > (par1[,1] - par1[,2]/par1[,3]) ) ){return(-1e300)}
  if( any(par2[,3] > 0) && any( thresh[2] < (par2[,1] - par2[,2]/par2[,3]) ) ){return(-1e300)}
  if( any(par2[,3] < 0) && any( thresh[2] > (par2[,1] - par2[,2]/par2[,3]) ) ){return(-1e300)}
  
  uv.data <- data$data.cens.uv
  
  # derive the betas coefficients
  beta <- net(coef, from='H')$beta
  k <- length(beta)-2
  
  # compute the difference of the coefficients:
  beta1 <- diff(beta); beta2 <- diff(beta1)
  
  indiv.cens.llik <- function(data, uv.data, bpb = NULL, bpb1=NULL, bpb2=NULL, par1, par2, thresh, kn){
    
    # data corresponds to the original censored data
    # the (marginal) observation above the thresholds are unit Frechet distributed
    # data is of length 2: data = (x, y)
    
    if( all(data<= thresh) ){ # x <= tx and y <= ty
      A.txty <- c(bpb %*% beta) # Pickands function A(t) at t = v(ty) / (u(tx) + v(ty))
      res <- - sum(uv.data) * A.txty # Returns -(u(tx) + v(ty)) * A(t)
    }
    
    
    if( data[1] > thresh[1] && data[2] <= thresh[2]){ # x > tx and y <= ty
      # t = v(ty) / (u(x) + v(ty))
      A.xty <- c(bpb %*% beta) # Pickands function at A(t)
      A1.xty <- c((k+1) * (bpb1 %*% beta1)) # 1st derivative Pickands function: A'(t)
      part1 <- -( A.xty - A1.xty * uv.data[2] / sum(uv.data)  ) * marg(data[1], par = par1, kn = kn[1], der = TRUE) # 1st part: - du(x)/dx * [ A(t) - t * A'(t) ]
      part2 <- - sum(uv.data) * A.xty # 2nd part: - (u(x) + v(ty)) * A(t)
      res <- log( part1 ) + part2
    }
    
    if( data[1] <= thresh[1] && data[2] > thresh[2]){ # x <= tx and y > ty
      # t = v(y) / (u(tx) + v(y))
      A.txy <- c(bpb %*% beta) # Pickands function at: A(t)
      A1.txy <- c((k+1) * (bpb1 %*% beta1)) # 1st derivative Pickands function: A'(t)
      part1 <- -( A.txy + A1.txy * uv.data[1] / sum(uv.data)  ) * marg(data[2], par = par2, kn = kn[2], der = TRUE) # 1st part: - dv(y)/dy * [ A(t) + (1 - t) * A'(t) ]
      part2 <- - sum(uv.data) * A.txy # 2nd part: - (u(tx) + v(y)) * A(t)
      res <- log( part1 ) + part2
    }
    
    if( all(data> thresh) ){ # x > tx and y > ty
      # t = v(y) / (u(x) + v(y))
      A.xy <- c(bpb %*% beta) # Pickands function at (u(x), v(y)): A(t)
      A1.xy <- c((k+1) * (bpb1 %*% beta1)) # 1st derivative Pickands function: A'(t)
      A2.xy <- c(k*(k+1) * (bpb2 %*% beta2)) # 2nd derivative Pickands function: A''(t)
      part1 <- marg(data[1], par = par1, kn = kn[1], der = TRUE) * marg(data[2], par = par2, kn = kn[2], der = TRUE) # 1st part: du(x)/dx * dv(y)/dy
      part2 <- ( A.xy - A1.xy * uv.data[2] / sum(uv.data)  ) # 2nd part: A( t ) - t * A'(t)
      part3 <- ( A.xy + A1.xy * uv.data[1] / sum(uv.data)  ) # 3rd part: A( t ) + (1-t) * A'(t)
      part4 <- - prod(uv.data) / sum(uv.data)^3 * A2.xy # 4th part: - t * (1-t) / ( u(x) + v(y) ) * A''(t)
      part5 <- - sum(uv.data) * A.xy # 5th part: - (u(x) + v(y)) * A(t)
      res <- log(part1) + log(part2 * part3 - part4) + part5
    }
    
    if(is.na(res)) return(-1e+300) else return(res)    
  }    
  
  Sum <- 0
  if(nrow(par1) == 1 && nrow(par2) == 1){ # There are no covariates, the parameters are the same for each observations
    if(ncol(data$data.cens)!=3){stop("data$data.cens should have 3 columns!")}
    for(i in 1:nrow(data$xy)){
      Sum <- Sum + data$data.cens[i,3] *indiv.cens.llik( data = data$xy[i,], uv.data = uv.data[i,], bpb = bpb[i,], bpb1 = bpb1[i,], bpb2 = bpb2[i,], par1 = par1, par2 = par2, thresh  = thresh, kn = kn )
    }
  }else{
    for(i in 1:nrow(data$xy)){ # There are covariates, the parameters change with the observations
      Sum <- Sum + indiv.cens.llik( data = data$xy[i,], uv.data = uv.data[i,], bpb = bpb[i,], bpb1 = bpb1[i,], bpb2 = bpb2[i,], par1 = par1[i,], par2 = par2[i,], thresh  = thresh, kn = kn )
    } 
  }
  
  return( Sum)  
  
}

summary.bbeed <- function(object, mcmc, burn, conf=0.95, plot=FALSE, ...) {
  
  w <- object
  nsim <- mcmc$nsim
  ww <- cbind(w, 1-w)
  
  # posterior k
  qk <- quantile(mcmc$k[burn:nsim], c(1-conf, 0.5, conf), type=3)
  k.median <- qk[2]
  k.up <- qk[3]
  k.low <- qk[1]
  
  bpg <- NULL
  for(i in 1:(max(mcmc$k, na.rm=TRUE)+1)) bpg[[i]] <- bp(ww, i)
  
  # h(w) = k sum_j^k-1 (eta_{j+1} - eta_j) bj(w,k-1)
  # A(w) = sum_j^k+1 beta_j bj(w,k+1)
  
  ngrid <- nrow(bpg[[1]]) 
  iters <- nsim-burn+1 
  eta.diff_post <- beta_post <- NULL
  h_post <- A_post <- matrix(NA,iters, ngrid)
  p0_post <- p1_post <- numeric(iters)
  if(length(mcmc$mar1) !=0 && length(mcmc$mar2) != 0){
    mar1_post <- mar2_post <- matrix(NA, iters, ncol(mcmc$mar1)) # Matrix for the marginal parameters
  }else{
    mar1_post <- mar2_post <- 0
  }
  
  for(i in burn:nsim){
    l <- i-burn+1
    ki <- mcmc$k[i]
    etai <- mcmc$eta[i,1:(ki+1)]
    eta.diff_post[[l]] <- diff(etai)
    h_post[l,] <- ki * c(bpg[[ki-1]] %*% eta.diff_post[[l]])
    p0_post[l] <- etai[1]
    p1_post[l] <- 1-etai[ki+1]
    beta_post[[l]] <- net(etai, from='H')$beta
    A_post[l,] <- c(bpg[[ki+1]] %*% beta_post[[l]])
    if(length(mcmc$mar1) != 0 && length(mcmc$mar2) != 0){ # If the marginal parameters have been fitted
      mar1_post[l,] <- mcmc$mar1[i,]
      mar2_post[l,] <- mcmc$mar2[i,]
    }
  }
  
  # Pointwise credibility bands and mean cuve of h(w)
  h.mean <- apply(h_post, 2, mean)#colMeans(h_post)
  h.up <- apply(h_post, 2, quantile, 1-conf, type=3)
  h.low <- apply(h_post, 2, quantile, conf, type=3)
  
  A.mean <- apply(A_post, 2, mean)#colMeans(A_post)
  A.up <- apply(A_post, 2, quantile, 1-conf, type=3)
  A.low <- apply(A_post, 2, quantile, conf, type=3)
  
  p0.mean <- mean(p0_post)
  p0.up <- quantile(p0_post, 1-conf, type=3)
  p0.low <- quantile(p0_post, conf, type=3)
  
  p1.mean <- mean(p1_post)
  p1.up <- quantile(p1_post, 1-conf, type=3)
  p1.low <- quantile(p1_post, conf, type=3)
  
  if(length(mcmc$mar1) != 0 && length(mcmc$mar2) != 0){
    mar1.mean <- apply(mar1_post, 2, mean) 
    mar1.up <- apply(mar1_post, 2, quantile, 1-conf, type=3) 
    mar1.low <- apply(mar1_post, 2, quantile, conf, type=3) 
    
    mar2.mean <- apply(mar2_post, 2, mean) 
    mar2.up <- apply(mar2_post, 2, quantile, 1-conf, type=3) 
    mar2.low <- apply(mar2_post, 2, quantile, conf, type=3)   
  }else{
    mar1.mean <- mar1.up <- mar1.low <- 0
    mar2.mean <- mar2.up <- mar2.low <- 0    
  }
  
  out <- list(k.median=k.median, k.up=k.up, k.low=k.low,
              h.mean=h.mean, h.up=h.up, h.low=h.low,
              A.mean=A.mean, A.up=A.up, A.low=A.low,
              p0.mean=p0.mean, p0.up=p0.up, p0.low=p0.low,
              p1.mean=p1.mean, p1.up=p1.up, p1.low=p1.low,
              mar1.mean=mar1.mean, mar1.up=mar1.up, mar1.low=mar1.low,
              mar2.mean=mar2.mean, mar2.up=mar2.up, mar2.low=mar2.low,
              A_post=A_post, h_post=h_post, 
              eta.diff_post=eta.diff_post, 
              beta_post=beta_post,
              mar1_post=mar1_post, mar2_post=mar2_post,
              p0_post=p0_post, p1_post=p1_post,
              w=w, burn=burn)
  
  if(plot)
    plot.bbeed(type = "summary", x=w, mcmc=mcmc, summary.mcmc=out, nsim=nsim, burn=burn, ...)
  
  return(out)
  
}

func <- function(y){
  y <- y[is.finite(y)]
  return(c(quantile(y, 0.05), mean(y), quantile(y, 0.95)))
}

###
### Hidden functions for frequentist estimation (EDhK estimator)
###

AngularMeasure <- function(data.x, data.y, data = NULL, k, method = "u", plot = TRUE) {
  
  # -------------------------------
  Weights <- function(theta) {
    k <- length(theta)
    f <- cos(theta) - sin(theta)
    MaxIter <- 20
    Tol <- 10^(-4)
    Converged <- FALSE
    iter <- 0
    lambda <- 0
    while(!Converged & iter < MaxIter) {
      t <- f / (lambda * f + k)
      dlambda <- sum(t) / sum(t^2)
      lambda <- lambda + dlambda
      iter <- iter + 1
      Converged <- abs(dlambda) < Tol
    }
    if (!Converged) warning("Algorithm to find Lagrange multiplier has not converged.", call. = FALSE)
    p <- 1 / (lambda * f + k)
    w <- p / sum(cos(theta) * p)
    return(w)
  }
  # -------------------------------
  
  if (!is.null(data)) {
    stopifnot(isTRUE(all.equal(length(dim(data)), 2)), isTRUE(all.equal(dim(data)[2], 2)))
    data.x <- data[,1]
    data.y <- data[,2]
    main <- paste("spectral measure -- data =", deparse(substitute(data)))
  }
  else {
    main <- paste("spectral measure -- data = ", deparse(substitute(data.x)), ", ", deparse(substitute(data.y)), sep = "")
  }
  n <- length(data.x)
  stopifnot(identical(length(data.x), length(data.y)))
  stopifnot(min(k) > 0, max(k) < n+1)
  r.x <- n / (n + 1 - rank(data.x))
  r.y <- n / (n + 1 - rank(data.y))
  t <- sqrt(r.x*r.x + r.y*r.y)
  if (length(method) < length(k)) method <- array(method, dim = length(k))
  if (length(k) < length(method)) k <- array(k, dim = length(method))
  
  if (plot) {
    plot(x = c(0, pi/2), y = c(0, 2), type = "n", main = main, xlab = "theta", ylab = "Phi(theta)", xaxt = "n")
    axis(side = 1, at = c((0:4)*pi/8), labels = c("0", "pi/8", "pi/4", "3*pi/8", "pi/2"))
    if (length(k) > 6) 
      lty <- rep(1, times = length(k)) 
    else {
      lty <- c("solid", "11", "1343", "22", "44", "2151")[1:length(k)]
      legend(x = "bottomright", legend = paste("k =", k), lty = lty, col = "black", lwd = 2)
    }
  }
  
  for (i in 1:length(k)) {
    j <- which(t > n/k[i])
    theta <- sort(atan(r.y[j]/r.x[j]))
    if (method[i] == "u") 
      w <- rep(1, times = length(j)) / k[i]
    else if (method[i] == "c")
      w <- Weights(theta)
    if (plot) lines(c(0, theta, pi/2), cumsum(c(0, w, 0)), type = "s", lty = lty[i], lwd = 2)
  }
  
  out <- list(angles = theta, weights = w, radii = t, indices = j)
  invisible(structure(out, class = "AngularMeasure"))
}
