fExtDep.np <- function(method, data, cov1=NULL, cov2=NULL, u, mar.fit=TRUE, mar.prelim=TRUE,
                       par10, par20, sig10, sig20, param0=NULL, k0=NULL,
                       pm0=NULL, prior.k="nbinom", prior.pm="unif", nk=70, lik=TRUE,
                       hyperparam = list(mu.nbinom = 3.2, var.nbinom = 4.48), nsim, warn=FALSE,
                       type="rawdata"){
  
  
  methods <- c("Bayesian", "Frequentist", "Empirical")
  if(!any(method ==  methods)){ stop("estimation method wrongly specified")}
  
  if(method=="Bayesian"){
    
    if(mar.fit){
      if(is.null(par10) || missing(par10)){stop("Need to specify par10 parameter")} 
      if(is.null(par20) || missing(par20)){stop("Need to specify par20 parameter")} 
      
      if(is.null(sig10) || missing(sig10)){stop("Missing initial values for sig10")}
      if(is.null(sig20) || missing(sig20)){stop("Missing initial values for sig20")}
      
    }  
    
    if(is.null(nsim)){stop("Missing number of replicates in algorithm")}
    
    if(missing(u)){ # Block maxima approach
      
      if(mar.fit){
        
        if(is.null(cov1)){cov1 <- as.matrix(rep(1,nrow(data)))}
        if(is.null(cov2)){cov2 <- as.matrix(rep(1,nrow(data)))}      
        
        out <- bbeed.mar(data=data, cov1=cov1, cov2=cov2, mar=mar.prelim, 
                         par10=par10, par20=par20, sig10=sig10, sig20=sig20, 
                         param0=param0, k0=k0, pm0=pm0, prior.k=prior.k, prior.pm=prior.pm, 
                         nk=nk, lik=lik, hyperparam=hyperparam, nsim=nsim, warn=warn)   
      }else{
        out <- bbeed(data=data, param0=param0, k0=k0, pm0=pm0, 
                     prior.k=prior.k, prior.pm=prior.pm, 
                     nk=nk, lik=lik, hyperparam=hyperparam, nsim=nsim, warn=warn)
      }
      
    }else{ # Threshold Exceedances
      
      if(!all( is.vector(u), is.numeric(u)) ){ 
        u <- apply(data,2,function(x) quantile(x, probs=0.9, type=3))
        message("u set to 90% quantile by default on both margins \n")
      }
      if(all( is.vector(u), is.numeric(u), length(u)!=2) ){
        u <- rep(u[1],2)
        message(paste("Warning, u incorrectly specified and set to u=(", u[1], ",", u[1],") \n",sep=""))
      }
      
      if(mar.fit){
        
        if(is.null(cov1)){cov1 <- as.matrix(rep(1,nrow(data)))}
        if(is.null(cov2)){cov2 <- as.matrix(rep(1,nrow(data)))}      
        
        out <- bbeed.thresh.mar(data=data, cov1=cov1, cov2=cov2, U=u, mar=mar.prelim, 
                                par10=par10, par20=par20, sig10=sig10, sig20=sig20, 
                                param0=param0, k0=k0, pm0=pm0, prior.k=prior.k, prior.pm=prior.pm, 
                                nk=nk, hyperparam=hyperparam, nsim=nsim, warn=warn)   
        
      }else{
        
        out <- bbeed.thresh(data=data, U=u, param0=param0, k0=k0, pm0=pm0, 
                            prior.k=prior.k, prior.pm=prior.pm, nk=nk, 
                            hyperparam=hyperparam, nsim=nsim, warn=warn) 
  
      }
      
    }
    
  }
  
  if(method == "Empirical"){
    out <- Fit.Empirical(data=data)
  }
  
  if(method == "Frequentist"){
    
    types <- c("maxima", "rawdata")
    if(!any(type == types)){ stop("type needs to be `maxima` or `rawdata` when method=`Frequentist`")}
    
    # Do we need to standardise to unit Frechet?
    if(type == "rawdata" && mar.fit){

      # Compute empirical distributions: 
      # y1 <- ecdf(data[,1])
      # y2 <- ecdf(data[,2])
      
      # Transform marginal distributions into unit-frechet:
      # fdata <- cbind(1/(-log(y1(data[,1]))), 1/(-log(y2(data[,2]))))
      
      # f1 <- fGEV(data=data[,1], method="Frequentist", optim.method="Nelder-Mead")
      # f2 <- fGEV(data=data[,2], method="Frequentist", optim.method="Nelder-Mead")
      # fdata <- cbind(sapply( data[,1], function(x) ExtremalDep:::trans.UFrech(c(x,f1$est)) ),
      #                sapply( data[,2], function(x) ExtremalDep:::trans.UFrech(c(x,f2$est)) )
      #               )
      
      fdata <- trans2UFrechet(data, type="Empirical")
      
    }else{
      fdata <- data
    }
    
    # Set the grid point for angles:
    nw <- 100
    w <- seq(0.00001, .99999, length=nw)
    
    if(type=="maxima"){
      
      q <- NULL
      r0 <- NULL
      extdata <- fdata
      
    }else if(type == "rawdata"){
      
      if(missing(u) || is.null(u)){
        q <- 0.9
        warning("u not provided, 90% quantile is considered.")
      }else if(is.vector(u) && length(u)==1){
        q <- sum(data[,1]<=u[1])/nrow(data)
      }else{
        stop("Error in providing u")
      }
      
      rdata <- rowSums(fdata) # Radial components
      wdata <- fdata[,1]/rdata # Angular components
      r0 <- quantile(rdata, probs=q) # Set high threshold
      extdata <- fdata[rdata>=r0,] # Extract data from the tail
     
    }
    
    # Compute starting value:
    Aest_proj <- beed(extdata, cbind(w, 1-w), 2, "md", "emp", k=k0) 
    
    # Optimize the angular density:
    Aest_mle <- constmle(fdata, k0, w, start=round(Aest_proj$beta,10), type=type, r0=r0)
    
    # Compute the angular density:
    hhat <- sapply(1:nw, function(i, w, beta) dh(w[i], beta), w, Aest_mle$beta)
    
    # Compute the angular measure
    Hhat <- sapply(1:nw, function(i, w, beta) ph(w[i], beta), w, Aest_mle$beta)
    
    # Compute the point masses on the corners
    p0 <- (1+k0*(Aest_mle$beta[2]-1))/2
    p1 <- (1-k0*(1-Aest_mle$beta[k0]))/2
    
    # if(type == "rawdata" && mar.fit){
    #   out <- list(method=method, type=type, hhat=hhat, Hhat=Hhat, p0=p0, p1=p1, 
    #               Ahat=Aest_mle, w=w, f1=f1, f2=f2, q=q, extdata=extdata )
    # }else{
    #   out <- list(method=method, type=type, hhat=hhat, Hhat=Hhat, p0=p0, p1=p1, 
    #               Ahat=Aest_mle, w=w, q=q, extdata=extdata )      
    # }

    out <- list(method=method, type=type, hhat=hhat, Hhat=Hhat, p0=p0, p1=p1, 
                Ahat=Aest_mle, w=w, q=q, extdata=extdata ) 
    
  }
  
  return(out)
}


###############################################################################
## Function to return a dataset on unit Frechet scale 
## either by fitting the GEV, transforming from GEV parameters or empirically
###############################################################################

trans2UFrechet <- function(data, pars, type="Empirical"){
  
  types <- c("Empirical", "GEV")
  if(!(type %in% types)){stop(" 'type' should be 'Empirical' or 'GEV'.")}
  
  orig.mat <- TRUE
  if(is.vector(data)){
    orig.mat <- FALSE
    data <- matrix(data, ncol=1)
    d <- 1
  }else if(is.matrix(data)){
    d <- ncol(data)  
  }else{
    stop("'data' must be a vector or a matrix.")
  }
  
  data_uf <- matrix(ncol=d, nrow=nrow(data))
  
  if(!missing(pars)){ type <- "GEV"}
  
  if(type == "Empirical"){
    
    for(i in 1:d){
      tmp_ec <- ecdf2(data[,i])
      data_uf[,i] <- 1/(-log(tmp_ec(data[,i])))
    }
    
  }
  
  if(type == "GEV"){
    
    if(!missing(pars)){
        if(is.vector(pars) && length(pars)==3){
            for(i in 1:d){
                data_uf[,i] <- sapply(data[,i], function(x) trans.UFrech(c(x,pars)) )
            }
        }else if(is.matrix(pars) && ncol(pars)==3 && nrow(pars)==d){
            for(i in 1:d){
                data_uf[,i] <- sapply(data[,i], function(x) trans.UFrech(c(x,pars[i,])) )
            }
        }
    }else{
        for(i in 1:d){
            pars <- fGEV(data[,i], method="Frequentist",optim.method="Nelder-Mead")$est
            data_uf[,i] <- sapply(data[,i], function(x) trans.UFrech(c(x,pars)) )
        }
    }
    
    
  }
  
  if(!orig.mat){
    return(as.vector(data_uf))
  }else{
    return(data_uf)    
  }

}

###############################################################################
## Function to return a dataset on GEV scale 
###############################################################################

trans2GEV <- function(data, pars){
  
  orig.mat <- TRUE
  if(is.vector(data)){
    orig.mat <- FALSE
    data <- matrix(data, ncol=1)
    d <- 1
  }else if(is.matrix(data)){
    d <- ncol(data)  
  }else{
    stop("'data' must be a vector or a matrix.")
  }
  
  data_gev <- matrix(ncol=d, nrow=nrow(data))
  
  if(missing(pars)){stop(" 'pars' must be specified.")}
  
  trans.GEV <- function(x){
    if(x[4]!=0){
      return( (x[1]^x[4]-1)*x[3]/x[4]+x[2])      
    }else if(x[4]==0){
      return( x[3]/x[1]+x[2] )
    }

  }
  
  if(is.vector(pars) && length(pars)==3){
    for(i in 1:d){
      data_gev[,i] <- sapply(data[,i], function(x) trans.GEV(c(x,pars)) )
    }
  }else if(is.matrix(pars) && ncol(pars)==3 && nrow(pars)==d){
    for(i in 1:d){
      data_gev[,i] <- sapply(data[,i], function(x) trans.GEV(c(x,pars[i,])) )
    }
  }
  
  if(!orig.mat){
    return(as.vector(data_gev))
  }else{
    return(data_gev)    
  }
}

###############################################################################
###############################################################################
## Hidden functions
###############################################################################
###############################################################################


###############################################################################
## bbeed function (Bayesian Bernstein Estimation of Extremal Dependence)     ##
##                                                                           ##
## INPUT:                                                                    ##
## data is a (n x 2)-matrix/dataframe                                        ##
## pm0, param0, k0 are the initial values of the parameters                  ##
## hyperparam is a list of the hyperparameters                               ##
## nsim is the number of the iterations of the chain                         ##
## prior.k and prior.pm specify the prior distributions                      ##
############################################################################### 

bbeed <- function(data, pm0, param0, k0, hyperparam, 
                  nsim, nk = 70, prior.k = c('pois','nbinom'), 
                  prior.pm = c('unif','beta', 'null'), lik = TRUE, warn=FALSE){
  
  # set info data:
  #      z = 1/x + 1/y (Frechet scale)                                         
  #      w = angular of data                                                   
  #      r = radius of the data                                                
  #      w2 = w^2                                                              
  #      r2 = r^2                                                              
  #      den = x^2, y^2                                                        
  z <-rowSums(1/data)
  r <- rowSums(data)
  w <- data/r
  den <- data^2
  data <- list(z=z, r=r, w=w, r2=r^2, w2=w^2, den=den)
  
  # Checks:
  Check <- check.bayes(prior.k = prior.k, prior.pm = prior.pm, hyperparam = hyperparam, pm0 = pm0, param0 = param0, k0 = k0, warn=warn)
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
  
  #param0 = eta.start
  n <- nrow(data$w)
  nkm <- nk-1
  nkp <- nk+2
  # set the chains
  spm <- array(0, dim=c(nsim,2))
  seta <- array(0, dim=c(nsim,nkp))
  sk <- rep(0,nsim) 
  
  accepted <- acc.vec <- rep(0, nsim)
  param <- param0
  pm <- pm0
  
  k <- k0
  if(k==3) q <- 0.5 else q <- 1
  
  # compute the polynomial bases:
  bpb_mat <- NULL 
  for(j in 2:nk) bpb_mat[[j]] <- bp(data$w, j)
  
  pb = txtProgressBar(min = 0, max = nsim, initial = 0, style=3) #creates a progress bar
  print.i <- seq(0, nsim, by=100)
  
  for(i in 1:nsim){
    # simulation from the proposal:
    # polynomial order
    k_new <- prior_k_sampler(k)
    
    if(k_new==nk){
      message('maximum value k reached')
      break
    }
    
    # point masses
    pm_new <- prior_p_sampler(a = a, b = b, prior.pm = prior.pm)
    while(check.p(pm_new,k_new+1)) pm_new <- prior_p_sampler(a = a, b = b, prior.pm = prior.pm)
    
    
    if(prior.k=='pois'){
      p <- dpois(k_new-3,mu.pois)/dpois(k-3,mu.pois)      
    }else if(prior.k=='nbinom'){
      p <- dnbinom(k_new-3,prob=pnb,size=rnb)/dnbinom(k-3,prob=pnb,size=rnb)
    }
    
    # polynomial coefficients
    param_new <- rcoef(k_new, pm=pm_new)
    
    # derive the polynomial bases:
    # compute the acceptance probability of   
    if(lik==TRUE){
      ratio <- exp(llik(data=data, pm=pm_new, coef=param_new$eta, k=k_new, bpb=bpb_mat[[k_new+1]], bpb1=bpb_mat[[k_new]], bpb2=bpb_mat[[k_new-1]]) + log(q) -
                     llik(data=data, pm=pm, coef=param$eta, k=k, bpb=bpb_mat[[k+1]], bpb1=bpb_mat[[k]], bpb2=bpb_mat[[k-1]]) + log(p) )
    }else{
      ratio <- exp(llik(coef=param_new$eta, k=k_new, bpb1=bpb_mat[[k_new-1]], approx=TRUE) + log(q) -
                     llik(coef=param$eta, k=k, bpb1=bpb_mat[[k-1]], approx=TRUE) + log(p) )
    } 
    
    acc.vec[i] <- ratio
    
    u <- runif(1)
    if(u<ratio){
      pm <- pm_new
      param <- param_new
      k <- k_new
      if(k==3) q <- 0.5 else q <- 1
      accepted[i] <- 1
    }
    spm[i,] <- c(pm$p0, pm$p1)
    seta[i, 1:(k+1)] <- param$eta
    sk[i] <- k
    
    if(i %in% print.i) setTxtProgressBar(pb, i) 
  }
  
  close(pb) 
  
  return(list(method="Bayesian", mar.fit=FALSE, type="maxima", data=data,
              pm=spm, eta=seta, k=sk, accepted=accepted, acc.vec=acc.vec,
              prior = list(hyperparam=hyperparam, k=prior.k, pm=prior.pm),
              nsim=nsim))
}

###############################################################################
## bbeed.mar function (bbeed with estimation of the margins jointly)         ##

bbeed.mar <- function(data, cov1=as.matrix(rep(1,nrow(data))), cov2=as.matrix(rep(1,nrow(data))), 
                       mar=TRUE, par10, par20, sig10, sig20, param0=NULL, k0=NULL,
                       pm0=NULL, prior.k="nbinom", prior.pm="unif", nk=70, lik=TRUE,
                       hyperparam = list(mu.nbinom = 3.2, var.nbinom = 4.48), nsim=NULL, warn=FALSE){
  
  if(nrow(data) != nrow(cov1) || nrow(data) != nrow(cov2)){stop("Different number of observations/covariates")}
  
  if(length(par10) != (ncol(cov1)+2)){stop("Wrong parameter length for margin 1")}
  if(length(par20) != (ncol(cov2)+2)){stop("Wrong parameter length for margin 2")} 
  
  if(mar){
    message("\n Preliminary on margin 1 \n")
    mar1.fit <- fGEV(data=data[,1], par.start = par10, 
                     method="Bayesian", cov=cov1, sig0 = sig10, nsim = 5e+4)
    par10 <- apply(mar1.fit$param_post[-c(1:30000),], 2, mean)
    sig10 <- tail(mar1.fit$sig.vec,1)
    
    message("\n Preliminary on margin 2 \n")
    mar2.fit <- fGEV(data=data[,2], par.start = par20, 
                     method="Bayesian", cov=cov2, sig0 = sig20, nsim = 5e+4)
    par20 <- apply(mar2.fit$param_post[-c(1:30000),], 2, mean)
    sig20 <- tail(mar2.fit$sig.vec,1)
  }

  ###############################################################
  # START Estimation of the extremal dependence (angular density)
  ###############################################################  
  
  message("\n Estimation of the extremal dependence and margins \n")
    
  Par10 <- cbind( cov1 %*% par10[1:(ncol(cov1))], rep(par10[ncol(cov1)+1], nrow(cov1) ), rep(par10[ncol(cov1)+2], nrow(cov1) )  )
  Par20 <- cbind( cov2 %*% par20[1:(ncol(cov2))], rep(par20[ncol(cov2)+1], nrow(cov2) ), rep(par20[ncol(cov2)+2], nrow(cov2) )  )
  
  data.uf <- cbind( apply(rbind(data[,1], t(Par10)),2,trans.UFrech),
                    apply(rbind(data[,2], t(Par20)),2,trans.UFrech) )
  
  # radial and angular components
  r <- r_new <- rowSums(data.uf)
  w <- w_new <- data.uf/r
  z <- rowSums(1/data.uf)
  den <- data.uf^2
  
  # Jacobian for change of variables (log scale)
  
  lder <- cbind( apply(rbind(data[,1], t(Par10)),2,function(x) log(trans.UFrech(x, der=TRUE))),
                 apply(rbind(data[,2], t(Par20)),2,function(x) log(trans.UFrech(x, der=TRUE)))
               )
  
  data0 <- list(z=z, r=r, w=w, r2=r^2, w2=w^2, den=den, data=data, data.uf=data.uf, lder=lder)
  
  # derive the polynomial bases:
  # bpb <- bp(cbind(w[,2], w[,1]), k0 + 1)
  # bpb1 <- bp(cbind(w[,2], w[,1]), k0)
  # bpb2 <- bp(cbind(w[,2], w[,1]), k0 - 1)
  bpb <- bp(w, k0 + 1)
  bpb1 <- bp(w, k0)
  bpb2 <- bp(w, k0 - 1)
  
  # Checks:
  Check <- check.bayes(prior.k = prior.k, prior.pm = prior.pm, hyperparam = hyperparam, pm0 = pm0, param0 = param0, k0 = k0, warn=warn)
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
  
  p.star <- 0.234
  n0 <- round(5/(p.star*(1-p.star)))
  iMax <- 100 # iMax is the max number of iterations before the last restart
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
    
    par1_new <- as.vector(rmvnorm(1, mean = par1, sigma = sig10^2 * sig1Mat))
    Par1_new <- cbind( cov1 %*% par1_new[1:(ncol(cov1))], rep(par1_new[ncol(cov1)+1], nrow(cov1) ), rep(par1_new[ncol(cov1)+2], nrow(cov1) )  )
    
    if( (any(data[,1] < (Par1_new[,1]-Par1_new[,2]/Par1_new[,3]) ) && any(Par1_new[,3]>0))
        || (any(data[,1] > (Par1_new[,1]-Par1_new[,2]/Par1_new[,3]) ) && any(Par1_new[,3]<0))
        || any(Par1_new[,2]<0) ){
      straight.reject1[i] <- 1
      acc.vec.mar1[i] <- 0
    }else{
      # Marginalised data
      data.uf_new <- cbind( apply(rbind(data0$data[,1], t(Par1_new)),2,trans.UFrech),
                            data0$data.uf[,2] )
      
      r_new <- rowSums(data.uf_new)
      w_new <- data.uf_new/r_new
      
      lder_new <- cbind( apply(rbind(data0$data[,1], t(Par1_new)),2,function(x) log(trans.UFrech(x, der=TRUE))),
                     data0$lder[,2] )
      
      data_new <- list(z=rowSums(1/data.uf_new), r=r_new, w=w_new, r2=r_new^2, 
                       w2=w_new^2, den=data.uf_new^2, 
                       data=data, data.uf=data.uf_new, lder=lder_new)
      
      # derive the polynomial bases:
      # bpb_new <- bp(cbind(w_new[,2], w_new[,1]), k + 1)
      # bpb1_new <- bp(cbind(w_new[,2], w_new[,1]), k)
      # bpb2_new <- bp(cbind(w_new[,2], w_new[,1]), k - 1)
      bpb_new <- bp(w_new, k + 1)
      bpb1_new <- bp(w_new, k)
      bpb2_new <- bp(w_new, k - 1)
      
      ratio.mar1 <- min( exp(llik(data=data_new, pm=pm, coef=param$eta, k=k, bpb=bpb_new, bpb1=bpb1_new, bpb2=bpb2_new) + sum(data_new$lder)
                             - llik(data=data0, pm=pm, coef=param$eta, k=k, bpb=bpb, bpb1=bpb1, bpb2=bpb2) - sum(data0$lder)
                             + log(Par1[1,2]) - log(Par1_new[1,2])), 1)
      
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
    
    par2_new <- as.vector(rmvnorm(1, mean = par2, sigma = sig20^2 * sig2Mat))
    Par2_new <- cbind( cov2 %*% par2_new[1:(ncol(cov2))], rep(par2_new[ncol(cov2)+1], nrow(cov2) ), rep(par2_new[ncol(cov2)+2], nrow(cov2) )  )
    
    if( (any(data[,2] < (Par2_new[,1]-Par2_new[,2]/Par2_new[,3]) ) && any(Par2_new[,3]>0))
        || (any(data[,2] > (Par2_new[,1]-Par2_new[,2]/Par2_new[,3]) ) && any(Par2_new[,3]<0))
        || any(Par2_new[,2]<0) ){
      straight.reject2[i] <- 1
      acc.vec.mar2[i] <- 0
    }else{
      # Marginalised data
      data.uf_new <- cbind( data0$data.uf[,1],
                            apply(rbind(data0$data[,2], t(Par2_new)),2,trans.UFrech) )

      r_new <- rowSums(data.uf_new)
      w_new <- data.uf_new/r_new
      
      lder_new <- cbind( data0$lder[,1],
                         apply(rbind(data0$data[,2], t(Par2_new)),2,function(x) log(trans.UFrech(x, der=TRUE))) )
      
      data_new <- list(z=rowSums(1/data.uf_new), r=r_new, w=w_new, r2=r_new^2, 
                       w2=w_new^2, den=data.uf_new^2, 
                       data=data, data.uf=data.uf_new, lder=lder_new)
      
      # derive the polynomial bases:
      # bpb_new <- bp(cbind(w_new[,2], w_new[,1]), k + 1)
      # bpb1_new <- bp(cbind(w_new[,2], w_new[,1]), k)
      # bpb2_new <- bp(cbind(w_new[,2], w_new[,1]), k - 1)
      bpb_new <- bp(w_new, k + 1)
      bpb1_new <- bp(w_new, k)
      bpb2_new <- bp(w_new, k - 1)
      
      ratio.mar2 <- min( exp(llik(data=data_new, pm=pm, coef=param$eta, k=k, bpb=bpb_new, bpb1=bpb1_new, bpb2=bpb2_new) + sum(data_new$lder)
                             - llik(data=data0, pm=pm, coef=param$eta, k=k, bpb=bpb, bpb1=bpb1, bpb2=bpb2) - sum(data0$lder) 
                             + log(Par2[1,2]) - log(Par2_new[1,2])), 1)
      
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
    
    if(k_new==nk){
      message('maximum value k reached')
      break
    }
    
    # derive the polynomial bases:
    # bpb_new <- bp(cbind(data0$w[,2], data0$w[,1]), k_new + 1)
    # bpb1_new <- bp(cbind(data0$w[,2], data0$w[,1]), k_new)
    # bpb2_new <- bp(cbind(data0$w[,2], data0$w[,1]), k_new - 1)
    bpb_new <- bp(data0$w, k_new + 1)
    bpb1_new <- bp(data0$w, k_new)
    bpb2_new <- bp(data0$w, k_new - 1)
    
    if(prior.k=="pois"){p <- dpois(k_new-3,mu.pois)/dpois(k-3,mu.pois)}
    if(prior.k=="nbinom"){p <- dnbinom(k_new-3,prob=pnb,size=rnb)/dnbinom(k-3,prob=pnb,size=rnb)}
    
    # point masses
    pm_new <- prior_p_sampler(a=a, b=b, prior.pm=prior.pm)
    while(check.p(pm_new,k_new)) pm_new <- prior_p_sampler(a=a, b=b, prior.pm=prior.pm)
    
    # polynomial coefficients
    param_new <- rcoef(k=k_new, pm=pm_new)
    
    # compute the acceptance probability of   
    if(lik==TRUE){
      ratio <- exp(llik(data=data0, pm=pm_new, coef=param_new$eta, k=k_new, bpb=bpb_new, bpb1=bpb1_new, bpb2=bpb2_new) + log(q) -
                   llik(data=data0, pm=pm, coef=param$eta, k=k, bpb=bpb, bpb1=bpb1, bpb2=bpb2) + log(p) )
    }else{
      ratio <- exp(llik(coef=param_new$eta, k=k_new, bpb1=bpb1_new, approx=TRUE) + log(q) -
                     llik(coef=param$eta, k=k, bpb1=bpb1, approx=TRUE) + log(p) )
    } 
    
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
          # restart the algorithm
          # message("\n restart the program at ", i, "th iteration", "\n")
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
          # restart the algorithm
          # message("\n restart the program at", i, "th iteration", "\n")
          sig2.restart <- c(sig2.restart, sig2)
          Numbig2 <- Numbig2 + Toobig2
          Numsmall2 <- Numsmall2 + Toosmall2
          i <- n0
          sig2.start <- sig2
        }
        
      } #end iMax mar 2
      
    }
    
  } # End loop i in 1:nsim  

  close(pb)
  
  return(list(method="Bayesian", mar.fit=TRUE, type="maxima", 
              data=data, cov1=cov1, cov2=cov2, 
              pm=spm, eta=seta, k=sk, mar1=smar1, mar2=smar2, 
              accepted.mar1=accepted.mar1, accepted.mar2=accepted.mar2, accepted=accepted,
              straight.reject1=straight.reject1, straight.reject2=straight.reject2,
              acc.vec=acc.vec, acc.vec.mar1=acc.vec.mar1, acc.vec.mar2=acc.vec.mar2, nsim=nsim, 
              sig1.vec=sig1.vec, sig2.vec=sig2.vec,
              prior = list(hyperparam=hyperparam, k=prior.k, pm=prior.pm) ))
 
}  
  

# Function to transform to unit Frechet margins

trans.UFrech <- function(x, der=FALSE){
  if(length(x) != 4){stop("trans.UFrech takes vectors of length 4")}
  
  if(x[4] !=0){
    
    q <- 1 + (x[1] - x[2]) * x[4] / x[3] 
    # max(q,0) is omitted since some checks are put so that q is always positive.
    
    if(!der){
      return( q^(1/x[4]) )
    }else{
      return( q^(1/x[4]-1) / x[3] )      
    }  
    
  }
  
  if(x[4]==0){
    if(!der){
      return(x[3] / (x[1] - x[2]))
    }
  }
  
  
}


###############################################################################
## bbeed.thresh.mar function (bbeed for peaks over threshold with estimation ##
##  of the margins jointly)                                                  ##

bbeed.thresh.mar <- function(data, cov1=as.matrix(rep(1,nrow(data))), cov2=as.matrix(rep(1,nrow(data))), 
                      U, mar=TRUE, par10, par20, sig10, sig20, param0=NULL, k0=NULL,
                      pm0=NULL, prior.k="nbinom", prior.pm="unif", nk=70,
                      hyperparam = list(mu.nbinom = 3.2, var.nbinom = 4.48), nsim=NULL, warn=FALSE){ # need d=5?

  if(nrow(cov1) != nrow(data) || nrow(cov2) != nrow(data)){stop("Wrong dimensions of the covariates")}

  if(length(par10) != (ncol(cov1)+2)){stop("Wrong length of initial first marginal parameter vector")}
  if(length(par20) != (ncol(cov2)+2)){stop("Wrong length of initial second marginal parameter vector")}
  
  kn <- c(sum(data[,1]>U[1]), sum(data[,2]>U[2])) / nrow(data)
  
  if(mar){
    message("\n Preliminary on margin 1 \n")
    mar1.fit <- fGEV(data=data[,1], par.start = par10, u=U[1],
                     method="Bayesian", cov=cov1, sig0 = sig10, nsim = 5e+4)
    par10 <- apply(mar1.fit$param_post[-c(1:30000),], 2, mean)
    sig10 <- tail(mar1.fit$sig.vec,1)
    
    message("\n Preliminary on margin 2 \n")
    mar2.fit <- fGEV(data=data[,2], par.start = par20, u=U[2],
                     method="Bayesian", cov=cov2, sig0 = sig20, nsim = 5e+4)
    par20 <- apply(mar2.fit$param_post[-c(1:30000),], 2, mean)
    sig20 <- tail(mar2.fit$sig.vec,1)
  }
  
  ###############################################################
  # START Estimation of the extremal dependence (angular density)
  ###############################################################
  
  message("\n Estimation of the extremal dependence and margins \n")
  mcmc <- cens.bbeed.mar(data=data, cov1=cov1, cov2=cov2, pm0=pm0, param0=param0, 
                         k0=k0, par10=par10, par20=par20, sig10=sig10, sig20=sig20, 
                         kn=kn, prior.k = prior.k, prior.pm=prior.pm, 
                         hyperparam=hyperparam, nsim=nsim, U=U, p.star=0.234, warn=warn)
  
  return(mcmc)
  
}  
  
cens.bbeed.mar <- function(data, cov1, cov2, pm0, param0, k0, par10, par20, sig10, 
                           sig20, kn, prior.k = c('pois','nbinom'), 
                           prior.pm = c('unif', 'beta'), hyperparam, nk=70, 
                           nsim, U=NULL, p.star, warn=FALSE){
  
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
  r.cens.uv <- r.cens.uv_new <- rowSums(data.cens.uv)
  w.cens.uv <- w.cens.uv_new <- data.cens.uv/r.cens.uv
  data0 <- list(r.cens.uv=r.cens.uv, w.cens.uv=w.cens.uv, xy = as.matrix(xy), data.cens = as.matrix(data.cens), data.cens.uv = as.matrix(data.cens.uv) )
  
  # derive the polynomial bases:
  bpb <- bp(cbind(w.cens.uv[,2], w.cens.uv[,1]), k0 + 1)
  bpb1 <- bp(cbind(w.cens.uv[,2], w.cens.uv[,1]), k0)
  bpb2 <- bp(cbind(w.cens.uv[,2], w.cens.uv[,1]), k0 - 1)
  
  # Checks:
  Check <- check.bayes(prior.k = prior.k, prior.pm = prior.pm, hyperparam = hyperparam, pm0 = pm0, param0 = param0, k0 = k0, warn=warn)
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
    
    if( any(U[1] < (Par1_new[,1]-Par1_new[,2]/Par1_new[,3]) ) || any(Par1_new[,2]<0) || any(Par1_new[,3]<0) ){
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
      
      ratio.mar1 <- min( exp(cens.llik.marg(data = data_new, coef = param$eta, bpb = bpb_new, bpb1 = bpb1_new, bpb2 = bpb2_new, thresh = U, par1 = Par1_new, par2 = Par2, kn = kn)
                             - cens.llik.marg(data = data0, coef = param$eta, bpb = bpb, bpb1 = bpb1, bpb2 = bpb2, thresh = U, par1 = Par1, par2 = Par2, kn = kn) 
                             + log(Par1[1,2]) - log(Par1_new[1,2])), 1)
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
    
    if( any(U[2] < (Par2_new[,1]-Par2_new[,2]/Par2_new[,3]) ) || any(Par2_new[,2]<0) || any(Par2_new[,3]<0) ){
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
      ratio.mar2 <- min( exp(cens.llik.marg(data = data_new, coef = param$eta, bpb = bpb_new, bpb1 = bpb1_new, bpb2 = bpb2_new, thresh = U, par1 = Par1, par2 = Par2_new, kn= kn) 
                             - cens.llik.marg(data = data0, coef = param$eta, bpb = bpb, bpb1 = bpb1, bpb2 = bpb2, thresh = U, par1 = Par1, par2 = Par2, kn = kn) 
                             + log(Par2[1,2]) - log(Par2_new[1,2])), 1)
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
    ratio <- min( exp(cens.llik.marg(data = data0, coef = param_new$eta, bpb = bpb_new, bpb1 = bpb1_new, bpb2 = bpb2_new, thresh = U, par1 = Par1, par2 = Par2, kn = kn) + log(q) -
                        cens.llik.marg(data = data0, coef = param$eta, bpb = bpb, bpb1 = bpb1, bpb2 = bpb2, thresh = U, par1 = Par1, par2 = Par2, kn = kn) + log(p)), 1)
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
          # restart the algorithm
          # message("\n restart the program at", i, "th iteration", "\n")
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
          # restart the algorithm
          # message("\n restart the program at", i, "th iteration", "\n")
          sig2.restart <- c(sig2.restart, sig2)
          Numbig2 <- Numbig2 + Toobig2
          Numsmall2 <- Numsmall2 + Toosmall2
          i <- n0
          sig2.start <- sig2
        }
        
      } #end iMax mar 2
      
    }
    
  }
  
  close(pb)
  
  return(list(method="Bayesian", mar.fit=TRUE, type="rawdata", 
              data=data, cov1=cov1, cov2=cov2, 
              pm=spm, eta=seta, k=sk, mar1=smar1, mar2=smar2, 
              accepted.mar1=accepted.mar1, accepted.mar2=accepted.mar2, accepted=accepted,
              straight.reject1=straight.reject1, straight.reject2=straight.reject2,
              acc.vec=acc.vec, acc.vec.mar1=acc.vec.mar1, acc.vec.mar2=acc.vec.mar2, nsim=nsim, 
              sig1.vec=sig1.vec, sig2.vec=sig2.vec, threshold=U, kn=kn, 
              prior = list(hyperparam=hyperparam, k=prior.k, pm=prior.pm) ))
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
cens.llik.marg <- function(coef, data = NULL, bpb = NULL, bpb1 = NULL, 
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

###############################################################################
## bbeed.thresh function (bbeed for peaks over threshold with estimation)    ##

bbeed.thresh <- function(data, U, param0=NULL, k0=NULL, pm0=NULL, 
                         prior.k="nbinom", prior.pm="unif", nk=70,
                         hyperparam = list(mu.nbinom = 3.2, var.nbinom = 4.48), nsim=NULL, warn=FALSE){
  
  kn <- c(sum(data[,1]>U[1]), sum(data[,2]>U[2])) / nrow(data)
  
  ###############################################################
  # START Estimation of the extremal dependence (angular density)
  ###############################################################
  
  message("\n Estimation of the extremal dependence \n")
  mcmc <- cens.bbeed(data=data, pm0=pm0, param0=param0, 
                     k0=k0, kn=kn, prior.k=prior.k, prior.pm=prior.pm, 
                     hyperparam=hyperparam, nsim=nsim, U=U, warn=warn)
  
  return(mcmc)
  
}  

cens.bbeed <- function(data, pm0, param0, k0, kn, prior.k = c('pois','nbinom'), 
                       prior.pm = c('unif', 'beta'), hyperparam, nk=70, 
                       nsim, U=NULL, warn=FALSE){
  
  if(is.null(U)){
    U <- apply(data,2,function(x) quantile(x, probs=0.9, type=3))
    message("U set to 90% quantile by default on both margins \n")
  }
  
  data.cens <- data
  data.cens[data.cens[,1]<=U[1],1] <- U[1]
  data.cens[data.cens[,2]<=U[2],2] <- U[2]
  
  data.censored <- apply(data.cens,1, function(x) ((x[1]==U[1]) && (x[2]==U[2])) )
  censored <- which(data.censored==1)[1]
  n.censored <- sum(data.censored)
  uncensored <- which(data.censored==0)
  data.cens <- cbind(data.cens[c(censored,uncensored),], c(n.censored, rep(1,length(uncensored))))
  xy <- as.matrix(data[c(censored,uncensored),])
    
  r.cens <- r.cens_new <- rowSums(data.cens)
  w.cens <- w.cens_new <- data.cens[,1:2]/r.cens
  data0 <- list(r.cens=r.cens, w.cens=w.cens, xy = as.matrix(xy), data.cens = as.matrix(data.cens) )
  
  # derive the polynomial bases:
  bpb_mat <- NULL 
  for(j in 2:nk) bpb_mat[[j]] <- bp(data0$w.cens, j)
  
  # Checks:
  Check <- check.bayes(prior.k = prior.k, prior.pm = prior.pm, hyperparam = hyperparam, pm0 = pm0, param0 = param0, k0 = k0, warn=warn)
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
  
  n <- nrow(data)
  
  acc.vec <- rep(NA, nsim) # vector of acceptances
  accepted <- rep(0, nsim)
  
  pm <- pm0
  param <- param0
  spm <- c(pm$p0, pm$p1) 
  sk <- k <- k_new <- k0
  seta <- array(0, dim=c(nsim+1,nk+2))
  seta[1,1:(k+1)] <- param$eta
  
  if(k==3) q <- 0.5 else q <- 1

  pb = txtProgressBar(min = 0, max = nsim, initial = 0, style=3) #creates a progress bar
  print.i <- seq(0, nsim, by=100)
  
  for(i in 1:nsim){
    
    ###
    ### Dependence structure
    ###
    
    # polynomial order
    k_new <- prior_k_sampler(k)
    
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
    ratio <- min( exp(cens.llik(data=data0, coef=param_new$eta, bpb=bpb_mat[[k_new+1]], bpb1=bpb_mat[[k_new]], bpb2=bpb_mat[[k_new-1]], thresh=U, kn=kn) + log(q) -
                        cens.llik(data=data0, coef=param$eta, bpb=bpb_mat[[k+1]], bpb1=bpb_mat[[k]], bpb2=bpb_mat[[k-1]], thresh=U, kn=kn) + log(p)), 1)
    acc.vec[i] <- ratio
    
    u <- runif(1)
    if(u < ratio){
      pm <- pm_new
      param <- param_new
      k <- k_new
      if(k==3) q <- 0.5 else q <- 1
      accepted[i] <- 1
    }
    
    sk <- c(sk,k)
    spm <- rbind(spm, c(pm$p0, pm$p1))
    seta[i+1,1:(k+1)] <- param$eta
    
    if(i %in% print.i) setTxtProgressBar(pb, i) 
    
  }
  
  close(pb)
  
  return(list(method="Bayesian", mar.fit=FALSE, type="rawdata", data=data, 
              pm=spm, eta=seta, k=sk, accepted=accepted, acc.vec=acc.vec, 
              nsim=nsim, threshold=U, kn=kn,
              prior = list(hyperparam=hyperparam, k=prior.k, pm=prior.pm) ))
}



# Considering exponential margins, exp(-u(x))
cens.llik <- function(coef, data = NULL, bpb = NULL, bpb1 = NULL, 
                      bpb2 = NULL, thresh, par1, kn){
  
  # coef:             A vector of the eta coefficients
  # data:             A list containing:
  #                      xy: the original data, 
  #                      data.cens: original censored data
  #                      w.cens: angular data of data.cens
  #                      r.cens: radius of data.cens
  # bpb, bpb1, bpb2:  Matrices of the Bernstein polynomial basis
  # thresh:           Some high thresholds for each variable
  # par1, par2:       Vectors of length 3 representing the marginal parameters as (mu, sigma, gamma)
  
  if( length(thresh) != 2){stop("Wrong length of threshold parameters")}
  
  # derive the betas coefficients
  beta <- net(coef, from='H')$beta
  k <- length(beta)-2
  
  # compute the difference of the coefficients:
  beta1 <- diff(beta); beta2 <- diff(beta1)
  
  indiv.cens.llik <- function(data, bpb = NULL, bpb1=NULL, bpb2=NULL, thresh, kn){
    
    # data corresponds to the original censored data
    # the (marginal) observation above the thresholds are unit Frechet distributed
    # data is of length 2: data = (x, y)
    
    if( all(data<= thresh) ){ # x <= tx and y <= ty
      A.txty <- c(bpb %*% beta) # Pickands function A(t) at t = y / (x + y)
      res <- - sum(data) * A.txty # Returns -(x + y) * A(t)
    }
    
    if( data[1] > thresh[1] && data[2] <= thresh[2]){ # x > tx and y <= ty
      # t = y / (x + y)
      A.xty <- c(bpb %*% beta) # Pickands function at A(t)
      A1.xty <- c((k+1) * (bpb1 %*% beta1)) # 1st derivative Pickands function: A'(t)
      part1 <- -( A.xty - A1.xty * data[2] / sum(data)  ) # 1st part: -[ A(t) - t * A'(t) ]
      part2 <- - sum(data) * A.xty # 2nd part: - (x + y) * A(t)
      res <- log( part1 ) + part2
    }
    
    if( data[1] <= thresh[1] && data[2] > thresh[2]){ # x <= tx and y > ty
      # t = y / (x + y)
      A.txy <- c(bpb %*% beta) # Pickands function at: A(t)
      A1.txy <- c((k+1) * (bpb1 %*% beta1)) # 1st derivative Pickands function: A'(t)
      part1 <- -( A.txy + A1.txy * data[1] / sum(data)  ) # 1st part: -[ A(t) + (1 - t) * A'(t) ]
      part2 <- - sum(data) * A.txy # 2nd part: - (x + y) * A(t)
      res <- log( part1 ) + part2
    }
    
    if( all(data> thresh) ){ # x > tx and y > ty
      # t = y / (x + y)
      A.xy <- c(bpb %*% beta) # Pickands function at (x, y): A(t)
      A1.xy <- c((k+1) * (bpb1 %*% beta1)) # 1st derivative Pickands function: A'(t)
      A2.xy <- c(k*(k+1) * (bpb2 %*% beta2)) # 2nd derivative Pickands function: A''(t)
      part1 <- ( A.xy - A1.xy * data[2] / sum(data)  ) # 1st part: A( t ) - t * A'(t)
      part2 <- ( A.xy + A1.xy * data[1] / sum(data)  ) # 2nd part: A( t ) + (1-t) * A'(t)
      part3 <- - prod(data) / sum(data)^3 * A2.xy # 3rd part: - t * (1-t) / ( x + y ) * A''(t)
      part4 <- - sum(data) * A.xy # 4th part: - (x + y) * A(t)
      res <- log(part1 * part2 - part3) + part4
    }
    
    if(is.na(res)) return(-1e+300) else return(res)    
  }    
  
  Sum <- 0
  
  if(ncol(data$data.cens)!=3){stop("data$data.cens should have 3 columns!")}
  for(i in 1:nrow(data$xy)){
    Sum <- Sum + data$data.cens[i,3] * indiv.cens.llik(data=data$xy[i,], bpb=bpb[i,], bpb1=bpb1[i,], bpb2=bpb2[i,], thresh=thresh, kn=kn )
  }
  
  return(Sum)  
  
}

check.p <- function(pm,k){
  p0 <- pm$p0
  p1 <- pm$p1
  up <- 1/(k-1)*(k/2-1+p1)
  low <- (2-k)/2 + (k-1)*p1
  (p0 < low | p0 > up)
} # if FALSE, the prior is ok.

# Sampler k
prior_k_sampler <- function(k){
  if(k>3) k_new <- k+sample(c(-1,1),1,prob=c(1/2,1/2)) else k_new <- 4
  return(k_new)
}

# Sampler p
prior_p_sampler <- function(a, b, prior.pm){
  
  if(all(prior.pm != c("unif", "beta"))){stop("Wrong prior for pm")}
  
  if(prior.pm=="unif"){
    p0 <- runif(1, a, b)
    p1 <- runif(1, a, b)
  }else if(prior.pm=="beta"){
    p0 <- 1/2 * rbeta(1, a[1], b[1])
    p1 <- 1/2 * rbeta(1, a[2], b[2])
  }else if(prior.pm=="NULL"){
    p0 <- p1 <- 0
  }
  
  return(list(p0=p0, p1=p1))
}

check.bayes <- function(prior.k = c('pois','nbinom'), prior.pm = c('unif', 'beta', 'null'), hyperparam, pm0, param0, k0, warn=TRUE){
  
  if(length(prior.k) != 1 || all(prior.k != c('pois','nbinom') )){
    prior.k <- 'nbinom'
    if(warn){
      warning('prior on k not specified: prior.k = "nbinom" by default.')
    }  
  }
  
  if(length(prior.pm) != 1 || all(prior.pm != c('unif', 'beta', 'null') )){
    prior.pm <- 'unif'
    if(warn){
      warning('prior on p not specified: prior.pm = "unif" by default.')
    }  
  }
  
  if(missing(hyperparam)){
    hyperparam <- NULL
    hyperparam$a.unif <- 0
    hyperparam$b.unif <- 0.5
    hyperparam$a.beta <- c(0.8,0.8)
    hyperparam$b.beta <- c(5,5)
    mu.pois <- hyperparam$mu.pois <- 4
    mu.nbinom <- hyperparam$mu.nbinom <- 4
    var.nbinom <- hyperparam$var.nbinom <-8
    pnb <- hyperparam$pnb <- mu.nbinom/var.nbinom
    rnb <- hyperparam$rnb <- mu.nbinom^2 / (var.nbinom-mu.nbinom)
    if(warn){
      warning('hyperparam missing and set by default.')
    }  
  }
  
  if(prior.k=='pois'){ 
    if(length(hyperparam$mu.pois) != 1){
      hyperparam$mu.pois <- 4
      if(warn){
        warning('hyperparam$mu.pois missing, set to 4 by default')
      }  
    }
    mu.pois <- hyperparam$mu.pois
  }
  if(!exists("mu.pois")){mu.pois <- 4}
  
  if(prior.k=='nbinom'){
    if(length(hyperparam$mu.nbinom) != 1){
      hyperparam$mu.nbinom <- 4
      if(warn){
        warning('hyperparam$mu.nbinom missing, set to 4 by default')
      }  
    }
    if(length(hyperparam$var.nbinom) != 1){
      hyperparam$var.nbinom <- 8
      if(warn){
        warning('hyperparam$var.nbinom missing, set to 8 by default')
      }  
    }
    mu.nbinom <- hyperparam$mu.nbinom
    var.nbinom <-hyperparam$var.nbinom
    pnb <- mu.nbinom/var.nbinom
    rnb <- mu.nbinom^2 / (var.nbinom-mu.nbinom)
  }
  if(!exists("pnb")){pnb <- 1 - 4/8}
  if(!exists("rnb")){pnb <- 4^2/(8-4)}
  
  if(prior.pm=='unif'){
    
    if(length(hyperparam$a.unif) != 1){
      hyperparam$a.unif <- 0
      if(warn){
        warning('hyperparam$a.unif missing, set to 0 by default')
      }  
    }
    if(length(hyperparam$b.unif) != 1){
      hyperparam$b.unif <- 0.5
      if(warn){
        warning('hyperparam$b.unif missing, set to 0.5 by default')
      }  
    }
    a <- hyperparam$a.unif
    b <- hyperparam$b.unif
  }
  
  if(prior.pm=='beta'){
    
    if(length(hyperparam$a.beta) != 1){
      hyperparam$a.beta <- c(.8,.8)
      if(warn){
        warning('hyperparam$a.beta missing, set to (0.8,0.8) by default')
      }  
    }
    if(length(hyperparam$b.beta) != 1){
      hyperparam$b.beta <- c(5,5)
      if(warn){
        warning('hyperparam$b.beta missing, set to (5,5) by default')
      }  
    }
    a <- hyperparam$a.beta
    b <- hyperparam$b.beta
  }
  
  if(missing(k0) || is.null(k0) || length(k0)!=1){
    if(prior.k=='pois'){
      k0 <- rpois(1, mu.pois)  
      if(warn){
        warning('k0 missing or length not equal to 1, random generated from a poisson distr.')
      }  
    }
    if(prior.k=='nbinom'){
      k0 <- rnbinom(1, size=rnb, prob=pnb)  
      if(warn){
        warning('k0 missing or length not equal to 1, random generated from a negative binomial distr.')
      }  
    }
  }
  
  if(missing(pm0) || is.null(pm0) || !is.list(pm0) || any(names(pm0) != c("p0", "p1")) ){
    pm0 <- prior_p_sampler(a=a, b=b, prior.pm=prior.pm )
    while(check.p(pm0,k0)) pm0 <- prior_p_sampler(a=a, b=b, prior.pm=prior.pm)
    if(warn){
      warning(paste("p0 missing or not correctly defined, random generated from a ", prior.pm ," distr.",sep=""))
    }  
  }
  
  if(missing(param0) || is.null(param0) || !is.list(param0) || any(names(param0) != c("eta", "beta")) || length(param0$eta) != (k0+1) || length(param0$beta) != (k0+2) ){
    param0 <- rcoef(k0, pm0)  
    if(warn){
      warning('param0 missing or not correctly defined, random generated from rcoef(k0, pm0)')
    }  
  }
  
  return(list(prior.k = prior.k, prior.pm = prior.pm, hyperparam = hyperparam, pm0 = pm0, param0 = param0, k0 = k0, a=a, b=b, mu.pois=mu.pois, pnb=pnb, rnb=rnb))
  
}


###############################################################################
## Fit.Empirical function (EKdH estimator of extremal dependence)            ##

Fit.Empirical <- function(data){
  
  nw <- 100
  fi <- AngularMeasure(data[,1], data[,2], k=nw, method="c", plot=FALSE)

  w <- fi$weights
  loc <- fi$angles
  lok <- pi/2 - loc
  m <- length(loc)
  h <- 0.45
  
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
  
  hhat <- function(w){
    0.5 * (w^2 + (1-w)^2)^(-3/2) * fi_hat(atan((1-w)/w))
  }
  
  theta.seq <- seq(from=0, to=pi/2, length=nw)
  psi_hat <- cbind(theta.seq, sapply(theta.seq, fi_hat))
  
  w.seq <- seq(from=0, to=1, length=nw)
  h_hat <- cbind(w.seq, sapply(w.seq, hhat))
  
  return(list(method="Empirical", fi=fi, psi_hat=psi_hat, h_hat=h_hat, theta.seq=theta.seq, data=data))
  
}

### Modification of the ecdf function to be divided by n+1 rather than n

ecdf2 <- function(x){ 
  x <- sort(x)
  n <- length(x)
  if (n < 1) 
    stop("'x' must have 1 or more non-missing values")
  vals <- unique(x)
  rval <- approxfun(vals, cumsum(tabulate(match(x, vals)))/(n+1), 
                    method = "constant", yleft = 0, yright = 1, f = 0, ties = "ordered")
  class(rval) <- c("ecdf", "stepfun", class(rval))
  assign("nobs", n, envir = environment(rval))
  attr(rval, "call") <- sys.call()
  rval
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



