#######################################################
### Authors: Giulia Marcon and Simone Padoan        ###
### Emails: giulia.marcon@phd.unibocconi.it,        ###
### simone.padoan@unibocconi.it                     ###
### Institution: Department of Decision Sciences,   ###
### University Bocconi of Milan                     ###
### File name: bbeed.r                              ###
### Description:                                    ###
### This file provides the posterior distribution   ###
### given the Bernstein representation of the       ###
### Pickands dependence function                    ###
### or the angular distribution.                    ###
### as proposed in Marcon et al. (2016)             ###
### Last change: 15/08/2016                         ###
#######################################################

###################### START METROPOLIS-HASTINGS ALGORITHM #####################

################################################################################
# INPUT:                                                                     ###
# data is a (n x 2)-matrix/dataframe                                         ###
# pm0, param0, k0 are the initial values of the parameters                   ###
# hyperparam is a list of the hyperparameters                                ###
# nsim is the number of the iterations of the chain                          ###
# prior.k and prior.pm specify the prior distributions                       ###
################################################################################ 


bbeed <- function(data, pm0, param0, k0, hyperparam, 
                  nsim, nk = 70, prior.k = c('pois','nbinom'), 
                  prior.pm = c('unif','beta', 'null'), lik = TRUE){
  
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
  
  #param0 = eta.start
  n <- nrow(data$w)
  nkm <- nk-1
  nkp <- nk+2
  # set the chains
  spm <- array(0, dim=c(nsim,2))
  seta <- array(0, dim=c(nsim,nkp))
  sk <- rep(0,nsim) 
  
  accepted <-0
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
      cat('maximum value k reached')
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
    
    u <- runif(1)
    if(u<ratio){
      pm <- pm_new
      param <- param_new
      k <- k_new
      if(k==3) q <- 0.5 else q <- 1
      accepted <- accepted + 1
    }
    spm[i,] <- c(pm$p0, pm$p1)
    seta[i, 1:(k+1)] <- param$eta
    sk[i] <- k
    
    if(i %in% print.i) setTxtProgressBar(pb, i) 
  }

  return(list(pm=spm, eta=seta, k=sk, accepted=accepted/nsim, nsim=nsim,
              prior = list(hyperparam=hyperparam, k=prior.k, pm=prior.pm)))
}

###################### END METROPOLIS-HASTINGS ALGORITHM ######################

####################### START RETURN VALUES ####################################

################################################################################
# P(Y_1 > y_1,Y_2 > y_2) = int_0^1 min(w/y_1, 1-w/y_2) h(w) dw               ###
# INPUT:                                                                     ###
# mcmc is the output of bbeed function                                       ###
# summary.mcmc is the output of summary.bbeed function                       ###
# y = c(y_1, y_2) vector/matrix of thresholds                                ###
################################################################################

returns <- function(mcmc, summary.mcmc, y, plot=FALSE){
  # P(X>x,Y>y) = int_0^1 min(w/x,1-w/y) h(w) dw
  # y = c(y_1, y_2) vector of thresholds
  
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
  
  id1 <- which(mcmc$k==summary.mcmc$k.median)
  idb1 <- id1[which(id1>=summary.mcmc$burn)]
  
  etalow <- apply(mcmc$eta[idb1,1:(summary.mcmc$k.median+1)],2,quantile,.05)
  etaup <- apply(mcmc$eta[idb1,1:(summary.mcmc$k.median+1)],2,quantile,.95)
  etamean <- colMeans(mcmc$eta[idb1,1:(summary.mcmc$k.median+1)])
  etamedian <- apply(mcmc$eta[idb1,1:(summary.mcmc$k.median+1)], 2, quantile,.5)
  
  uno <- due <- tre <- quattro <- nrow(y)
  for(i in 1:nrow(y)){
    uno[i] <- returns.int(y=y[i,], etalow)
    due[i] <- returns.int(y=y[i,], etaup)
    tre[i] <- returns.int(y=y[i,], etamean)
    quattro[i] <- returns.int(y=y[i,], etamedian)
  }

  rrmcmc <- rowMeans(cbind(uno,due,tre,quattro))
  
  if(plot) plot.bbeed(type = "returns", mcmc=mcmc, summary.mcmc=summary.mcmc, 
                      y=y, probs=rrmcmc)
    
  return(rrmcmc)
}

#################### END RETURN VALUES #########################################


### Subroutines:

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

check.bayes <- function(prior.k = c('pois','nbinom'), prior.pm = c('unif', 'beta', 'null'), hyperparam, pm0, param0, k0){
  
  if(length(prior.k) != 1 || all(prior.k != c('pois','nbinom') )){
    prior.k <- 'nbinom'
    warning('prior on k not specified: prior.k = "nbinom" by default.')
  }
  
  if(length(prior.pm) != 1 || all(prior.pm != c('unif', 'beta', 'null') )){
    prior.pm <- 'unif'
    warning('prior on p not specified: prior.pm = "unif" by default.')
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
    warning('hyperparam missing and set by default.')
  }
  
  if(prior.k=='pois'){ 
    if(length(hyperparam$mu.pois) != 1){
      hyperparam$mu.pois <- 4
      warning('hyperparam$mu.pois missing, set to 4 by default')
    }
    mu.pois <- hyperparam$mu.pois
  }
  if(!exists("mu.pois")){mu.pois <- 4}
  
  if(prior.k=='nbinom'){
    if(length(hyperparam$mu.nbinom) != 1){
      hyperparam$mu.nbinom <- 4
      warning('hyperparam$mu.nbinom missing, set to 4 by default')
    }
    if(length(hyperparam$var.nbinom) != 1){
      hyperparam$var.nbinom <- 8
      warning('hyperparam$var.nbinom missing, set to 8 by default')
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
      warning('hyperparam$a.unif missing, set to 0 by default')
    }
    if(length(hyperparam$b.unif) != 1){
      hyperparam$b.unif <- 0.5
      warning('hyperparam$b.unif missing, set to 0.5 by default')
    }
    a <- hyperparam$a.unif
    b <- hyperparam$b.unif
  }
  
  if(prior.pm=='beta'){
    
    if(length(hyperparam$a.beta) != 1){
      hyperparam$a.beta <- c(.8,.8)
      warning('hyperparam$a.beta missing, set to (0.8,0.8) by default')
    }
    if(length(hyperparam$b.beta) != 1){
      hyperparam$b.beta <- c(5,5)
      warning('hyperparam$b.beta missing, set to (5,5) by default')
    }
    a <- hyperparam$a.beta
    b <- hyperparam$b.beta
  }
  
  if(missing(k0) || is.null(k0) || length(k0)!=1){
    if(prior.k=='pois'){
      k0 <- rpois(1, mu.pois)  
      warning('k0 missing or length not equal to 1, random generated from a poisson distr.')
    }
    if(prior.k=='nbinom'){
      k0 <- rnbinom(1, size=rnb, prob=pnb)  
      warning('k0 missing or length not equal to 1, random generated from a negative binomial distr.')
    }
  }
  
  if(missing(pm0) || is.null(pm0) || !is.list(pm0) || any(names(pm0) != c("p0", "p1")) ){
    pm0 <- prior_p_sampler(a=a, b=b, prior.pm=prior.pm )
    while(check.p(pm0,k0)) pm0 <- prior_p_sampler(a=a, b=b, prior.pm=prior.pm)
    warning(paste("p0 missing or not correctly defined, random generated from a ", prior.pm ," distr.",sep=""))
  }
  
  if(missing(param0) || is.null(param0) || !is.list(param0) || any(names(param0) != c("eta", "beta")) || length(param0$eta) != (k0+1) || length(param0$beta) != (k0+2) ){
    param0 <- rcoef(k0, pm0)  
    warning('param0 missing or not correctly defined, random generated from rcoef(k0, pm0)')
  }
  
  return(list(prior.k = prior.k, prior.pm = prior.pm, hyperparam = hyperparam, pm0 = pm0, param0 = param0, k0 = k0, a=a, b=b, mu.pois=mu.pois, pnb=pnb, rnb=rnb))
  
}

