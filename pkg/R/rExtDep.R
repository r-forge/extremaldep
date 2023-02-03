rExtDep <- function(n, model, par, angular=FALSE, mar=c(1,1,1), num,
                    threshold, exceed.type){
  
  models <- c("semi.bvevd", "semi.bvexceed", "HR", "ET", "EST")  
  if(!any(model ==  models)){ stop("model wrongly specified")}
  
  if(missing(num)){
    if( model %in% c("HR", "ET", "EST")){ num <- 5e+5}
    if(model == "semi.bvevd"){ num <- 100}
  }
  
  if(model == "semi.bvexceed"){
    if(missing(threshold)){stop(" 'threshold' must be provided.")}
    if(length(threshold)!=2){stop("'threshold' should be a bivariate vector")}
    if(!(exceed.type) %in% c("and", "or")){stop(" 'exceed.type' must be 'and' or 'or'.")}
  }
  
  if(model=="HR"){
    data <- rhusler(n=n, lambda=par, num=num)
  }
  if(model=="ET"){
    data <- rextremalt(n=n, param=par, num=num)
  }
  if(model=="EST"){
    data <- rmextst(n=n, param=par, num=num)
  }
  
  if(angular){
    
    if(model %in% c("semi.bvevd", "semi.bvexceed")){
      data <- rh(n, beta=par)
    }
    
    if( model %in% c("HR", "ET", "EST")){
      data <- data / rowSums(data)
    }
  }
  
  if(!angular){ # modify the marginal distributions
    
    if(model == "semi.bvevd"){
      data <- rsemibevd(n=n, beta=par, num=num)
    }
    
    if(model == "semi.bvexceed"){
      
      u.star <- function(threshold, mar){
        if(mar[3] == 0){
          return( exp((threshold - mar[1])/mar[2]) )
        }else{
          return( (1 + mar[3]*(threshold-mar[1])/mar[2])^(1/mar[3]) )
        }
      }
      
      if(is.vector(mar)){
        ustar1 <- u.star(threshold[1], mar) 
        ustar2 <- u.star(threshold[2], mar) 
      }else if(is.matrix(mar)){
        ustar1 <- u.star(threshold[1], mar[1,]) 
        ustar2 <- u.star(threshold[2], mar[2,])         
      }
    
      data <- rsemiexceed(n=n, beta=par, ustar=c(ustar1,ustar2), type=exceed.type)
      
    }
    
    if(!all(mar==1)){
      
      gev.mar <- function(data, mar){
        
        if(length(mar) != 3){ stop("Wrong length of parameter vector")}
        
        if(all(mar==1)){
          out <- data
        }else{
          if(mar[3] !=0){
            out <- (data^mar[3]-1) * mar[2] / mar[3] + mar[1]  
          }else{
            out <- apply(data, 2, function(x) qGEV(pGEV(q=x, loc=1, scale=1, shape=1),loc=mar[1], scale=mar[2], shape=mar[3]) )
          }
        }
        
        return(out)
      }
      
      d <- ncol(data)
      
      if(is.vector(mar) && length(mar)==3){
        mar <- matrix(rep(mar,d), ncol=3)
      }
      if(is.matrix(mar)){
        if(ncol(mar)!=3){stop("mar must have 3 columns")}
        if(nrow(mar)!=d){stop(cat("mar must have ", d, " rows"))}
      }
      
      for(i in 1:d){
        data[,i] <- gev.mar(data=data[,i], mar=mar[i,])
      }
        
    }

  }
  
  return(data)
}


###
### Husler-Reiss
###

rhusler <- function(n, lambda, num=5e+5){
  dim <- dim_ExtDep(model="HR", par=lambda)
  randU <- matrix(ncol=dim);
  rho <- 1 - lambda^2/log(num)
  Sigma <- diag(dim)
  Sigma[lower.tri(Sigma)] = rho
  Sigma[row(Sigma) < col(Sigma)] = t(Sigma)[row(Sigma) < col(Sigma)]
  bn <- sqrt( 2*log(num) - log(log(num))-log(4*pi) )
  an <- 1/bn
  
  while(nrow(randU)<n+1){
    sim <- rmvnorm(num, sigma = Sigma ) 
    sim <- ( sim - bn ) /an 	# Normalization
    randU <- rbind(randU,apply(sim,2,max))
    
  }
  randU <- exp(randU) # Transform from Frechet with parameter df to Standard Frechet
  return(randU[c(2:(n+1)),])		
}

###
### Extremal-t model
###

rextremalt <- function(n, param, num=5e+5){
  
  dim <- dim_ExtDep(model="ET", par=param)
  
  Sigma <- diag(dim)
  Sigma[lower.tri(Sigma)] <- param[1:choose(dim,2)]
  Sigma[upper.tri(Sigma)] <- param[1:choose(dim,2)]			
  
  df <- param[-(1:choose(dim,2))]
  
  res <- matrix(double(n*dim), ncol=dim, nrow=n)
  cst <- gamma((df+1)/2) / gamma(df/2) * (df*pi)^(-1/2) * df^((df-1)/2);
  an <- (num*cst)^(1/df)
  
  for(i in 1:n){
    sim <- rmst(n=num, Omega=Sigma, alpha=rep(0,dim), nu=df) 
    maxima <- apply(sim, 2, max)
    res[i,] <- maxima/an
  }
  return(res^df) # Transform from Frechet with parameter df to Standard Frechet
}

###
### Extremal Skew-t model
###

rmextst <- function(n, param, num=5e+5){
  
  dim <- dim_ExtDep(model="EST", par=param)
  
  Sigma <- diag(dim)
  Sigma[lower.tri(Sigma)] <- param[1:choose(dim,2)]
  Sigma[upper.tri(Sigma)] <- param[1:choose(dim,2)]
  
  shape <- param[choose(dim,2)+(1:dim)]
  df <- param[-(1:(choose(dim,2)+dim))]
  
  res <- matrix(double(n*dim), ncol=dim, nrow=n)
  df1 <- df+1
  df2 <- df+2
  an <- (gamma(0.5*df1)*df^(df/2)*pt(shape*sqrt(df1),df))/(gamma(0.5*df2)*sqrt(pi))
  an <- (num*an)^(1/df)
  for(i in 1:n){
    sim <- rmst(n=num, Omega=Sigma, alpha=shape, nu=df)
    maxima <- apply(sim, 2, max)
    res[i,] <- maxima/an
  }
  return(res^df) # Transform from Frechet with parameter df to Standard Frechet
}

###
### Semi-parametric
###

rsemibevd <- function(n, beta, num){

  y <- matrix(nrow=n, ncol=2)
  for(i in 1:n){
    # simulate observations from the angular density
    w <- rh(num, beta)
    # simulate the radial component
    # r <- rfrechet(num)
    #r <- rexp(num)^(-1)
    r <- cumsum(rexp(num))^(-1)
    # compute the componentwise maxima
    #y[i,] <- 2*apply(cbind(r*w, r*(1-w))/num,2,max)
    y[i,] <- 2*apply(cbind(r*w, r*(1-w)),2,max)
  }
  return(y)
}

rsemiexceed <- function(n, beta, ustar, type){
  
  y <- matrix(nrow=n, ncol=2)
  i <- 1
  
  while(i <= n){
    # simulate observations from the angular density
    w <- rh(1, beta)
    while(w==0 | w==1){w <- rh(1, beta)}
    # simulate the radial component
    r <- runif(1)^(-1)
    # Unit Frechet marginal distributions
    y_temp <- cbind(2*r*w, 2*r*(1-w))
    
    if((type == "and" & all(y_temp>ustar)) | (type == "or" & any(y_temp>ustar))  ){
      y[i,] <- y_temp
      i <- i+1
    }
    
  }
  return(y)
}

### Simulate from a mixture of beta densities 
# 
#
#
rh.mixt <- function(N=1000, beta){
  k <- length(beta) - 1
  # Sample N random uniforms U
  U = runif(N)
  
  #Variable to store the samples from the mixture distribution                                             
  rand.samples = rep(NA,N)
  
  #probabilities mixture
  prob <- diff(diff(beta)) / (2-beta[2]-beta[k])
  
  #Sampling from the mixture
  rand.samples <- sapply(1:N, function(i, prob, k, U){ 
    j <- min(which((U[i]<cumsum(prob)))); return(rbeta(1,j,k-j))},
    prob, k, U)
  
  return(rand.samples)
}

### Simulate from an angular density
#
#
#
rh <- function(N=1000, beta){
  k <- length(beta) - 1
  #Sample N random uniforms U
  U = runif(N)
  
  p0 <- (1+k*(beta[2]-1))/2
  p1 <- (1-k*(1-beta[k]))/2
  prob <- c(p0, 1-p0-p1, p1)
  
  rand.samples <- sapply(1:N, function(i, prob, U){
    j <- min(which((U[i]<cumsum(prob))));
    return(switch(j, 0, rh.mixt(1,beta),1))},
    prob,U)
}


