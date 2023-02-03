#######################################################
### Authors: Giulia Marcon and Simone Padoan        ###
### Emails: giulia.marcon@phd.unibocconi.it,        ###
### simone.padoan@unibocconi.it                     ###
### Institution: Department of Decision Sciences,   ###
### University Bocconi of Milan                     ###
### File name: beed.r                               ###
### Description:                                    ###
### This file provides an estimate of the uni- and  ###
### multi-variate Pickands dependence function,     ###
### based on an Bernstein polynomials approximation ###
### as proposed in Marcon et al. (2016)             ###
### Last change: 15/08/2016                         ###
#######################################################

############################################################################
# Bernstein Estimation of Extremal Dependence.                           ###
# Pickands Dependence Function Estimate is provided.                     ###
# NB. This relies on the representation of the                           ###
# exponent measure function assuming exponential marginals if using the  ###
# Madogram, whereas use frechet marginals when using the other           ###
# estimators for the pilot estimates                                     ###
############################################################################
################################################################################
# INPUT:                                                                     ###
# x is the d-dimensional simplex                                             ###
# y is the pilot estimate. If missing, it is estimated                       ###
# est = c('md','cfg','ht','pick') choses the method for the initial estimate ###
# k is the order of polynomials                                              ###
# d is the dimension of data                                                 ###
# beta is the vector of coefficient for the Bernstein representation.        ###
#   If it's missing, it is estimated                                         ###
################################################################################ 

beed <- function(data, x, d=3, est = c("ht","cfg","md","pick"), 
                 margin=c("emp", "est", "exp", "frechet", "gumbel"),
                 k = 13, y = NULL, beta = NULL, 
                 plot = FALSE){
  
  datadim <- dim(data)
  xdim <- dim(x)
  
  # Check conditions
  if(d!=datadim[2]) 
    stop("`data' must be a matrix with `d' columns")
  if (!is.matrix(x) || xdim[2] != d) 
    stop("`x' must be a matrix with `d' columns")
  if (any(x < 0, na.rm = TRUE)) 
    stop("`x' must be non-negative")
  rs <- rowSums(x)
  nx <- xdim[1]
  
  if (any(rs < 0, na.rm = TRUE)) 
    warning("row(s) of `x' which do not lay in the simplex are removed")
  
  # Subroutine that provides index beta coefficients
  # p(k,d) = choose(k+d-1,d-1) number of coefficients 
  index <- function(k, d){
    beta.index <- expand.grid(rep(list(0:k),d))
    beta.index.sum <- rowSums(beta.index)
    restr <- which(beta.index.sum<=k)
    v <- as.matrix(beta.index[restr,ncol(beta.index):1])
    return(v)
  }
  
  # Bivariate Case
  if(d==2){
    # Create design matrix:
    xx <- as.matrix(x[,-d])
    vb <- index(k=k,d-1)
    q <- nrow(vb) # q = choose(k+d-1,d-1)
    Z <- bp(x=x,k=k,v=vb)
    
    # Preliminary Estimation if missing
    if(is.null(y)){
      
      if(length(est)>1){
        est='md'
        warning("invalid argument for `est', madogram by default")
      }
      if(length(margin)>1){
        margin='emp'
        warning("invalid argument for `margin', empirical transformation dy default")
      }
      
      if(margin=='emp')
        y <- switch(est, ht = abvnonpar(x=xx, data=data, method="pickands", d=d, madj=2, epmar = TRUE), 
                    md = madogram(w = x, data = data, margin = "emp"),
                    pick = An.biv(x=data, w=xx, estimator = "Pickands"),
                    cfg = An.biv(x=data, w=xx, estimator = "CFG")
        )
      if(margin=='est')
        y <- switch(est, ht = abvnonpar(x=xx, data=data, method="pickands", d=d, madj=2, epmar = FALSE), 
                    md = madogram(w = x, data = data, margin = "est"),
                    pick = An.biv(x=data, w=xx, estimator = "Pickands"),
                    cfg = An.biv(x=data, w=xx, estimator = "CFG")
        )
      if((margin!="emp")&(margin!="est")){
        y <- switch(est, ht = abvnonpar(x=xx, data=data, method="pickands", d=d, madj=2, epmar = FALSE), 
                    md = madogram(w = x, data = data, margin = margin),
                    pick = An.biv(x=data, w=xx, estimator = "Pickands"),
                    cfg = An.biv(x=data, w=xx, estimator = "CFG")
        )
        if((est=="ht")|(est=="pick")|(est=="cfg")) 
          warning("Marginal distributions estimated")
      }
    }
    
    n <- length(y)
    
    # Estimation procedure, else assign beta vector input
    if(is.null(beta)){
      
      # Create the A, b0 matrix for constraints
      # R1) Convexity
      R1 <- matrix(0,k-1,k+1)
      for (i in 1:(k-1))
        R1[i,i:(i+2)] <- c(1,-2,1)
      r1 <- rep(0,nrow(R1))
      
      # R2) upper bound
      # A(ei) = 1
      R2 <- matrix(0,2*d,q)
      vertix <- which(apply(vb, 1, function(x) any(x==k)))
      R2[1,1] <- 1
      R2[(d+1),1]<- -1
      
      for(i in 1:(d-1)) {
        R2[i+1,vertix[i]] <- 1
        R2[i+1+d,vertix[i]] <- -1
      }
      r2 <- c(rep(1,d),rep(-1,d))
      
      # R3) lower bound
      # A(w) >= max(w)
      nr <- d*(d-1)
      R3 <- matrix(0,nr,q)
      
      vertix_10 <- which(vb==1) 
      vertix_k_10 <- which(vb==k-1)
      vertix <- c(vertix_10,vertix_k_10)
      
      for(i in 1:nr) 
        R3[i,vertix[i]] <- 1
      r3 <- rep(1-1/k,nr)
      
      R = rbind(R1, R2, R3)
      r0 = c(r1, r2, r3) 
      
      Dmat  <- t(Z)%*%Z          
      dvec  <- t(y)%*%Z
      Rmat  <- t(R)     
      
      fit <- solve.QP(Dmat,dvec,Rmat,r0)
      beta.tilde <- as.matrix(fit$solution)
      
    }
    else beta.tilde <- beta
    
    y.tilde <- Z%*%beta.tilde
    yy <- y
    
    if (plot == TRUE){
      plot(x[,1],x[,2],type='n',xlab='t',ylab='A(t)',ylim=c(.5,1))
      polygon(c(0, 0.5, 1), c(1, 0.5, 1), lty = 1, lwd = 1, border = 'grey')
      lines(x[,1],y.tilde,lty=1,col=1)
    }
  }
  
  # multivariate case
  else{
    zero1 <- which(rs>1)
    zero0 <- which(rs==0)
    zeroNa <- which(is.na(rs))
    if((length(zero1)==0)&(length(zeroNa)==0)&(length(zero0)==0)) 
      X <- x
    else{ 
      X <- x[-c(zero1,zeroNa,zero0),]
    }
    
    # Create design matrix:
    xx <- as.matrix(X[,-d])
    vb <- index(k,d-1)
    q <- nrow(vb) 
    Z <- bp(x=X,k=k,v=vb)
    # Preliminary Estimation if missing
    if(is.null(y)){
      
      if(length(est)>1){
        est='md'
        warning("invalid argument for `est', madogram by default")
      }
      if(length(margin)>1){
        margin='emp'
        warning("invalid argument for `margin', empirical transformation dy default")
      }
      
      if(margin=='emp')
        y <- switch(est, ht = amvnonpar(x=X, data=data, d=d, madj=2, epmar = TRUE),
                    md = madogram(w = X, data = data, margin="emp"),
                    cfg = An(x = data, w = X)$CFG,
                    pick = An(x = data, w = X)$P
        )
      if(margin=='est')
        y <- switch(est, ht = amvnonpar(x=X, data = data, d=d, madj=2, epmar = FALSE),
                    md = madogram(w = X, data = data, margin = "est"),
                    cfg = An(x = data, w = X)$CFG,
                    pick = An(x = data, w = X)$P
        )
      if((margin!="emp")&(margin!="est")){
        y <- switch(est, ht = amvnonpar(x=X, data = data, d=d, madj=2, epmar = FALSE),
                    md = madogram(w = X, data = data, margin = margin),
                    cfg = An(x = data, w = X)$CFG,
                    pick = An(x = data, w = X)$P
        )
        if((est=="ht")|(est=="pick")|(est=="cfg")) 
          warning("Marginal distributions estimated")
      }
    }
    
    n <- length(y)
    
    yy <- rep(NA,nx)
    if((length(zero1)==0)&(length(zeroNa)==0)&(length(zero0)==0)) yy <- y
    else yy[-c(zero1,zeroNa,zero0)] <- y
    
    
    # Estimation procedure, else assign beta vector input
    if(is.null(beta)){
      
      # Create the A, b0 matrix for constraints
      # R1) Convexity
      #R1 <- convexity(d = d, k = k)
      #r1 <- rep(0,nrow(R1))
      R1 <- convexity_old(v = vb, d = d)
      r1 <- rep(0,nrow(R1))
      
      # R2) upper bound
      # A(ei) = 1
      R2 <- matrix(0,2*d,q)
      vertix <- which(apply(vb, 1, function(x) any(x==k)))
      R2[1,1] <- 1
      R2[(d+1),1]<- -1
      
      for(i in 1:(d-1)) {
        R2[i+1,vertix[i]] <- 1
        R2[i+1+d,vertix[i]] <- -1
      }
      r2 <- c(rep(1,d),rep(-1,d))
      
      # R3) lower bound
      # A(w) >= max(w)
      if(d==3){
        nr <- d*(d-1)
        R3 <- matrix(0,nr,q)
        
        v_10 <- which(apply(vb, 1, function(x) any(x==1) & 
                              any(x==0) & sum(x)==1)) 
        v_k_10 <- which(apply(vb, 1, function(x) any(x==k-1) & 
                                any(x==0)& sum(x)==k-1))
        v_k_11 <- which(apply(vb, 1, function(x) any(x==k-1) &
                                any(x==1)& sum(x)==k)) 
        vertix <- c(v_10,v_k_10,v_k_11)
        
        for(i in 1:nr) 
          R3[i,vertix[i]] <- 1
        r3 <- rep(1-1/k,nr)
      }
      if(d>3){
        R3 <- Z
        r3 <- apply(X,1,max)
      }
      
      R = rbind(R1, R2, R3)
      r0 = c(r1, r2, r3) 
      
      Dmat  <- t(Z)%*%Z          
      dvec  <- t(y)%*%Z
      Rmat  <- t(R)     
      
      fit <- solve.QP(Dmat,dvec,Rmat,r0)
      beta.tilde <- as.matrix(fit$solution)
      
    }
    else beta.tilde <- beta
    
    y.tilde <- rep(NA,nx)
    if((length(zero1)==0)&(length(zeroNa)==0)&(length(zero0)==0)) y.tilde <- Z%*%beta.tilde
    else y.tilde[-c(zero1,zeroNa,zero0)] <- Z%*%beta.tilde
    
    if (plot == TRUE){
      if(d == 3){
        if((length(zero1)==0)&(length(zeroNa)==0)&(length(zero0)==0)){
          numg <- round((-1 + sqrt(1 + 8*xdim[1]))/2)
          mat <- matrix(NA, numg, numg)
          mat[lower.tri(mat, diag = T)] <- y.tilde
        }
        else{
          numg <- sqrt(nx)
          mat <- matrix(y.tilde, numg, numg)
        }
        xy <- seq(0,1,length=numg)
        plot(xy, xy, type='n', xlab=expression(t[1]), ylab=expression(t[2]))
        image(x=xy, y=xy, z=mat, col=heat.colors(numg),add=TRUE)
        contour(x=xy, y=xy, z=mat, add=T,col='black',labcex=.6,lty=1, cex=2)
      }
      if (d >= 4) 
        warning("can not plot in high dimensions")
    }
    
    
  }
  
  med <- matrix(rep(1/d,d),ncol=d)
  Zmed <- bp(med,v=vb,k=k)
  extind <- as.numeric( d * (Zmed%*%beta.tilde) )
  
  out <- list(beta=beta.tilde,A=y.tilde, Anonconvex=yy, extind=extind)
  return(out)
}


########################################################################
# Define a multivariate Bernstein Polynomials on the Simplex         ###
# X %*% beta = Bernstein Representation                              ###
# Input: x is the d-dimensional simplex, v index set, k polynomial   ###
# order                                                              ###
# bp returns a design matrix X                                       ###
########################################################################

bp <- function(x,k,v=NULL)
{
  
  x <- as.data.frame(x)
  dd <- ncol(x)-1
  xx <- x[,(1:dd)]
  if (dd == 1){
    #x <- as.matrix(sort(xx))
    q <- nrow(x) # number of points of the grid
    basis <- array(0, dim=c(q, k+1))
    
    for(j in 0:k){
      basis[,j+1] <- choose(k,j) * x[,1]^j * x[,2]^(k-j)}
    
    return(basis)
  }
  
  else x <- as.matrix(xx)
  q <- nrow(v)
  n <- nrow(x)
  b <- matrix(NA,n,q)
  xv <- matrix(NA,n,dd)
  
  for (i in 1:q) {
    for (j in 1:dd) {
      xv[,j] <- x[,j]^v[i,j]
    }
    b[,i] <- factorial(k)/( prod(factorial(v[i,])) * factorial(k-sum(v[i,]))) * 
      apply(xv,1,prod) * (1-rowSums(x))^(k-sum(v[i,]))
  }
  return(b)
}

########################################################################
########################## BIVARIATE CASE ##############################
bp2d <- function(w, k)
{
  #ind <- index(k) # compute the index
  q <- nrow(w) # number of points of the grid
  basis <- array(0, dim=c(q, k+1))
  
  for(j in 0:k){
    basis[,j+1] <- choose(k,j) * w[,1]^j * w[,2]^(k-j)}
  
  return(basis)
}

########################## BIVARIATE CASE ##############################
########################################################################

########################################################################
# Define a multivariate Unit-Simplex                                 ###
# Input: n is the length of the one-dimension simplex                ###
# simplex returns a matrix of the d-dimensional Simplex              ###
########################################################################

simplex <- function(d, n=50, a=0, b=1){
  if(d<2){stop("'d' must be an integer taking minimum value 2")}
  if(n<2){stop("'n' must be an integer taking minimum value 2")}
  if(a >=b){stop("'a' must be smaller than 'b'")}
  if(a <0){stop("'a' must be positive")}
  w <- seq(a, b, length = n)
  if(d==2) return(cbind(w,1-w))
  else{
    ww <- expand.grid(rep(list(w), d-1))
    os <- which(rowSums(ww)>1)
    x <- cbind(ww, 1 - rowSums(ww))
    x <- x[-os,]
    colnames(x) = paste("w", 1:d, sep="")
    return(as.matrix(x))
  }
}