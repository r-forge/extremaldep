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
### as proposed in Marcon et al. (2014)             ###
### Last change: 01/12/2014                         ###
#######################################################

####################################################################
# Bernstein Estimation of Extremal Dependence.
# Pickands Dependence Function Estimate is provided.
# NB. This relies on the representation of the 
# exponent measure function assuming exponential marginals if using the 
# Madogram, whereas use frechet marginals when using the other 
# estimators for the pilot estimates
####################################################################
####################################################################
# INPUT:
# x is the d-dimensional simplex
# y is the pilot estimate. If missing, it is estimated
# est =c('md','cfg','ht') choses the method for the initial estimate
# k is the order of polynomials 
# d is the dimension of data
# beta is the vector of coefficient for the Bernstei representation.
#   If it's missing, it is estimated 
####################################################################

beed <- function(data, x, d=3, est = c("ht","cfg","md"), 
                 margin = c("emp","Gev"), k = 13, y = NULL, beta = NULL, 
                 matrix = FALSE, plot = FALSE){
  ddata <- ncol(data)
  
  # Check conditions
  if(d!=ddata) 
    stop("`data' must be a matrix with `d' columns")
  if (is.vector(x)) 
    x <- as.matrix(t(x))
  if (!is.matrix(x) || ncol(x) != d) 
    stop("`x' must be a vector/matrix with `d' elements/columns")
  if (any(x < 0, na.rm = TRUE)) 
    stop("`x' must be non-negative")
  rs <- rowSums(x)
  nx <- nrow(x)

  if (any(rs < 0, na.rm = TRUE)) 
    warning("row(s) of `x' which do not lay in the simplex are removed")
  zero1 <- which(rs>1)
  zero0 <- which(rs==0)
  zeroNa <- which(is.na(rs))
  if((length(zero1)==0)&(length(zeroNa)==0)&(length(zero0)==0)) 
    X <- x
  else{ 
    X <- x[-c(zero1,zeroNa,zero0),]
  }
  
  # Subroutine that provides index beta coefficients
  # nrow(v) = p(k,d) number of coefficients 
  index <- function(k,d){
    #vec <- matrix(rep(0:k,d),d,k+1,byrow=T)
    if(d==1)
      beta.index <- expand.grid(0:k)
    if(d==2)
      beta.index <- expand.grid(0:k,0:k)
    if(d==3)
      beta.index <- expand.grid(0:k,0:k,0:k)
    if(d==4)
      beta.index <- expand.grid(0:k,0:k,0:k,0:k)
    if(d==5)
      beta.index <- expand.grid(0:k,0:k,0:k,0:k,0:k)
    if(d==6)
      beta.index <- expand.grid(0:k,0:k,0:k,0:k,0:k,0:k)
    if(d==7)
      beta.index <- expand.grid(0:k,0:k,0:k,0:k,0:k,0:k,0:k)
    beta.index.sum <- apply(beta.index,1,sum)
    restr <- which(beta.index.sum<=k)
    v <- as.matrix(beta.index[restr,ncol(beta.index):1])
    return(v)
  }
  
  # Create design matrix:
  xx <- as.matrix(X[,-2])
  vb <- index(k,d-1)
  q <- nrow(vb) # q= nrow(combn_eqk)
  Z <- bp(X,vb,k)
  
  # Estimation procedure, else assign beta vector input
  if(is.null(beta)){
    
    # Preliminary Estimation if missing
    if(is.null(y)){
      
      if(length(est)>1) stop("invalid argument for `est'")
      if(length(margin)>1) stop("invalid argument for `margin'")
      
      if(margin=='emp'){
        # as x is dimension d with NaN
        data_emp <- apply(data, 2, rank, na.last = "keep")
        nasm <- apply(data_emp, 2, function(x) sum(!is.na(x)))
        data_emp <- data_emp/rep(nasm + 1, each = nrow(data_emp))
        data <- data_emp
      }
      if(margin=='Gev'){
        data <- Dist2Dist(data,from='Gev',to='sFrechet')
      }
      
      if(d==2)
        y <- switch(est, ht = abvnonpar(x=xx, data=data, method="pickands", d=d, madj=2, epmar = TRUE), 
                  md = madogram(w = X, data = data),
                  cfg = abvnonpar(x=xx, data=data, method="cfg", epmar = TRUE)
      )
      else
        y <- switch(est, ht = amvnonpar(x=X, data=data, d=d, madj=2, epmar = TRUE),
                    md = madogram(w = X, data = data),
                    cfg = An(x = data, w = X)$CFG
        )
    }
    n <- length(y)
    
    # Create the A, b0 matrix for constraints
    # A(ei) = 1
    A2 <- matrix(0,d,q)
    vertix <- which(apply(vb, 1, function(x) any(x==k)))
    A2[1,1] <- 1
    for(i in 1:(d-1)) A2[i+1,vertix[i]] <- 1
    a2 <- rep(1,d)
    
    A4 <- matrix(0,d,q)
    A4[1,1] <- -1
    for(i in 1:(d-1)) A4[i+1,vertix[i]] <- -1
    a4 <- rep(-1,d)
    
    # A(w) >= max(w)
    A3 <- Z
    a3 <- apply(X,1,max)
    
    # Convexity
      
      Aconvx <- convexity(v = vb, d = d)
      aconvx <- rep(0,nrow(Aconvx))
      A = rbind(A2,A4,A3,Aconvx)
      b0 = c(a2,a4,a3,aconvx) 
    
    
    Dmat  <- t(Z)%*%Z          
    dvec  <- t(y)%*%Z
    Amat  <- t(A)     # transpose(Amat) = Am
    
    fit <- solve.QP(Dmat,dvec,Amat,b0)
    beta.tilde <- as.matrix(fit$solution)
  }
  else beta.tilde <- beta
  
  y.tilde <- rep(NA,nx)
  if((length(zero1)==0)&(length(zeroNa)==0)&(length(zero0)==0)) y.tilde <- Z%*%beta.tilde
    else y.tilde[-c(zero1,zeroNa,zero0)] <- Z%*%beta.tilde
  
  if (matrix == TRUE){
    if(d!=3) 
      warning("Not possible to replace in matrix form")
    numg <- abs(-1/2 + sqrt(1/4 + 2*nx))
    xy <- seq(0,1,length=numg)
    math <- matrix(NA, numg,numg)
    math[1,1:numg] <- y.tilde[1:numg]
    for(i in 0:(numg-2))
      math[i+2,1:(numg-(i+1))] <- y.tilde[((i+1)*numg-i*(i-1)/2+1-i):((i+2)*numg-i*(i+1)/2-(i+1))]
    y.tilde <- math
  } 
  
  if (plot == TRUE){
    if(d == 2){
      plot(x[,1],x[,2],type='n',xlab='w',ylab='A(w)',ylim=c(.5,1))
      polygon(c(0, 0.5, 1), c(1, 0.5, 1), lty = 1, lwd = 1, border = 'grey')
      lines(x[,1],y.tilde,lty=1,col=1)
    }
    if(d == 3){
      numg <- sqrt(nx)
      xy <- seq(0,1,length=numg)
      mat <- matrix(y.tilde, numg, numg)
      plot(xy, xy, type='n', xlab=expression(w[1]), ylab=expression(w[2]))
      image(x=xy, y=xy, z=mat, col=heat.colors(numg),add=TRUE)
      contour(x=xy, y=xy, z=mat, add=TRUE,col='black',labcex=.6,lty=1)
    }
    if (d >= 4) 
      stop("cannot plot in high dimensions")
  }
  
  med <- matrix(rep(1/d,d),ncol=d)
  Zmed <- bp(med,vb,k)
  extind <- as.numeric( d * (Zmed%*%beta.tilde) )
  
  out <- list(beta=beta.tilde,A=y.tilde, Anonconvex=y, extind=extind)
  return(out)
}




# Define a multivariate Bernstein Polynomials on the Simplex
# X %*% beta = Bernstein Representation
# Input: x is the d-dimensional simplex, v index set, k polynomial
# order
# bp returns a design matrix X 
bp <- function(x,v,k)
{
  
  x <- as.data.frame(x)
  dd <- ncol(x)-1
  xx <- x[,(1:dd)]
  if (dd == 1) x <- as.matrix(sort(xx))
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
      apply(xv,1,prod) * (1-apply(x,1,sum))^(k-sum(v[i,]))
  }
  return(b)
}

