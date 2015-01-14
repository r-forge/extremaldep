####################################################
### Authors: Giulia Marcon and Simone Padoan
### Emails: giulia.marcon@phd.unibocconi.it,
### simone.padoan@unibocconi.it
### Institution: Department of Decision Sciences,
### University Bocconi of Milan
### File name: ConfidenceBands.r
### Description:
### This file contains routines that compute
### (1-a)% bootstrap confident bands for the estimated 
### Pickands dependence function
### Last change: 05/10/2014.
####################################################

# Subroutine that computes resampled estimates
beed.boot <- function(data, x, d=3, est = c("ht", "md", "cfg"),
                      margin = c("emp","Gev"),
                      k = 13, nboot = 200, y = NULL, matrix = FALSE,
                      print = FALSE){
  nx <- nrow(x)
  ndata <- nrow(data)
  ddata <- ncol(data)
  if(d!=ddata) 
    stop("`data' must be a matrix/data.frame with `d' columns")
  bootsamp <- matrix(0, nrow=ndata, ncol=d)
  if(is.null(y))
    fit <- beed(data = data, x = x, d = d, est = est, k = k, 
                matrix = matrix, plot = FALSE, margin = margin)
  else fit <- y
  A.tilde <- fit$A
  beta.tilde <- fit$beta
  
  beta <- bootA <- NULL
  for(i in 1 : nboot)
  {
    indx <- sample(1:ndata, ndata, replace = TRUE)
    bootsamp <- data[indx, ]
    mod <- beed(data = bootsamp, x = x, d = d, est = est, k = k,
                matrix = FALSE, plot = FALSE, margin = margin)
    b <- mod$beta
    beta <- cbind(beta,b)
    bootA <- cbind(bootA, mod$A)
    if(print) print(i)
  }
  out=list(A=A.tilde,bootA=bootA,beta=beta)
  invisible(out)
}

# Main routine that computes the confidence bands
beed.confband <- function(data, x, d=3, est = c("ht", "md", "cfg"), 
                          margin = c("emp","Gev"),
                          k = 13, nboot = 200, y = NULL, conf = 0.95,
                          matrix = FALSE, plot = FALSE, print = FALSE){
  nx <- nrow(x)
  ndata <- nrow(data)
  ddata <- ncol(data)
  if(d!=ddata) 
    stop("`data' must be a matrix/data.frame with `d' columns")
  
  # Boostrap
  if(is.null(y))
    fit <- beed.boot(data = data, x = x, d = d, est = est, k = k, 
                     nboot = nboot, matrix = matrix, print = print,
                     margin = margin)
  else 
    fit <- beed.boot(data = data, x = x, d = d, est = est, k = k, 
                     y = y, nboot = nboot, matrix = matrix, 
                     print = print, margin = margin)
  A.tilde <- fit$A
  bootA <- fit$bootA
  beta <- fit$beta
  
  alpha <- 1-conf
  # Confidence bands pointwise on x
  A.low.pointwise <- A.up.pointwise <- numeric(nx)
  for (i in 1:nx){
    ord <- sort(bootA[i,])
    A.low.pointwise[i] <- ord[round(nboot*(alpha/2))]
    A.up.pointwise[i] <- ord[round(nboot*(1-alpha/2))]
  }
  
  # Nnumber of coefficients
  p <- nrow(beta)
  
  # Confidence bands on each beta
  low.beta <- up.beta <- numeric(p)
  for (i in 1:p){
    ord <- sort(beta[i,])
    low.beta[i] <- ord[round(nboot*(alpha/2))]
    up.beta[i] <- ord[round(nboot*(1-alpha/2))]
  }
  A.up.beta <- beed(data = data, x = x, d = d, est = est, k = k, 
                    beta = up.beta, matrix = matrix, 
                    margin = margin)$A
  A.low.beta <- beed(data = data, x = x, d = d, est = est, k = k, 
                     beta = low.beta, matrix = matrix,
                     margin = margin)$A
  
  if (plot == TRUE){
    if(d == 2){
      plot(x[,2],x[,1],type='n',xlab='w1',ylab='w2',ylim=c(.5,1))
      polygon(c(0, 0.5, 1), c(1, 0.5, 1), lty = 1, lwd = 1, border = 'grey')
      lines(x[,2],A.tilde,lty=1,col=1)
      lines(x[,2],A.low.beta,lty=2,col=1)
      lines(x[,2],A.up.beta,lty=2,col=1)
    }
    if(d == 3){
      if(matrix == TRUE){
        mat <- A.tilde
        matup <- A.up.beta
        matlow <- A.low.beta
      }
      else{
        numg <- abs(-1/2 + sqrt(1/4 + 2*nx))
        xy <- seq(0,1,length=numg)
        mat <- matup <- matlow <- matrix(NA, numg,numg)
        
        mat[1,1:numg] <- A.tilde[1:numg]
        for(i in 0:(numg-2))
          mat[i+2,1:(numg-(i+1))] <- A.tilde[((i+1)*numg-i*(i-1)/2+1-i):((i+2)*numg-i*(i+1)/2-(i+1))]
        
        matup[1,1:numg] <- A.up.beta[1:numg]
        for(i in 0:(numg-2))
          matup[i+2,1:(numg-(i+1))] <- A.up.beta[((i+1)*numg-i*(i-1)/2+1-i):((i+2)*numg-i*(i+1)/2-(i+1))]
        
        matlow[1,1:numg] <- A.low.beta[1:numg]
        for(i in 0:(numg-2))
          matlow[i+2,1:(numg-(i+1))] <- A.low.beta[((i+1)*numg-i*(i-1)/2+1-i):((i+2)*numg-i*(i+1)/2-(i+1))]
      }
      plot(xy, xy, type='n', xlab='w1', ylab='w2')
      image(x=xy, y=xy, z=mat, col=heat.colors(numg),add=TRUE)
      contour(x=xy, y=xy, z=mat, add=T,col='black',labcex=.6,lty=1)
      contour(x=xy, y=xy, z=matup, add=T,col='black',labcex=.6,lty=2)
      contour(x=xy, y=xy, z=matlow, add=T,col='black',labcex=.6,lty=2)
    }
    if (d >= 4) 
      stop("cannot plot in high dimensions")
  }

  out=list(A=A.tilde,A.up.beta=A.up.beta,A.low.beta=A.low.beta,
           bootA=bootA,A.up.pointwise=A.up.pointwise,
           A.low.pointwise=A.low.pointwise,up.beta=up.beta,
           low.beta=low.beta)
  invisible(out)
}
