#############################################################
### Authors: Giulia Marcon and Simone Padoan              ###
### Emails: giulia.marcon@phd.unibocconi.it,              ###
### simone.padoan@unibocconi.it                           ###
### Institution: Department of Decision Sciences,         ###
### University Bocconi of Milan                           ###
### File name: ConfidenceBands.r                          ### 
### Description:                                          ### 
### This file contains routines that compute              ###
### (1-a)% bootstrap confident bands for the estimated    ###
### Pickands dependence function                          ###
### Last change: 15/08/2016                               ###
#############################################################  

# Subroutine that computes resampled estimates
beed.boot <- function(data, x, d=3, est = c("ht", "md", "cfg", "pick"),
                      margin=c("emp", "est", "exp", "frechet", "gumbel"),
                      k = 13, nboot = 500, y = NULL, 
                      print = FALSE){
  nx <- nrow(x)
  ndata <- nrow(data)
  ddata <- ncol(data)
  if(d!=ddata) 
    stop("`data' must be a matrix/data.frame with `d' columns")
  bootsamp <- matrix(0, nrow=ndata, ncol=d)
  if(is.null(y))
    fit <- beed(data = data, x = x, d = d, est = est, k = k, 
                plot = FALSE, margin = margin)
  else fit <- y
  A.tilde <- fit$A
  beta.tilde <- fit$beta
  
  beta <- bootA <- NULL
  for(i in 1 : nboot)
  {
    indx <- sample(1:ndata, ndata, replace = TRUE)
    bootsamp <- data[indx, ]
    mod <- beed(data = bootsamp, x = x, d = d, est = est, k = k,
                plot = FALSE, margin = margin)
    b <- mod$beta
    beta <- cbind(beta,b)
    bootA <- cbind(bootA, mod$A)
    if(print){
      print.i <- seq(0, nboot, by=100)
      if(i %in% print.i) message(paste(c('iteration', i, 'out of', nboot,'\n')))
    } 
  }
  out=list(A=A.tilde,bootA=bootA,beta=beta)
  invisible(out)
}

# Main routine that computes the confidence bands
beed.confband <- function(data, x, d=3, est = c("ht", "md", "cfg", "pick"), 
                          margin=c("emp", "est", "exp", "frechet", "gumbel"),
                          k = 13, nboot = 500, y = NULL, conf = 0.95,
                          plot = FALSE, print = FALSE){
  nx <- nrow(x)
  ndata <- nrow(data)
  ddata <- ncol(data)
  if(d!=ddata) 
    stop("`data' must be a matrix/data.frame with `d' columns")
  
  # Boostrap
  if(is.null(y))
    fit <- beed.boot(data = data, x = x, d = d, est = est, k = k, 
                     nboot = nboot, print = print,
                     margin = margin)
  else 
    fit <- beed.boot(data = data, x = x, d = d, est = est, k = k, 
                     y = y, nboot = nboot, print = print, margin = margin)
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
                    beta = up.beta, margin = margin)$A
  A.low.beta <- beed(data = data, x = x, d = d, est = est, k = k, 
                     beta = low.beta, margin = margin)$A
  
  if (plot == TRUE){
    if(d == 2){
      plot(x[,1],x[,2],type='n',xlab='w',ylab='A(w)',ylim=c(.5,1))
      polygon(c(0, 0.5, 1), c(1, 0.5, 1), lty = 1, lwd = 1, border = 'grey')
      lines(x[,1],A.tilde,lty=1,col=1)
      lines(x[,1],A.low.beta,lty=2,col=1)
      lines(x[,1],A.up.beta,lty=2,col=1)
    }
    if(d == 3){
      numg <- sqrt(nx)
      xy <- seq(0,1,length=numg)
      mat <- matrix(A.tilde, numg, numg)
      matup <- matrix(A.up.beta, numg, numg)
      matlow <- matrix(A.low.beta, numg, numg)
      
      plot(xy, xy, type='n', xlab=expression(w[1]), ylab=expression(w[2]))
      image(x=xy, y=xy, z=mat, col=heat.colors(numg),add=TRUE)
      contour(x=xy, y=xy, z=mat, add=T,col='black',labcex=.6,lty=1)
      contour(x=xy, y=xy, z=matup, add=T,col='black',labcex=.6,lty=2)
      contour(x=xy, y=xy, z=matlow, add=T,col='black',labcex=.6,lty=2)
    }
    if (d >= 4) 
      warning("cannot plot in high dimensions")
  }
  
  out=list(A=A.tilde,A.up.beta=A.up.beta,A.low.beta=A.low.beta,
           bootA=bootA,A.up.pointwise=A.up.pointwise,
           A.low.pointwise=A.low.pointwise,up.beta=up.beta,
           low.beta=low.beta)
  invisible(out)
}
