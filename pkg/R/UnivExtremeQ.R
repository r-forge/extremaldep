#######################################################
### Authors: Boris Beranger and Simone Padoan       ###
### Emails:                                         ###
### (BB) borisberanger@gmail.com,                   ###
### (SP) simone.padoan@unibocconi.it                ###
### Institutions:                                   ###
### (BB) School of Mathematics and Statistics       ###
###      UNSW Sydney, Australia                     ###
### (SP) Department of Decision Sciences            ###
###      University Bocconi of Milan                ###
### File name: UnivExtremeQ.R                       ###
### Description:                                    ###
### This file enables to compute univariate extreme ###
### quantiles as proposed in Beranger et al. (2019) ###
### Last change: 19/07/2019                         ###
#######################################################

ExtQ <- function(P=NULL, method="Frequentist", pU=NULL, 
                 cov=NULL, param=NULL,  param_post=NULL){
  
  methods <- c("Frequentist", "Bayesian")
  if(!any(method ==  methods)){ stop("method wrongly specified")}
  
  if(is.null(P)){stop("Must specify the probability associated with the quantile(s)")}
  if(is.null(cov)){cov <- as.matrix(c(1))}

  if(method=="Frequentist"){
    
    if(length(param) != (ncol(cov)+2)){stop("Wrong length of parameter vector")}
    
    mu <- param[1:(ncol(cov))] %*% t(cov)
    sigma <- param[ncol(cov)+1] 
    gamma <- param[ncol(cov)+2]
    
    if(!is.null(pU)){ # Threshold exceedances
      
      if(ncol(cov)==1){
        Q.est <- mu[1] + sigma * ((pU/P)^gamma-1) / gamma
      }else{
        Q.est <- sapply(P, function(y) mu + sigma * ((pU/y)^gamma-1) / gamma )
      }
    }else{ # Block maxima
      if(ncol(cov)==1){
        Q.est <-  qGEV(p=1-P, loc=mu[1], scale=sigma, shape=gamma)
      }else{
        Q.est <- t(sapply(P,function(y) qGEV(p=1-y, loc=mu, scale=sigma, shape=gamma) ) )
      }
    }
    
  }else if(method=="Bayesian"){
    
    if(missing(param_post)){stop("Matrix of posterior parameters needs to be provided")}
    if(ncol(param_post) != (ncol(cov)+2)){stop("Wrong dimension of param_post")}
    
    mu.ps <-  param_post[,1:(ncol(cov))] %*% t(cov)
    sig.ps <- param_post[, ncol(cov)+1]
    gam.ps <- param_post[, ncol(cov)+2]
    
    if(!is.null(pU)){ # Threshold exceedances
      
      if(ncol(cov)==1){
        Q.est <- matrix(ncol=length(P), nrow=nrow(param_post))
        for(j in 1:length(P)){
          Q.est[,j] <- mu.ps + sig.ps * ((pU/P[j])^gam.ps-1) / gam.ps
        }  
      }else{
        # Q.est <- list(length=nrow(cov))
        # for(j in 1:nrow(cov)){
        #   mat.temp <- matrix(ncol=length(P), nrow=nrow(param_post))
        #   for(l in 1:length(P)){
        #     mat.temp[,l] <- mu.ps[,j] + sig.ps * ((pU/P[l])^gam.ps-1) / gam.ps
        #   }
        #   Q.est[[j]] <- mat.temp  
        # } 
        
        Q.est <- list(length=length(P))
        for(j in 1:length(P)){
          mat.temp <- matrix(ncol=nrow(cov), nrow=nrow(param_post))
          for(l in 1:nrow(cov)){
            mat.temp[,l] <- mu.ps[,l] + sig.ps * ((pU/P[j])^gam.ps-1) / gam.ps
          }
          Q.est[[j]] <- mat.temp  
        } 
        
      }
    }else{ # Block maxima
      if(ncol(cov)==1){
        Q.est <- matrix(ncol=length(P), nrow=nrow(param_post))
        for(j in 1:length(P)){
          Q.est[,j] <- qGEV(p=1-P[j], loc=mu.ps, scale=sig.ps, shape=gam.ps)
        }  
      }else{
        Q.est <- list(length=nrow(cov))
        for(j in 1:nrow(cov)){
          mat.temp <- matrix(ncol=length(P), nrow=nrow(param_post))
          for(l in 1:length(P)){
            mat.temp[,l] <- qGEV(p=1-P[l], loc=mu.ps[,j], scale=sig.ps, shape=gam.ps)
          }
          Q.est[[j]] <- mat.temp  
        }  
      }
    }  
    
  }
      
  return(Q.est)    
  
}

