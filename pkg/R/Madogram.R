#######################################################
### Authors: Giulia Marcon and Simone Padoan        ###
### Emails: giulia.marcon@phd.unibocconi.it,        ###
### simone.padoan@unibocconi.it                     ###
### Institution: Department of Decision Sciences,   ###
### University Bocconi of Milan                     ###
### File name: Madogram.r                           ###
### Description:                                    ###
### This file enables to compute the multivariate   ###
### madogram as proposed in Marcon et al. (2014)    ###
### Last change: 03/12/2014                         ###
#######################################################

# The routine estimates the Pickands function using
# the multivariate madogram
# Multivariate case

madogram <- function(w,data){
  
  lmadogram <- function(x,data){
    sumdata <- dim(data)
    if(!is.matrix(w)) w <- matrix(w, ncol= sumdata[2])
    sumw <- dim(w)
    if (sumw[2] != sumdata[2]) 
      stop("`x' must be a vector/matrix with `d' elements/columns")
    
    ans <- numeric(sumw[1])
    u <- matrix(data,sumdata[1],sumdata[2])

    for (i in 1:sumw[1]){
      for (j in 1:sumdata[2]){
        e <- ecdf(data[,j])
        u[,j] <- e(data[,j])
        u[,j] <- u[,j]^(1/w[i,j])  
      }
      ma <- apply(u,1,max)
      me <- apply(u,1,mean)
      ans[i] <- mean(ma-me)
    }
    
    return(ans)
  }
  
  ans <- lmadogram(x=w,data)
  sumdata <- dim(data)
  if(!is.matrix(w)) w <- matrix(w, ncol= sumdata[2])
  W <- w/(1+w)
  cst <- apply(W,1,mean)
  A <- ((cst + ans) / (1 - cst - ans))
  return(A)
}
