#######################################################
### Authors: Boris Beranger and Simone Padoan        ###
### Emails: borisberanger@gmail.com,        ###
### simone.padoan@unibocconi.it                     ###
### Institutions: Department of Decision Sciences,   ###
### University Bocconi of Milan, School of Mathematics and  ###
### Statistics, University of New South Wales ###
### File name: Distribution.r                               ###
### Description:                                    ###
### This file provides the probability distribution functions ###
### of some family of skew-symmetric models     ###
### Last change: 27/11/2014                         ###
#######################################################


# univariate extended skew-normal probability distribution function
pesn <- function(x, loc=0, scale=1, shape=0, ext=0)
{
  par <- c(x,loc,scale,shape,ext)
  res <- .C("pesn", as.double(-1), as.double(0), as.double(par),
     double(1), double(1), PACKAGE="ExtremalDep", NAOK=TRUE)
  return(res[[4]])
}
# bivariate and trivariate extended skew-normal probability distribution function
pmesn <- function(x, loc=rep(0,length(x)), scale=diag(length(x)), shape=rep(0,length(x)), ext=0)
{
  if(any(eigen(scale)$value<0)){stop("Scale matrix is not positive definite")}		
  par <- c(x,loc,scale,shape,ext)
  if(length(x)==2){
    res <- .C("pmesn", as.double(c(-1,-1)), as.double(c(0,0)), as.double(par),
    double(1), double(1), PACKAGE="ExtremalDep", NAOK=TRUE)
  }else if(length(x)==3){
    res <- .C("pmesn3", as.double(c(-1,-1,-1)), as.double(c(0,0,0)), as.double(par),
    double(1), double(1), PACKAGE="ExtremalDep", NAOK=TRUE)		
  }else{
  	stop("x must be of length 2 or 3")
  }	
  return(res[[4]])
}
# univariate extended skew-normal probability density function
desn <- function(x, loc=0, scale=1, shape=0, ext=0)
{
  res <- .C("desn", as.double(x), as.double(loc), as.double(scale), as.double(shape),
     as.double(ext), double(1), PACKAGE="ExtremalDep", NAOK=TRUE)
  return(res[[6]])
}
# Bivariate and trivariate extended skew-normal probability density function
dmesn <- function(x, loc=rep(0,length(x)), scale=diag(length(x)), shape=rep(0,length(x)), ext=0)
{
  if(any(eigen(scale)$value<0)){stop("Scale matrix is not positive definite")}	
  if(length(x)==2){
  	res <- .C("dmesn", as.double(x), as.double(loc), as.double(scale), as.double(shape),
     as.double(ext), double(1), PACKAGE="ExtremalDep", NAOK=TRUE)
  }else if(length(x)==3){
  	res <- .C("dmesn3", as.double(x), as.double(loc), as.double(scale), as.double(shape),
     as.double(ext), double(1), PACKAGE="ExtremalDep", NAOK=TRUE)  	
  }else{
  	stop("x must be of length 2 or 3")	
  }		   
  return(res[[6]])
}
# univariate extended skew-t probability distribution function
pest <- function(x, loc=0, scale=1, shape=0, ext=0, df=1)
{
  if(df==Inf)
    return( pesn(x, loc,scale,shape,ext) )
  else{
    par <- c(x,loc,scale,df,shape,ext)
    res <- .C("pest", as.double(-1), as.double(0), as.double(par),
       double(1), double(1), PACKAGE="ExtremalDep", NAOK=TRUE)
	return(res[[4]])
  }	
}
# bivariate and trivariate extended skew-t probability distribution function
pmest <- function(x, loc=rep(0,length(x)), scale=diag(length(x)), shape=rep(0,length(x)), ext=0, df=1)
{	
  if(any(eigen(scale)$value<0)){stop("Scale matrix is not positive definite")}	
  if(df==Inf)
    return( pmesn(x, loc,scale,shape,ext) )
  else{
    par <- c(x,loc,scale,df,shape,ext)
    if(length(x)==2){
	    res <- .C("pmest", as.double(c(-1,-1)), as.double(c(0,0)), as.double(par),
       double(1), double(1), PACKAGE="ExtremalDep", NAOK=TRUE)
	}else if(length(x)==3){
	    res <- .C("pmest3", as.double(c(-1,-1,-1)), as.double(c(0,0,0)), as.double(par),
       double(1), double(1), PACKAGE="ExtremalDep", NAOK=TRUE)		
	}else{
		stop("x must be of length 2 or 3")
	}	
	return(res[[4]])
  }  
}
# univariate extended skew-t probability density function
dest <- function(x, loc=0, scale=1, shape=0, ext=0, df=1)
{
  if(df==Inf)
    return( desn(x, loc,scale,shape,ext) )
  else{
    res <- .C("dest", as.double(x), as.double(loc), as.double(scale), as.double(df),
       as.double(shape), as.double(ext), double(1), PACKAGE="ExtremalDep", NAOK=TRUE)
  return(res[[7]])
  }
}
# multivariate extended skew-t probability density function
dmest <- function(x, loc=rep(0,length(x)), scale=diag(length(x)), shape=rep(0,length(x)), ext=0, df=1)
{
  if(any(eigen(scale)$value<0)){stop("Scale matrix is not positive definite")}	
  if(df==Inf)
    return( dmesn(x, loc,scale,shape,ext) )
  else{
  	if(length(x)==2){
    		res <- .C("dmest", as.double(x), as.double(loc), as.double(scale), as.double(df),
    		   as.double(shape), as.double(ext), double(1), PACKAGE="ExtremalDep", NAOK=TRUE)
    }else if(length(x)==3){
    		res <- .C("dmest3", as.double(x), as.double(loc), as.double(scale), as.double(df),
   	 	   as.double(shape), as.double(ext), double(1), PACKAGE="ExtremalDep", NAOK=TRUE)   	
    }else{
    		stop("x must be of length 2 or 3")
    }  	
    return(res[[7]])	 
  }     
  
}
