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
pesn <- function(x=0, location=0, scale=1, shape=0, extended=0)
{
  if(!is.atomic(x) || !is.numeric(x) || length(x)!=1)
    stop("x argument has incorrect format")
  if(!is.atomic(location) || length(location) != length(x))
    stop("Location argument has incorrect length")
  if(length(scale) != length(x))
    stop("Scale argument has incorrect dimensions")
  if(!is.atomic(shape) || length(shape) != length(x))
    stop("Shape argument has incorrect length")  
  if(length(extended)!=1)
    stop("Extended argument is not of length 1")
  
  par <- c(x,location,scale,shape,extended)
  .C("pesn", as.double(-1), as.double(0), as.double(par),
     out=double(1), double(1), NAOK=TRUE)$out
}
# bivariate and trivariate extended skew-normal probability distribution function
pmesn <- function(x=c(0,0), location=rep(0,length(x)), scale=diag(length(x)), shape=rep(0,length(x)), extended=0)
{
  if(!is.atomic(x) || !is.numeric(x))
    stop("x argument has incorrect format")
  if(!is.atomic(location) || length(location) != length(x))
    stop("Location argument has incorrect length")
  if(nrow(scale)!=ncol(scale) && nrow(scale) != length(x))
    stop("Scale argument has incorrect dimensions")
  if(!is.atomic(shape) || length(shape) != length(x))
    stop("Shape argument has incorrect length")  
  if(length(extended)!=1)
    stop("Extended argument is not of length 1")
  
  if(any(eigen(scale)$value<0)){stop("Scale matrix is not positive definite")}		
  par <- c(x,location,scale,shape,extended)
  if(length(x)==2){
    res = .C("pmesn", as.double(c(-1,-1)), as.double(c(0,0)), as.double(par),
            out=double(1), double(1), NAOK=TRUE)$out
  }else if(length(x)==3){
    res = .C("pmesn3", as.double(c(-1,-1,-1)), as.double(c(0,0,0)), as.double(par),
            out=double(1), double(1), NAOK=TRUE)$out		
  }else{
  	stop("x must be of length 2 or 3")
  }	
  return(res)
}
# univariate extended skew-normal probability density function
desn <- function(x=0, location=0, scale=1, shape=0, extended=0)
{
  if(!is.atomic(x) || !is.numeric(x) || length(x)!=1)
    stop("x argument has incorrect format")
  if(!is.atomic(location) || length(location) != length(x))
    stop("Location argument has incorrect length")
  if(length(scale) != length(x))
    stop("Scale argument has incorrect dimensions")
  if(!is.atomic(shape) || length(shape) != length(x))
    stop("Shape argument has incorrect length")  
  if(length(extended)!=1)
    stop("Extended argument is not of length 1")
  
  .C("desn", as.double(x), as.double(location), as.double(scale), as.double(shape),
     as.double(extended), out = double(1), NAOK=TRUE)$out
}
# Bivariate and trivariate extended skew-normal probability density function
dmesn <- function(x=c(0,0), location=rep(0,length(x)), scale=diag(length(x)), shape=rep(0,length(x)), extended=0)
{
  if(!is.atomic(x) || !is.numeric(x))
    stop("x argument has incorrect format")
  if(!is.atomic(location) || length(location) != length(x))
    stop("Location argument has incorrect length")
  if(nrow(scale)!=ncol(scale) && nrow(scale) != length(x))
    stop("Scale argument has incorrect dimensions")
  if(!is.atomic(shape) || length(shape) != length(x))
    stop("Shape argument has incorrect length")  
  if(length(extended)!=1)
    stop("Extended argument is not of length 1")
    
  if(any(eigen(scale)$value<0)){stop("Scale matrix is not positive definite")}	
  if(length(x)==2){
  	res = .C("dmesn", as.double(x), as.double(location), as.double(scale), as.double(shape),
     as.double(extended), out=double(1), NAOK=TRUE)$out
  }else if(length(x)==3){
  	res = .C("dmesn3", as.double(x), as.double(location), as.double(scale), as.double(shape),
     as.double(extended), out=double(1), NAOK=TRUE)$out  	
  }else{
  	stop("x must be of length 2 or 3")	
  }		   
  return(res)
}
# univariate extended skew-t probability distribution function
pest <- function(x=0, location=0, scale=1, shape=0, extended=0, df=Inf)
{
  if(!is.atomic(x) || !is.numeric(x) || length(x)!=1)
    stop("x argument has incorrect format")
  if(!is.atomic(location) || length(location) != length(x))
    stop("Location argument has incorrect length")
  if(length(scale) != length(x))
    stop("Scale argument has incorrect dimensions")
  if(!is.atomic(shape) || length(shape) != length(x))
    stop("Shape argument has incorrect length")  
  if(length(extended)!=1)
    stop("Extended argument is not of length 1")
  
  if(df==Inf)
    res = pesn(x, location,scale,shape,extended)
  else{
    par <- c(x,location,scale,df,shape,extended)
    res = .C("pest", as.double(-1), as.double(0), as.double(par),
       out=double(1), double(1), NAOK=TRUE)$out
  }
  return(res)
}
# bivariate and trivariate extended skew-t probability distribution function
pmest <- function(x=c(0,0), location=rep(0,length(x)), scale=diag(length(x)), shape=rep(0,length(x)), extended=0, df=Inf)
{	
  if(!is.atomic(x) || !is.numeric(x))
    stop("x argument has incorrect format")
  if(!is.atomic(location) || length(location) != length(x))
    stop("Location argument has incorrect length")
  if(nrow(scale)!=ncol(scale) && nrow(scale) != length(x))
    stop("Scale argument has incorrect dimensions")
  if(!is.atomic(shape) || length(shape) != length(x))
    stop("Shape argument has incorrect length")  
  if(length(extended)!=1)
    stop("Extended argument is not of length 1")
  
  if(any(eigen(scale)$value<0)){stop("Scale matrix is not positive definite")}	
  if(df==Inf)
    res = pmesn(x, location,scale,shape,extended)
  else{
    par <- c(x,location,scale,df,shape,extended)
    if(length(x)==2){
	    res = .C("pmest", as.double(c(-1,-1)), as.double(c(0,0)), as.double(par),
       out=double(1), double(1), NAOK=TRUE)$out
	}else if(length(x)==3){
	    res = .C("pmest3", as.double(c(-1,-1,-1)), as.double(c(0,0,0)), as.double(par),
       out=double(1), double(1), NAOK=TRUE)$out		
	}else{
		stop("x must be of length 2 or 3")
	}	
  }
  return(res)
}
# univariate extended skew-t probability density function
dest <- function(x=0, location=0,scale=1,shape=0,extended=0,df=Inf)
{
  if(!is.atomic(x) || !is.numeric(x) || length(x)!=1)
    stop("x argument has incorrect format")
  if(!is.atomic(location) || length(location) != length(x))
    stop("Location argument has incorrect length")
  if(length(scale) != length(x))
    stop("Scale argument has incorrect dimensions")
  if(!is.atomic(shape) || length(shape) != length(x))
    stop("Shape argument has incorrect length")  
  if(length(extended)!=1)
    stop("Extended argument is not of length 1")
  
  if(df==Inf)
    res = desn(x, location,scale,shape,extended)
  else
    res = .C("dest", as.double(x), as.double(location), as.double(scale), as.double(df),
       as.double(shape), as.double(extended), out=double(1), NAOK=TRUE)$out
  return(res)
}
# multivariate extended skew-t probability density function
dmest <- function(x=c(0,0), location=rep(0,length(x)), scale=diag(length(x)), shape=rep(0,length(x)), extended=0, df=Inf)
{
  if(!is.atomic(x) || !is.numeric(x))
    stop("x argument has incorrect format")
  if(!is.atomic(location) || length(location) != length(x))
    stop("Location argument has incorrect length")
  if(nrow(scale)!=ncol(scale) && nrow(scale) != length(x))
    stop("Scale argument has incorrect dimensions")
  if(!is.atomic(shape) || length(shape) != length(x))
    stop("Shape argument has incorrect length")  
  if(length(extended)!=1)
    stop("Extended argument is not of length 1")
  
  if(any(eigen(scale)$value<0)){stop("Scale matrix is not positive definite")}	
  if(df==Inf)
    res = dmesn(x, location,scale,shape,extended)
  else{
  	if(length(x)==2){
    		res = .C("dmest", as.double(x), as.double(location), as.double(scale), as.double(df),
    		   as.double(shape), as.double(extended), out=double(1), NAOK=TRUE)$out
    }else if(length(x)==3){
    		res = .C("dmest3", as.double(x), as.double(location), as.double(scale), as.double(df),
   	 	   as.double(shape), as.double(extended), out=double(1), NAOK=TRUE)$out   	
    }else{
    		stop("x must be of length 2 or 3")
    }  		 
  }     
  return(res)
}
