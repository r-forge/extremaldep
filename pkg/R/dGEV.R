
dGEV <- function(x, loc, scale, shape, log=FALSE)
{
  
  if(is.vector(x) || is.vector(loc) || is.vector(scale) || is.vector(shape)){
    
    n <- max(c(length(x), length(loc), length(scale), length(shape) ))
    vals <- cbind(rep(x,length=n), rep(loc,length=n), rep(scale,length=n), rep(shape,length=n))
    result <- apply(vals,1, function(y) .C('dGEV', as.double(y[1]), as.double(y[2]), as.double(y[3]), as.double(y[4]),
                                           res=double(1), DUP=TRUE, NAOK=TRUE)$res )
    
  }else{
    result <- .C('dGEV', as.double(x), as.double(loc), as.double(scale), as.double(shape),
                 res=double(1), DUP=TRUE, NAOK=TRUE)$res  
  }
  
  if(log){
    return(log(result))
  }else{
    return(result)    
  }
}


pGEV <- function(q, loc, scale, shape, lower.tail=TRUE)
{
  
  if(is.vector(q) || is.vector(loc) || is.vector(scale) || is.vector(shape)){
    
    n <- max(c(length(q), length(loc), length(scale), length(shape) ))
    vals <- cbind(rep(q,length=n), rep(loc,length=n), rep(scale,length=n), rep(shape,length=n))
    result <- apply(vals,1, function(y) .C('pGEV', as.double(y[1]), as.double(y[2]), as.double(y[3]), as.double(y[4]),
                                           res=double(1), DUP=TRUE, NAOK=TRUE)$res )
    
  }else{
    result <- .C('pGEV', as.double(q), as.double(loc), as.double(scale), as.double(shape),
                 res=double(1), DUP=TRUE, NAOK=TRUE)$res  
  }
  
  if(lower.tail){
    return(result)
  }else{
    return(1-result)    
  }
}

qGEV <- function(p, loc, scale, shape)
{
  
  if(is.vector(p) || is.vector(loc) || is.vector(scale) || is.vector(shape)){
    
    n <- max(c(length(p), length(loc), length(scale), length(shape) ))
    vals <- cbind(rep(p,length=n), rep(loc,length=n), rep(scale,length=n), rep(shape,length=n))
    result <- apply(vals,1, function(y) .C('qGEV', as.double(y[1]), as.double(y[2]), as.double(y[3]), as.double(y[4]),
                                           res=double(1), DUP=TRUE, NAOK=TRUE)$res )
    
  }else{
    result <- .C('qGEV', as.double(p), as.double(loc), as.double(scale), as.double(shape),
                 res=double(1), DUP=TRUE, NAOK=TRUE)$res  
  }
  
  return(result)
}