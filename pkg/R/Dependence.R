index.ExtDep <- function(object, model, par, x, u){
  
  objects <- c("extremal", "pickands", "upper.tail", "lower.tail")
  if(!any(object ==  objects)){ stop("Wrong object specified")}
  
  if(object %in% objects[-4]){
    models <- c("HR", "ET", "EST")
    if(!any(model ==  models)){ stop("Wrong model specified")}
    d <- dim_ExtDep(model=model, par=par)
  }
  
  # Etxremal coefficient
  
  if(object == "extremal"){
    return(-log(pExtDep(q=rep(1,d), type="lower", method="Parametric", model=model, par=par)))
  }
  
  # Pickands dependence function
  
  if(object == "pickands"){
    
    if(missing(x)){stop("'x' should be provided")}
    if(!is.vector(x)){stop("'x' should be a vector")}
    if(round(sum(x),10) != 1){
      message("'x' should be defined on the unit simplex")
      return(NA)
    }
    if(d != length(x)){stop("wrong length of 'x'")}
    if(!any(d == c(2,3) )){stop("Pickands dependence function available in 2 or 3 dimensions only.")}
    
    dim.check <- dim_ExtDep(model=model, par=par, dim=d)
    if(!dim.check){stop('Wrong length of parameters')}
    
    if(any(x==0)){
      return(NA)
    }
    
    if(model == "HR"){
      return(-log(pExtDep(q=1/x, type="lower", method="Parametric", model="HR", par=par)))
    }
    
    if(model == "EST"){
      return(pk.extst(x=x, param=par))
    }
    
    if(model == "ET"){
      if(length(x)==2 && length(par)==2){
        par2 <- c(par[1], rep(0,2), par[2])
      }
      if(length(x)==3 && length(par)==7){
        par2 <- c(par[1:3], rep(0,3), par[4])
      }
      return(pk.extst(x=x, param=par2))
    }
    
  }
  
  # Coefficient of upper tail dependence
  
  if(object == "upper.tail"){
    if(d!=2){ stop("upper tail dependence coefficient only available in dimension 2")}
    return(2+log(pExtDep(q=rep(1,d), type="lower", method="Parametric", model=model, par=par)))
  }
  
  # Coefficient of lower tail dependence
  
  if(object == "lower.tail"){
    models <- c("ET","EST", "SN")
    if(!any(model ==  models)){ stop("Wrong model specified")}
    
    if(model == "EST"){
      
      if(length(par) !=4){ stop("Length of 'par' should be 4.")}
      return(chi.extst(corr=par[1], shape=par[2:3], df=par[4], tail="lower"))
    }

    if(model == "ET"){
      
      if(length(par) !=2){ stop("Length of 'par' should be 2.")}
      return(chi.extst(corr=par[1], shape=rep(0,2), df=par[2], tail="lower"))
    }
    
    if(model == "SN"){
      
      if(length(par) !=3){ stop("Length of 'par' should be 3.")}
      corr <- par[1]
      shape <- par[2:3]
      if(missing(u)){stop("'u' needs to be provided.")}
      
      return( chi.bsn(u=u, corr=corr, shape=shape, tail="lower") )
    }
    
  }  
  
}


###############################################################################
###############################################################################
## Hidden functions for tailDep function
###############################################################################
###############################################################################

chi.extst <- function(corr=0, shape=rep(0,2), df=1, tail="upper"){

	chistup <- function(scale=1, shape=rep(0,2), df=1){
	    .C("chistup", as.double(scale), as.double(df), as.double(shape), out=double(1), NAOK=TRUE)$out
	}

	chistlo <- function(scale=1, shape=rep(0,2), df=1){
	    .C("chistlo", as.double(scale), as.double(df), as.double(shape), out=double(1), NAOK=TRUE)$out
	}

	if(tail=="upper"){
		return(chistup(corr, shape, df))
	}else if(tail=="lower"){
		return(chistlo(corr, shape, df))
	}else{
	  stop("tail should be 'lower' or 'upper'.")
	}
}

chi.bsn <- function(u, corr=0, shape=rep(0,2), tail="upper"){

	chibsnup <- function(u, scale=diag(2), shape=rep(0,2)){
	    res <- double(1)
	    x <- c(qsn(u, omega=scale[1,1], alpha=shape[1]),
	           qsn(u, omega=scale[2,2], alpha=shape[2]))
	    pbsn <- pmesn(x=x, scale=scale, shape=shape)
	    res <- (2*log(1-u))/log(1-2*u+pbsn)-1
	    return(res)
	}
	
	chibsnlo <- function(u, scale=diag(2), shape=rep(0,2)){
	    res <- double(1)
	    x <- c(qsn(u, omega=scale[1,1], alpha=shape[1]),
	           qsn(u, omega=scale[2,2], alpha=shape[2]))
	    pbsn <- pmesn(x=x, scale=scale, shape=shape)
	    res <- (2*log(u))/log(pbsn)-1
	    return(res)
	}

	scale <- matrix(c(1,corr,corr,1),ncol=2)

	if(tail=="upper"){
		return(chibsnup(u, scale, shape))
	}else if(tail=="lower"){
		return(chibsnlo(u, scale, shape))
	}else{
	  stop("tail should be 'lower' or 'upper'.")
	}	

}

###############################################################################
###############################################################################
## Hidden functions for pickands function
###############################################################################
###############################################################################

pk.extst <- function(x, param=c(rep(0,choose(length(x),2)),rep(0,length(x)),1)){
	
	bivpkst <- function(x,scale, shape, df){
	    if(any(is.na(x)))
	      return(NA)
	    .C("bivpkst", as.double(x), as.double(scale), as.double(df), as.double(shape), out=double(1), NAOK=TRUE)$out
	}

	trivpkst <- function(x, scale, shape, df){
	    if(any(is.na(x)))
	      return(NA)
	    .C("trivpkst", as.double(x), as.double(scale), as.double(df), as.double(shape), out=double(1), NAOK=TRUE)$out
	}

	if(length(x)==2 && length(param)==4){
		return(bivpkst(x,scale=param[1], shape=param[2:3], df=param[4]))
	}
	if(length(x)==3 && length(param)==7){
		Sigma <- diag(3)
		Sigma[lower.tri(Sigma)] = param[1:3]
		Sigma[upper.tri(Sigma)] = param[1:3]
		return(trivpkst(x,scale=Sigma, shape=param[4:5], df=param[7]))
	}
}


