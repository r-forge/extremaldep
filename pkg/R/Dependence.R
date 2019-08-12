chi.extst <- function(corr=0, shape=rep(0,2), df=1, tail="upper"){

	chistup <- function(scale=1, shape=rep(0,2), df=1){
	    .C("chistup", as.double(scale), as.double(df), as.double(shape), out=double(1), NAOK=TRUE)$out
	}

	chistlo <- function(scale=1, shape=rep(0,2), df=1){
	    .C("chistlo", as.double(scale), as.double(df), as.double(shape), out=double(1), NAOK=TRUE)$out
	}

	if(tail=="upper"){
		return(chistup(corr, shape, df))
	}else{
		return(chistlo(corr, shape, df))
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
	}else{
		return(chibsnlo(u, scale, shape))
	}	

}
#
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

#