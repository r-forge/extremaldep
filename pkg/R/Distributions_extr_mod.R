################################################################################################
### Authors: Boris Beranger and Simone Padoan        	 									 ###
### 																							 ###	
###	Emails: borisberanger@gmail.com, simone.padoan@unibocconi.it								 ###
### 																							 ###
###	Institutions: Department of Decision Sciences, University Bocconi of Milan				 ###
### 			  	School of Mathematics and Statistics, University of New South Wales 			 ###
### 																							 ###
### File name: dens_extremal_models.r                 						             	 ###
### 																							 ###
### Description:                                  							      		     ###
### This file provides:        									 							 ###
### 		1- 	The density of the bivariate	 Extremal-t and Extremal skew-t models.				 ###
### 		2-	The exponent function for the Husler-Reiss model,					 			 ###
### 			the Extremal-t model and Extremal Skew-t in 2,3 dimensions					     ###
### 		3-  The probability of excess for the Husler-Reiss, Extremal-t 						 ###
### 			and Extremal Skew-t models													     ###
### 		4-  Pairwise composite likelihood estimation procedure for up to quadrivarite		 ###
### 			data, for the Husler-Reiss, Extremal-t and Extremal Skew-t models				 ###
### 	5-	Random generation for the Extremal-t and Extremal-Skew-t models					 ###
### 																							 ###
### Last change: 06/09/2016                           		  								 ###
### 																							 ###
################################################################################################


dens_extr_mod <- function(model, z, param, log=TRUE){
	
	dmextst <- function(x, scale=1, shape=rep(0,length(x)), df=1){
		if(any(is.na(x)))
    		return(NA)
    	res <- .C("dmextst", as.double(x), as.double(scale), as.double(df),
       		as.double(shape), double(1), PACKAGE="ExtremalDep", NAOK=TRUE)
    	return(res[[5]])
	}
	
	if(model=="Extremalt"){
		if(length(param) !=2){stop("Wrong length of parameters")}
		if(log){
			return( log(dmextst(x=z, scale=param[1], df=param[2])) )
		}else{
			return( dmextst(x=z, scale=param[1], df=param[2]) )
		}	
	}

	if(model=="Skewt"){
		if(length(param) !=4){stop("Wrong length of parameters")}
		if(log){
			return( log(dmextst(x=z, scale=param[1], shape=param[2:3], df=param[4])) )
		}else{
			return(dmextst(x=z, scale=param[1], shape=param[2:3], df=param[4]))
		}	
	}

}

exponent_extr_mod <- function(model, z, param, dist){
	
	
	## Extremal-t model
	
	exponent.triv.ext <- function(z,rho,mu){
		z1 <- z[1]; z2<-z[2]; z3<-z[3];	
		rho12 <- rho[1]; rho13 <- rho[2]; rho23 <- rho[3];
		if( (any(rho<=-1) || any(rho >=1) || (mu <= 0))==1){return(1e-50)}

		R1 <- (rho23-rho12*rho13)/sqrt((1-rho12^2)*(1-rho13^2))
		R2 <- (rho13-rho12*rho23)/sqrt((1-rho12^2)*(1-rho23^2))
		R3 <- (rho12-rho13*rho23)/sqrt((1-rho13^2)*(1-rho23^2))
		if(sum(eigen(matrix(c(1,R1,R1,1),ncol=2))$values<0)>=1){return(1e-50)}
		if(sum(eigen(matrix(c(1,R2,R2,1),ncol=2))$values<0)>=1){return(1e-50)}
		if(sum(eigen(matrix(c(1,R3,R3,1),ncol=2))$values<0)>=1){return(1e-50)}
	
		x1 <- sqrt((mu+1)/(1-rho12^2))*((z2/z1)^(1/mu)-rho12)
		y1 <- sqrt((mu+1)/(1-rho13^2))*((z3/z1)^(1/mu)-rho13)
		x2 <- sqrt((mu+1)/(1-rho12^2))*((z1/z2)^(1/mu)-rho12)
		y2 <- sqrt((mu+1)/(1-rho23^2))*((z3/z2)^(1/mu)-rho23)
		x3 <- sqrt((mu+1)/(1-rho13^2))*((z1/z3)^(1/mu)-rho13) 
		y3 <- sqrt((mu+1)/(1-rho23^2))*((z2/z3)^(1/mu)-rho23)	

		I1 <- pmest(x=c(x1,y1), scale=matrix(c(1,R1,R1,1),ncol=2), df=mu+1)
		I2 <- pmest(x=c(x2,y2), scale=matrix(c(1,R2,R2,1),ncol=2), df=mu+1)
		I3 <- pmest(x=c(x3,y3), scale=matrix(c(1,R3,R3,1),ncol=2), df=mu+1)

		return( I1/z1 + I2/z2 + I3/z3 )
	}

	## Extremal Skew-t
	
	pmextst <- function(x, scale=1, shape=rep(0,length(x)), df=1){
    	if(any(is.na(x)))
    	  return(NA)
    	res <- .C("pmextst", as.double(x), as.double(scale), as.double(df),
    	   as.double(shape), double(1), PACKAGE="ExtremalDep", NAOK=TRUE)
    	return(res[[5]])
	}	

	exponent.triv.skewt <- function(x, rho, alpha, mu){
		z1 <- z[1]; z2<-z[2]; z3<-z[3];				
		p <- round(uniroot(function(x){length(rho)-choose(x,2)},lower=1, upper=10)$root) 
		Sig <- diag(p)
		Sig[lower.tri(Sig)] <- rho
		Sig[upper.tri(Sig)] <- rho

		Sigma_1 <- Sig[-1,-1] - Sig[-1,1] %*% t(Sig[1,-1])
		Sigma_2 <- Sig[-2,-2] - Sig[-2,2] %*% t(Sig[2,-2])
		Sigma_3 <- Sig[-3,-3] - Sig[-3,3] %*% t(Sig[3,-3])	
	
		s_1 <- s_2 <- s_3 <- diag(p-1)
		diag(s_1) <- sqrt(diag(Sigma_1))
		diag(s_2) <- sqrt(diag(Sigma_2))
		diag(s_3) <- sqrt(diag(Sigma_3))
	
		Sigma_bar_1 <- solve(s_1) %*% Sigma_1 %*% solve(s_1)
		Sigma_bar_2 <- solve(s_2) %*% Sigma_2 %*% solve(s_2)
		Sigma_bar_3 <- solve(s_3) %*% Sigma_3 %*% solve(s_3)
			
		alpha_tilde_1 <- t(alpha[-1])
		alpha_tilde_2 <- t(alpha[-2])
		alpha_tilde_3 <- t(alpha[-3])
	
		num1 <- alpha[1] + Sig[1,-1] %*% t(alpha_tilde_1)
		denom1 <- sqrt( 1 + alpha_tilde_1 %*% Sigma_1 %*% t(alpha_tilde_1)  )
		alpha_bar_1 <- num1/denom1

		num2 <- alpha[2] + Sig[2,-2] %*% t(alpha_tilde_2)
		denom2 <- sqrt( 1 + alpha_tilde_2 %*% Sigma_2 %*% t(alpha_tilde_2)  )
		alpha_bar_2 <- num2/denom2
	
		num3 <- alpha[3] + Sig[3,-3] %*% t(alpha_tilde_3)
		denom3 <- sqrt( 1 + alpha_tilde_3 %*% Sigma_3 %*% t(alpha_tilde_3)  )
		alpha_bar_3 <- num3/denom3

		nu_1 <- pt(sqrt(mu+1) * alpha_bar_1,df=mu+1) * 2^(mu/2-1) * gamma((mu+1)/2) / sqrt(pi)
		nu_2 <- pt(sqrt(mu+1) * alpha_bar_2,df=mu+1) * 2^(mu/2-1) * gamma((mu+1)/2) / sqrt(pi)
		nu_3 <- pt(sqrt(mu+1) * alpha_bar_3,df=mu+1) * 2^(mu/2-1) * gamma((mu+1)/2) / sqrt(pi)	
	
		x1_bar <- z1 *nu_1
		x2_bar <- z2 *nu_2
		x3_bar <- z3 *nu_3	
	
		comp1_1 <- sqrt((mu+1)/(1-rho[1]^2))*((x2_bar/x1_bar)^(1/mu)-rho[1])
		comp2_1 <- sqrt((mu+1)/(1-rho[2]^2))*((x3_bar/x1_bar)^(1/mu)-rho[2])
		comp1_2 <- sqrt((mu+1)/(1-rho[1]^2))*((x1_bar/x2_bar)^(1/mu)-rho[1])
		comp2_2 <- sqrt((mu+1)/(1-rho[3]^2))*((x3_bar/x2_bar)^(1/mu)-rho[3])
		comp1_3 <- sqrt((mu+1)/(1-rho[2]^2))*((x1_bar/x3_bar)^(1/mu)-rho[2])
		comp2_3 <- sqrt((mu+1)/(1-rho[3]^2))*((x2_bar/x3_bar)^(1/mu)-rho[3])	

		alpha_star_1 <- alpha_tilde_1 %*% s_1
		alpha_star_2 <- alpha_tilde_2 %*% s_2
		alpha_star_3 <- alpha_tilde_3 %*% s_3	
	
		tau_star_1 <- sqrt(mu+1) * (alpha[1] + Sig[-1,1] %*% t(alpha_tilde_1) )
		tau_star_2 <- sqrt(mu+1) * (alpha[2] + Sig[-2,2] %*% t(alpha_tilde_2) )
		tau_star_3 <- sqrt(mu+1) * (alpha[3] + Sig[-3,3] %*% t(alpha_tilde_3) )	
				
		I1 <- pmest(x=c(comp1_1,comp2_1), scale=Sigma_bar_1, shape=alpha_star_1, ext=tau_star_1, df=mu+1)
		I2 <- pmest(x=c(comp1_2,comp2_2), scale=Sigma_bar_2, shape=alpha_star_2, ext=tau_star_2, df=mu+1)
		I3 <- pmest(x=c(comp1_3,comp2_3), scale=Sigma_bar_3, shape=alpha_star_3, ext=tau_star_3, df=mu+1)
			
		return( I1/z1 + I2/z2 + I3/z3 )
	}
	
	## Husler-Reiss model
	
	exponent.hr <- function(X, par){

		expo.hr.2d <- function(X,par){
			x=X[1]; y=X[2];
			lambda12 <- par;
			pt1 <- 1/x * pnorm(lambda12+log(y/x)/(2*lambda12))
			pt2 <- 1/y * pnorm(lambda12+log(x/y)/(2*lambda12))
			return( pt1 + pt2 )
		}

		expo.hr.3d <- function(X,par){
			x=X[1]; y=X[2]; z=X[3];
			lambda12 <- par[1]; lambda13 <- par[2]; lambda23 <- par[3];
			S1 <- matrix(c(4*lambda12^2, rep(2*(lambda12^2 + lambda13^2 - lambda23^2),2), 4*lambda13^2),nrow=2)
			S2 <- matrix(c(4*lambda12^2, rep(2*(lambda12^2 + lambda23^2 - lambda13^2),2), 4*lambda23^2),nrow=2)
			S3 <- matrix(c(4*lambda13^2, rep(2*(lambda13^2 + lambda23^2 - lambda12^2),2), 4*lambda23^2),nrow=2)
			pt1 <- 1/x * pmesn(x=c(log(y/x)+2*lambda12^2, log(z/x)+2*lambda13^2),scale = S1 )
			pt2 <- 1/y * pmesn(x=c(log(x/y)+2*lambda12^2, log(z/y)+2*lambda23^2),scale = S2 )
			pt3 <- 1/z * pmesn(x=c(log(x/z)+2*lambda13^2, log(y/z)+2*lambda23^2),scale = S3 )
			return(pt1 + pt2 + pt3 )
		}

		if(any(par<=0)){break('wrong value of parameters')}	
		if(length(X)==2 && length(par)==1){return( expo.hr.2d(X,par) )} # in 2d there is only 1 parameter
		if(length(X)==3 && length(par)==3){return( expo.hr.3d(X,par) )} # in 3d there are 3 parameters
	}
	
	res <- 0

	if(model=="Extremalt"){
		if(length(z)==2){
			if(length(param) !=2){stop("Wrong length of parameters")}
			res <- -log(pmextst(x=z, scale=param[1], df=param[2]))
		}
		if(length(z)==3){
			if(length(param) !=4){stop("Wrong length of parameters")}
			res <- exponent.triv.ext(z, param[1:3],param[4])
		}
	}	
	
	if(model=="Skewt"){
		if(length(z)==2){
			if(length(param) !=4){stop("Wrong length of parameters")}
			res <- -log(pmextst(x=z, scale=param[1], shape=param[2:3], df=param[4]))
		}
		if(length(z)==3){
			if(length(param) !=7){stop("Wrong length of parameters")}
			res <- exponent.triv.skewt(z, param[1:3],param[4:6],param[7])
		}
	}	

	if(model=="hr"){
		res <- exponent.hr(X=z, par=param)
	}	

	if(dist==TRUE){return(exp(-res))}else{
		return(res)
	}		
}

excess_pr_extr_mod <- function(model, z, param ){
	
	if(length(z)==1){
		return(1-exp(-1/z))
	}
	if(length(z)==2){
		return( 1 - exp(-1/z[1]) - exp(-1/z[2]) + exponent_extr_mod(model,z,param,dist=TRUE) )
	}
	if(length(z)==3){

		if(model=="Extremalt" && length(param)==4){
			return( 1 - exp(-1/z[1]) - exp(-1/z[2]) - exp(-1/z[3]) + exponent_extr_mod(model,z[-3],param[c(1,4)],dist=TRUE) + exponent_extr_mod(model,z[-2],param[c(2,4)],dist=TRUE) + exponent_extr_mod(model,z[-1],param[c(3,4)],dist=TRUE) - exponent_extr_mod(model,z,param,dist=TRUE))
		}		
		if(model=="Skewt" && length(param)==7){
			return( 1 - exp(-1/z[1]) - exp(-1/z[2]) - exp(-1/z[3]) + exponent_extr_mod(model,z[1:2],param[c(1,4,5,7)],dist=TRUE) + exponent_extr_mod(model,z[c(1,3)],param[c(2,4,6,7)],dist=TRUE) + exponent_extr_mod(model,z[2:3],param[c(3,5:7)],dist=TRUE) - exponent_extr_mod(model,z,param,dist=TRUE))
		}
		if(model=="hr" && length(param)==3){
			return( 1 - exp(-1/z[1]) - exp(-1/z[2]) - exp(-1/z[3]) + exponent_extr_mod(model,z[1:2],param[1],dist=TRUE) + exponent_extr_mod(model,z[c(1,3)],param[2],dist=TRUE) + exponent_extr_mod(model,z[2:3],param[3],dist=TRUE) - exponent_extr_mod(model,z,param,dist=TRUE))
		}		
	}	
	
}

fit_pclik_extr_mod <- function(model, data, parastart, trace){

	## Husler-Reiss model

	llHR <- function(x, lambda=1){
    	if(is.matrix(x)){
    		n <- nrow(x)
    	}else{ n <- 1}
    res <- .C("llHRmax", as.double(x), as.double(lambda), as.integer(n), double(1), 
    		PACKAGE="ExtremalDep", NAOK=TRUE)
    	return(res[[4]])
	}

	pair.cllik.hr <- function(z,par){
		pllik <- 0
		if( length(z)==2 && length(par)==1){
			pllik <- llHR(z, lambda=par)	
		}
		if( length(z)==3 && length(par)==3 ){
			pllik <-  llHR(z[-3], lambda=par[1]) + llHR(z[-2], lambda=par[2]) + llHR(z[-1], lambda=par[3])
		}
		if( length(z)==4 && length(par)==6 ){
			pllik <- llHR(z[c(1,2)], lambda=par[1]) + llHR(z[c(1,3)], lambda=par[2]) + llHR(z[c(1,4)], lambda=par[3]) + llHR(z[c(2,3)], lambda=par[4]) + llHR(z[c(2,4)], lambda=par[5]) + llHR(z[c(3,4)], lambda=par[6])	
		}
		return(pllik)
	}


	## Extremal-t model

	llET <- function(x, scale=1, df=1){
    	if(is.matrix(x)){
    		n <- nrow(x)
    	}else{ n <- 1}
    	par <- c(df,scale)
    	res <- .C("llETmax", as.double(x), as.double(par), as.integer(n), double(1), 
    		PACKAGE="ExtremalDep", NAOK=TRUE)
    	return(res[[4]])
	}

	pair.cllik.ext <- function(z,par){
		pllik <- 0
		if( length(z)==2 && length(par)==2){
			pllik <- llET(z, scale=par[1], df=par[2])	
		}
		if( length(z)==3 && length(par)==4 ){
			pllik <-  llET(z[-3], scale=par[1], df=par[4]) + llET(z[-2], scale=par[2], df=par[4]) + llET(z[-1], scale=par[3], df=par[4])
		}
		if( length(z)==4 && length(par)==7 ){
			pllik <- llET(z[c(1,2)], scale=par[1], df=par[7]) + llET(z[c(1,3)], scale=par[2], df=par[7]) + llET(z[c(1,4)], scale=par[3], df=par[7]) + llET(z[c(2,3)], scale=par[4], df=par[7]) + llET(z[c(2,4)], scale=par[5], df=par[7]) + llET(z[c(3,4)], scale=par[6], df=par[7])	
		}
		return(pllik)
	}

	## Extremal Skew-t model

	llextst <- function(x, scale=1, shape=rep(0,length(x)), df=1){
    	if(is.matrix(x)){
    		n <- nrow(x)
    	}else{ n <- 1}
    	
    	res <- .C("llextst", as.double(x), as.integer(n), as.double(scale), as.double(df),
    	   as.double(shape), double(1), PACKAGE="ExtremalDep", NAOK=TRUE)
    	return(res[[6]])
	}

	pair.cllik.skewt <- function(z,par){
		pllik <- 0
		if( length(z)==2 && length(par)==4){
			pllik <- llextst(x=z, scale=par[1], shape=par[2:3], df=par[4])
		}
		if( length(z)==3 && length(par)==7 ){
			pllik <- llextst(x=z[-3], scale=par[1], shape=par[4:5], df=par[7]) + llextst(x=z[-2], scale=par[2], shape=par[c(4,6)], df=par[7]) + llextst(x=z[-1], scale=par[3], shape=par[5:6], df=par[7]) 
		}
		if( length(z)==4 && length(par)==11 ){
			pllik <- llextst(x=z[c(1,2)], scale=par[1], shape=par[c(7,8)], df=par[11]) + llextst(x=z[c(1,3)], scale=par[2], shape=par[c(7,9)], df=par[11]) + llextst(x=z[c(1,4)], scale=par[3], shape=par[c(7,10)], df=par[11]) + llextst(x=z[c(2,3)], scale=par[4], shape=par[c(8,9)], df=par[11]) + llextst(x=z[c(2,4)], scale=par[5], shape=par[c(8,10)], df=par[11]) + llextst(x=z[c(3,4)], scale=par[6], shape=par[c(9,10)], df=par[11])	
		}
		return(pllik)
	}

	biv.cllik <- function(par){ # model and data are extrenal parameters
		if(model=="hr"){
			return( sum( apply(data,1,function(x) { pair.cllik.hr(x,par) } ) ) )
		}
		if(model=="Extremalt"){
			return( sum( apply(data,1,function(x) { pair.cllik.ext(x,par) } ) ) )
		}
		if(model=="Skewt"){
			return( sum( apply(data,1,function(x) { pair.cllik.skewt(x,par) } ) ) )
		}		
	}
	
	est=optim(parastart, biv.cllik, method="Nelder-Mead",control = list(maxit = 1e8, trace = trace, fnscale=-1),hessian=FALSE) 	
	return(list(par=est$par,LL=est$value))
}


r_extr_mod <- function(model, n, param){

	## Extremal-t model
	
	rextremalt <- function(n, ndim, scale=diag(ndim), df=1, nblock=500){
		res <- matrix(double(n*ndim),ncol=ndim,nrow=n)
		cst <- gamma((df+1)/2) / gamma(df/2) * (df*pi)^(-1/2) * df^((df-1)/2);
		an <- (nblock*cst)^(1/df)
 
		for(i in 1:n){
			sim <- rmst(n=nblock, Omega=scale, alpha=rep(0,ndim), nu=df) 
			maxima <- apply(sim, 2, max)
    		res[i,] <- maxima/an
		}
		return(res^df) # Transform from Frechet with parameter df to Standard Frechet
	}

	## Extremal Skew-t model
	
	rmextst <- function(n, ndim, scale=diag(ndim), shape=rep(0,ndim), df=1, nblock=500){
		res <- matrix(double(n*ndim),ncol=ndim,nrow=n)
		df1 <- df+1
		df2 <- df+2
		an <- (gamma(0.5*df1)*df^(df/2)*pt(shape*sqrt(df1),df))/(gamma(0.5*df2)*sqrt(pi))
		an <- (nblock*an)^(1/df)
		for(i in 1:n){
    		sim <- rmst(n=nblock, Omega=scale, alpha=shape, nu=df)
    		maxima <- apply(sim, 2, max)
    		res[i,] <- maxima/an
    		}
  		return(res^df) # Transform from Frechet with parameter df to Standard Frechet
	}

	dim <- nrow(scale)
	num <- 5e+5
	
	if(model=="Extremalt"){
		dim <- round(uniroot(function(x) choose(x,2)+1-length(param), interval = c(1,10))$root)	
		Sigma <- diag(dim)
		Sigma[lower.tri(Sigma)] = param[1:choose(dim,2)]
		Sigma[upper.tri(Sigma)] = param[1:choose(dim,2)]			
		return(rextremalt(n=n, ndim=dim, scale=Sigma, df=param[-(1:choose(dim,2))], nblock=num))
	}
	if(model=="Skewt"){
		dim <- round(uniroot(function(x) choose(x,2)+x+1-length(param), interval = c(1,10))$root)	
		Sigma <- diag(dim)
		Sigma[lower.tri(Sigma)] = param[1:choose(dim,2)]
		Sigma[upper.tri(Sigma)] = param[1:choose(dim,2)]		
		return(rmextst(n=n, ndim=dim, scale=Sigma, shape=param[choose(dim,2)+(1:dim)], df=param[-(1:(choose(dim,2)+dim))], nblock=num))
	}

}




