################################################################################################
### Authors: Boris Beranger and Simone Padoan        	 		                    							 ###
### 																							                                           ###	
###	Emails: borisberanger@gmail.com, simone.padoan@unibocconi.it				            				 ###
### 																						                                          	 ###
###	Institutions: Department of Decision Sciences, University Bocconi of Milan		      		 ###
### School of Mathematics and Statistics, University of New South Wales 			        			 ###
### 																							                                           ###
### File name: Densities.r	                 							                            	     ###
### 																						                                          	 ###
### Description:                                  		              					      		     ###
### This file provides the Angular densities of extremal dependence models: 					       ###
### 1) Pairwise Beta, Dirichlet and Husler-Reiss (mass only on the interior of the simplex)  ###
### 2) Asymmetric logistic and Extremal-t (mass on all subsets of the simplex)				       ###
### It also provides the likelihood and log-likelihood functions								             ###
### 																							                                           ###
### Last change: 06/08/2019                         		  									                 ###
### 																					                                          		 ###
################################################################################################


dens <- function(x=rbind(c(0.1,0.3,0.6),c(0.1,0.2,0.7)), model="Pairwise", par=c(2,2,2,4), c=NULL, log=FALSE, vectorial=TRUE){

### Functions needed for the Pairwise Beta model 

# Model defined for a minimum of three dimensions 

dens_pb <- function (x, b, alpha, log, vectorial){
	if(any(b<=0) || b<=0){return(1e-300)}
	hij<-function(wi,wj,bij,alpha,p){
		wij=wi+wj;
		return(wij^(2*alpha-1)*(1-wij)^(alpha*(p-2)-p+2)*gamma(2*bij)/gamma(bij)^2*(wi/wij)^(bij-1)*(wj/wij)^(bij-1))
	}

	dens_pairb<-function(w,b,alpha){
		p<-length(w);
		dens=0;	
		K= 2*factorial(p-3)*gamma(alpha*p+1)/(p*(p-1)*gamma(2*alpha+1)*gamma(alpha*(p-2)))
		for(i in 1:(p-1)){
			for(j in (i+1):p){
				d= hij(w[i],w[j],b[(i-1)*p-sum(1:i)+j],alpha,p) 
				dens=dens+d
			}
		}
		dens=dens*K
		return(dens)
	}
	
	xvect = as.double(as.vector(t(x)))
    if (is.vector(x)) {
        dim = as.integer(length(x))
        n = as.integer(1)
        if(round(sum(x),7) != 1){ stop("Data is not angular")}
        result <- dens_pairb(x,b,alpha)
    }
    else {
        dim = as.integer(ncol(x))
        n = as.integer(nrow(x))
        if (sum(apply(x,1,sum)) != n){ stop("Data is not angular")}

	    if (vectorial) {
	        result = double(n)
	        result <- apply(x,1,function(y){dens_pairb(y,b,alpha)}) 
	    } else { # vectorial=FALSE mean we return the likelihood function
    	    result = as.double(1)
       		result <- prod(apply(x,1,function(y){dens_pairb(y,b,alpha)})) 
    	}
    }	
    if(log)
    	return(log(result))
    else return(result)
}

### Functions needed for the Husler-Reiss model

dens_hr <- function (x, lambda, log, vectorial){
	
	dens_husler<-function(w,lambda){
		p<-length(w);
		k=p-1
		w.tilde<-rep(0,k);
 		Sigma<-diag(k);
  
		for(i in 1:k){
 	 		if(w[1]==0 | w[i+1]==0){return(1e-50)} else{
 	    		w.tilde[i] <- log(w[i+1]/w[1])+2*lambda[i]^2;
      
  	    		for(j in i:k){
  	    			if(i==j){Sigma[i,j]=Sigma[j,i]=2*(lambda[i]^2+lambda[j]^2)}else{
  	        		Sigma[i,j]=Sigma[j,i]=2*(lambda[i]^2+lambda[j]^2-lambda[sum(k:(k-i+1))+j-i]^2)
 	       			}
 	     		}
 	   		}
		}
  
		if(sum(eigen(Sigma)$values<0)>=1){return(1e-50)} # Check if matrix is positive definite
		if(any(is.na(w.tilde))){return (1e-50)}

  		part1 = w[1]^2*prod(w[2:p])*(2*pi)^(k/2)* abs(det(Sigma))^(1/2)
 	 	part2 = exp(-0.5*t(w.tilde) %*% solve(Sigma) %*% t(t(w.tilde)))
   
  		return(part2/part1) 
		}
	
	xvect = as.double(as.vector(t(x)))
	if (is.vector(x)) {
    	dim = as.integer(length(x))
    	n = as.integer(1)
    	if(round(sum(x),7) != 1){ stop("Data is not angular")}
    	if(log){
      		result = log(dens_husler(x,lambda)) 
    	}else{
      		result = dens_husler(x,lambda)
    	}	
	}
	else {
    	dim = as.integer(ncol(x))
    	n = as.integer(nrow(x))
    	if (sum(apply(x,1,sum)) != n){ stop("Data is not angular")}
    
    if (vectorial) {
    	result = double(n)
    	if (log){
        	result <- apply(x,1,function(y){log(dens_husler(y,lambda))}) 
    	}else{
        	result <- apply(x,1,function(y){dens_husler(y,lambda)}) 
    	}
    } else { # vectorial=FALSE mean we return the likelihood function
    	result = as.double(1)
    	if (log){
       	 result <- sum(apply(x,1,function(y){log(dens_husler(y,lambda))})) 
    	}else{
       	 result <- prod(apply(x,1,function(y){dens_husler(y,lambda)})) 
    	}   
    }
  	}  
  	return(result)
}

### Functions needed for the Dirichlet model

dens_di <- function (x, para, log, vectorial){
	
	dens_diri <- function(w,para){
		d <- length(para)
	
		if(any(para <=0)){return(1e-50)}
		if(length(para) != d){stop("Wrong length of parameter")}
	
		part1 <- prod(para/gamma(para))
		part2 <- gamma(sum(para)+1)/(sum(para*w)^(d+1))
		part3 <- prod((para*w/sum(para*w))^(para-1))
		return(part1*part2*part3)
	}
	
	xvect = as.double(as.vector(t(x)))
    if (is.vector(x)) {
        dim = as.integer(length(x))
        n = as.integer(1)
        if(round(sum(x),7) != 1){ stop("Data is not angular")}
        result = dens_diri(x,para)
    }
    else {
        dim = as.integer(ncol(x))
        n = as.integer(nrow(x))
        if (sum(apply(x,1,sum)) != n){ stop("Data is not angular")}

	    if (vectorial) {
	        result = double(n)
	        result <- apply(x,1,function(y){dens_diri(y,para)}) 
	    } else { # vectorial=FALSE mean we return the likelihood function
    	    result = as.double(1)
       		result <- prod(apply(x,1,function(y){dens_diri(y,para)})) 
    	}
    }
    if(log)
    	return(log(result))
    else return(result)
}

### Functions needed for the Asymmetric logistic model

dens_al <- function (x, alpha, beta, c, log, vectorial){

	# The 2d case: 

	# in the following functions:
	#
	# alpha is a vector of size 1: for the subset {1,2}
	# beta is a vector of size 2: for [1,{1,2}] and [2,{1,2}]
	# beta for [1,{1}], [2,{2}] and [3,{3}] are omitted as obtained as [1,{1}] = 1 - [1,{1,2}] and [2,{2}] = 1 - [2,{1,2}] 

	interior_alm_2d <- function(w,alpha,beta){
		part1 <- (alpha-1) * (beta[1]*beta[2])^alpha * (w[1]*w[2])^(alpha-2)
		part2 <- ( (beta[1]*w[2])^alpha + (beta[2]*w[2])^alpha )^(1/alpha-2)
	
		return(part1*part2)	
	}	

	corners_alm_2d <- function(w,alpha,beta,s){ # mass on the corner of the s-th component
		if(s==1){return(1-beta[1])}
		if(s==2){return(1-beta[2])}
	}

	dens_alm_2d <- function(w,alpha,beta,c){
		if(length(alpha)!=1){return(stop("Wrong length of parameter alpha"))}
		if(length(beta)!=2){return(stop("Wrong length of parameter beta"))}
		if( (alpha < 1) || any(beta < 0) || any(beta > 1)){return(1e-50)}
  
		if(sum(w > (1-c) ) == 1){ # then we are in a corner 
			ind <- which(w > (1-c))
			return(corners_alm_2d(w,alpha,beta,ind))
		} else {
			return(interior_alm_2d(w,alpha,beta))	
		}
	}

	# The 3d case:

	# in the following functions:
	#
	# alpha is a vector of size 4: for the subsets {1,2}, {1,3}, {2,3} and {1,2,3}
	# beta is a vector of size 9: for [1,{1,2}], [2,{1,2}], [1,{1,3}], [3,{1,3}], [2,{2,3}], [3,{2,3}], [1,{1,2,3}], [2,{1,2,3}] and [3,{1,2,3}]
	# beta for [1,{1}], [2,{2}] and [3,{3}] are omitted as obtained as [1,{1}] = 1 - [1,{1,2}]+[1,{1,3}]+[1,{1,2,3}] etc... 

	interior_alm_3d <- function(w,alpha,beta){
	
		k <- c(1,2,3)
	
		part1 <- (alpha[4]-1) * (2*alpha[4]-1)
		part2 <- prod( beta[6+k]^alpha[4] * w[k]^(-alpha[4]-1) )
		part3 <- sum( (beta[6+k]/w[k])^alpha[4] )^(1/alpha[4]-3)
	
		return(part1*part2*part3)
	}

	edges_alm_3d <- function(w,alpha,beta,s,t){ # mass on the edge linking s-th and t-th components (in increasing order)
	
		if(t<s){return('t cannot be less than s')}
	
		if(s==1 && t==2){a=alpha[1];b1=beta[1];b2=beta[2]} ## s=1, t=2 or s=2, t=1 
		if(s==1 && t==3){a=alpha[2];b1=beta[3];b2=beta[4]} ## s=1, t=3 or s=3, t=1 
		if(s==2 && t==3){a=alpha[3];b1=beta[5];b2=beta[6]} ## s=3, t=2 or s=2, t=3 			
		w1=w[s];w2=w[t];
	
		part1 <- (a-1) * (b1*b2)^a * (w1*w2)^(-a-1)
		part2 <- ( (b1/w1)^a + (b2/w2)^a )^(1/a-2)
	
		return(part1*part2)
	}	

	corners_alm_3d <- function(w,alpha,beta,s){ # mass on the corner of the s-th component
		if(s==1){return(1-beta[1]-beta[2]-beta[4])}
		if(s==2){return(1-beta[2]-beta[5]-beta[8])}
		if(s==3){return(1-beta[4]-beta[6]-beta[9])}
	}

	dens_alm_3d <- function(w,alpha,beta,c){
		if(length(alpha)!=4){return(stop("Wrong length of parameter alpha"))}
		if(length(beta)!=9){return(stop("Wrong length of parameter beta"))}
	  	if(beta[1]+beta[3]+beta[7]>1){return(1e-50)} # see conditions on the beta parameters
	  	if(beta[2]+beta[5]+beta[8]>1){return(1e-50)}
	  	if(beta[4]+beta[6]+beta[9]>1){return(1e-50)}
		if(any(alpha < 1) || any(beta < 0) || any(beta > 1)){return(1e-50)}
 
   		if(c==0){return(interior_alm_3d(w,alpha,beta))}
  
		if(sum(w<c) == 2){ # then we are in a corner 
			ind <- which(w > c)
			return(corners_alm_3d(w,alpha,beta,ind)/c^2)
		} else if(sum(w<c) == 1){
			ind <- which(w >= c)
			w2 <- w[ind]/sum(w[ind])
			edge_surface <- c*sqrt(3)*(1-2*c)/2

			if(w[1]<=1-c && w[2]<=1-c && w[3]<=c && w[1]>=(1-w[2])/2 && w[1]>=1-2*w[2]){ #EDGE {1,2} 			
				edg01 <- integrate(Vectorize(function(x){edges_alm_3d(w=c(x,1-x,0), alpha=alpha, beta=beta, s=1, t=2)}), lower=0, upper=1)$value
				return(edges_alm_3d(c(w2,w[3]),alpha,beta,1,2) * edg01 / edge_surface)
			}

			if(w[1]<=1-c && w[2]<=c && w[3]<=1-c && w[2]<=w[1] && w[1]<=1-2*w[2]){ #EDGE {1,3} 			
				edg01 <- integrate(Vectorize(function(x){edges_alm_3d(w=c(x,0,1-x), alpha=alpha, beta=beta, s=1, t=3)}), lower=0, upper=1)$value
				return(edges_alm_3d(c(w2[1],w[2],w2[2]),alpha,beta,1,3) * edg01 / edge_surface)
			}

			if(w[1]<=c && w[2]<=1-c && w[3]<=1-c && w[1]<=w[2] && w[1]<=(1-w[2])/2){ #EDGE {2,3} 			
				edg01 <- integrate(Vectorize(function(x){edges_alm_3d(w=c(0,x,1-x), alpha=alpha, beta=beta, s=2, t=3)}), lower=0, upper=1)$value
				return(edges_alm_3d(c(w[1],w2),alpha,beta,2,3) * edg01 / edge_surface)
			}
		} else {
			int01 <- integrate(Vectorize(function(y) integrate(Vectorize(function(x) interior_alm_3d(c(x,y,1-x-y),alpha=alpha, beta=beta)), lower=0, upper=1-y )$value), lower=0, upper=1)$value
			intc <- integrate(Vectorize(function(y) integrate(Vectorize(function(x) interior_alm_3d(c(x,y,1-x-y),alpha=alpha, beta=beta)), lower=c, upper=1-y-c )$value), lower=c, upper=1-2*c)$value
			
			return(interior_alm_3d(w,alpha,beta)*int01/intc)	
		}
	}


	### Angular density for the Asymmetric Logistic model on the 2 and 3 dimensional simplex

		xvect = as.double(as.vector(t(x)))
   	if (is.vector(x)) {
   	     dim = as.integer(length(x))
   	     n = as.integer(1)
   	     if(round(sum(x),7) !=1){ stop("Data is not angular") }
   	     if(dim==2){ result <- dens_alm_2d(x,alpha,beta,c)}
   	     if(dim==3){ result <- dens_alm_3d(x,alpha,beta,c)}
   	}
   	else {
   		dim = as.integer(ncol(x))
   		n = as.integer(nrow(x))
	   	if (sum(apply(x,1,sum)) != n){ stop("Data is not angular") }
		if (vectorial) {
	        result = double(n)
	        if(dim==2){ result <- apply(x,1,function(y){dens_alm_2d(y,alpha,beta,c)}) }
	        if(dim==3){ result <- apply(x,1,function(y){dens_alm_3d(y,alpha,beta,c)}) }
	    } else { # vectorial=FALSE mean we return the likelihood function
    	    result = as.double(1)
        	if(dim==2){ result <- prod(apply(x,1,function(y){dens_alm_2d(y,alpha,beta,c)})) }
        	if(dim==3){ result <- prod(apply(x,1,function(y){dens_alm_3d(y,alpha,beta,c)})) }
    	}
    }	
    if(log)
    	return(log(result))
    else return(result)
}

### Functions needed for the Extremal-t model

dens_et <- function (x, rho, mu, c, log, vectorial){

	# in the following functions:
	#
	# rho is a vector of size choose(dim,2) which mean choose 2 from dim
	# mu is of size 1 (degree of freedom)


	interior_et_d <- function(w,rho,mu){ # mass on the interior of the d-dim simplex
		p<-length(w);
		k=p-1;
		w.tilde<-rep(0,k);
		Sigma<-diag(k);
	
		if(any(w==0) || any(w==1)){return(1e-300)}	
		
		for(i in 1:k){
			w.tilde[i] <- ((w[i+1]/w[1])^(1/mu)-rho[i])*sqrt((mu+1)/(1-rho[i]^2));
			for(j in i:k){
				if(i==j){Sigma[i,j]=Sigma[j,i]=1}else{
					Sigma[i,j]=Sigma[j,i]=(rho[sum(k:(k-i+1))+j-i]-rho[i]*rho[j])/sqrt((1-rho[i]^2)*(1-rho[j]^2))
				}
			}
		}
		
		if(sum(eigen(Sigma)$values<0)>=1){return(1e-50)} #Check if matrix is positive definite
		
		deriv = (w[-1]/w[1])^(1/mu-1)/mu*sqrt((mu+1)/(1-rho[1:k]^2))
		if(k==1){
			return(dest(x=w.tilde, scale=Sigma, df=mu+1)*w[1]^(-p-1)*prod(deriv) )	
		}else{
			return(dmest(x=w.tilde, scale=Sigma, df=mu+1)*w[1]^(-p-1)*prod(deriv) )	
		}	
	}

	# Mass on the other subsets of the 2d case: 

	corners_et_2d <- function(w,rho,mu,s){ # mass on the corner of the s-th component

		part1 <- pt( -rho * sqrt((mu+1)/(1-rho^2)) ,df=mu+1)	
	
		return(part1)
	}

	dens_et_2d <- function(w,rho,mu,c){
		if(length(rho)!=1){return(stop("Wrong length of parameter rho"))}
		if( (abs(rho) > 1) || (mu<1) ){return(1e-50)}
  
		if(sum(w<1-c) == 1){ # then we are in a corner 
			ind <- which(w > c)
			return(corners_et_2d(w,rho,mu,ind))
		} else {
			return(interior_et_d(w,rho,mu))	
		}
	}

	# Mass on the other subsets of the 3d case:

	corners_et_3d <- function(w,rho,mu,s){ # mass on the corner of the s-th component
		
		if(s==1){rho1=rho[1];rho2=rho[2];rho3=rho[3]}
		if(s==2){rho1=rho[1];rho2=rho[3];rho3=rho[2]}	
		if(s==3){rho1=rho[2];rho2=rho[3];rho3=rho[1]}
		
		R <- (rho3 - rho1*rho2)/sqrt((1-rho1^2) * (1-rho2^2))
		S <- matrix( c(1,R,R,1), ncol=2)
		if(sum(eigen(S)$values<0)>=1){return(1e-50)} # Check if matrix R1 is positive definite
		
		a = -rho1 * sqrt((mu+1)/(1-rho1^2))
		b = -rho2 * sqrt((mu+1)/(1-rho2^2))
		
		return( pmest(x=c(a,b), scale=S, df=mu+1) )	
	}

	# mass on the edge linking the s-th and t-th components (in inceasing order)

	edges_et_3d <- function(w,rho,mu,s,t){ 
		
		if(s==1 && t==2){rho1=rho[1];rho2=rho[2];rho3=rho[3]}
		if(s==1 && t==3){rho1=rho[2];rho2=rho[1];rho3=rho[3]}
		if(s==2 && t==3){rho1=rho[3];rho2=rho[1];rho3=rho[2]}
		x=w[s];y=w[t];
		
		A1 <- sqrt((mu+1)/(1-rho1^2)) * ( (y/x)^(1/mu) -rho1 )
		A2 <- sqrt((mu+1)/(1-rho1^2)) * ( (x/y)^(1/mu) -rho1 )
		B1 <- - rho2 * sqrt((mu+1)/(1-rho2^2))
		B2 <- - rho3 * sqrt((mu+1)/(1-rho3^2))
		d1 <- sqrt((mu+1)/(1-rho1^2)) / mu * y^(1/mu-1) /x^(1/mu)
		d2 <- sqrt((mu+1)/(1-rho1^2)) / mu * x^(1/mu-1) /y^(1/mu)
		R1 <- (rho3 - rho1*rho2)/sqrt((1-rho1^2) * (1-rho2^2))
		R2 <- (rho2 - rho1*rho3)/sqrt((1-rho1^2) * (1-rho3^2))
		
		S1 <- matrix(c(1,R1,R1,1),ncol=2)
		S2 <- matrix(c(1,R2,R2,1),ncol=2)
		
		if(sum(eigen(S1)$values<0)>=1){return(1e-50)} # Check if matrix R1 is positive definite
		if(sum(eigen(S2)$values<0)>=1){return(1e-50)} # Check if matrix R2 is positive definite
		if(any( abs(c(R1,R2)) >1 ) ){return(1e-50)}
		
		part11 <- dt(x=A1,df=mu+1) * pt(sqrt((mu+2)/(mu+1+A1^2))*(B1-R1*A1)/sqrt(1-R1^2),df=mu+2) * d1 / x^2
		part12 <- -1 - 1/mu + y * d1 * (mu+2) * A1 / (mu+1+A1^2)
		part13 <- dt(x=A1,df=mu+1) * dt(x=sqrt((mu+2)/(mu+1+A1^2))*(B1-R1*A1)/sqrt(1-R1^2),df=mu+2) * d1^2 * y / x^2
		part14 <- sqrt((mu+2)/(1-R1^2)) * (R1*(mu+1)+B1*A1) / (mu+1+A1^2)^(3/2)
		
		part21 <- dt(x=A2,df=mu+1) * pt(sqrt((mu+2)/(mu+1+A2^2))*(B2-R2*A2)/sqrt(1-R2^2),df=mu+2) * d2 / y^2
		part22 <- -1 - 1/mu + x * d2 * (mu+2) * A2 / (mu+1+A2^2)
		part23 <- dt(x=A2,df=mu+1) * dt(x=sqrt((mu+2)/(mu+1+A2^2))*(B2-R2*A2)/sqrt(1-R2^2),df=mu+2) * d2^2 * x / y^2
		part24 <- sqrt((mu+2)/(1-R2^2)) * (R2*(mu+1)+B2*A2) / (mu+1+A2^2)^(3/2)
		
		return( -(x+y)^3 * (part11*part12 + part13*part14 + part21*part22 + part23*part24))	
		
	}
	
	dens_et_3d <- function(w,rho,mu,c){
		
		if(length(rho)!=3){return(stop("Wrong length of parameter rho"))}
		if(any(abs(rho)>=1) || mu<=0 ){return(1e-50)}

		if(c==0){return(interior_et_d(w,rho,mu))}

		if(sum(w<c) == 2){ # then we are in a corner 
			ind <- which(w > c)
			return(corners_et_3d(w,rho,mu,ind)/c^2)
		}else if(sum(w<c) == 1){
			ind <- which(w >= c)
			w2 <- w[ind]/sum(w[ind])
			edge_surface <- c*sqrt(3)*(1-2*c)/2
		 
			if(w[1]<=1-c && w[2]<=1-c && w[3]<=c && w[1]>=(1-w[2])/2 && w[1]>=1-2*w[2]){ #EDGE {1,2} 			
				edg01 <- integrate(Vectorize(function(x){edges_et_3d(w=c(x,1-x,0), rho=rho, mu=mu, s=1, t=2)}), lower=0, upper=1)$value
				return(edges_et_3d(c(w2,w[3]),rho,mu,1,2) * edg01 / edge_surface)
			}

			if(w[1]<=1-c && w[2]<=c && w[3]<=1-c && w[2]<=w[1] && w[2]<=(1-w[1])/2){ #EDGE {1,3} 			
				edg01 <- integrate(Vectorize(function(x){edges_et_3d(w=c(x,0,1-x), rho=rho, mu=mu, s=1, t=3)}), lower=0, upper=1)$value
				return(edges_et_3d(c(w2[1],w[2],w2[2]),rho,mu,1,3) * edg01 / edge_surface)
			}

			if(w[1]<=c && w[2]<=1-c && w[3]<=1-c && w[1]<=w[2] && w[1]<=(1-w[2])/2){ #EDGE {2,3} 			
				edg01 <- integrate(Vectorize(function(x){edges_et_3d(w=c(0,x,1-x), rho=rho, mu=mu, s=2, t=3)}), lower=0, upper=1)$value
				return(edges_et_3d(c(w[1],w2),rho,mu,2,3) * edg01 / edge_surface)
			}
		}	else {
			int01 <- integrate(Vectorize(function(y) integrate(Vectorize(function(x) interior_et_d(c(x,y,1-x-y),rho=rho, mu=mu)), lower=0, upper=1-y )$value), lower=0, upper=1)$value
			intc <- integrate(Vectorize(function(y) integrate(Vectorize(function(x) interior_et_d(c(x,y,1-x-y),rho=rho, mu=mu)), lower=c, upper=1-y-c )$value), lower=c, upper=1-2*c)$value
			return(interior_et_d(w,rho,mu)*int01/intc)
		}
	}

	### Angular density for the Extremal-t model on the 2 and 3 dimensional simplex

	xvect = as.double(as.vector(t(x)))
    if (is.vector(x)) {
        dim = as.integer(length(x))
        n = as.integer(1)
        if(round(sum(x),7) != 1){ stop("Data is not angular")}
        if(dim==2){ result <- dens_et_2d(x,rho,mu,c)} 
    		if(dim==3){ result <- dens_et_3d(x,rho,mu,c)}  
    }
    else {
        dim = as.integer(ncol(x))
        n = as.integer(nrow(x))
    	
		if (sum(apply(x,1,sum)) != n){ stop("Data is not angular") }
    
    	if (vectorial) {
    	    result = double(n)
    	    if(dim==2){ result <- apply(x,1,function(y){dens_et_2d(y,rho,mu,c)}) }
    	    if(dim==3){ result <- apply(x,1,function(y){dens_et_3d(y,rho,mu,c)}) }
    	} else { # vectorial=FALSE mean we return the likelihood function
    	    result = as.double(1)
    	    if(dim==2){ result <- prod(apply(x,1,function(y){dens_et_2d(y,rho,mu,c)})) }
    	    if(dim==3){ result <- prod(apply(x,1,function(y){dens_et_3d(y,rho,mu,c)})) }
   		}
   	}	
    if(log){
    		if(any(result==0)){
    			result[which(result==0)] <- log(-1e+50)
    		result[which(result!=0)] <- log(result[which(result!=0)]) 
    		return(result)   			
    		}else{
    			return(log(result))
    		}
    	}else{ 
    		return(result)
	}
}

### Functions needed for the Extremal Skew-t model

dens_est <- function (x, rho, alpha, mu, c, log, vectorial){
	
	## Preliminary functions
	
	Sigma <- function(rho){
			p <- round(uniroot(function(x){length(rho)-choose(x,2)},lower=1, upper=10)$root) 
			Sig <- diag(p)
			Sig[lower.tri(Sig)] = rho
			Sig[row(Sig) < col(Sig)] = t(Sig)[row(Sig) < col(Sig)]
			return(Sig)
		}

	Sigma_j <-function(rho,j){
			Sig <- Sigma(rho)	
			return(Sig[-j,-j] - Sig[-j,j] %*% t(Sig[j,-j]) )
		}
		
	s_j <- function(rho,j){
			p <- round(uniroot(function(x){length(rho)-choose(x,2)},lower=1, upper=10)$root)
			k=p-1;

			sigma_j <- Sigma_j(rho,j)
			M <- diag(k)
			diag(M) = sqrt(diag(sigma_j)) 
			return( M )
		}
		
	Sigma_bar_j <- function(rho,j){
			sigma_j <- Sigma_j(rho,j)
			sj <- s_j(rho,j)
			return( solve(sj) %*% t(sigma_j) %*% solve(sj) )
		}
		
	alpha_tilde <- function(alpha,j){
			return(t(alpha[-j]))
		}
		
	alpha_bar_j <- function(rho,alpha,j){
			Sig <- Sigma(rho) 
			sigma_j <- Sigma_j(rho,j)
			Alpha_tilde <- alpha_tilde(alpha,j)
	
			num <- alpha[j] + Sig[j,-j] %*% t(Alpha_tilde)
			denom <- sqrt( 1 + Alpha_tilde %*% t(sigma_j) %*% t(Alpha_tilde)  )
			return(num/denom)
		}		
		
	alpha_star_j <- function(alpha,j){
			Alpha_tilde <- alpha_tilde(alpha,j)
			sj <- s_j(rho,j)
			return( Alpha_tilde %*% sj )
		}

	tau_star_j <- function(rho,alpha,j){
			Sig <- Sigma(rho)
			Alpha_tilde <- alpha_tilde(alpha,j)
			return( sqrt(mu+1) * (alpha[j] + Sig[-j,j] %*% t(Alpha_tilde) )	 )
		}	

	nu_p <- function(rho,alpha,mu,j){
			Alpha_bar <- alpha_bar_j(rho,alpha,j)
			return( pt(sqrt(mu+1)*Alpha_bar,df=mu+1) * 2^(mu/2-1) * gamma((mu+1)/2) / sqrt(pi) )
		}

	x_bar <- function(x,rho,alpha,mu,j){
			return( x * nu_p(rho,alpha,mu,j) )
		}

	comp <- function(z1,z2,z3,rho,alpha,mu,j,k){

			#function corresponding to the j-th component of I_k
			# Gives x1,y1, x2,y2, x3,y3
			# a1,b1, a2,b2, a3,b3
	

			#j = 1 or 2
			#k = 1,2 or 3
	
			if(j==1 & k==1){corr=rho[1];x=x_bar(z2,rho,alpha,mu,2); y=x_bar(z1,rho,alpha,mu,1)}
			if(j==1 & k==2){corr=rho[1];x=x_bar(z1,rho,alpha,mu,1); y=x_bar(z2,rho,alpha,mu,2)}	
			if(j==1 & k==3){corr=rho[2];x=x_bar(z1,rho,alpha,mu,1); y=x_bar(z3,rho,alpha,mu,3)}

			if(j==2 & k==1){corr=rho[2];x=x_bar(z3,rho,alpha,mu,3); y=x_bar(z1,rho,alpha,mu,1)}
			if(j==2 & k==2){corr=rho[3];x=x_bar(z3,rho,alpha,mu,3); y=x_bar(z2,rho,alpha,mu,2)}	
			if(j==2 & k==3){corr=rho[3];x=x_bar(z2,rho,alpha,mu,2); y=x_bar(z3,rho,alpha,mu,3)}
	
			if(z1==0 && z2==0 && z3==0){
				return( -corr * sqrt((mu+1)/(1-corr^2)) )
			}else{
				return( sqrt((mu+1)/(1-corr^2))*((x/y)^(1/mu)-corr) )
			}
		}

	R_j <- function(rho,alpha,j,matrix=TRUE){
			sigma_bar <- Sigma_bar_j(rho,j)
			d <- ncol(sigma_bar)
			delta <- delta_j(rho,alpha,j)

			S <- diag(d+1)
			S[(1:d),(1:d)] <- sigma_bar
			S[(1:d),d+1] <- S[d+1,(1:d)] <- -delta
			if (matrix==TRUE){return(S)}else{
				seq <- vector("numeric")
				for(i in 1:d){
					seq <- c(seq,S[(i+1):(d+1),i])
				}
				return(seq)
			}	
		}

	delta_j <- function(rho,alpha,j){
	
			sigma_bar <- Sigma_bar_j(rho,j)
			alpha_star <- alpha_star_j(alpha,j)
			return( (alpha_star %*% sigma_bar )/ as.numeric(sqrt( 1 + alpha_star %*% sigma_bar %*% t(alpha_star) )) )
		}

	tau_bar_j <- function(rho,alpha,j){
			tau_star <- tau_star_j(rho,alpha,j)
			alpha_star <- alpha_star_j(alpha,j)
			sigma_bar <- Sigma_bar_j(rho,j)	
			return( tau_star / sqrt( 1 + alpha_star %*% sigma_bar %*% t(alpha_star) ) )
		}

	# in the following functions:
	#
	# rho is a vector of size choose(dim,2) which mean choose 2 from dim (correlation coefficients)
	# alpha is a vector of size dim (skewness parameters)
	# mu is of size 1 (degree of freedom)

	## Mass on the interior of the simplex
	
	interior_skewt_d <- function(w,rho,alpha,mu){ # mass on the interior of the d-dim simplex
		p<-length(w);
		k=p-1;
		w.tilde<-rep(0,k);
				
		cond <- sapply(1:p,function(j){alpha_tilde(alpha,j) %*% t(Sigma_j(rho,j)) %*% t(alpha_tilde(alpha,j))})
		
		if( any(cond < -1)){return(1e-50)}else{		
			sigma_bar1 <- Sigma_bar_j(rho,1)
			alpha_star1 <- alpha_star_j(alpha,1)
			tau_star1 <- tau_star_j(rho,alpha,1)
	
			w_bar <- sapply(1:p,function(j){w[j]*nu_p(rho,alpha,mu,j)})
			
			for(i in 1:k){
				w.tilde[i] <- ((w_bar[i+1]/w_bar[1])^(1/mu)-rho[i])*sqrt((mu+1)/(1-rho[i]^2));
			}	
			
			# if( alpha_star1 %*% sigma_bar1 %*% t(alpha_star1) < -1){return(1e-50)}
			# if( t(w.tilde) %*% solve(sigma_bar1) %*% w.tilde < -1){return(1e-50)}
			
			deriv = (w_bar[-1]/w_bar[1])^(1/mu-1) / mu * sqrt((mu+1)/(1-rho[1:k]^2)) * sapply(2:p,function(x){nu_p(rho,alpha,mu,x)}) / nu_p(rho,alpha,mu,1)

		if(k==1){
			return(dest(x=w.tilde, scale=sigma_bar1, shape=t(alpha_star1), extended=tau_star1, df=mu+1)*w[1]^(-p-1)*prod(deriv) )	
		}else{
			return(dmest(x=w.tilde, location=rep(0,2), scale=sigma_bar1, shape=t(alpha_star1), extended=tau_star1, df=mu+1) * w[1]^(-p-1) * prod(deriv) )	
		}	
										
		}	
	}
	
	## mass on the corner of the s-th component of the 2-d simplex
	
	corners_skewt_2d <- function(w,rho,alpha,mu,s){ 

		alpha_star <- alpha_star_j(alpha,s)
		tau_star <- tau_star_j(rho,alpha,s)
		part1 <- pest(-rho * sqrt((mu+1)/(1-rho^2)),shape=alpha_star, extended=tau_star ,df=mu+1)	
	
		return(part1)
	}
	
	## density on 2-d simplex
	
	dens_skewt_2d <- function(w,rho,alpha,mu,c){
		if(length(rho)!=1){return(stop("Wrong length of parameter rho"))}
		if( (abs(rho) > 1) || (mu<1) ){return(1e-50)}
  
		if(sum(w<c) == 1){ # then we are in a corner 
			ind <- which(w > c)
			return(corners_skewt_2d(w,rho,alpha,mu,ind))
		} else {
			return(interior_skewt_d(w,rho,alpha,mu))	
		}
	}

	## mass on the corner of the s-th component of the 3-d simplex

	corners_skewt_3d <- function(w,rho,alpha,mu,s){
 
		s_bar <- Sigma_bar_j(rho,s)
		al=t(alpha_star_j(alpha,s))
		if(t(al) %*% s_bar %*% al < -1 ){return(1e-50)}
		
		up <- c(comp(0,0,0,rho,alpha,mu,1,s), comp(0,0,0,rho,alpha,mu,2,s))
		return(pmest(x=up, location=rep(0,2), scale=s_bar, shape=al, extended=tau_star_j(rho,alpha,s), df=mu+1 ))	
	}
			
	## mass on the edge between the s-th and t-th components of the 3-d simplex

	edges_skewt_3d <- function(w,rho,alpha,mu,s,t){
	
		sigma_j1 <- Sigma_j(rho,s)
		Alpha_tilde1 <- alpha_tilde(alpha,s)		
		sigma_j2 <- Sigma_j(rho,t)
		Alpha_tilde2 <- alpha_tilde(alpha,t)		
		alpha_star1 <- alpha_star_j(alpha,s)
		sigma_bar1 <- Sigma_bar_j(rho,s)	
		alpha_star2 <- alpha_star_j(alpha,t)
		sigma_bar2 <- Sigma_bar_j(rho,t)	
		if( Alpha_tilde1 %*% t(sigma_j1) %*% t(Alpha_tilde1)< -1){return(1e-50)}
		if( Alpha_tilde2 %*% t(sigma_j2) %*% t(Alpha_tilde2)< -1){return(1e-50)}
		if(	alpha_star1 %*% sigma_bar1 %*% t(alpha_star1)< -1){return(1e-50)} 	
		if(	alpha_star2 %*% sigma_bar2 %*% t(alpha_star2)< -1){return(1e-50)}
				
		ind= c(1,2,3)
		w[which((ind != s)  & (ind != t))]=0
	
		if(s==1 && t==2){
			a1_p <- comp(w[1],w[2],w[3],rho,alpha,mu,1,s)
			b1_p <- comp(w[1],w[2],w[3],rho,alpha,mu,2,s)		
			R1 <- R_j(rho,alpha,s,FALSE)
		
			a2_p <- comp(w[1],w[2],w[3],rho,alpha,mu,1,t)
			b2_p <- comp(w[1],w[2],w[3],rho,alpha,mu,2,t)	
			R2 <- R_j(rho,alpha,t,FALSE)
		}
		if(s==1 && t==3){
			a1_p <- comp(w[1],w[2],w[3],rho,alpha,mu,2,s)
			b1_p <- comp(w[1],w[2],w[3],rho,alpha,mu,1,s)		
			R1 <- R_j(rho,alpha,s,FALSE)
			r1 <- R1[2];
			r2 <- R1[3];
			R1[2] <- r2;
			R1[3] <- r1;	
		
			a2_p <- comp(w[1],w[2],w[3],rho,alpha,mu,1,t)
			b2_p <- comp(w[1],w[2],w[3],rho,alpha,mu,2,t)	
			R2 <- R_j(rho,alpha,t,FALSE)	
		}
		if(s==2 && t==3){
			a1_p <- comp(w[1],w[2],w[3],rho,alpha,mu,2,s)
			b1_p <- comp(w[1],w[2],w[3],rho,alpha,mu,1,s)		
			R1 <- R_j(rho,alpha,s,FALSE)
			r1 <- R1[2];
			r2 <- R1[3];
			R1[2] <- r2;
			R1[3] <- r1;	
		
			a2_p <- comp(w[1],w[2],w[3],rho,alpha,mu,2,t)
			b2_p <- comp(w[1],w[2],w[3],rho,alpha,mu,1,t)	
			R2 <- R_j(rho,alpha,t,FALSE)
			rr1 <- R2[2];
			rr2 <- R2[3];
			R2[2] <- rr2;
			R2[3] <- rr1;	
		}	

		if( R1[1]^2 >1 || R1[2]^2>1 || R2[1]^2 >1 || R2[2]^2>1 ){return(1e-50)}		
		
		c1 <- tau_bar_j(rho,alpha,s)
		u1_p <- c(a1_p,b1_p,c1)
	
		R1_p <- diag(2)
		R1_p[1,2] <- R1_p[2,1] <- (R1[3]-R1[1]*R1[2])/sqrt((1-R1[1]^2)*(1-R1[2]^2))
			
		c2 <- tau_bar_j(rho,alpha,t)
		u2_p <- c(a2_p,b2_p,c2)

		R2_p <- diag(2)
		R2_p[1,2] <- R2_p[2,1] <- (R2[3]-R2[1]*R2[2])/sqrt((1-R2[1]^2)*(1-R2[2]^2))
	
		x1_bar <- x_bar(w[s],rho,alpha,mu,s)
		x2_bar <- x_bar(w[t],rho,alpha,mu,t)
	
		v1_1 <- (b1_p-R1[1]*a1_p)/sqrt(1-R1[1]^2) * sqrt((mu+2)/(mu+1+a1_p^2))
		v2_1 <- (c1-R1[2]*a1_p)/sqrt(1-R1[2]^2) * sqrt((mu+2)/(mu+1+a1_p^2))
		
		v1_1_p <- (a1_p-R1[1]*b1_p)/sqrt(1-R1[1]^2) * sqrt((mu+2)/(mu+1+b1_p^2))
			
		v1_2 <- (b2_p-R2[1]*a2_p)/sqrt(1-R2[1]^2) * sqrt((mu+2)/(mu+1+a2_p^2))
		v2_2 <- (c2-R2[2]*a2_p)/sqrt(1-R2[2]^2) * sqrt((mu+2)/(mu+1+a2_p^2))
	
		v1_2_p <- (a2_p-R2[1]*b2_p)/sqrt(1-R2[1]^2) * sqrt((mu+2)/(mu+1+b2_p^2))
	
		if(s==1 & t==2){rho12=rho[1];rho13=rho[2];rho23=rho[3]}
		if(s==1 & t==3){rho12=rho[2];rho13=rho[1];rho23=rho[3]}
		if(s==2 & t==3){rho12=rho[3];rho13=rho[1];rho23=rho[2]}
		
		part1 <- pt(c1,df=mu+1)
		part11 <- -dt(a1_p,df=mu+1) * pmest(x=c( v1_1, v2_1), scale=R1_p, df=mu+2) * sqrt((mu+1)/(1-rho12^2)) * ( x2_bar/x1_bar )^(1/mu-1)
		part12 <- nu_p(rho,alpha,mu,s)^2 * nu_p(rho,alpha,mu,t) / (mu*x1_bar^3) * ( 1+ 1/mu )
	
			p1 <- dt(a1_p,df=mu+1)/(mu+1+a1_p^2)
			p2 <- (mu+2) * a1_p * pmest(x=c( v1_1, v2_1), scale=R1_p, df=mu+2 )	

			p3 <- dt(v1_1,df=mu+2) * sqrt((mu+2)/(1-R1[1]^2)) * (b1_p*a1_p+R1[1]*(mu+1))/sqrt(mu+1+a1_p^2)
			p4num <- sqrt(mu+3) * ( (c1 -R1[2]*a1_p)*(1-R1[1]^2) - (R1[3]-R1[1]*R1[2])*(b1_p-R1[1]*a1_p) )
			p4denom <- ((1-R1[1]^2)*(mu+1+a1_p^2)+(b1_p-R1[1]*a1_p)^2) * ((1-R1[1]^2)*(1-R1[2]^2)-(R1[3]-R1[1]*R1[2])^2)
			p4 <- pt(p4num/sqrt(p4denom),df=mu+2)

			p5 <- dt(v2_1,df=mu+2) * sqrt((mu+2)/(1-R1[2]^2)) * (c1*a1_p+R1[2]*(mu+1))/sqrt(mu+1+a1_p^2)
			p6num <- sqrt(mu+3) * ( (b1_p -R1[1]*a1_p)*(1-R1[2]^2) - (R1[3]-R1[1]*R1[2])*(c1-R1[2]*a1_p) )
			p6denom <- ((1-R1[2]^2)*(mu+1+a1_p^2)+(c1-R1[2]*a1_p)^2) * ((1-R1[1]^2)*(1-R1[2]^2)-(R1[3]-R1[1]*R1[2])^2)
			p6 <- pt(p6num/sqrt(p6denom),df=mu+2)

		part13 <- (p1 * (p2+ p3*p4 + p5*p6)) * (mu+1)/(1-rho12^2) * ( x2_bar/x1_bar )^(2/mu-1)
		part14 <- nu_p(rho,alpha,mu,s)^2 * nu_p(rho,alpha,mu,t) / (mu^2 * x1_bar^3)
		
		part2 <- pt(c2,df=mu+1)
		part21 <- -dt(a2_p,df=mu+1) * pmest(x=c( v1_2, v2_2), scale=R2_p, df=mu+2 ) * sqrt((mu+1)/(1-rho12^2)) * ( x1_bar/x2_bar )^(1/mu-1)
		part22 <- nu_p(rho,alpha,mu,s) * nu_p(rho,alpha,mu,t)^2 / (mu*x2_bar^3) * ( 1+ 1/mu )
	
			P1 <- dt(a2_p,df=mu+1)/(mu+1+a2_p^2)
			P2 <- (mu+2) * a2_p * pmest(x=c( v1_2, v2_2), scale=R2_p, df=mu+2 )	

			P3 <- dt(v1_2,df=mu+2) * sqrt((mu+2)/(1-R2[1]^2)) * (b2_p*a2_p+R2[1]*(mu+1))/sqrt(mu+1+a2_p^2)
			P4num <- sqrt(mu+3) * ( (c2 -R2[2]*a2_p)*(1-R2[1]^2) - (R2[3]-R2[1]*R2[2])*(b2_p-R2[1]*a2_p) )
			P4denom <- ((1-R2[1]^2)*(mu+1+a2_p^2)+(b2_p-R2[1]*a2_p)^2) * ((1-R2[1]^2)*(1-R2[2]^2)-(R2[3]-R2[1]*R2[2])^2)
			P4 <- pt(P4num/sqrt(P4denom),df=mu+2)

			P5 <- dt(v2_2,df=mu+2) * sqrt((mu+2)/(1-R2[2]^2)) * (c2*a2_p+R2[2]*(mu+1))/sqrt(mu+1+a2_p^2)
			P6num <- sqrt(mu+3) * ( (b2_p -R2[1]*a2_p)*(1-R2[2]^2) - (R2[3]-R2[1]*R2[2])*(c2-R2[2]*a2_p) )	
			P6denom <- ((1-R2[2]^2)*(mu+1+a2_p^2)+(c2-R2[2]*a2_p)^2) * ((1-R2[1]^2)*(1-R2[2]^2)-(R2[3]-R2[1]*R2[2])^2)
			P6 <- pt(P6num/sqrt(P6denom),df=mu+2)
	
		part23 <- (P1 * (P2+ P3*P4 + P5*P6)) * (mu+1)/(1-rho12^2) * ( x1_bar/x2_bar )^(2/mu-1)
		part24 <- nu_p(rho,alpha,mu,s) * nu_p(rho,alpha,mu,t)^2 / (mu^2 * x2_bar^3)
			
		return( -(w[s]+w[t])^3 * ((part11*part12+part13*part14)/part1 + (part21*part22+part23*part24)/part2)) 
	
	}
	
	## density on 3-d simplex
	
	dens_skewt_3d <- function(w,rho,alpha,mu,c){

		if(length(rho)!=3){return(stop("Wrong length of parameter rho"))}
		if(length(alpha)!=3){return(stop("Wrong length of parameter rho"))}
		if(any(abs(rho)>=1) || mu<=0){return(1e-50)}
  
  		if(c==0){return(interior_skewt_d(w,rho,alpha,mu))}	
  
		if(sum(w<c) == 2){ # then we are in a corner 
			ind <- which(w > c)
			return(corners_skewt_3d(w,rho,alpha,mu,ind)/c^2)
		} else 
		if(sum(w<c) == 1){
			ind <- which(w >= c)
			w2 <- w[ind]/sum(w[ind])
			edge_surface <- c*sqrt(3)*(1-2*c)/2
		 
			if(w[1]<=1-c && w[2]<=1-c && w[3]<=c && w[1]>=(1-w[2])/2 && w[1]>=1-2*w[2]){ #EDGE {1,2} 			
				edg01 <- integrate(Vectorize(function(x){edges_skewt_3d(w=c(x,1-x,0), rho=rho, alpha=alpha, mu=mu, s=1, t=2)}), lower=0, upper=1)$value
				return(edges_skewt_3d(c(w2,w[3]),rho,alpha,mu,1,2) * edg01 / edge_surface)
			}

			if(w[1]<=1-c && w[2]<=c && w[3]<=1-c && w[2]<=w[1] && w[1]<=1-2*w[2]){ #EDGE {1,3} 			
				edg01 <- integrate(Vectorize(function(x){edges_skewt_3d(w=c(x,0,1-x), rho=rho, alpha=alpha, mu=mu, s=1, t=3)}), lower=0, upper=1)$value
				return(edges_skewt_3d(c(w2[1],w[2],w2[2]),rho,alpha,mu,1,3) * edg01 / edge_surface)
			}

			if(w[1]<=c && w[2]<=1-c && w[3]<=1-c && w[1]<=w[2] && w[1]<=(1-w[2])/2){ #EDGE {2,3} 			
				edg01 <- integrate(Vectorize(function(x){edges_skewt_3d(w=c(0,x,1-x), rho=rho, alpha=alpha, mu=mu, s=2, t=3)}), lower=0, upper=1)$value
				return(edges_skewt_3d(c(w[1],w2),rho,alpha,mu,2,3) * edg01 / edge_surface)
			}
		}	else {
			int01 <- integrate(Vectorize(function(y) integrate(Vectorize(function(x) interior_skewt_d(c(x,y,1-x-y),rho=rho, alpha=alpha, mu=mu)), lower=0, upper=1-y )$value), lower=0, upper=1)$value
			intc <- integrate(Vectorize(function(y) integrate(Vectorize(function(x) interior_skewt_d(c(x,y,1-x-y),rho=rho, alpha=alpha, mu=mu)), lower=c, upper=1-y-c )$value), lower=c, upper=1-2*c)$value
			
			return(interior_skewt_d(w,rho,alpha,mu)*int01/intc)
		}
	}
	
	### Angular density for the Extremal Skew-t model on the 2 and 3 dimensional simplex

	xvect = as.double(as.vector(t(x)))
    if (is.vector(x)) {
        dim = as.integer(length(x))
        n = as.integer(1)
        if(round(sum(x),7) != 1){ stop("Data is not angular")}
        if(dim==2){ result <- dens_skewt_2d(x,rho,alpha,mu,c)} 
    	if(dim==3){ result <- dens_skewt_3d(x,rho,alpha,mu,c)} 
    }
    else {
        dim = as.integer(ncol(x))
        n = as.integer(nrow(x))
    	
		if (sum(apply(x,1,sum)) != n){ stop("Data is not angular") }
    
    	if (vectorial) {
    	    result = double(n)
    	    if(dim==2){ result <- apply(x,1,function(y){dens_skewt_2d(y,rho,alpha,mu,c)}) }
    	    if(dim==3){ result <- apply(x,1,function(y){dens_skewt_3d(y,rho,alpha,mu,c)}) }
    	} else { # vectorial=FALSE mean we return the likelihood function
    	    result = as.double(1)
    	    if(dim==2){ result <- prod(apply(x,1,function(y){dens_skewt_2d(y,rho,alpha,mu,c)})) }
    	    if(dim==3){ result <- prod(apply(x,1,function(y){dens_skewt_3d(y,rho,alpha,mu,c)})) }
   		}
   	}	
    if(log){
		if(result==0){return(log(1e-50))}else{
			return(log(result))
		}
	}else return(result)	
}	

# Dimension of the density
if (is.vector(x)) { d <- as.integer(length(x))}else{d <- as.integer(ncol(x)) }	

if(model=='Pairwise'){
	if(length(par)!= (choose(d,2)+1) ){stop('Wrong length of parameters')}
	return(dens_pb(x=x, b=par[1:choose(d,2)], alpha=par[choose(d,2)+1], log=log,  vectorial=vectorial))
}
if(model=='Husler'){
	if(length(par)!= choose(d,2) ){stop('Wrong length of parameters')}
	return(dens_hr(x=x, lambda=par, log=log,  vectorial=vectorial))
}
if(model=='Dirichlet'){
	if(length(par)!= d ){stop('Wrong length of parameters')}
	return(dens_di(x=x, para=par, log=log,  vectorial=vectorial))
}
if(model=='Extremalt'){
	if(length(par)!= (choose(d,2)+1) ){stop('Wrong length of parameters')}
	if(is.null(c)){stop('c needs to be specified')}
	return(dens_et(x=x, rho=par[1:choose(d,2)], mu=par[choose(d,2)+1], c=c, log=log,  vectorial=vectorial))
}
if(model=='Skewt'){
	if(length(par)!= (choose(d,2)+d+1) ){stop('Wrong length of parameters')}
	if(is.null(c)){stop('c needs to be specified')}
	return(dens_est(x=x, rho=par[1:choose(d,2)], alpha=par[choose(d,2)+1:d],mu=par[choose(d,2)+d+1], c=c, log=log,  vectorial=vectorial))
}
if(model=='Asymmetric'){
	if(is.null(c)){stop('c needs to be specified')}	
	if(d==2){
		if(length(par)!=3 ){
			stop('Wrong length of parameters')
		}else{
			return(dens_al(x=x, alpha=par[1], beta=par[2:3], c=c, log=log,  vectorial=vectorial))
		}
	}
	if(d==3){
		if(length(par)!=13 ){
			stop('Wrong length of parameters')
		}else{
			return(dens_al(x=x, alpha=par[1:4], beta=par[5:13], c=c, log=log,  vectorial=vectorial))
		}	
	}
}
	
	
}

### Definition of the angular density function and its rescaled version, i.e.
#
#  h*(w) := A''(w)/(A'(1) - A'(0))
#
dh <- function(w, beta, mixture = FALSE){
  k <- length(beta) - 1  
  j <- 1:(k-1)
  const <- 2/(k * (2-beta[2]-beta[k]))
  res <- diff(diff(beta)) * dbeta(w, j, k-j) 
  if(mixture) return(k/2 * const * sum(res))
  else return(k/2 * sum(res))
}

### Definition of the angular probability measure

ph <- function(w, beta)
{
  k <- length(beta) - 1  
  j <- 1:k
  #  res <- diff(beta) * dbeta(w, j, k-j+1)
  res <- 0.5 * (diff(beta) + 1/k) * dbeta(w, j, k-j+1)
  return(sum(res))
}

### Simulate from a mixture of beta densities 
# 
#
#
rh.mixt <- function(N=1000, beta){
  k <- length(beta) - 1
  # Sample N random uniforms U
  U = runif(N)
  
  #Variable to store the samples from the mixture distribution                                             
  rand.samples = rep(NA,N)
  
  #probabilities mixture
  prob <- diff(diff(beta)) / (2-beta[2]-beta[k])
  
  #Sampling from the mixture
  rand.samples <- sapply(1:N, function(i, prob, k, U){ 
    j <- min(which((U[i]<cumsum(prob)))); return(rbeta(1,j,k-j))},
    prob, k, U)
  
  return(rand.samples)
}

### Simulate from an angular density
#
#
#
rh <- function(N=1000, beta){
  k <- length(beta) - 1
  #Sample N random uniforms U
  U = runif(N)
  
  p0 <- (1+k*(beta[2]-1))/2
  p1 <- (1-k*(1-beta[k]))/2
  prob <- c(p0, 1-p0-p1, p1)
  
  rand.samples <- sapply(1:N, function(i, prob, U){
    j <- min(which((U[i]<cumsum(prob))));
    return(switch(j, 0, rh.mixt(1,beta),1))},
    prob,U)
}

### Estimation and generation from angular density

angular <- function(data, model, n, dep, asy, alpha, beta, df, seed, k, nsim, plot=TRUE, nw=100){
  
  w <- seq(0.00001, .99999, length=nw)
  
  # Simulation of synthetic data
  if(missing(data)){
    if(missing(model) || missing(n) || missing(dep) || missing(seed) ){stop("model, n. dep and missing must be specified when data is not provided")}
    # warning("Data not provided and generated according to the parameters provided")
    
    if(any(model==c("log", "alog", "hr", "neglog", "aneglog", "bilog", "negbilog", "ct", "amix"))){
      set.seed(seed)
      data <- rbvevd(n=n, dep = dep, asy, alpha, beta, model = model, mar1 = c(1,1,1))
      Atrue <- abvevd(w, dep=dep, asy, alpha, beta, model=model) # True pickands dependence function
      htrue <- hbvevd(1-w, dep=dep, asy, alpha, beta, model=model, half=TRUE) # True angular density
      if(any(model==c("alog","aneglog"))){Htrue <- (1-asy)/2}
      if(model=="amix"){Htrue <- c(1-alpha-beta, 1-alpha-2*beta)/2}
    }else if(model=="Extremalt"){
      set.seed(seed)
      data <- r_extr_mod(model, n=n, param=c(dep, df))
      Atrue <- rep(NA,nw)
      for(i in 1:nw){Atrue[i] <- pk.extst(c(w[i],1-w[i]), c(dep,0,0,df))}        
      htrue <- dens(x=cbind(w,1-w), model=model, par=c(dep,df), c=0, log=FALSE, vectorial=TRUE)/2
      Htrue <- dens(x=cbind(c(1,0),c(0,1)), model=model, par=c(dep,df), c=0, log=FALSE, vectorial=TRUE)/2
    }
  }else{
    if(!is.matrix(data) || ncol(data)!=2){stop("data must be a matrix with 2 columns")}
    n <- nrow(data)
    model <- dep <- seed <- NULL
    Atrue <- htrue <- 0
    warning("True Pickands function and angular density not computed as data is provided and true model is unsure")
  }
  
  if(missing(k)){stop("k, the polynomial order must be specified")}
  if(missing(nsim)){stop("nsim, the number of of generation from the estimated angular density must be specified")}  
  
  # Compute the angles:
  wdata <- data[,1]/rowSums(data)
  
  # Estimate the pickands function:
  Aest <- beed(data, cbind(w, 1-w), 2, 'md', 'emp', k=k, plot=FALSE)
  beta <- Aest$beta
  
  # Compute the angular density in the interior
  hest <- sapply(1:nw, function(i, w, beta) dh(w[i], beta), w, beta)
  
  # Compute the angular measure
  Hest <- sapply(1:nw, function(i, w, beta) ph(w[i], beta), w, beta)
  
  # Compute the point masses on the corners
  p0 <- (1+k*(beta[2]-1))/2
  p1 <- (1-k*(1-beta[k]))/2
  
  # simulate from the semiparametric angular density
  wsim <- rh(nsim, beta)
  
  if(plot){
    par(mai = c(0.85, 0.85, 0.1, 0.1), mgp = c(2.8, 1, 0), cex.axis=2, cex.lab=2)
    hist(wsim[wsim!=0 & wsim!=1], freq=FALSE, ylim=c(0,3.5), xlim=c(0,1), xlab='w', ylab='h(w)', main="",col="gray")
    lines(w, hest, col=1, lty=2, lwd=2.5)
    if(is.vector(htrue)){ lines(w[2:99], htrue[2:99],col=2, lty=1, lwd=1.5)}
    points(1,sum(wsim==1)/nsim, pch=19, cex=2)
    points(0,sum(wsim==0)/nsim, pch=19, cex=2)
    if(any(model==c('alog','aneglog','amix','Extremalt'))){
      points(0, Htrue[1] , pch=21, cex=2,col=2, lwd=2)
      points(1, Htrue[2], pch=21, cex=2,col=2, lwd=2)
    }
  }
  
  return(list(model=model, n=n, dep=dep, data=data, Aest=Aest, hest=hest, Hest=Hest, p0=p0, p1=p1, wsim=wsim, Atrue=Atrue, htrue=htrue))
  
}
