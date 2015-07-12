################################################################################################
### Authors: Boris Beranger and Simone Padoan        	 									 ###
### 																						 ###	
###	Emails: borisberanger@gmail.com, simone.padoan@unibocconi.it							 ###
### 																						 ###
###	Institutions: Department of Decision Sciences, University Bocconi of Milan				 ###
### School of Mathematics and Statistics, University of New South Wales 					 ###
### 																						 ###
### File name: Densities.r	                 							             	     ###
### 																						 ###
### Description:                                  							      		     ###
### This file provides the Angular densities of extremal dependence models: 				 ###
### 1) Pairwise Beta, Dirichlet and H??sler-Reiss (mass only on the interior of the simplex)  ###
### 2) Asymmetric logistic and Extremal-t (mass on all subsets of the simplex)				 ###
### It also provides the likelihood and log-likelihood functions							 ###
### 																						 ###
### Last change: 11/07/2015                         		  								 ###
### 																						 ###
################################################################################################


dens <- function(x=rbind(c(0.1,0.3,0.6),c(0.1,0.2,0.7)), model="Pairwise", par=c(2,2,2,4), log=FALSE, vectorial=TRUE){

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


### Functions neeed for the H??sler-Reiss model


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




### Functions neeed for the Asymmetric logistic model

dens_al <- function (x, alpha, beta, log, vectorial){

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

	dens_alm_2d <- function(w,alpha,beta){
		if(length(alpha)!=1){return(stop("Wrong length of parameter alpha"))}
		if(length(beta)!=2){return(stop("Wrong length of parameter beta"))}
		if( (alpha < 1) || any(beta < 0) || any(beta > 1)){return(1e-50)}
  
		if(sum(w<0.01) == 1){ # then we are in a corner 
			ind <- which(w > 0.99)
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

	dens_alm_3d <- function(w,alpha,beta){
		if(length(alpha)!=4){return(stop("Wrong length of parameter alpha"))}
		if(length(beta)!=9){return(stop("Wrong length of parameter beta"))}
		if(any(alpha < 1) || any(beta < 0) || any(beta > 1)){return(1e-50)}
	  	if(beta[1]+beta[3]+beta[7]>1){return(1e-50)} # see conditions on the beta parameters
	  	if(beta[2]+beta[5]+beta[8]>1){return(1e-50)}
	  	if(beta[4]+beta[6]+beta[9]>1){return(1e-50)}
  
		if(sum(w<0.01) == 2){ # then we are in a corner 
			ind <- which(w > 0.98)
			return(corners_alm_3d(w,alpha,beta,ind))
		} else if(sum(w<0.01) == 1){
			ind <- which(w >= 0.01)
			return(edges_alm_3d(w,alpha,beta,ind[1],ind[2]))		
		} else {
			return(interior_alm_3d(w,alpha,beta))	
		}
	}

	# The 4d case:

	# in the following functions:
	#
	# alpha is a vector of size 11: for the subsets {1,2}, {1,3}, {1,4}, {2,3}, {2,4}, {3,4}, {1,2,3}, etc...
	# beta is a vector of size 28: for [1,{1,2}], [2,{1,2}], [1,{1,3}], [3,{1,3}], [1,{1,4}], [4,{1,4}], [2,{2,3}], [3,{2,3}], etc...
	# beta for [1,{1}], [2,{2}] and [3,{3}] are omitted as obtained as :
	# 	[1,{1}] = 1 - [1,{1,2}]+[1,{1,3}]+[1,{1,4}]+[1,{1,2,3}]+[1,{1,2,4}]+[1,{1,3,4}]+[1,{1,2,3,4}] etc... 

	interior_alm_4d <- function(w,alpha,beta){
	
		k <- c(1,2,3,4)
	
		part1 <- (alpha[11]-1) * (2*alpha[11]-1) * (3*alpha[11]-1)
		part2 <- prod(beta[24+k]^alpha[11]*w[k]^(-alpha[11]-1))
		part3 <- sum( (beta[24+k]/w[k])^alpha[11] )^(1/alpha[11]-4)

		return(part1*part2*part3)
	}

	faces_alm_4d <- function(w,alpha,beta,s,t,u){ # mass put on a face of the 4-dim simplex between the s-th, t-th and u-th components (in increasing order)
	
		if(t>s || u>s || u>t){return('s, t and u should be in increasing order')}	
	
		if(s==1 && t==2 && u==3){a=alpha[7];b1=beta[13];b2=beta[14];b3=beta[15]}
		if(s==1 && t==2 && u==4){a=alpha[8];b1=beta[16];b2=beta[17];b3=beta[18]}
		if(s==1 && t==3 && u==4){a=alpha[9];b1=beta[19];b2=beta[20];b3=beta[21]}
		if(s==2 && t==3 && u==4){a=alpha[10];b1=beta[22];b2=beta[23];b3=beta[24]}
		w1=w[s];w2=w[t];w3=w[u];
	
		part1 <- (a-1) * (2*a-1) * (b1*b2*b3)^a * (w1*w2*w3)^(-a-1)
		part2 <- ( (b1/w1)^a + (b2/w2)^a + (b3/w3)^a )^(1/a-3) 			

		return(part1*part2)
	}

	edges_alm_4d <- function(w,alpha,beta,s,t){ # mass put on an edge of the 4-dim simplex between the s-th and t-th components (in increasing order)

		if(t>s){return('s and t should be in increasing order')}	
	
		if(s==1 && t==2){a=alpha[1];b1=beta[1];b2=beta[2]}
		if(s==1 && t==3){a=alpha[2];b1=beta[3];b2=beta[4]}
		if(s==1 && t==4){a=alpha[3];b1=beta[5];b2=beta[6]}
		if(s==2 && t==3){a=alpha[4];b1=beta[7];b2=beta[8]}
		if(s==2 && t==4){a=alpha[5];b1=beta[9];b2=beta[10]}
		if(s==3 && t==4){a=alpha[6];b1=beta[11];b2=beta[12]}
		w1=w[s];w2=w[t];
		
		part1 <- (a-1)  * (b1*b2)^a * (w1*w2)^(-a-1)
		part2 <- ( (b1/w1)^a + (b2/w2)^a )^(1/a-2) 			

		return(part1*part2)
	
	}	

	corners_alm_4d <- function(w,alpha,beta,s){ # mass put on a corner of the 4-dim simplex of th s-th component
		if(s==1){return( 1-beta[1]-beta[3]-beta[5]-beta[13]-beta[16]-beta[19]-beta[25] )}
		if(s==2){return( 1-beta[2]-beta[7]-beta[9]-beta[14]-beta[17]-beta[22]-beta[26] )}
		if(s==3){return( 1-beta[4]-beta[8]-beta[11]-beta[15]-beta[20]-beta[23]-beta[27] )}
		if(s==4){return( 1-beta[6]-beta[10]-beta[12]-beta[18]-beta[21]-beta[24]-beta[28] )}			
	}

	dens_alm_4d <- function(w,alpha,beta){
	
		if(length(alpha)!=11){return(stop("Wrong length of parameter alpha"))}
		if(length(beta)!=28){return(stop("Wrong length of parameter beta"))}
		if(any(alpha < 1) || any(beta < 0) || any(beta > 1)){return(1e-50)}
  		if(beta[1]+beta[3]+beta[5]+beta[13]+beta[16]+beta[19]+beta[25]>1){return(1e-50)} # see conditions on the beta parameters
  		if(beta[2]+beta[7]+beta[9]+beta[14]+beta[17]+beta[22]+beta[26]>1){return(1e-50)}
  		if(beta[4]+beta[8]+beta[11]+beta[15]+beta[20]+beta[23]+beta[27]>1){return(1e-50)}
  		if(beta[6]+beta[10]+beta[12]+beta[18]+beta[21]+beta[24]+beta[28]>1){return(1e-50)}
			
		if(sum(w < 0.03) == 3){ # in one of the corners of the 4-d simplex
			ind <- c(1:4)[-which(w<0.03)]
			return(corners_alm_4d(w,alpha,beta,ind))			
		} else if(sum(w < 0.02) == 2){ # in one of the edges of the 4-d simplex
		    ind <- c(1:4)[-which(w<0.02)]
		    return(edges_alm_4d(w,alpha,beta,ind[1],ind[2]))
		}else	if(sum(w < 0.01) == 1){ # in one of the faces of the 4-d simplex
		    ind <- c(1:4)[-which(w<0.01)]
			return(faces_alm_4d(w,alpha,beta,ind[1],ind[2],ind[3]))
		}else{ # in the interior of the 4-d simplex
			return(interior_alm_4d(w,alpha,beta))	
		}	
	}


	### Angular density for the Asymmetric Logistic model on the 2, 3 and 4 dimensional simplex

		xvect = as.double(as.vector(t(x)))
   	if (is.vector(x)) {
   	     dim = as.integer(length(x))
   	     n = as.integer(1)
   	     if(round(sum(x),7) !=1){ stop("Data is not angular") }
   	     if(dim==2){ result <- dens_alm_2d(x,alpha,beta)}
   	     if(dim==3){ result <- dens_alm_3d(x,alpha,beta)}
   	     if(dim==4){ result <- dens_alm_4d(x,alpha,beta)}
   	}
   	else {
   		dim = as.integer(ncol(x))
   		n = as.integer(nrow(x))
	   	if (sum(apply(x,1,sum)) != n){ stop("Data is not angular") }
		if (vectorial) {
	        result = double(n)
	        if(dim==2){ result <- apply(x,1,function(y){dens_alm_2d(y,alpha,beta)}) }
	        if(dim==3){ result <- apply(x,1,function(y){dens_alm_3d(y,alpha,beta)}) }
	        if(dim==4){ result <- apply(x,1,function(y){dens_alm_4d(y,alpha,beta)}) }
	    } else { # vectorial=FALSE mean we return the likelihood function
    	    result = as.double(1)
        	if(dim==2){ result <- prod(apply(x,1,function(y){dens_alm_2d(y,alpha,beta)})) }
        	if(dim==3){ result <- prod(apply(x,1,function(y){dens_alm_3d(y,alpha,beta)})) }
        	if(dim==4){ result <- prod(apply(x,1,function(y){dens_alm_4d(y,alpha,beta)})) }
    	}
    }	
    if(log)
    	return(log(result))
    else return(result)
}


### Functions neeed for the Extremal-t model

dens_et <- function (x, rho, mu, log, vectorial){

	# in the following functions:
	#
	# rho is a vector of size choose(dim,2) which mean choose 2 from dim
	# mu is of size 1 (degree of freedom)


	interior_et_d <- function(w,rho,mu){ # mass on the interior of the d-dim simplex
		p<-length(w);
		k=p-1;
		w.tilde<-rep(0,k);
		Sigma<-diag(k);
		
		for(i in 1:k){
			w.tilde[i] <- ((w[i+1]/w[1])^(1/mu)-rho[i])*sqrt(mu+1);
			for(j in i:k){
				if(i==j){Sigma[i,j]=Sigma[j,i]=1-rho[i]^2}else{
					Sigma[i,j]=Sigma[j,i]=(rho[sum(k:(k-i+1))+j-i]-rho[i]*rho[j])
				}
			}
		}
	
		if(sum(eigen(Sigma)$values<0)>=1){return(1e-50)} #Check if matrix is positive definite
	
		deriv = (w[-1]/w[1])^(1/mu-1)/mu*sqrt(mu+1)
		return(dmvt(w.tilde,rep(0,k),sigma=Sigma,df=mu+1,log=FALSE)*w[1]^(-p-1)*prod(deriv) )	
	
	}

	# Mass on the other subsets of the 2d case: 

	corners_et_2d <- function(w,rho,mu,s){ # mass on the corner of the s-th component

		part1 <- pt( -rho * sqrt((mu+1)/(1-rho^2)) ,df=mu+1)	
	
		return(part1)
	}

	dens_et_2d <- function(w,rho,mu){ # mass on the interior of the 2-dim simplex
		if(length(rho)!=1){return(stop("Wrong length of parameter rho"))}
		if( (abs(rho) > 1) || (mu<1) ){return(1e-50)}
  
		if(sum(w<0.01) == 1){ # then we are in a corner 
			ind <- which(w > 0.99)
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

		a = -rho1 * sqrt((mu+1)/(1-rho1^2))
		b = -rho2 * sqrt((mu+1)/(1-rho2^2))

		return( pmvt(lower=rep(-Inf,2), upper=c(a,b), sigma=S, df=round(mu)+1)[1] )
	
	}

	edges_et_3d <- function(w,rho,mu,s,t){ # mass on the edge linking the s-th and t-th components (in inceasing order)

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
		
		part11 <- dt(A1,df=mu+1) * pt(sqrt((mu+2)/(mu+1+A1^2))*(B1-R1*A1)/sqrt(1-R1^2),df=mu+2) * d1 / x^2
		part12 <- -1 - 1/mu + y * d1 * (mu+2) * A1 / (mu+1+A1^2)
		part13 <- dt(A1,df=mu+1) * dt(sqrt((mu+2)/(mu+1+A1^2))*(B1-R1*A1)/sqrt(1-R1^2),df=mu+2) * d1^2 * y / x^2
		part14 <- sqrt((mu+2)/(1-R1^2)) * (R1*(mu+1)+B1*A1) / (mu+1+A1^2)^(3/2)
	
		part21 <- dt(A2,df=mu+1) * pt(sqrt((mu+2)/(mu+1+A2^2))*(B2-R2*A2)/sqrt(1-R2^2),df=mu+2) * d2 / y^2
		part22 <- -1 - 1/mu + x * d2 * (mu+2) * A2 / (mu+1+A2^2)
		part23 <- dt(A2,df=mu+1) * dt(sqrt((mu+2)/(mu+1+A2^2))*(B2-R2*A2)/sqrt(1-R2^2),df=mu+2) * d2^2 * x / y^2
		part24 <- sqrt((mu+2)/(1-R2^2)) * (R2*(mu+1)+B2*A2) / (mu+1+A2^2)^(3/2)
	
		return( -(x+y)^3 * (part11*part12 + part13*part14 + part21*part22 + part23*part24))	
	
	}

	dens_et_3d <- function(w,rho,mu){

		if(length(rho)!=3){return(stop("Wrong length of parameter rho"))}
		if(any(abs(rho)>=1) || mu<=0 ){return(1e-50)}
  
		if(sum(w<0.01) == 2){ # then we are in a corner 
			ind <- which(w > 0.98)
			return(corners_et_3d(w,rho,mu,ind))
		} else if(sum(w<0.01) == 1){
			ind <- which(w >= 0.01)
			return(edges_et_3d(w,rho,mu,ind[1],ind[2]))
		}	else {
			return(interior_et_d(w,rho,mu))
		}
	}


	# Mass on the other subsets of the 4d case:

		###########################################
		#										  #	
		#   In all of the following functions :   #
		#										  #	
		#	v is of size 3						  #
		#	rho is a 3x3 correlation matrix		  #
		#	mu is the degree of freedom			  #
		#										  #	
		###########################################



	## First derivative of the trivariate t CDF w.r.t the k-th component

	dT_dx <- function(v,rho,mu,k){
		if(k==1){x=v[1];y=v[2];z=v[3];rho12=rho[1];rho13=rho[2];rho23=rho[3]}
		if(k==2){x=v[2];y=v[1];z=v[3];rho12=rho[1];rho13=rho[3];rho23=rho[2]}
		if(k==3){x=v[3];y=v[1];z=v[2];rho12=rho[2];rho13=rho[3];rho23=rho[1]}

		up1 <- (y-rho12*x)/sqrt(1-rho12^2) * sqrt((mu+1)/(mu+x^2))	
		up2 <- (z-rho13*x)/sqrt(1-rho13^2) * sqrt((mu+1)/(mu+x^2))
		rho <- (rho23-rho12*rho13)/sqrt((1-rho12^2)*(1-rho13^2))
	
		part1 <- dt(x,df=mu)
		part2 <- pmvt(lower=rep(-Inf,2),upper=c(up1,up2),sigma=matrix(c(1,rho,rho,1),ncol=2),df=round(mu)+1)
		return((part1*part2)[1])		
	}

	## Second derivative of the trivariate t CDF w.r.t the k-th component

	dT2_dx2 <- function(v,rho,mu,k){
		if(k==1){x=v[1];y=v[2];z=v[3];rho12=rho[1];rho13=rho[2];rho23=rho[3]}
		if(k==2){x=v[2];y=v[1];z=v[3];rho12=rho[1];rho13=rho[3];rho23=rho[2]}
		if(k==3){x=v[3];y=v[1];z=v[2];rho12=rho[2];rho13=rho[3];rho23=rho[1]}

		up1 <- (y-rho12*x)/sqrt(1-rho12^2) * sqrt((mu+1)/(mu+x^2))
		up2 <- (z-rho13*x)/sqrt(1-rho13^2) * sqrt((mu+1)/(mu+x^2))
		rho <- (rho23-rho12*rho13)/sqrt((1-rho12^2)*(1-rho13^2))

		part1 <- dt(x,df=mu)/(mu+x^2)
		part2 <- pmvt(lower=rep(-Inf,2),upper=c(up1,up2),sigma=matrix(c(1,rho,rho,1),ncol=2),df=round(mu)+1) * (-x*(mu+1))
		part31 <- dt(up1,df=mu+1) * sqrt((mu+1)/(1-rho12^2)) * (-rho12*(mu+x^2)-(y-rho12*x)*x)/sqrt(mu+x^2)
		part32 <- pt(sqrt((mu+2)/(mu+1+up1^2))*(up2-rho*up1)/sqrt(1-rho^2),df=mu+2)
		part41 <- dt(up2,df=mu+1) * sqrt((mu+1)/(1-rho13^2)) * (-rho13*(mu+x^2)-(z-rho13*x)*x)/sqrt(mu+x^2)
		part42 <- pt(sqrt((mu+2)/(mu+1+up2^2))*(up1-rho*up2)/sqrt(1-rho^2),df=mu+2)
	
		return((part1*(part2+part31*part32+part41*part42))[1]) 	
	}

	## Second derivative of the trivariate t CDF w.r.t the k and j-th components

	dT2_dxdy <- function(v,rho,mu,k,j){
		s <- sort(c(k,j))
		k=s[1]; j=s[2];
		if(k==1 && j==2){x=v[1];y=v[2];z=v[3];rho12=rho[1];rho13=rho[2];rho23=rho[3]}
		if(k==1 && j==3){x=v[1];y=v[3];z=v[2];rho12=rho[2];rho13=rho[1];rho23=rho[3]}
		#if(k==2 && j==1){x=v[2];y=v[1];z=v[3];rho12=rho[1];rho13=rho[3];rho23=rho[2]}
		if(k==2 && j==3){x=v[2];y=v[3];z=v[1];rho12=rho[3];rho13=rho[1];rho23=rho[2]}
		#if(k==3 && j==1){x=v[3];y=v[1];z=v[2];rho12=rho[2];rho13=rho[3];rho23=rho[1]}
		#if(k==3 && j==2){x=v[3];y=v[2];z=v[1];rho12=rho[3];rho13=rho[2];rho23=rho[1]}	

		up1 <- (x-rho12*y)/sqrt(1-rho12^2) * sqrt((mu+1)/(mu+y^2))
		up2 <- (z-rho23*y)/sqrt(1-rho23^2) * sqrt((mu+1)/(mu+y^2))
		rho <- (rho13-rho12*rho23)/sqrt((1-rho12^2)*(1-rho23^2))

		part1 <- dt(y,df=mu) * dt(up1,df=mu+1) * sqrt((mu+1)/(1-rho12^2)/(mu+y^2))
		part2 <- pt(sqrt((mu+2)/(mu+1+up1^2))*(up2-rho*up1)/sqrt(1-rho^2),df=mu+2)
		return(part1*part2)
	}

	## Third derivative of the trivariate t CDF w.r.t the k-th component

	dT3_dx3 <- function(v,rho,mu,k){
		if(k==1){x=v[1];y=v[2];z=v[3];rho12=rho[1];rho13=rho[2];rho23=rho[3]}
		if(k==2){x=v[2];y=v[1];z=v[3];rho12=rho[1];rho13=rho[3];rho23=rho[2]}
		if(k==3){x=v[3];y=v[1];z=v[2];rho12=rho[2];rho13=rho[3];rho23=rho[1]}

		up1 <- (y-rho12*x)/sqrt(1-rho12^2) * sqrt((mu+1)/(mu+x^2))
		up2 <- (z-rho13*x)/sqrt(1-rho13^2) * sqrt((mu+1)/(mu+x^2))
		rho <- (rho23-rho12*rho13)/sqrt((1-rho12^2)*(1-rho13^2))

		part1 <- dt(x,df=mu)/(mu+x^2)
	
		part21 <- pmvt(lower=rep(-Inf,2),upper=c(up1,up2),sigma=matrix(c(1,rho,rho,1),ncol=2),df=round(mu)+1) 
		part22 <- ((mu+1)+2)/(mu+x^2)*(mu+1)*x^2 - (mu+1)
		part2 <- part21*part22
	
		part31 <- dt(up1,df=mu+1) * pt(sqrt((mu+2)/(mu+1+up1^2))*(up2-rho*up1)/sqrt(1-rho^2),df=mu+2)
		part32 <- 2*(mu+2)*x/(mu+x^2) * (rho12*sqrt(mu+x^2)+(y-rho12*x)*x/sqrt(mu+x^2)) * sqrt((mu+1)/(1-rho12^2))
		part33 <- (-(mu+2)*up1/((mu+1)+up1^2)) * (mu+1)/(1-rho12^2) / (mu+x^2) * (rho12*sqrt(mu+x^2)+(y-rho12*x)*x/sqrt(mu+x^2))^2
		part34 <- (y-rho12*x)*(-mu)/(mu+x^2)^(3/2) * sqrt((mu+1)/(1-rho12^2))
		part3 <- part31*(part32+part33+part34)
	
		part41 <- dt(up2,df=mu+1) * pt(sqrt((mu+2)/(mu+1+up2^2))*(up1-rho*up2)/sqrt(1-rho^2),df=mu+2)
		part42 <- 2*(mu+2)*x/(mu+x^2) * (rho13*sqrt(mu+x^2)+(z-rho13*x)*x/sqrt(mu+x^2)) * sqrt((mu+1)/(1-rho13^2))
		part43 <- (-(mu+2)*up2/((mu+1)+up2^2)) * (mu+1)/(1-rho13^2) / (mu+x^2) * (rho13*sqrt(mu+x^2)+(z-rho13*x)*x/sqrt(mu+x^2))^2
		part44 <- (z-rho13*x)*(-mu)/(mu+x^2)^(3/2) * sqrt((mu+1)/(1-rho13^2))
		part4 <- part41*(part42+part43+part44)
	
		part50 <- sqrt((mu+1)/(1-rho12^2))*(-rho12*sqrt(mu+x^2)-(y-rho12*x)*x/sqrt(mu+x^2))
		part51 <- dt(up1,df=mu+1) * dt(sqrt((mu+2)/(mu+1+up1^2))*(up2-rho*up1)/sqrt(1-rho^2),df=mu+2)
		part52 <- (1-rho12^2)*sqrt(mu+2)/sqrt((1-rho12^2)*(1-rho13^2)-(rho23-rho12*rho13)^2)/((1-rho12^2)*(mu+x^2)+(y-rho12*x)^2)^(3/2)
		part53 <- (-rho13+rho12*(rho23-rho12*rho13)/(1-rho12^2)) * ((1-rho12^2)*(mu+x^2)+(y-rho12*x)^2)
		part54 <- ((z-rho13*x)-(rho23-rho12*rho13)/(1-rho12^2)*(y-rho12*x)) * (x*(1-rho12^2)-rho12*(y-rho12*x))	
		part5 <- part50*part51*part52*(part53-part54)
	
		part60 <- sqrt((mu+1)/(1-rho13^2))*(-rho13*sqrt(mu+x^2)-(z-rho13*x)*x/sqrt(mu+x^2))	
		part61 <- dt(up2,df=mu+1) * dt(sqrt((mu+2)/(mu+1+up2^2))*(up1-rho*up2)/sqrt(1-rho^2),df=mu+2)
		part62 <- (1-rho13^2)*sqrt(mu+2)/sqrt((1-rho12^2)*(1-rho13^2)-(rho23-rho12*rho13)^2)/((1-rho13^2)*(mu+x^2)+(z-rho13*x)^2)^(3/2)
		part63 <- (-rho12+rho13*(rho23-rho12*rho13)/(1-rho13^2)) * ((1-rho13^2)*(mu+x^2)+(z-rho13*x)^2)
		part64 <- ((y-rho12*x)-(rho23-rho12*rho13)/(1-rho13^2)*(z-rho13*x)) * (x*(1-rho13^2)-rho13*(z-rho13*x))	
		part6 <- part60*part61*part62*(part63-part64)
	
		
		return(part1*(part2+part3+part4+part5+part6)) 	
	}

	## Third derivative of the trivariate t CDF w.r.t the k (twice) and j-th (once) components

	dT3_dx2dy <- function(v,rho,mu,k,j){
		if(k==1 && j==2){x=v[1];y=v[2];z=v[3];rho12=rho[1];rho13=rho[2];rho23=rho[3]}
		if(k==1 && j==3){x=v[1];y=v[3];z=v[2];rho12=rho[2];rho13=rho[1];rho23=rho[3]}
		if(k==2 && j==1){x=v[2];y=v[1];z=v[3];rho12=rho[1];rho13=rho[3];rho23=rho[2]}
		if(k==2 && j==3){x=v[2];y=v[3];z=v[1];rho12=rho[3];rho13=rho[1];rho23=rho[2]}
		if(k==3 && j==1){x=v[3];y=v[1];z=v[2];rho12=rho[2];rho13=rho[3];rho23=rho[1]}
		if(k==3 && j==2){x=v[3];y=v[2];z=v[1];rho12=rho[3];rho13=rho[2];rho23=rho[1]}	

		up1 <- (x-rho12*y)/sqrt(1-rho12^2) * sqrt((mu+1)/(mu+y^2))
		up2 <- (z-rho23*y)/sqrt(1-rho23^2) * sqrt((mu+1)/(mu+y^2))
		rho <- (rho13-rho12*rho23)/sqrt((1-rho12^2)*(1-rho23^2))

		part1 <- dt(y,df=mu) * sqrt((mu+1)/(1-rho12^2)/(mu+y^2)) 
		part21 <- dt(up1,df=mu+1) * pt(sqrt((mu+2)/(mu+1+up1^2))*(up2-rho*up1)/sqrt(1-rho^2),df=mu+2)
		part22 <- -((mu+2)*up1)/(mu+1+up1^2) * sqrt((mu+1)/(1-rho12^2)/(mu+y^2))
		part31 <- dt(up1,df=mu+1) * dt(sqrt((mu+2)/(mu+1+up1^2))*(up2-rho*up1)/sqrt(1-rho^2),df=mu+2)
		part32 <- sqrt(mu+2)*(1-rho12^2)/sqrt((1-rho12^2)*(1-rho23^2)-(rho13-rho12*rho23)^2)
		part33 <- -((rho13-rho12*rho23)*(mu+y^2)+(z-rho23*y)*(x-rho12*y))/((1-rho12^2)*(mu+y^2)+(x-rho12*y)^2)^(3/2)
		return(part1*(part21*part22+part31*part32*part33))
	}


	######################################################################################
	# 	Definition of the functions x,y,w and their derivatives
	######################################################################################


		###########################################
		#										  #	
		#   In all of the following functions :   #
		#										  #	
		#	z is of size 4						  #
		#	rho is of size 6 (4dim = 6corr coeff) #
		#	mu is the degree of freedom			  #
		#										  #	
		###########################################




	#function corresponding to the j-th component of I_k
	# Gives x1,y1,w1, x2,y2,w2, ..., x4,y4,w4


	comp <- function(z1,z2,z3,z4,rho,mu,j,k){
		#j = 1,2 or 3
		#k = 1,2,3 or 4
	
		if(j==1 & k==1){corr=rho[1];x=z2;y=z1}
		if(j==1 & k==2){corr=rho[1];x=z1;y=z2}	
		if(j==1 & k==3){corr=rho[2];x=z1;y=z3}
		if(j==1 & k==4){corr=rho[3];x=z1;y=z4}	

		if(j==2 & k==1){corr=rho[2];x=z3;y=z1}
		if(j==2 & k==2){corr=rho[4];x=z3;y=z2}	
		if(j==2 & k==3){corr=rho[4];x=z2;y=z3}
		if(j==2 & k==4){corr=rho[6];x=z2;y=z4}	

		if(j==3 & k==1){corr=rho[3];x=z4;y=z1}
		if(j==3 & k==2){corr=rho[5];x=z4;y=z2}	
		if(j==3 & k==3){corr=rho[6];x=z4;y=z3}
		if(j==3 & k==4){corr=rho[6];x=z3;y=z4}	

	
		return( sqrt((mu+1)/(1-corr^2))*((x/y)^(1/mu)-corr) )
	}

	# first derivative of the j-th component of I_k w.r.t the s-th component of z (z_s)

	d_comp <- function(z1,z2,z3,z4,rho,mu,j,k,s){
		#j = 1,2 or 3
		#k = 1,2,3 or 4
		#ind = 1,2,3 or 4	
	
		if(s==1){return(jacobian(function(i){comp(i,z2,z3,z4,rho,mu,j,k)},z1))}
		if(s==2){return(jacobian(function(i){comp(z1,i,z3,z4,rho,mu,j,k)},z2))}
		if(s==3){return(jacobian(function(i){comp(z1,z2,i,z4,rho,mu,j,k)},z3))}
		if(s==4){return(jacobian(function(i){comp(z1,z2,z3,i,rho,mu,j,k)},z4))}	
	}


	# Second derivative of the j-th component of I_k w.r.t the s-th and t-th component of z (z_s and z_t)

	d2_comp_dsdt <- function(z1,z2,z3,z4,rho,mu,j,k,s,t){
		# j = 1,2 or 3
		# k = 1,2,3 or 4
		# s = 1,2,3 or 4	
		# t = 1,2,3 or 4
		# s != t	
	
		so <- sort(c(s,t))
		s<- so[1]; t<- so[2];
		
		if(s==1 & t==2){return(hessian(function(i){comp(i[1],i[2],z3,z4,rho,mu,j,k)},c(z1,z2))[1,2])}
		if(s==1 & t==3){return(hessian(function(i){comp(i[1],z2,i[2],z4,rho,mu,j,k)},c(z1,z3))[1,2])}
		if(s==1 & t==4){return(hessian(function(i){comp(i[1],z2,z3,i[2],rho,mu,j,k)},c(z1,z4))[1,2])}
	
		if(s==2 & t==3){return(hessian(function(i){comp(z1,i[1],i[2],z4,rho,mu,j,k)},c(z2,z3))[1,2])}
		if(s==2 & t==4){return(hessian(function(i){comp(z1,i[1],z3,i[2],rho,mu,j,k)},c(z2,z4))[1,2])}
		if(s==3 & t==4){return(hessian(function(i){comp(z1,z2,i[1],i[2],rho,mu,j,k)},c(z3,z4))[1,2])}
	
	}



	######################################################################################
	# 	Definition of the I_k and their derivatives
	######################################################################################


		###########################################
		#										  #	
		#   In all of the following functions :   #
		#										  #	
		#	z is of size 4						  #
		#	rho is of size 6 (4dim = 6corr coeff) #
		#	mu is the degree of freedom			  #
		#										  #	
		###########################################


	# Correlation matrix of I_k

	R_k <- function(rho,k,matrix){
	
		rho12=rho[1];rho13=rho[2];rho14=rho[3];rho23=rho[4];rho24=rho[5];rho34=rho[6];
	
		if(k==1){
			R12=(rho23-rho12*rho13)/sqrt((1-rho12^2)*(1-rho13^2));
			R13=(rho24-rho12*rho14)/sqrt((1-rho12^2)*(1-rho14^2));
			R23=(rho34-rho13*rho14)/sqrt((1-rho13^2)*(1-rho14^2));
		}
	
		if(k==2){
			R12=(rho13-rho12*rho23)/sqrt((1-rho12^2)*(1-rho23^2));
			R13=(rho14-rho12*rho24)/sqrt((1-rho12^2)*(1-rho24^2));
			R23=(rho34-rho23*rho24)/sqrt((1-rho23^2)*(1-rho24^2));
		}
	
		if(k==3){
			R12=(rho12-rho13*rho23)/sqrt((1-rho13^2)*(1-rho23^2));
			R13=(rho14-rho13*rho34)/sqrt((1-rho13^2)*(1-rho34^2));
			R23=(rho24-rho23*rho34)/sqrt((1-rho23^2)*(1-rho34^2));
		}
	
		if(k==4){
			R12=(rho12-rho14*rho24)/sqrt((1-rho14^2)*(1-rho24^2));
			R13=(rho13-rho14*rho34)/sqrt((1-rho14^2)*(1-rho34^2));
			R23=(rho23-rho24*rho34)/sqrt((1-rho24^2)*(1-rho34^2));
		}	
		if(matrix==TRUE){return(matrix(c(1,R12,R13,R12,1,R23,R13,R23,1),ncol=3))} # return under matrix form
		if(matrix==FALSE){return(c(R12,R13,R23))}	
	}


	# Function I_k

	I_k <- function(z,rho,mu,k){
	
		up <- c(comp(z[1],z[2],z[3],z[4],rho,mu,1,k),comp(z[1],z[2],z[3],z[4],rho,mu,2,k),comp(z[1],z[2],z[3],z[4],rho,mu,3,k))
		return(pmvt(lower=rep(-Inf,3),upper=up,sigma=R_k(rho,k,TRUE),df=round(mu)+1,algorithm=TVPACK)[1])
	
	}


	# First derivative of I_k w.r.t its s-th component (z_s)

	dI_k_ds <- function(z,rho,mu,k,s){
		# j = 1,2 or 3
		# k = 1,2,3 or 4
		# s = 1,2,3 or 4
	
		xk <- comp(z[1],z[2],z[3],z[4],rho,mu,1,k)
		yk <- comp(z[1],z[2],z[3],z[4],rho,mu,2,k)
		zk <- comp(z[1],z[2],z[3],z[4],rho,mu,3,k)

		part1 <- dT_dx(c(xk,yk,zk),R_k(rho,k,FALSE),mu+1,1) * d_comp(z[1],z[2],z[3],z[4],rho,mu,1,k,s)
		part2 <- dT_dx(c(xk,yk,zk),R_k(rho,k,FALSE),mu+1,2) * d_comp(z[1],z[2],z[3],z[4],rho,mu,2,k,s)
		part3 <- dT_dx(c(xk,yk,zk),R_k(rho,k,FALSE),mu+1,3) * d_comp(z[1],z[2],z[3],z[4],rho,mu,3,k,s)		

		return(part1+part2+part3)
	}

	# Second derivative of I_k w.r.t its s-th and t-th components (z_s and z_t)

	d2I_k_dsdt <- function(z,rho,mu,k,s,t){
		# j = 1,2 or 3
		# k = 1,2,3 or 4
		# s = 1,2,3 or 4
		# t = 1,2,3 or 4
		# s != t

		xk <- comp(z[1],z[2],z[3],z[4],rho,mu,1,k)
		yk <- comp(z[1],z[2],z[3],z[4],rho,mu,2,k)
		zk <- comp(z[1],z[2],z[3],z[4],rho,mu,3,k)
	
		part1 <- dT2_dx2(c(xk,yk,zk),R_k(rho,k,FALSE),mu+1,1) * d_comp(z[1],z[2],z[3],z[4],rho,mu,1,k,s) * d_comp(z[1],z[2],z[3],z[4],rho,mu,1,k,t)
		part2 <- dT2_dx2(c(xk,yk,zk),R_k(rho,k,FALSE),mu+1,2) * d_comp(z[1],z[2],z[3],z[4],rho,mu,2,k,s) * d_comp(z[1],z[2],z[3],z[4],rho,mu,2,k,t)
		part3 <- dT2_dx2(c(xk,yk,zk),R_k(rho,k,FALSE),mu+1,3) * d_comp(z[1],z[2],z[3],z[4],rho,mu,3,k,s) * d_comp(z[1],z[2],z[3],z[4],rho,mu,3,k,t)		
		part4 <- dT_dx(c(xk,yk,zk),R_k(rho,k,FALSE),mu+1,1) * d2_comp_dsdt(z[1],z[2],z[3],z[4],rho,mu,1,k,s,t) 
		part5 <- dT_dx(c(xk,yk,zk),R_k(rho,k,FALSE),mu+1,2) * d2_comp_dsdt(z[1],z[2],z[3],z[4],rho,mu,2,k,s,t) 
		part6 <- dT_dx(c(xk,yk,zk),R_k(rho,k,FALSE),mu+1,3) * d2_comp_dsdt(z[1],z[2],z[3],z[4],rho,mu,3,k,s,t)

		part71 <- dT2_dxdy(c(xk,yk,zk),R_k(rho,k,FALSE),mu+1,1,2)
		part72 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,1,k,t) * d_comp(z[1],z[2],z[3],z[4],rho,mu,2,k,s)
		part73 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,1,k,s) * d_comp(z[1],z[2],z[3],z[4],rho,mu,2,k,t)

		part81 <- dT2_dxdy(c(xk,yk,zk),R_k(rho,k,FALSE),mu+1,1,3)
		part82 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,1,k,t) * d_comp(z[1],z[2],z[3],z[4],rho,mu,3,k,s)
		part83 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,1,k,s) * d_comp(z[1],z[2],z[3],z[4],rho,mu,3,k,t)
	
		part91 <- dT2_dxdy(c(xk,yk,zk),R_k(rho,k,FALSE),mu+1,2,3)
		part92 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,2,k,t) * d_comp(z[1],z[2],z[3],z[4],rho,mu,3,k,s)
		part93 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,2,k,s) * d_comp(z[1],z[2],z[3],z[4],rho,mu,3,k,t)
	
	
		return(part1+part2+part3+part4+part5+part6+part71*(part72+part73)+part81*(part82+part83)+part91*(part92+part93))
	}	

	# Third derivative of I_k w.r.t its s-th, t-th and u-th components (z_s, z_t and z_u)

	d3I_k_dsdtdu <- function(z,rho,mu,k,s,t,u){
		# j = 1,2 or 3
		# k = 1,2,3 or 4
		# s = 1,2,3 or 4
		# t = 1,2,3 or 4
		# u = 1,2,3 or 4
		# s != t != u
	
		xk <- comp(z[1],z[2],z[3],z[4],rho,mu,1,k)
		yk <- comp(z[1],z[2],z[3],z[4],rho,mu,2,k)
		zk <- comp(z[1],z[2],z[3],z[4],rho,mu,3,k)

		part10 <- dT2_dxdy(c(xk,yk,zk),R_k(rho,k,FALSE),mu+1,1,2)
		part11 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,2,k,s) * d2_comp_dsdt(z[1],z[2],z[3],z[4],rho,mu,1,k,u,t)
		part12 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,1,k,s) * d2_comp_dsdt(z[1],z[2],z[3],z[4],rho,mu,2,k,u,t)
		part13 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,2,k,t) * d2_comp_dsdt(z[1],z[2],z[3],z[4],rho,mu,1,k,u,s)
		part14 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,1,k,t) * d2_comp_dsdt(z[1],z[2],z[3],z[4],rho,mu,2,k,u,s)
		part15 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,2,k,u) * d2_comp_dsdt(z[1],z[2],z[3],z[4],rho,mu,1,k,s,t)
		part16 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,1,k,u) * d2_comp_dsdt(z[1],z[2],z[3],z[4],rho,mu,2,k,s,t)
		part1 <- part10*(part11+part12+part13+part14+part15+part16)
	
		part20 <- dT2_dxdy(c(xk,yk,zk),R_k(rho,k,FALSE),mu+1,1,3)
		part21 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,3,k,s) * d2_comp_dsdt(z[1],z[2],z[3],z[4],rho,mu,1,k,u,t)
		part22 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,1,k,s) * d2_comp_dsdt(z[1],z[2],z[3],z[4],rho,mu,3,k,u,t)
		part23 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,3,k,t) * d2_comp_dsdt(z[1],z[2],z[3],z[4],rho,mu,1,k,u,s)
		part24 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,1,k,t) * d2_comp_dsdt(z[1],z[2],z[3],z[4],rho,mu,3,k,u,s)
		part25 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,3,k,u) * d2_comp_dsdt(z[1],z[2],z[3],z[4],rho,mu,1,k,s,t)
		part26 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,1,k,u) * d2_comp_dsdt(z[1],z[2],z[3],z[4],rho,mu,3,k,s,t)
		part2 <- part20*(part21+part22+part23+part24+part25+part26)
	 
		part30 <- dT2_dxdy(c(xk,yk,zk),R_k(rho,k,FALSE),mu+1,2,3)
		part31 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,3,k,s) * d2_comp_dsdt(z[1],z[2],z[3],z[4],rho,mu,2,k,u,t)
		part32 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,2,k,s) * d2_comp_dsdt(z[1],z[2],z[3],z[4],rho,mu,3,k,u,t)
		part33 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,3,k,t) * d2_comp_dsdt(z[1],z[2],z[3],z[4],rho,mu,2,k,u,s)
		part34 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,2,k,t) * d2_comp_dsdt(z[1],z[2],z[3],z[4],rho,mu,3,k,u,s)
		part35 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,3,k,u) * d2_comp_dsdt(z[1],z[2],z[3],z[4],rho,mu,2,k,s,t)
		part36 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,2,k,u) * d2_comp_dsdt(z[1],z[2],z[3],z[4],rho,mu,3,k,s,t)
		part3 <- part30*(part31+part32+part33+part34+part35+part36)	
	
		part40 <- dT3_dx2dy(c(xk,yk,zk),R_k(rho,k,FALSE),mu+1,1,2)
		part41 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,2,k,s) * d_comp(z[1],z[2],z[3],z[4],rho,mu,1,k,t) * d_comp(z[1],z[2],z[3],z[4],rho,mu,1,k,u)
		part42 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,1,k,s) * d_comp(z[1],z[2],z[3],z[4],rho,mu,2,k,t) * d_comp(z[1],z[2],z[3],z[4],rho,mu,1,k,u)
		part43 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,1,k,s) * d_comp(z[1],z[2],z[3],z[4],rho,mu,1,k,t) * d_comp(z[1],z[2],z[3],z[4],rho,mu,2,k,u)
		part4 <- part40*(part41+part42+part43)

		part50 <- dT3_dx2dy(c(xk,yk,zk),R_k(rho,k,FALSE),mu+1,1,3)
		part51 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,3,k,s) * d_comp(z[1],z[2],z[3],z[4],rho,mu,1,k,t) * d_comp(z[1],z[2],z[3],z[4],rho,mu,1,k,u)
		part52 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,1,k,s) * d_comp(z[1],z[2],z[3],z[4],rho,mu,3,k,t) * d_comp(z[1],z[2],z[3],z[4],rho,mu,1,k,u)
		part53 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,1,k,s) * d_comp(z[1],z[2],z[3],z[4],rho,mu,1,k,t) * d_comp(z[1],z[2],z[3],z[4],rho,mu,3,k,u)
		part5 <- part50*(part51+part52+part53)
	
		part60 <- dT3_dx2dy(c(xk,yk,zk),R_k(rho,k,FALSE),mu+1,2,1)
		part61 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,1,k,s) * d_comp(z[1],z[2],z[3],z[4],rho,mu,2,k,t) * d_comp(z[1],z[2],z[3],z[4],rho,mu,2,k,u)
		part62 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,2,k,s) * d_comp(z[1],z[2],z[3],z[4],rho,mu,1,k,t) * d_comp(z[1],z[2],z[3],z[4],rho,mu,2,k,u)
		part63 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,2,k,s) * d_comp(z[1],z[2],z[3],z[4],rho,mu,2,k,t) * d_comp(z[1],z[2],z[3],z[4],rho,mu,1,k,u)
		part6 <- part60*(part61+part62+part63)

		part70 <- dT3_dx2dy(c(xk,yk,zk),R_k(rho,k,FALSE),mu+1,2,3)
		part71 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,3,k,s) * d_comp(z[1],z[2],z[3],z[4],rho,mu,2,k,t) * d_comp(z[1],z[2],z[3],z[4],rho,mu,2,k,u)
		part72 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,2,k,s) * d_comp(z[1],z[2],z[3],z[4],rho,mu,3,k,t) * d_comp(z[1],z[2],z[3],z[4],rho,mu,2,k,u)
		part73 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,2,k,s) * d_comp(z[1],z[2],z[3],z[4],rho,mu,2,k,t) * d_comp(z[1],z[2],z[3],z[4],rho,mu,3,k,u)
		part7 <- part70*(part71+part72+part73)
	
		part80 <- dT3_dx2dy(c(xk,yk,zk),R_k(rho,k,FALSE),mu+1,3,1)
		part81 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,1,k,s) * d_comp(z[1],z[2],z[3],z[4],rho,mu,3,k,t) * d_comp(z[1],z[2],z[3],z[4],rho,mu,3,k,u)
		part82 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,3,k,s) * d_comp(z[1],z[2],z[3],z[4],rho,mu,1,k,t) * d_comp(z[1],z[2],z[3],z[4],rho,mu,3,k,u)
		part83 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,3,k,s) * d_comp(z[1],z[2],z[3],z[4],rho,mu,3,k,t) * d_comp(z[1],z[2],z[3],z[4],rho,mu,1,k,u)
		part8 <- part80*(part81+part82+part83)

		part90 <- dT3_dx2dy(c(xk,yk,zk),R_k(rho,k,FALSE),mu+1,3,2)
		part91 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,2,k,s) * d_comp(z[1],z[2],z[3],z[4],rho,mu,3,k,t) * d_comp(z[1],z[2],z[3],z[4],rho,mu,3,k,u)
		part92 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,3,k,s) * d_comp(z[1],z[2],z[3],z[4],rho,mu,2,k,t) * d_comp(z[1],z[2],z[3],z[4],rho,mu,3,k,u)
		part93 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,3,k,s) * d_comp(z[1],z[2],z[3],z[4],rho,mu,3,k,t) * d_comp(z[1],z[2],z[3],z[4],rho,mu,2,k,u)
		part9 <- part90*(part91+part92+part93)	
	
		part100 <- dmvt(c(xk,yk,zk),sigma=R_k(rho,k,TRUE),df=mu+1,log=FALSE)
		part1001 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,1,k,s) * d_comp(z[1],z[2],z[3],z[4],rho,mu,2,k,t) * d_comp(z[1],z[2],z[3],z[4],rho,mu,3,k,u)			
		part1002 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,1,k,s) * d_comp(z[1],z[2],z[3],z[4],rho,mu,3,k,t) * d_comp(z[1],z[2],z[3],z[4],rho,mu,2,k,u)			
		part1003 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,2,k,s) * d_comp(z[1],z[2],z[3],z[4],rho,mu,1,k,t) * d_comp(z[1],z[2],z[3],z[4],rho,mu,3,k,u)			
		part1004 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,2,k,s) * d_comp(z[1],z[2],z[3],z[4],rho,mu,3,k,t) * d_comp(z[1],z[2],z[3],z[4],rho,mu,1,k,u)			
		part1005 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,3,k,s) * d_comp(z[1],z[2],z[3],z[4],rho,mu,1,k,t) * d_comp(z[1],z[2],z[3],z[4],rho,mu,2,k,u)			
		part1006 <- d_comp(z[1],z[2],z[3],z[4],rho,mu,3,k,s) * d_comp(z[1],z[2],z[3],z[4],rho,mu,2,k,t) * d_comp(z[1],z[2],z[3],z[4],rho,mu,1,k,u)			
		part <- part100*(part1001+part1002+part1003+part1004+part1005+part1006)
	
		return(part1+part2+part3+part4+part5+part6+part7+part8+part9+part)
	}

	######################################################################################
	# 	Definition of the V and its derivatives
	######################################################################################


		###########################################
		#										  #	
		#   In all of the following functions :   #
		#										  #	
		#	z is of size 4						  #
		#	rho is of size 6 (4dim = 6corr coeff) #
		#	mu is the degree of freedom			  #
		#										  #	
		###########################################


	# V function (4 dimensions)


	V <- function(z,rho,mu){
  
	  # check if the 4 correlation matrices are all positive definite (and correlation coefficient between -1 and 1)
	  M1 <- R_k(rho,1,TRUE); M2 <- R_k(rho,2,TRUE); M3 <- R_k(rho,3,TRUE); M4 <- R_k(rho,4,TRUE);
	  if(sum(eigen(M1)$values<0)>=1 || sum(eigen(M2)$values<0)>=1 || sum(eigen(M3)$values<0)>=1 || sum(eigen(M4)$values<0)>=1 ){return(1e-50)}
 	 if(any(abs(M1)>1) || any(abs(M2)>1) || any(abs(M3)>1) || any(abs(M4)>1)){return(1e-50)}
  
		return(I_k(z,rho,mu,1)/z[1] + I_k(z,rho,mu,2)/z[2] + I_k(z,rho,mu,3)/z[3] + I_k(z,rho,mu,4)/z[4])
	}


	# first derivative of V with respect to its s-th component (z_s)

	dV_ds <- function(z,rho,mu,s){
	
		# s= 1,2,3 or 4

 	 # check if the 4 correlation matrices are all positive definite (and correlation coefficient between -1 and 1)
	  M1 <- R_k(rho,1,TRUE); M2 <- R_k(rho,2,TRUE); M3 <- R_k(rho,3,TRUE); M4 <- R_k(rho,4,TRUE);
	  if(sum(eigen(M1)$values<0)>=1 || sum(eigen(M2)$values<0)>=1 || sum(eigen(M3)$values<0)>=1 || sum(eigen(M4)$values<0)>=1 ){return(1e-50)}
	  if(any(abs(M1)>1) || any(abs(M2)>1) || any(abs(M3)>1) || any(abs(M4)>1)){return(1e-50)}
    
		part1 <- -I_k(z,rho,mu,s)/z[s]^2 
		part2 <- dI_k_ds(z,rho,mu,1,s)/z[1] + dI_k_ds(z,rho,mu,2,s)/z[2] + dI_k_ds(z,rho,mu,3,s)/z[3] + dI_k_ds(z,rho,mu,4,s)/z[4]
		return(part1+part2)
	}

	# second dereivative of V with respect to its s-th and t-th components (z_s and z_t)

	d2V_dsdt <- function(z,rho,mu,s,t){
	
		# s= 1,2,3 or 4
		# t= 1,2,3 or 4
		# s != t

 	 # check if the 4 correlation matrices are all positive definite (and correlation coefficient between -1 and 1)
	  M1 <- R_k(rho,1,TRUE); M2 <- R_k(rho,2,TRUE); M3 <- R_k(rho,3,TRUE); M4 <- R_k(rho,4,TRUE);
	  if(sum(eigen(M1)$values<0)>=1 || sum(eigen(M2)$values<0)>=1 || sum(eigen(M3)$values<0)>=1 || sum(eigen(M4)$values<0)>=1 ){return(1e-50)}
	  if(any(abs(M1)>1) || any(abs(M2)>1) || any(abs(M3)>1) || any(abs(M4)>1)){return(1e-50)}
  
		part1 <- -dI_k_ds(z,rho,mu,s,t)/z[s]^2 -dI_k_ds(z,rho,mu,t,s)/z[t]^2 
		part2 <- d2I_k_dsdt(z,rho,mu,1,s,t)/z[1] + d2I_k_dsdt(z,rho,mu,2,s,t)/z[2] + d2I_k_dsdt(z,rho,mu,3,s,t)/z[3] + d2I_k_dsdt(z,rho,mu,4,s,t)/z[4]

		return(part1+part2)	
		
	}	

	# third derivative of V with respect to its s-th, t-th and u-th components (z_s, z_t and z_u)

	d3V_dsdtdu <- function(z,rho,mu,s,t,u){
	
		# s= 1,2,3 or 4
		# t= 1,2,3 or 4
		# u= 1,2,3 or 4
		# s != t != u
	
	  # check if the 4 correlation matrices are all positive definite (and correlation coefficient between -1 and 1)
	  M1 <- R_k(rho,1,TRUE); M2 <- R_k(rho,2,TRUE); M3 <- R_k(rho,3,TRUE); M4 <- R_k(rho,4,TRUE);
	  if(sum(eigen(M1)$values<0)>=1 || sum(eigen(M2)$values<0)>=1 || sum(eigen(M3)$values<0)>=1 || sum(eigen(M4)$values<0)>=1 ){return(1e-50)}
	  if(any(abs(M1)>1) || any(abs(M2)>1) || any(abs(M3)>1) || any(abs(M4)>1)){return(1e-50)}
  
		part1 <- -d2I_k_dsdt(z,rho,mu,s,t,u)/z[s]^2 -d2I_k_dsdt(z,rho,mu,t,s,u)/z[t]^2 -d2I_k_dsdt(z,rho,mu,u,s,t)/z[u]^2
		part2 <- d3I_k_dsdtdu(z,rho,mu,1,s,t,u)/z[1] + d3I_k_dsdtdu(z,rho,mu,2,s,t,u)/z[2] + d3I_k_dsdtdu(z,rho,mu,3,s,t,u)/z[3] + d3I_k_dsdtdu(z,rho,mu,4,s,t,u)/z[4]
	
		return(part1+part2)
	
	}	



	dens_et_4d <- function(w,rho,mu){
		p <- 4; # 4-variate data
		k=p-1;
		if(sum(abs(rho)>=1)){return(1e-50)}
		w.tilde<-rep(0,k);
		Sigma<-diag(k);
			
		if(sum(w < 0.03) == 3){ # in one of the corners of the 4-d simplex
		  ind.lim <-which(w <0.03)
		  ind <- c(1:4)[-which(w<0.03)]
		  z <- vector(length=4)
		  z[ind.lim] <- 1e-20 # In order to consider the limit when these components are going to zero
		  z[ind] <- w[ind] 
			return( -sum(z[ind])^(length(ind)+1) * dV_ds(z,rho,mu,ind) )		
		} else if(sum(w < 0.02) == 2){ # in one of the edges of the 4-d simplex
		    ind.lim <-which(w <0.02)
		    ind <- c(1:4)[-which(w<0.02)]
		    z <- vector(length=4)
		    z[ind.lim] <- 1e-20 # In order to consider the limit when these components are going to zero
		    z[ind] <- w[ind] 
	      return( -sum(z[ind])^(length(ind)+1) * d2V_dsdt(z,rho,mu,ind[1],ind[2]) )
		}else	if(sum(w < 0.01) == 1){ # in one of the faces of the 4-d simplex
		    ind.lim <-which(w <0.01)
		    ind <- c(1:4)[-which(w<0.01)]
		    z <- vector(length=4)
		    z[ind.lim] <- 1e-20 # In order to consider the limit when these components are going to zero
		    z[ind] <- w[ind] 
	      return( -sum(z[ind])^(length(ind)+1) * d3V_dsdtdu(z,rho,mu,ind[1],ind[2],ind[3]) )		
		}else{ # in the interior of the 4-d simplex
			return(interior_et_d(w,rho,mu))
		}	
	}


### Angular density for the Extremal-t model on the 2, 3 and 4 dimensional simplex

	xvect = as.double(as.vector(t(x)))
    if (is.vector(x)) {
        dim = as.integer(length(x))
        n = as.integer(1)
        if(round(sum(x),7) != 1){ stop("Data is not angular")}
        if(dim==2){ result <- dens_et_2d(x,rho,mu)} 
    	if(dim==3){ result <- dens_et_3d(x,rho,mu)} 
    	if(dim==4){ result <- dens_et_4d(x,rho,mu)} 
    }
    else {
        dim = as.integer(ncol(x))
        n = as.integer(nrow(x))
    	
		if (sum(apply(x,1,sum)) != n){ stop("Data is not angular") }
    
    	if (vectorial) {
    	    result = double(n)
    	    if(dim==2){ result <- apply(x,1,function(y){dens_et_2d(y,rho,mu)}) }
    	    if(dim==3){ result <- apply(x,1,function(y){dens_et_3d(y,rho,mu)}) }
    	    if(dim==4){ result <- apply(x,1,function(y){dens_et_4d(y,rho,mu)}) }
    	} else { # vectorial=FALSE mean we return the likelihood function
    	    result = as.double(1)
    	    if(dim==2){ result <- prod(apply(x,1,function(y){dens_et_2d(y,rho,mu)})) }
    	    if(dim==3){ result <- prod(apply(x,1,function(y){dens_et_3d(y,rho,mu)})) }
    	    if(dim==4){ result <- prod(apply(x,1,function(y){dens_et_4d(y,rho,mu)})) }
   		}
   	}	
    if(log)
    	return(log(result))
    else return(result)
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
	return(dens_et(x=x, rho=par[1:choose(d,2)], mu=par[choose(d,2)+1], log=log,  vectorial=vectorial))
}
if(model=='Asymmetric'){
	if(d==2){
		if(length(par)!=3 ){
			stop('Wrong length of parameters')
		}else{
			return(dens_al(x=x, alpha=par[1], beta=par[2:3], log=log,  vectorial=vectorial))
		}
	}
	if(d==3){
		if(length(par)!=13 ){
			stop('Wrong length of parameters')
		}else{
			return(dens_al(x=x, alpha=par[1:4], beta=par[5:13], log=log,  vectorial=vectorial))
		}
	}
	if(d==4){
		if(length(par)!=39 ){
			stop('Wrong length of parameters')
		}else{
			return(dens_al(x=x, alpha=par[1:11], beta=par[12:39], log=log,  vectorial=vectorial))
		}
	}
}
	
	
}



