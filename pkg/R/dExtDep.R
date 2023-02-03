dExtDep <- function(x, method="Parametric", model, par, angular=TRUE, log=FALSE, c=NULL, vectorial=TRUE, 
                    mixture=FALSE){
  
  # Checking method
  methods <- c("Parametric", "NonParametric")
  if(!any(method ==  methods)){ stop("Wrong method specified")}
  
  if(method == "Parametric"){
    
    if(angular){
      models <- c("PB", "HR", "ET", "EST", "TD", "AL")
      if(!any(model ==  models)){ stop("model wrongly specified")}
      
      # Dimension of the model
      if (is.vector(x)) { 
        d <- as.integer(length(x))
        if(round(sum(x),7)!= 1){stop("'x' should be a vector in the unit simplex")}
      }else if(is.matrix(x)){
        d <- as.integer(ncol(x)) 
        if(!all(round(rowSums(x),7)==1)){stop("'x' should be a matrix with row vectors belonging to the unit simplex")}
      }else{
        stop(" 'x' should be a vector or a matrix")
      }
      
      if(model %in% models[-c(1,2,5)]){
        if(d>3){stop("The angular density is only available for dimension 2 or 3.")}
      }
      
      dim.check <- dim_ExtDep(model=model, par=par, dim=d)
      if(!dim.check){stop('Wrong length of parameters')}
      
      if(model=='PB'){
        return(dens_pb(x=x, b=par[1:choose(d,2)], alpha=par[choose(d,2)+1], log=log,  vectorial=vectorial))
      }
      if(model=='HR'){
        return(dens_hr(x=x, lambda=par, log=log,  vectorial=vectorial))
      }
      if(model=='TD'){
        return(dens_di(x=x, para=par, log=log,  vectorial=vectorial))
      }
      if(model=='ET'){
        if(is.null(c)){stop('c needs to be specified')}
        return(dens_et(x=x, rho=par[1:choose(d,2)], mu=par[choose(d,2)+1], c=c, log=log,  vectorial=vectorial))
      }
      if(model=='EST'){
        if(is.null(c)){stop('c needs to be specified')}
        return(dens_est(x=x, rho=par[1:choose(d,2)], alpha=par[choose(d,2)+1:d],mu=par[choose(d,2)+d+1], c=c, log=log,  vectorial=vectorial))
      }
      if(model=='AL'){
        if(is.null(c)){stop('c needs to be specified')}
        if(d==2){
          return(dens_al(x=x, alpha=par[1], beta=par[2:3], c=c, log=log,  vectorial=vectorial))
        }
        if(d==3){
          return(dens_al(x=x, alpha=par[1:4], beta=par[5:13], c=c, log=log,  vectorial=vectorial))
        }
      }
      
    }else{ # IF NOT ANGULAR
      models <- c("HR", "ET", "EST")
      if(!any(model ==  models)){ stop("model wrongly specified")}
      
      # Dimension of the model
      if (is.vector(x)) { 
        d <- as.integer(length(x))
      }else if(is.matrix(x)){
        d <- as.integer(ncol(x)) 
      }else{
        stop(" 'x' should be a vector or a matrix")
      }
      
      if(d!=2){stop("Density of parametric models only available for d=2")}
      
      dim.check <- dim_ExtDep(model=model, par=par, dim=d)
      if(!dim.check){stop('Wrong length of parameters')}
      
      if(model=="HR"){ # Husler-Reiss model
        
        if(is.vector(x)){
          out <- dHuslerReiss(x=x, lambda=par)
        }else if(is.matrix(x)){
          out <- apply(x, 1, function(z) dHuslerReiss(x=z, lambda=par))
        }
        
      } #END HR Model
      
      if(model=="ET"){ # Extremal-t model
        
        if(is.vector(x)){
          out <- dmextst(x=x, scale=par[1], df=par[2])
        }else if(is.matrix(x)){
          out <- apply(x, 1, function(z) dmextst(x=z, scale=par[1], df=par[2]))
        }
        
      } #END ET Model
      
      if(model=="EST"){ # Extremal skew-t model
        
        if(is.vector(x)){
          out <- dmextst(x=x, scale=par[1], shape=par[2:3], df=par[4])
        }else if(is.matrix(x)){
          out <- apply(x, 1, function(z) dmextst(x=z, scale=par[1], shape=par[2:3], df=par[4]))
        }
        
      } #END EST Model
      
      if(log){out <- log(out)}
      if(!vectorial){
        if(log){
          out <- sum(out)
        }else{
          out <- prod(out)
        }
      }
      
      return(out)
      
    } # END IF NOT ANGULAR
    
  }else if(method == "NonParametric"){
    
    if(angular){
      
      if (is.vector(x)) { 
        if(length(x)!=2){stop("'x' should of length 2")}
        if(sum(x)!= 1){stop("'x' should be a vector in the unit simplex")}
        dens <- dh(w=x[1], beta=par, mixture = mixture)
      }else if(is.matrix(x)){
        if(ncol(x)!=2){stop("'x' should have 2 columns")}
        if(!all(rowSums(x)==1)){stop("'x' should be a matrix with row vectors belonging to the unit simplex")}
        dens <- apply(x,1,function(y){dh(w=y[1], beta=par, mixture=mixture)})
      }else{
        stop(" 'x' should be a vector or a matrix")
      }
      
      if(log){dens <- log(dens)}
      if(!vectorial){
        if(log){
          dens <- sum(dens)
        }else{
          dens <- prod(dens)
        }
      }
      
      return(dens)
      
    }else{
      stop("Only the angular density available in non-parametric form")
    }
    
  }
  
  
}




####################################################################
####################################################################
####################################################################
### Internal functions for PARAMETRIC and ANGULAR
####################################################################
####################################################################
####################################################################

###
### Pairwise Beta model
###

dens_pb <- function (x, b, alpha, log, vectorial){
  
  if(any(b<=0) || alpha<=0){return(1e-300)}
  hij<-function(wi,wj,bij,alpha,p){
    wij = wi + wj;
    return(wij^(2*alpha-1)*(1-wij)^(alpha*(p-2)-p+2)*gamma(2*bij)/gamma(bij)^2*(wi/wij)^(bij-1)*(wj/wij)^(bij-1))
  }
  
  dens_pairb<-function(w,b,alpha){
    p<-length(w);
    if(p==2){
      dens = gamma(2*b) / gamma(b)^2 * prod(w)^(b-1)
    }else{
      dens = 0;
      K = 2 * factorial(p-3) * gamma(alpha*p+1) / (p * (p-1) * gamma(2*alpha+1) * gamma(alpha*(p-2)) )
      for(i in 1:(p-1)){
        for(j in (i+1):p){
          d = hij(w[i],w[j],b[(i-1)*p-sum(1:i)+j],alpha,p)
          dens = dens + d
        }
      }
      dens = dens * K
    }
    
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
    if (any(round(apply(x,1,sum),7) != 1)){ stop("Data is not angular")}
    
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

###
### Husler-Reiss model
###

dens_hr <- function (x, lambda, log, vectorial){
  
  dens_husler<-function(w,lambda){
    p<-length(w);
    k=p-1
    w.tilde<-rep(0,k);
    Sigma<-diag(k);
    
    for(i in 1:k){
      if(w[1]==0 | w[i+1]==0){return(1e-300)} else{
        w.tilde[i] <- log(w[i+1]/w[1])+2*lambda[i]^2;
        
        for(j in i:k){
          if(i==j){Sigma[i,j]=Sigma[j,i]=2*(lambda[i]^2+lambda[j]^2)}else{
            Sigma[i,j]=Sigma[j,i]=2*(lambda[i]^2+lambda[j]^2-lambda[sum(k:(k-i+1))+j-i]^2)
          }
        }
      }
    }
    
    if(sum(eigen(Sigma)$values<0)>=1){return(1e-300)} # Check if matrix is positive definite
    if(any(is.na(w.tilde))){return (1e-300)}
    
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
    if (any(round(apply(x,1,sum),7) != 1)){ stop("Data is not angular")}
    
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

###
### Tilted Dirichlet model
###

dens_di <- function (x, para, log, vectorial){
  
  dens_diri <- function(w,para){
    d <- length(para)
    
    if(any(para <=0)){return(1e-300)}
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
    if (any(round(apply(x,1,sum),7) != 1)){ stop("Data is not angular")}
    
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

###
### Asymmetric logistic model
###

dens_al <- function (x, alpha, beta, c, log, vectorial){
  
  # The 2d case:
  
  # in the following functions:
  #
  # alpha is a vector of size 1: for the subset {1,2}
  # beta is a vector of size 2: for [1,{1,2}] and [2,{1,2}]
  # beta for [1,{1}], [2,{2}] and [3,{3}] are omitted as obtained as [1,{1}] = 1 - [1,{1,2}] and [2,{2}] = 1 - [2,{1,2}]
  
  interior_alm_2d <- function(w,alpha,beta){
    part1 <- (alpha-1) * (beta[1]*beta[2])^alpha * (w[1]*w[2])^(alpha-2)
    part2 <- ( (beta[1]*w[2])^alpha + (beta[2]*w[1])^alpha )^(1/alpha-2)
    
    return(part1*part2)
  }
  
  corners_alm_2d <- function(alpha,beta,s){ # mass on the corner of the s-th component
    if(s==1){return(1-beta[1])}
    if(s==2){return(1-beta[2])}
  }
  
  dens_alm_2d <- function(data,alpha,beta,c){
    if(length(alpha)!=1){return(stop("Wrong length of parameter alpha"))}
    if(length(beta)!=2){return(stop("Wrong length of parameter beta"))}
    if( (alpha < 1) || any(beta < 0) || any(beta > 1)){return(1e-300)}
    
    ###
    ### Only one observation
    ###
    
    if(is.vector(data) && length(data)==2){
      
      if(c==0){
        hdens <- interior_alm_2d(w=data, alpha=alpha, beta=beta)
      }else if(c>0){
        subs <- subset.c(w=data, c=c)
        if(length(subs) == 1){ #corner
          K <- 1 / c
          hdens <- K * corners_alm_2d(alpha=alpha, beta=beta, s=subs)
        }else if (length(subs)==2){
          int01 <- 2 - corners_alm_2d(alpha=alpha, beta=beta, s=1) - corners_alm_2d(alpha=alpha, beta=beta, s=2)
          intc <- tryCatch(integrate(Vectorize(function(x) interior_alm_2d(c(x,1-x), alpha=alpha, beta=beta)), lower=c, upper=1-c )$value, error=function(e) -1)
          if(intc==-1){
            hdens <- 0
          }else{
            K.inter <- int01/intc
            hdens <- K.inter * interior_alm_2d(w=data, alpha=alpha, beta=beta)            
          }
        }
      }
      return(hdens)
    }
    
    ###
    ### Multiple observation in matrix form
    ###
    
    if(is.matrix(data) && ncol(data)==2){
      
      hdens <- vector(length=nrow(data))
      
      if(c == 0){
        hdens <- apply(data,1, function(x) interior_alm_3d(w=x, alpha=alpha, beta=beta))
      }else if(c>0){
        
        subs <- apply(data, 1, function(x) subset.c(x,c=c))
        
        ## Corners
        
        c1 <- which(lapply(subs, function(x) (length(x)==1 && any(x==1) )) == TRUE)
        c2 <- which(lapply(subs, function(x) (length(x)==1 && any(x==2) )) == TRUE)
        hdens[c1] <- corners_alm_2d(alpha=alpha, beta=beta, s=1) / c
        hdens[c2] <- corners_alm_2d(alpha=alpha, beta=beta, s=2) / c
        
        ## interior
        
        inter <- which(lapply(subs, function(x) (length(x)==2 )) == TRUE)
        int01 <- 2 - corners_alm_2d(alpha=alpha, beta=beta, s=1) - corners_alm_2d(alpha=alpha, beta=beta, s=2)
        intc <- tryCatch(integrate(Vectorize(function(x) interior_alm_2d(c(x,1-x), alpha=alpha, beta=beta)), lower=c, upper=1-c )$value, error=function(e) -1)
        if(intc==-1){
          hdens[inter] <- rep(0, length(inter))
        }else{
          K.inter <- int01/intc
          
          if(length(inter)==1){
            hdens[inter] <- K.inter * interior_alm_2d(w=data[inter,], alpha=alpha, beta=beta)
          }else if(length(inter)>1){
            hdens[inter] <- K.inter * apply(data[inter,],1, function(x) interior_alm_2d(w=x, alpha=alpha, beta=beta))
          }
        }
      }
      return(hdens)
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
    w1=w[s]/sum(w[c(s,t)]);w2=w[t]/sum(w[c(s,t)]);
    
    part1 <- (a-1) * (b1*b2)^a * (w1*w2)^(-a-1)
    part2 <- ( (b1/w1)^a + (b2/w2)^a )^(1/a-2)
    
    return(part1*part2)
  }
  
  corners_alm_3d <- function(alpha,beta,s){ # mass on the corner of the s-th component
    if(s==1){return(1-beta[1]-beta[2]-beta[4])}
    if(s==2){return(1-beta[2]-beta[5]-beta[8])}
    if(s==3){return(1-beta[4]-beta[6]-beta[9])}
  }
  
  dens_alm_3d <- function(data,alpha,beta,c){
    if(length(alpha)!=4){return(stop("Wrong length of parameter alpha"))}
    if(length(beta)!=9){return(stop("Wrong length of parameter beta"))}
    if(beta[1]+beta[3]+beta[7]>1){return(1e-300)} # see conditions on the beta parameters
    if(beta[2]+beta[5]+beta[8]>1){return(1e-300)}
    if(beta[4]+beta[6]+beta[9]>1){return(1e-300)}
    if(any(alpha < 1) || any(beta < 0) || any(beta > 1)){return(1e-300)}
    
    ###
    ### Only one observation
    ###
    
    if(is.vector(data) && length(data)==3){
      
      if(c==0){
        hdens <- interior_alm_3d(w=data, alpha=alpha, beta=beta)
      }else if(c>0){
        subs <- subset.c(w=data, c=c)
        if(length(subs) == 1){ #corner
          K = sqrt(3) / c^2
          hdens <- K * corners_alm_3d(alpha=alpha, beta=beta, s=subs)
        }else if (length(subs)==2){ #edge
          if(subs[1]==1 && subs[2]==2){
            edg01 <- tryCatch(integrate(Vectorize(function(x){edges_alm_3d(w=c(x,1-x,0), alpha=alpha, beta=beta, s=1, t=2)}), lower=0, upper=1)$value, error = function(e) -1)
            edge_surface <- tryCatch(integrate(Vectorize(function(x){edges_alm_3d(w=c(x,1-x,0), alpha=alpha, beta=beta, s=1, t=2)}), lower=c/(2*(1-c/2)), upper=(1-c)/(1-c/2) )$value, error = function(e) -1)
          }else if(subs[1]==1 && subs[2]==3){
            edg01 <- tryCatch(integrate(Vectorize(function(x){edges_alm_3d(w=c(x,0,1-x), alpha=alpha, beta=beta, s=1, t=3)}), lower=0, upper=1)$value, error = function(e) -1)
            edge_surface <- tryCatch(integrate(Vectorize(function(x){edges_alm_3d(w=c(x,0,1-x), alpha=alpha, beta=beta, s=1, t=3)}), lower=c/(2*(1-c/2)), upper=(1-c)/(1-c/2))$value, error = function(e) -1)
          }else if(subs[1]==2 && subs[2]==3){
            edg01 <- tryCatch(integrate(Vectorize(function(x){edges_alm_3d(w=c(0,x,1-x), alpha=alpha, beta=beta, s=2, t=3)}), lower=0, upper=1)$value, error = function(e) -1)
            edge_surface <- tryCatch(integrate(Vectorize(function(x){edges_alm_3d(w=c(0,x,1-x), alpha=alpha, beta=beta, s=2, t=3)}), lower=c/(2*(1-c/2)), upper=(1-c)/(1-c/2))$value, error = function(e) -1)
          }
          #edge_surface <- c - 4/sqrt(3)*c^2
          if(any(c(edg01, edge_surface)==-1)){
            hdens <- 0
          }else{
            K <- edg01 / edge_surface
            hdens <- K * edges_alm_3d(w=data, alpha=alpha, beta=beta, s=subs[1], t=subs[2])            
          }

        }else{ #interior
          edg12 <- tryCatch(integrate(Vectorize(function(x){edges_alm_3d(w=c(x,1-x,0), alpha=alpha, beta=beta, s=1, t=2)}), lower=0, upper=1)$value, error = function(e) -1)
          edg13 <- tryCatch(integrate(Vectorize(function(x){edges_alm_3d(w=c(x,0,1-x), alpha=alpha, beta=beta, s=1, t=3)}), lower=0, upper=1)$value, error = function(e) -1)
          edg23 <- tryCatch(integrate(Vectorize(function(x){edges_alm_3d(w=c(0,x,1-x), alpha=alpha, beta=beta, s=2, t=3)}), lower=0, upper=1)$value, error = function(e) -1)
          
          if(any(c(edg12, edg13, edg23)==-1)){
            hdens <- 0
          }else{
            int01 <- 3 - corners_alm_3d(alpha=alpha, beta=beta, s=1) - corners_alm_3d(alpha=alpha, beta=beta, s=2) - corners_alm_3d(alpha=alpha, beta=beta, s=3) - edg12 - edg13 - edg23
            intc <- tryCatch(integrate(Vectorize(function(y) integrate(Vectorize(function(x) interior_alm_3d(c(x,y,1-x-y),alpha=alpha, beta=beta)), lower=c, upper=1-y-c )$value), lower=c, upper=1-2*c)$value, error = function(e) -1)
            if(intc==-1){
              hdens <- 0
            }else{
              K <- int01/intc
              hdens <- K * interior_alm_3d(w=data, alpha=alpha, beta=beta)
            }            
          }
        }
      }
      return(hdens)
    }
    
    ###
    ### Multiple observation in matrix form
    ###
    
    if(is.matrix(data) && ncol(data)==3){
      
      hdens <- vector(length=nrow(data))
      
      if(c == 0){
        hdens <- apply(data,1, function(x) interior_alm_3d(w=x, alpha=alpha, beta=beta))
      }else if(c>0){
        subs <- apply(data, 1, function(x) subset.c(x,c=c))
        
        if(is.matrix(subs) && nrow(subs)==3){ # all points in the interior (c is too small)
          hdens <- apply(data,1, function(x) interior_alm_3d(w=x, alpha=alpha, beta=beta))
        }
        
        if(is.list(subs)){
          
          ## Corners
          
          c1 <- which(lapply(subs, function(x) (length(x)==1 && any(x ==1) )) == TRUE)
          c2 <- which(lapply(subs, function(x) (length(x)==1 && any(x ==2) )) == TRUE)
          c3 <- which(lapply(subs, function(x) (length(x)==1 && any(x ==3) )) == TRUE)
          K.c1 <- K.c2 <- K.c3 <- sqrt(3) / c^2
          
          hdens[c1] <-  K.c1 * corners_alm_3d(alpha=alpha, beta=beta, s=1)
          hdens[c2] <-  K.c2 * corners_alm_3d(alpha=alpha, beta=beta, s=2)
          hdens[c3] <-  K.c3 * corners_alm_3d(alpha=alpha, beta=beta, s=3)
          
          ## Edges
          
          #edge_surface <- c - 4/sqrt(3)*c^2
          
          ## Edge {1,2}
          
          e12 <- which(lapply(subs, function(x) (length(x)==2 && any(x ==1) && any(x ==2) )) == TRUE)
          edg12 <- tryCatch(integrate(Vectorize(function(x){edges_alm_3d(w=c(x,1-x,0), alpha=alpha, beta=beta, s=1, t=2)}), lower=0, upper=1)$value, error = function(e) -1)
          edge_surface12 <- tryCatch(integrate(Vectorize(function(x){edges_alm_3d(w=c(x,1-x,0), alpha=alpha, beta=beta, s=1, t=2)}), lower=c/(2*(1-c/2)), upper=(1-c)/(1-c/2))$value, error = function(e) -1)
          if(any(c(edg12, edge_surface12)==-1)){
            hdens[e12] <- rep(0, length(e12))
          }else{
            K.e12 <- edg12 / edge_surface12
            if(length(e12)==1){
              hdens[e12] <- K.e12 * edges_alm_3d(w=data[e12,], alpha=alpha, beta=beta, s=1, t=2)
            }else if(length(e12)>1){
              hdens[e12] <- K.e12 * apply(data[e12,],1, function(x) edges_alm_3d(w=x, alpha=alpha, beta=beta, s=1, t=2))
            }
          }
          
          ## Edge {1,3}
          
          e13 <- which(lapply(subs, function(x) (length(x)==2 && any(x ==1) && any(x ==3) )) == TRUE)
          edg13 <- tryCatch(integrate(Vectorize(function(x){edges_alm_3d(w=c(x,0,1-x), alpha=alpha, beta=beta, s=1, t=3)}), lower=0, upper=1)$value, error = function(e) -1)
          edge_surface13 <- tryCatch(integrate(Vectorize(function(x){edges_alm_3d(w=c(x,0,1-x), alpha=alpha, beta=beta, s=1, t=3)}), lower=c/(2*(1-c/2)), upper=(1-c)/(1-c/2))$value, error = function(e) -1)
          if(any(c(edg13, edge_surface13)==-1)){
            hdens[e13] <- rep(0, length(e13))
          }else{
            K.e13 <- edg13 / edge_surface13
            if(length(e13)==1){
              hdens[e13] <- K.e13 * edges_alm_3d(w=data[e13,], alpha=alpha, beta=beta, s=1, t=3)
            }else if(length(e13)>1){
              hdens[e13] <- K.e13 * apply(data[e13,],1, function(x) edges_alm_3d(w=x, alpha=alpha, beta=beta, s=1, t=3))
            }
          }
          
          ## Edge {2,3}
          
          e23 <- which(lapply(subs, function(x) (length(x)==2 && any(x ==2) && any(x ==3) )) == TRUE)
          edg23 <- tryCatch(integrate(Vectorize(function(x){edges_alm_3d(w=c(0,x,1-x), alpha=alpha, beta=beta, s=2, t=3)}), lower=0, upper=1)$value, error = function(e) -1)
          edge_surface23 <- tryCatch(integrate(Vectorize(function(x){edges_alm_3d(w=c(0,x,1-x), alpha=alpha, beta=beta, s=2, t=3)}), lower=c/(2*(1-c/2)), upper=(1-c)/(1-c/2))$value, error = function(e) -1)
          if(any(c(edg23, edge_surface23)==-1)){
            hdens[e23] <- rep(0, length(e23))
          }else{
            K.e23 <- edg23 / edge_surface23
            if(length(e23)==1){
              hdens[e23] <- K.e23 * edges_alm_3d(w=data[e23,], alpha=alpha, beta=beta, s=2, t=3)
            }else if(length(e23)>1){
              hdens[e23] <- K.e23 * apply(data[e23,],1, function(x) edges_alm_3d(w=x, alpha=alpha, beta=beta, s=2, t=3))
            }
          }
          
          ## Interior
          
          inter <- which(lapply(subs, function(x) (length(x)==3 )) == TRUE)
          if(any(c(edg12, edg13, edg23)==-1)){
            hdens[inter] <- rep(0, length(inter))
          }else{
            int01 <- 3 - corners_alm_3d(alpha=alpha, beta=beta, s=1) - corners_alm_3d(alpha=alpha, beta=beta, s=2) - corners_alm_3d(alpha=alpha, beta=beta, s=3) - edg12 - edg13 - edg23
            intc <- tryCatch(integrate(Vectorize(function(y) integrate(Vectorize(function(x) interior_alm_3d(c(x,y,1-x-y), alpha=alpha, beta=beta)), lower=c, upper=1-y-c )$value), lower=c, upper=1-2*c)$value, error = function(e) -1)
            if(intc==-1){
              hdens[inter] <- rep(0, length(inter))
            }else{
              K.inter <- int01/intc
              if(length(inter)==1){
                hdens[inter] <- K.inter * interior_alm_3d(w=data[inter,], alpha=alpha, beta=beta)
              }else if(length(inter)>1){
                hdens[inter] <- K.inter * apply(data[inter,],1, function(x) interior_alm_3d(w=x, alpha=alpha, beta=beta))
              }  
            }
          }
          
        }
        
      }
      
      return(hdens)
      
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
    
    if (any(round(apply(x,1,sum),7) != 1)){ stop("Data is not angular")}
    
    if(dim==2){ result <- dens_alm_2d(data=x, alpha=alpha, beta=beta, c=c)}
    if(dim==3){ result <- dens_alm_3d(data=x, alpha=alpha, beta=beta, c=c)}
    
    if(!vectorial){
      result <- prod(result)
    }
  }
  if(any(is.na(result))){
    result[is.na(result)] <- 0
  }
  
  if(log){
    if(any(result==0)){
      ind0 <- which(result==0)
      ind <- which(result!=0)
      result[ind0] <- -1e+300
      result[ind] <- log(result[ind])
      return(result)
    }else{
      return(log(result))
    }
  }else{
    return(result)
  }
}

###
### Extremal-t model
###

subset.c <- function(w, c){
  
  d <- length(w)
  
  if(d==2){
    if(any(w > 1-c)){ # Corner
      return(which(w > 1-c))
    }else{
      return(c(1,2))
    }
  }else if(d==3){
    
    if(sum(w > 1-c)==1){ # Corner
      return(which(w > 1-c))
    }else if( all(w > c) && all(w < 1-2*c) ){ # interior
      return(c(1,2,3))
    }else{ #edges
      if(w[1]<=1-c && w[2]<=1-c && w[3]<=c && w[1]>=(1-w[2])/2 && w[1]>=1-2*w[2]){ #EDGE {1,2}
        return(c(1,2))
      }else if(w[1]<=1-c && w[2]<=c && w[3]<=1-c && w[2]<=w[1] && w[2]<=(1-w[1])/2){ #EDGE {1,3}
        return(c(1,3))
      }else if(w[1]<=c && w[2]<=1-c && w[3]<=1-c && w[1]<=w[2] && w[1]<=(1-w[2])/2){ #EDGE {2,3}
        return(c(2,3))
      }
    }
  }
  
}

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
    
    if(sum(eigen(Sigma)$values<0)>=1){return(1e-300)} #Check if matrix is positive definite
    
    deriv = (w[-1]/w[1])^(1/mu-1)/mu*sqrt((mu+1)/(1-rho[1:k]^2))
    if(k==1){
      return(dest(x=w.tilde, scale=Sigma, df=mu+1)*w[1]^(-p-1)*prod(deriv) )
    }else{
      return(dmest(x=w.tilde, scale=Sigma, df=mu+1)*w[1]^(-p-1)*prod(deriv) )
    }
  }
  
  # Mass on the other subsets of the 2d case:
  
  corners_et_2d <- function(rho,mu){ # mass on the corner of the s-th component
    
    part1 <- pt( -rho * sqrt((mu+1)/(1-rho^2)) ,df=mu+1)
    
    return(part1)
  }
  
  dens_et_2d <- function(data, rho, mu, c){
    
    if(length(rho)!=1){return(stop("Wrong length of parameter rho"))}
    if( (abs(rho) > 1) || (mu<1) ){return(1e-300)}
    
    ###
    ### Only one observation
    ###
    
    if(is.vector(data) && length(data)==2){
      
      if(c==0){
        hdens <- interior_et_d(w=data,rho,mu)
      }else if(c>0){
        subs <- subset.c(w=data, c=c)
        if(length(subs) == 1){ #corner
          K <- 1 / c
          hdens <- K * corners_et_2d(rho=rho, mu=mu)
        }else if (length(subs)==2){
          int01 <- 2 - 2*corners_et_2d(rho=rho, mu=mu)
          intc <- tryCatch(integrate(Vectorize(function(x) interior_et_d(c(x,1-x),rho=rho, mu=mu)), lower=c, upper=1-c )$value, error=function(e) -1)
          if(intc==-1){
            hdens <- 0
          }else{
            K.inter <- int01/intc
            hdens <- K.inter * interior_et_d(w=data, rho=rho, mu=mu)          
          }
  
        }
      }
      return(hdens)
    }
    
    ###
    ### Multiple observation in matrix form
    ###
    
    if(is.matrix(data) && ncol(data)==2){
      
      hdens <- vector(length=nrow(data))
      
      if(c == 0){
        hdens <- apply(data,1, function(x) interior_et_d(w=x, rho=rho, mu=mu))
      }else if(c>0){
        
        subs <- apply(data, 1, function(x) subset.c(x,c=c))
        
        ## Corners
        
        corners <- which(lapply(subs, function(x) (length(x)==1 )) == TRUE)
        hdens[corners] <-  corners_et_2d(rho=rho, mu=mu) / c
        
        ## interior
        
        inter <- which(lapply(subs, function(x) (length(x)==2 )) == TRUE)
        int01 <- 2 - 2*corners_et_2d(rho=rho, mu=mu)
        intc <- tryCatch(integrate(Vectorize(function(x) interior_et_d(c(x,1-x),rho=rho, mu=mu)), lower=c, upper=1-c )$value, error=function(e) -1)
        if(intc==-1){
          hdens[inter] <- rep(0, length(inter))
        }else{
          K.inter <- int01/intc
          
          if(length(inter)==1){
            hdens[inter] <- K.inter * interior_et_d(w=data[inter,], rho=rho, mu=mu)
          }else if(length(inter)>1){
            hdens[inter] <- K.inter * apply(data[inter,],1, function(x) interior_et_d(w=x, rho=rho, mu=mu))
          }  
        }
        
      }
      return(hdens)
    }
    
  }
  
  #dens_et_2d <- function(w,rho,mu,c){
  #    if(length(rho)!=1){return(stop("Wrong length of parameter rho"))}
  #    if( (abs(rho) > 1) || (mu<1) ){return(1e-300)}
  #
  #    if(sum(w >= (1-c)) == 1){ # then we are in a corner
  #        ind <- which(w >= (1-c))
  #        return(corners_et_2d(w,rho,mu,ind))
  #    } else {
  #        return(interior_et_d(w,rho,mu))
  #    }
  #}
  
  # Mass on the other subsets of the 3d case:
  
  corners_et_3d <- function(w,rho,mu,s){ # mass on the corner of the s-th component
    
    if(s==1){rho1=rho[1];rho2=rho[2];rho3=rho[3]}
    if(s==2){rho1=rho[1];rho2=rho[3];rho3=rho[2]}
    if(s==3){rho1=rho[2];rho2=rho[3];rho3=rho[1]}
    
    R <- (rho3 - rho1*rho2)/sqrt((1-rho1^2) * (1-rho2^2))
    S <- matrix( c(1,R,R,1), ncol=2)
    if(sum(eigen(S)$values<0)>=1){return(1e-300)} # Check if matrix R1 is positive definite
    
    a = -rho1 * sqrt((mu+1)/(1-rho1^2))
    b = -rho2 * sqrt((mu+1)/(1-rho2^2))
    
    return( pmest(x=c(a,b), scale=S, df=mu+1) )
  }
  
  # mass on the edge linking the s-th and t-th components (in inceasing order)
  
  edges_et_3d <- function(w,rho,mu,s,t){
    
    if(s==1 && t==2){rho1=rho[1];rho2=rho[2];rho3=rho[3]}
    if(s==1 && t==3){rho1=rho[2];rho2=rho[1];rho3=rho[3]}
    if(s==2 && t==3){rho1=rho[3];rho2=rho[1];rho3=rho[2]}
    x=w[s]/sum(w[c(s,t)]);
    y=w[t]/sum(w[c(s,t)]);
    
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
    
    if(sum(eigen(S1)$values<0)>=1){return(1e-300)} # Check if matrix R1 is positive definite
    if(sum(eigen(S2)$values<0)>=1){return(1e-300)} # Check if matrix R2 is positive definite
    if(any( abs(c(R1,R2)) >1 ) ){return(1e-300)}
    
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
  
  #dens_et_3d <- function(w,rho,mu,c){
  #
  #    if(length(rho)!=3){return(stop("Wrong length of parameter rho"))}
  #    if(any(abs(rho)>=1) || mu<=0 ){return(1e-300)}
  #
  #    if(c==0){return(interior_et_d(w,rho,mu))}
  #
  #    if(sum(w<=c) == 2){ # then we are in a corner
  #        ind <- which(w > c)
  #        return(corners_et_3d(w,rho,mu,ind)/c^2)
  #    }else if(sum(w<=c) == 1){
  #        ind <- which(w > c)
  #        w2 <- w[ind]/sum(w[ind])
  #        edge_surface <- c*sqrt(3)*(1-2*c)/2
  #
  #        if(w[1]<=1-c && w[2]<=1-c && w[3]<=c && w[1]>=(1-w[2])/2 && w[1]>=1-2*w[2]){ #EDGE {1,2}
  #            edg01 <- integrate(Vectorize(function(x){edges_et_3d(w=c(x,1-x,0), rho=rho, mu=mu, s=1, t=2)}), lower=0, upper=1)$value
  #            return(edges_et_3d(c(w2,w[3]),rho,mu,1,2) * edg01 / edge_surface)
  #        }
  #
  #        if(w[1]<=1-c && w[2]<=c && w[3]<=1-c && w[2]<=w[1] && w[2]<=(1-w[1])/2){ #EDGE {1,3}
  #            edg01 <- integrate(Vectorize(function(x){edges_et_3d(w=c(x,0,1-x), rho=rho, mu=mu, s=1, t=3)}), lower=0, upper=1)$value
  #            return(edges_et_3d(c(w2[1],w[2],w2[2]),rho,mu,1,3) * edg01 / edge_surface)
  #        }
  #
  #        if(w[1]<=c && w[2]<=1-c && w[3]<=1-c && w[1]<=w[2] && w[1]<=(1-w[2])/2){ #EDGE {2,3}
  #            edg01 <- integrate(Vectorize(function(x){edges_et_3d(w=c(0,x,1-x), rho=rho, mu=mu, s=2, t=3)}), lower=0, upper=1)$value
  #            return(edges_et_3d(c(w[1],w2),rho,mu,2,3) * edg01 / edge_surface)
  #        }
  #    }else {
  #        int01 <- integrate(Vectorize(function(y) integrate(Vectorize(function(x) interior_et_d(c(x,y,1-x-y),rho=rho, mu=mu)), lower=0, upper=1-y )$value), lower=0, upper=1)$value
  #         intc <- integrate(Vectorize(function(y) integrate(Vectorize(function(x) interior_et_d(c(x,y,1-x-y),rho=rho, mu=mu)), lower=c, upper=1-y-c )$value), lower=c, upper=1-2*c)$value
  #        return(interior_et_d(w,rho,mu)*int01/intc)
  #    }
  #}
  
  dens_et_3d <- function(data, rho, mu, c){
    
    if(length(rho)!=3){return(stop("Wrong length of parameter rho"))}
    if(any(abs(rho)>=1) || mu<=0 ){return(1e-300)}
    
    ###
    ### Only one observation
    ###
    
    if(is.vector(data) && length(data)==3){
      
      if(c==0){
        hdens <- interior_et_d(w=data, rho=rho, mu=mu)
      }else if(c>0){
        subs <- subset.c(w=data, c=c)
        if(length(subs) == 1){ #corner
          K = sqrt(3) / c^2
          hdens <- K * corners_et_3d(rho=rho, mu=mu, s=subs)
        }else if (length(subs)==2){ #edge
          if(subs[1]==1 && subs[2]==2){
            edg01 <- tryCatch(integrate(Vectorize(function(x){edges_et_3d(w=c(x,1-x,0), rho=rho, mu=mu, s=1, t=2)}), lower=0, upper=1)$value, error=function(e) -1)
            edge_surface <- tryCatch(integrate(Vectorize(function(x){edges_et_3d(w=c(x,1-x,0), rho=rho, mu=mu, s=1, t=2)}), lower=c/(2*(1-c/2)), upper=(1-c)/(1-c/2) )$value, error=function(e) -1)
          }else if(subs[1]==1 && subs[2]==3){
            edg01 <- tryCatch(integrate(Vectorize(function(x){edges_et_3d(w=c(x,0,1-x), rho=rho, mu=mu, s=1, t=3)}), lower=0, upper=1)$value, error=function(e) -1)
            edge_surface <- tryCatch(integrate(Vectorize(function(x){edges_et_3d(w=c(x,0,1-x), rho=rho, mu=mu, s=1, t=3)}), lower=c/(2*(1-c/2)), upper=(1-c)/(1-c/2))$value, error=function(e) -1)
          }else if(subs[1]==2 && subs[2]==3){
            edg01 <- tryCatch(integrate(Vectorize(function(x){edges_et_3d(w=c(0,x,1-x), rho=rho, mu=mu, s=2, t=3)}), lower=0, upper=1)$value, error=function(e) -1)
            edge_surface <- tryCatch(integrate(Vectorize(function(x){edges_et_3d(w=c(0,x,1-x), rho=rho, mu=mu, s=2, t=3)}), lower=c/(2*(1-c/2)), upper=(1-c)/(1-c/2))$value, error=function(e) -1)
          }
          #edge_surface <- c - 4/sqrt(3)*c^2
          if(any(c(edg01, edge_surface)==-1)){
            hdens <- 0
          }else{
            K <- edg01 / edge_surface
            hdens <- K * edges_et_3d(w=data, rho=rho, mu=mu, s=subs[1], t=subs[2])            
          }

        }else{ #interior
          edg12 <- tryCatch(integrate(Vectorize(function(x){edges_et_3d(w=c(x,1-x,0), rho=rho, mu=mu, s=1, t=2)}), lower=0, upper=1)$value, error=function(e) -1)
          edg13 <- tryCatch(integrate(Vectorize(function(x){edges_et_3d(w=c(x,0,1-x), rho=rho, mu=mu, s=1, t=3)}), lower=0, upper=1)$value, error=function(e) -1)
          edg23 <- tryCatch(integrate(Vectorize(function(x){edges_et_3d(w=c(0,x,1-x), rho=rho, mu=mu, s=2, t=3)}), lower=0, upper=1)$value, error=function(e) -1)
          
          if(any(c(edg12, edg13, edg23)==-1)){
            hdens <- 0
          }else{
            int01 <- 3 - corners_et_3d(rho=rho, mu=mu, s=1) - corners_et_3d(rho=rho, mu=mu, s=2) - corners_et_3d(rho=rho, mu=mu, s=3) - edg12 - edg13 - edg23
            intc <- tryCatch(integrate(Vectorize(function(y) integrate(Vectorize(function(x) interior_et_d(c(x,y,1-x-y),rho=rho, mu=mu)), lower=c, upper=1-y-c )$value), lower=c, upper=1-2*c)$value, error=function(e) -1)
            if(intc==-1){
              hdens <- 0
            }else{
              K <- int01/intc
              hdens <- K * interior_et_d(w=data, rho=rho, mu=mu)
            }
            
          }
          
        }
      }
      return(hdens)
    }
    
    ###
    ### Multiple observation in matrix form
    ###
    
    if(is.matrix(data) && ncol(data)==3){
      
      hdens <- vector(length=nrow(data))
      
      if(c == 0){
        hdens <- apply(data,1, function(x) interior_et_d(w=x, rho=rho, mu=mu))
      }else if(c>0){
        subs <- apply(data, 1, function(x) subset.c(x,c=c))
        
        if(is.matrix(subs) && nrow(subs)==3){ # all points in the interior (c is too small)
          hdens <- apply(data,1, function(x) interior_et_d(w=x, rho=rho, mu=mu))
        }
        
        if(is.list(subs)){
          
          ## Corners
          
          c1 <- which(lapply(subs, function(x) (length(x)==1 && any(x ==1) )) == TRUE)
          c2 <- which(lapply(subs, function(x) (length(x)==1 && any(x ==2) )) == TRUE)
          c3 <- which(lapply(subs, function(x) (length(x)==1 && any(x ==3) )) == TRUE)
          K.c1 <- K.c2 <- K.c3 <- sqrt(3) / c^2
          
          hdens[c1] <-  K.c1 * corners_et_3d(rho=rho, mu=mu, s=1)
          hdens[c2] <-  K.c2 * corners_et_3d(rho=rho, mu=mu, s=2)
          hdens[c3] <-  K.c3 * corners_et_3d(rho=rho, mu=mu, s=3)
          
          ## Edges
          
          #edge_surface <- c - 4/sqrt(3)*c^2
          
          ## Edge {1,2}
          
          e12 <- which(lapply(subs, function(x) (length(x)==2 && any(x ==1) && any(x ==2) )) == TRUE)
          edg12 <- tryCatch(integrate(Vectorize(function(x){edges_et_3d(w=c(x,1-x,0), rho=rho, mu=mu, s=1, t=2)}), lower=0, upper=1)$value, error=function(e) -1)
          edge_surface12 <- tryCatch(integrate(Vectorize(function(x){edges_et_3d(w=c(x,1-x,0), rho=rho, mu=mu, s=1, t=2)}), lower=c/(2*(1-c/2)), upper=(1-c)/(1-c/2))$value, error=function(e) -1)
          if(any(c(edg12, edge_surface12)==-1)){
            hdesn[e12] <- rep(0, length(e12))
          }else{
            K.e12 <- edg12 / edge_surface12
            if(length(e12)==1){
              hdens[e12] <- K.e12 * edges_et_3d(w=data[e12,], rho=rho, mu=mu, s=1, t=2)
            }else if(length(e12)>1){
              hdens[e12] <- K.e12 * apply(data[e12,],1, function(x) edges_et_3d(w=x, rho=rho, mu=mu, s=1, t=2))
            }
          }
          
          ## Edge {1,3}
          
          e13 <- which(lapply(subs, function(x) (length(x)==2 && any(x ==1) && any(x ==3) )) == TRUE)
          edg13 <- tryCatch(integrate(Vectorize(function(x){edges_et_3d(w=c(x,0,1-x), rho=rho, mu=mu, s=1, t=3)}), lower=0, upper=1)$value, error=function(e) -1)
          edge_surface13 <- tryCatch(integrate(Vectorize(function(x){edges_et_3d(w=c(x,0,1-x), rho=rho, mu=mu, s=1, t=3)}), lower=c/(2*(1-c/2)), upper=(1-c)/(1-c/2))$value, error=function(e) -1)
          if(any(c(edg13, edge_surface13)==-1)){
            hdesn[e13] <- rep(0, length(e13))
          }else{
            K.e13 <- edg13 / edge_surface13
            if(length(e13)==1){
              hdens[e13] <- K.e13 * edges_et_3d(w=data[e13,], rho=rho, mu=mu, s=1, t=3)
            }else if(length(e13)>1){
              hdens[e13] <- K.e13 * apply(data[e13,],1, function(x) edges_et_3d(w=x, rho=rho, mu=mu, s=1, t=3))
            }
          }
          
          ## Edge {2,3}
          
          e23 <- which(lapply(subs, function(x) (length(x)==2 && any(x ==2) && any(x ==3) )) == TRUE)
          edg23 <- tryCatch(integrate(Vectorize(function(x){edges_et_3d(w=c(0,x,1-x), rho=rho, mu=mu, s=2, t=3)}), lower=0, upper=1)$value, error=function(e) -1)
          edge_surface23 <- tryCatch(integrate(Vectorize(function(x){edges_et_3d(w=c(0,x,1-x), rho=rho, mu=mu, s=2, t=3)}), lower=c/(2*(1-c/2)), upper=(1-c)/(1-c/2))$value, error=function(e) -1)
          if(any(c(edg23, edge_surface23)==-1)){
            hdesn[e23] <- rep(0, length(e23))
          }else{
            K.e23 <- edg23 / edge_surface23
            if(length(e23)==1){
              hdens[e23] <- K.e23 * edges_et_3d(w=data[e23,], rho=rho, mu=mu, s=2, t=3)
            }else if(length(e23)>1){
              hdens[e23] <- K.e23 * apply(data[e23,],1, function(x) edges_et_3d(w=x, rho=rho, mu=mu, s=2, t=3))
            }
          }
          
          ## Interior
          
          inter <- which(lapply(subs, function(x) (length(x)==3 )) == TRUE)
          if(any(c(edg12, edg13, edg23)==-1)){
            hdens[inter] <- rep(0, length(inter))
          }else{
            int01 <- 3 - corners_et_3d(rho=rho, mu=mu, s=1) - corners_et_3d(rho=rho, mu=mu, s=2) - corners_et_3d(rho=rho, mu=mu, s=3) - edg12 - edg13 - edg23
            intc <- tryCatch(integrate(Vectorize(function(y) integrate(Vectorize(function(x) interior_et_d(c(x,y,1-x-y),rho=rho, mu=mu)), lower=c, upper=1-y-c )$value), lower=c, upper=1-2*c)$value, error=function(e) -1)
            if(intc==-1){
              hdens[inter] <- rep(0, length(inter))
            }else{
              K.inter <- int01/intc
              if(length(inter)==1){
                hdens[inter] <- K.inter * interior_et_d(w=data[inter,], rho=rho, mu=mu)
              }else if(length(inter)>1){
                hdens[inter] <- K.inter * apply(data[inter,],1, function(x) interior_et_d(w=x, rho=rho, mu=mu))
              }
            }
 
          }
          
        }
        
      }
      
      return(hdens)
      
    }
    
  }
  
  ### Angular density for the Extremal-t model on the 2 and 3 dimensional simplex
  
  xvect = as.double(as.vector(t(x)))
  if (is.vector(x)) {
    dim = as.integer(length(x))
    n = as.integer(1)
    if(round(sum(x),7) != 1){ stop("Data is not angular")}
    if(dim==2){ result <- dens_et_2d(data=x, rho=rho, mu=mu, c=c)}
    if(dim==3){ result <- dens_et_3d(data=x, rho=rho, mu=mu, c=c)}
  }
  else {
    dim = as.integer(ncol(x))
    n = as.integer(nrow(x))
    
    if(any(round(apply(x,1,sum),7) != 1)){ stop("Data is not angular") }
    
    if(dim==2){ result <- dens_et_2d(data=x, rho=rho, mu=mu, c=c)}
    if(dim==3){ result <- dens_et_3d(data=x, rho=rho, mu=mu, c=c)}
    
    if(!vectorial){
      result <- prod(result)
    }
    
    #if (vectorial) {
    #    result = double(n)
    #    if(dim==2){ result <- apply(x,1,function(y){dens_et_2d(y,rho,mu,c)}) }
    #    if(dim==3){ result <- apply(x,1,function(y){dens_et_3d(y,rho,mu,c)}) }
    #} else { # vectorial=FALSE mean we return the likelihood function
    #    result = as.double(1)
    #    if(dim==2){ result <- prod(apply(x,1,function(y){dens_et_2d(y,rho,mu,c)})) }
    #    if(dim==3){ result <- prod(apply(x,1,function(y){dens_et_3d(y,rho,mu,c)})) }
    #   }
  }
  if(any(is.na(result))){
    result[is.na(result)] <- 0
  }
  if(log){
    if(any(result==0)){
      ind0 <- which(result==0)
      ind <- which(result!=0)
      result[ind0] <- -1e+300
      result[ind] <- log(result[ind])
      return(result)
    }else{
      return(log(result))
    }
  }else{
    return(result)
  }
}

###
### Extremal Skew-t model
###

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
  
  alpha_star_j <- function(rho, alpha,j){
    Alpha_tilde <- alpha_tilde(alpha,j)
    sj <- s_j(rho,j)
    return( Alpha_tilde %*% sj )
  }
  
  tau_star_j <- function(rho,alpha,mu,j){
    Sig <- Sigma(rho)
    Alpha_tilde <- alpha_tilde(alpha,j)
    return( sqrt(mu+1) * (alpha[j] + Sig[-j,j] %*% t(Alpha_tilde) )     )
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
    alpha_star <- alpha_star_j(rho,alpha,j)
    return( (alpha_star %*% sigma_bar )/ as.numeric(sqrt( 1 + alpha_star %*% sigma_bar %*% t(alpha_star) )) )
  }
  
  tau_bar_j <- function(rho,alpha,mu,j){
    tau_star <- tau_star_j(rho,alpha,mu,j)
    alpha_star <- alpha_star_j(rho, alpha,j)
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
    
    if( any(cond < -1)){return(1e-300)}else{
      sigma_bar1 <- Sigma_bar_j(rho,1)
      alpha_star1 <- alpha_star_j(rho,alpha,1)
      tau_star1 <- tau_star_j(rho,alpha,mu,1)
      
      w_bar <- sapply(1:p,function(j){w[j]*nu_p(rho,alpha,mu,j)})
      
      for(i in 1:k){
        w.tilde[i] <- ((w_bar[i+1]/w_bar[1])^(1/mu)-rho[i])*sqrt((mu+1)/(1-rho[i]^2));
      }
      
      # if( alpha_star1 %*% sigma_bar1 %*% t(alpha_star1) < -1){return(1e-300)}
      # if( t(w.tilde) %*% solve(sigma_bar1) %*% w.tilde < -1){return(1e-300)}
      if(any(eigen(sigma_bar1)$value<=0)){return(1e-300)}
      
      deriv = (w_bar[-1]/w_bar[1])^(1/mu-1) / mu * sqrt((mu+1)/(1-rho[1:k]^2)) * sapply(2:p,function(x){nu_p(rho,alpha,mu,x)}) / as.vector(nu_p(rho,alpha,mu,1))
      
      if(k==1){
        return(dest(x=w.tilde, scale=sigma_bar1, shape=t(alpha_star1), extended=tau_star1, df=mu+1)*w[1]^(-p-1)*prod(deriv) )
      }else{
        return(dmest(x=w.tilde, location=rep(0,2), scale=sigma_bar1, shape=t(alpha_star1), extended=tau_star1, df=mu+1) * w[1]^(-p-1) * prod(deriv) )
      }
      
    }
  }
  
  ## mass on the corner of the s-th component of the 2-d simplex
  
  corners_skewt_2d <- function(w,rho,alpha,mu,s){
    
    alpha_star <- alpha_star_j(rho,alpha,s)
    tau_star <- tau_star_j(rho,alpha,mu,s)
    part1 <- pest(-rho * sqrt((mu+1)/(1-rho^2)),shape=alpha_star, extended=tau_star ,df=mu+1)
    
    return(part1)
  }
  
  ## density on 2-d simplex
  
  dens_skewt_2d <- function(data, rho, alpha, mu, c){
    if(length(rho)!=1){stop("Wrong length of parameter rho")}
    if(length(alpha)!=2){stop("Wrong length of parameter alpha")}
    if( (abs(rho) > 1) || (mu<1) ){return(1e-300)}
    
    ###
    ### Only one observation
    ###
    
    if(is.vector(data) && length(data)==2){
      
      if(c==0){
        hdens <- interior_skewt_d(w=data, rho=rho, alpha=alpha, mu=mu)
      }else if(c>0){
        subs <- subset.c(w=data, c=c)
        if(length(subs) == 1){ #corner
          K <- 1 / c
          hdens <- K * corners_skewt_2d(rho=rho, alpha=alpha, mu=mu, s=subs)
        }else if (length(subs)==2){
          int01 <- 2 - corners_skewt_2d(rho=rho, alpha=alpha, mu=mu, s=1) - corners_skewt_2d(rho=rho, alpha=alpha, mu=mu, s=2)
          intc <- tryCatch(integrate(Vectorize(function(x) interior_skewt_d(c(x,1-x),rho=rho, alpha=alpha, mu=mu)), lower=c, upper=1-c )$value, error=function(e) -1)
          if(intc==-1){
            hdens <- 0
          }else{
            K.inter <- int01/intc
            hdens <- K.inter * interior_skewt_d(w=data, rho=rho, alpha=alpha, mu=mu)            
          }
        }
      }
      return(hdens)
    }
    
    ###
    ### Multiple observation in matrix form
    ###
    
    if(is.matrix(data) && ncol(data)==2){
      
      hdens <- vector(length=nrow(data))
      
      if(c == 0){
        hdens <- apply(data,1, function(x) interior_skewt_d(w=x, rho=rho, alpha=alpha, mu=mu))
      }else if(c>0){
        
        subs <- apply(data, 1, function(x) subset.c(x,c=c))
        
        ## Corners
        
        corners1 <- which(lapply(subs, function(x) (length(x)==1 && x==1 )) == TRUE)
        corners2 <- which(lapply(subs, function(x) (length(x)==1 && x==2 )) == TRUE)
        hdens[corners1] <- corners_skewt_2d(rho=rho, alpha=alpha, mu=mu, s=1) / c
        hdens[corners2] <- corners_skewt_2d(rho=rho, alpha=alpha, mu=mu, s=2) / c
        
        ## interior
        
        inter <- which(lapply(subs, function(x) (length(x)==2 )) == TRUE)
        int01 <- 2 - corners_skewt_2d(rho=rho, alpha=alpha, mu=mu, s=1) - corners_skewt_2d(rho=rho, alpha=alpha, mu=mu, s=2)
        intc <- tryCatch(integrate(Vectorize(function(x) interior_skewt_d(c(x,1-x),rho=rho, alpha=alpha, mu=mu)), lower=c, upper=1-c )$value, error=function(e) -1)
        if(intc==-1){
          hdens[inter] <- rep(0, length(inter))
        }else{
          K.inter <- int01/intc
          
          if(length(inter)==1){
            hdens[inter] <- K.inter * interior_skewt_d(w=data[inter,], rho=rho, alpha=alpha, mu=mu)
          }else if(length(inter)>1){
            hdens[inter] <- K.inter * apply(data[inter,],1, function(x) interior_skewt_d(w=x, rho=rho, alpha=alpha, mu=mu))
          }
        }
        
      }
      return(hdens)
    }
    
  }
  
  ## mass on the corner of the s-th component of the 3-d simplex
  
  corners_skewt_3d <- function(w,rho,alpha,mu,s){
    
    s_bar <- Sigma_bar_j(rho,s)
    al=t(alpha_star_j(rho,alpha,s))
    if(t(al) %*% s_bar %*% al < -1 ){return(1e-300)}
    if(sum(eigen(s_bar)$values<0)>=1){return(1e-300)} #Check if matrix is positive definite
    
    up <- c(comp(0,0,0,rho,alpha,mu,1,s), comp(0,0,0,rho,alpha,mu,2,s))
    return(pmest(x=up, location=rep(0,2), scale=s_bar, shape=al, extended=tau_star_j(rho,alpha,mu,s), df=mu+1 ))
  }
  
  ## mass on the edge between the s-th and t-th components of the 3-d simplex
  
  edges_skewt_3d <- function(w,rho,alpha,mu,s,t){
    
    sigma_j1 <- Sigma_j(rho,s)
    Alpha_tilde1 <- alpha_tilde(alpha,s)
    sigma_j2 <- Sigma_j(rho,t)
    Alpha_tilde2 <- alpha_tilde(alpha,t)
    alpha_star1 <- alpha_star_j(rho,alpha,s)
    sigma_bar1 <- Sigma_bar_j(rho,s)
    alpha_star2 <- alpha_star_j(rho,alpha,t)
    sigma_bar2 <- Sigma_bar_j(rho,t)
    if( Alpha_tilde1 %*% t(sigma_j1) %*% t(Alpha_tilde1)< -1){return(1e-300)}
    if( Alpha_tilde2 %*% t(sigma_j2) %*% t(Alpha_tilde2)< -1){return(1e-300)}
    if(    alpha_star1 %*% sigma_bar1 %*% t(alpha_star1)< -1){return(1e-300)}
    if(    alpha_star2 %*% sigma_bar2 %*% t(alpha_star2)< -1){return(1e-300)}
    
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
    
    if( R1[1]^2 >1 || R1[2]^2>1 || R2[1]^2 >1 || R2[2]^2>1 ){return(1e-300)}
    
    c1 <- tau_bar_j(rho,alpha,mu,s)
    u1_p <- c(a1_p,b1_p,c1)
    
    R1_p <- diag(2)
    R1_p[1,2] <- R1_p[2,1] <- (R1[3]-R1[1]*R1[2])/sqrt((1-R1[1]^2)*(1-R1[2]^2))
    if(sum(eigen(R1_p)$values<0)>=1){return(1e-300)} #Check if matrix is positive definite
    
    c2 <- tau_bar_j(rho,alpha,mu,t)
    u2_p <- c(a2_p,b2_p,c2)
    
    R2_p <- diag(2)
    R2_p[1,2] <- R2_p[2,1] <- (R2[3]-R2[1]*R2[2])/sqrt((1-R2[1]^2)*(1-R2[2]^2))
    if(sum(eigen(R2_p)$values<0)>=1){return(1e-300)} #Check if matrix is positive definite
    
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
  
  dens_skewt_3d <- function(data, rho, alpha, mu, c){
    
    if(length(rho)!=3){return(stop("Wrong length of parameter rho"))}
    if(length(alpha)!=3){return(stop("Wrong length of parameter rho"))}
    if(any(abs(rho)>=1) || mu<=0){return(1e-300)}
    
    ###
    ### Only one observation
    ###
    
    if(is.vector(data) && length(data)==3){
      
      if(c==0){
        hdens <- interior_skewt_d(w=data, rho=rho, alpha=alpha, mu=mu)
      }else if(c>0){
        subs <- subset.c(w=data, c=c)
        if(length(subs) == 1){ #corner
          K = sqrt(3) / c^2
          hdens <- K * corners_skewt_3d(rho=rho, alpha=alpha, mu=mu, s=subs)
        }else if (length(subs)==2){ #edge
          if(subs[1]==1 && subs[2]==2){
            edg01 <- tryCatch(integrate(Vectorize(function(x){edges_skewt_3d(w=c(x,1-x,0), rho=rho, alpha=alpha, mu=mu, s=1, t=2)}), lower=0, upper=1)$value, error=function(e) -1)
            edge_surface <- tryCatch(integrate(Vectorize(function(x){edges_skewt_3d(w=c(x,1-x,0), rho=rho, alpha=alpha, mu=mu, s=1, t=2)}), lower=c/(2*(1-c/2)), upper=(1-c)/(1-c/2) )$value, error=function(e) -1)
          }else if(subs[1]==1 && subs[2]==3){
            edg01 <- tryCatch(integrate(Vectorize(function(x){edges_skewt_3d(w=c(x,0,1-x), rho=rho, alpha=alpha, mu=mu, s=1, t=3)}), lower=0, upper=1)$value, error=function(e) -1)
            edge_surface <- tryCatch(integrate(Vectorize(function(x){edges_skewt_3d(w=c(x,0,1-x), rho=rho, alpha=alpha, mu=mu, s=1, t=3)}), lower=c/(2*(1-c/2)), upper=(1-c)/(1-c/2))$value, error=function(e) -1)
          }else if(subs[1]==2 && subs[2]==3){
            edg01 <- tryCatch(integrate(Vectorize(function(x){edges_skewt_3d(w=c(0,x,1-x), rho=rho, alpha=alpha, mu=mu, s=2, t=3)}), lower=0, upper=1)$value, error=function(e) -1)
            edge_surface <- tryCatch(integrate(Vectorize(function(x){edges_skewt_3d(w=c(0,x,1-x), rho=rho, alpha=alpha, mu=mu, s=2, t=3)}), lower=c/(2*(1-c/2)), upper=(1-c)/(1-c/2))$value, error=function(e) -1)
          }
          #edge_surface <- c - 4/sqrt(3)*c^2
          if(any(c(edg01, edge_surface)==-1)){
            hdens <- 0
          }else{
            K <- edg01 / edge_surface
            hdens <- K * edges_skewt_3d(w=data, rho=rho, alpha=alpha, mu=mu, s=subs[1], t=subs[2])
          }
          
        }else{ #interior
          edg12 <- tryCatch(integrate(Vectorize(function(x){edges_skewt_3d(w=c(x,1-x,0), rho=rho, alpha=alpha, mu=mu, s=1, t=2)}), lower=0, upper=1)$value, error=function(e) -1)
          edg13 <- tryCatch(integrate(Vectorize(function(x){edges_skewt_3d(w=c(x,0,1-x), rho=rho, alpha=alpha, mu=mu, s=1, t=3)}), lower=0, upper=1)$value, error=function(e) -1)
          edg23 <- tryCatch(integrate(Vectorize(function(x){edges_skewt_3d(w=c(0,x,1-x), rho=rho, alpha=alpha, mu=mu, s=2, t=3)}), lower=0, upper=1)$value, error=function(e) -1)
          
          if(any(c(edg12, edg13, edg23)==-1)){
            hdens <- 0
          }else{
            int01 <- 3 - corners_skewt_3d(rho=rho, alpha=alpha, mu=mu, s=1) - corners_skewt_3d(rho=rho, alpha=alpha, mu=mu, s=2) - corners_skewt_3d(rho=rho, alpha=alpha, mu=mu, s=3) - edg12 - edg13 - edg23
            intc <- tryCatch(integrate(Vectorize(function(y) integrate(Vectorize(function(x) interior_skewt_d(c(x,y,1-x-y),rho=rho, alpha=alpha, mu=mu)), lower=c, upper=1-y-c )$value), lower=c, upper=1-2*c)$value, error=function(e) -1)
            if(intc==-1){
              hdens <- 0
            }else{
              K <- int01/intc
              hdens <- K * interior_skewt_d(w=data, rho=rho, alpha=alpha, mu=mu)
            }
          }
          
        }
      }
      return(hdens)
    }
    
    ###
    ### Multiple observation in matrix form
    ###
    
    if(is.matrix(data) && ncol(data)==3){
      
      hdens <- vector(length=nrow(data))
      
      if(c == 0){
        hdens <- apply(data,1, function(x) interior_skewt_d(w=x, rho=rho, alpha=alpha, mu=mu))
      }else if(c>0){
        subs <- apply(data, 1, function(x) subset.c(x,c=c))
        
        if(is.matrix(subs) && nrow(subs)==3){ # all points in the interior (c is too small)
          hdens <- apply(data,1, function(x) interior_skewt_d(w=x, rho=rho, alpha=alpha, mu=mu))
        }
        
        if(is.list(subs)){
          
          ## Corners
          
          c1 <- which(lapply(subs, function(x) (length(x)==1 && any(x ==1) )) == TRUE)
          c2 <- which(lapply(subs, function(x) (length(x)==1 && any(x ==2) )) == TRUE)
          c3 <- which(lapply(subs, function(x) (length(x)==1 && any(x ==3) )) == TRUE)
          K.c1 <- K.c2 <- K.c3 <- sqrt(3) / c^2
          
          hdens[c1] <-  K.c1 * corners_skewt_3d(rho=rho, alpha=alpha, mu=mu, s=1)
          hdens[c2] <-  K.c2 * corners_skewt_3d(rho=rho, alpha=alpha, mu=mu, s=2)
          hdens[c3] <-  K.c3 * corners_skewt_3d(rho=rho, alpha=alpha, mu=mu, s=3)
          
          ## Edges
          
          #edge_surface <- c - 4/sqrt(3)*c^2
          
          ## Edge {1,2}
          
          e12 <- which(lapply(subs, function(x) (length(x)==2 && any(x ==1) && any(x ==2) )) == TRUE)
          edg12 <- tryCatch(integrate(Vectorize(function(x){edges_skewt_3d(w=c(x,1-x,0), rho=rho, alpha=alpha, mu=mu, s=1, t=2)}), lower=0, upper=1)$value, error=function(e) -1)
          edge_surface12 <- tryCatch(integrate(Vectorize(function(x){edges_skewt_3d(w=c(x,1-x,0), rho=rho, alpha=alpha, mu=mu, s=1, t=2)}), lower=c/(2*(1-c/2)), upper=(1-c)/(1-c/2))$value, error=function(e) -1)
          if(any(c(edg12, edge_surface12)==-1)){
            hdens[e12] <- rep(0, length(e12))
          }else{
            K.e12 <- edg12 / edge_surface12
            if(length(e12)==1){
              hdens[e12] <- K.e12 * edges_skewt_3d(w=data[e12,], rho=rho, alpha=alpha, mu=mu, s=1, t=2)
            }else if(length(e12)>1){
              hdens[e12] <- K.e12 * apply(data[e12,],1, function(x) edges_skewt_3d(w=x, rho=rho, alpha=alpha, mu=mu, s=1, t=2))
            }
          }
          
          ## Edge {1,3}
          
          e13 <- which(lapply(subs, function(x) (length(x)==2 && any(x ==1) && any(x ==3) )) == TRUE)
          edg13 <- tryCatch(integrate(Vectorize(function(x){edges_skewt_3d(w=c(x,0,1-x), rho=rho, alpha=alpha, mu=mu, s=1, t=3)}), lower=0, upper=1)$value, error=function(e) -1)
          edge_surface13 <- tryCatch(integrate(Vectorize(function(x){edges_skewt_3d(w=c(x,0,1-x), rho=rho, alpha=alpha, mu=mu, s=1, t=3)}), lower=c/(2*(1-c/2)), upper=(1-c)/(1-c/2))$value, error=function(e) -1)
          if(any(c(edg13, edge_surface13)==-1)){
            hdens[e13] <- rep(0, length(e13))
          }else{
            K.e13 <- edg13 / edge_surface13
            if(length(e13)==1){
              hdens[e13] <- K.e13 * edges_skewt_3d(w=data[e13,], rho=rho, alpha=alpha, mu=mu, s=1, t=3)
            }else if(length(e13)>1){
              hdens[e13] <- K.e13 * apply(data[e13,],1, function(x) edges_skewt_3d(w=x, rho=rho, alpha=alpha, mu=mu, s=1, t=3))
            }
          }
          
          
          ## Edge {2,3}
          
          e23 <- which(lapply(subs, function(x) (length(x)==2 && any(x ==2) && any(x ==3) )) == TRUE)
          edg23 <- tryCatch(integrate(Vectorize(function(x){edges_skewt_3d(w=c(0,x,1-x), rho=rho, alpha=alpha, mu=mu, s=2, t=3)}), lower=0, upper=1)$value, error=function(e) -1)
          edge_surface23 <- tryCatch(integrate(Vectorize(function(x){edges_skewt_3d(w=c(0,x,1-x), rho=rho, alpha=alpha, mu=mu, s=2, t=3)}), lower=c/(2*(1-c/2)), upper=(1-c)/(1-c/2))$value, error=function(e) -1)
          if(any(c(edg23, edge_surface23)==-1)){
            hdens[e23] <- rep(0, length(e23))
          }else{
            K.e23 <- edg23 / edge_surface23
            if(length(e23)==1){
              hdens[e23] <- K.e23 * edges_skewt_3d(w=data[e23,], rho=rho, alpha=alpha, mu=mu, s=2, t=3)
            }else if(length(e23)>1){
              hdens[e23] <- K.e23 * apply(data[e23,],1, function(x) edges_skewt_3d(w=x, rho=rho, alpha=alpha, mu=mu, s=2, t=3))
            }
          }
          
          ## Interior
          
          inter <- which(lapply(subs, function(x) (length(x)==3 )) == TRUE)
          if(any(c(edg12, edg13, edg23)==-1)){
            hdens[inter] <- rep(0, length(inter))
          }else{
            int01 <- 3 - corners_skewt_3d(rho=rho, alpha=alpha, mu=mu, s=1) - corners_skewt_3d(rho=rho, alpha=alpha, mu=mu, s=2) - corners_skewt_3d(rho=rho, alpha=alpha, mu=mu, s=3) - edg12 - edg13 - edg23
            intc <- tryCatch(integrate(Vectorize(function(y) integrate(Vectorize(function(x) interior_skewt_d(c(x,y,1-x-y),rho=rho, alpha=alpha, mu=mu)), lower=c, upper=1-y-c )$value), lower=c, upper=1-2*c)$value, error=function(e) -1)
            if(intc==-1){
              hdens[inter] <- rep(0, length(inter))
            }else{
              K.inter <- int01/intc
              if(length(inter)==1){
                hdens[inter] <- K.inter * interior_skewt_d(w=data[inter,], rho=rho,  alpha=alpha, mu=mu)
              }else if(length(inter)>1){
                hdens[inter] <- K.inter * apply(data[inter,],1, function(x) interior_skewt_d(w=x, rho=rho, alpha=alpha, mu=mu))
              }
            }
            
          }
          
        }
        
      }
      
      return(hdens)
      
    }
    
  }
  
  ### Angular density for the Extremal-t model on the 2 and 3 dimensional simplex
  
  xvect = as.double(as.vector(t(x)))
  if (is.vector(x)) {
    dim = as.integer(length(x))
    n = as.integer(1)
    if(round(sum(x),7) != 1){ stop("Data is not angular")}
    if(dim==2){ result <- dens_skewt_2d(data=x, rho=rho, alpha=alpha, mu=mu, c=c)}
    if(dim==3){ result <- dens_skewt_3d(data=x, rho=rho, alpha=alpha, mu=mu, c=c)}
  }
  else {
    dim = as.integer(ncol(x))
    n = as.integer(nrow(x))
    
    if(any(round(apply(x,1,sum),7) != 1)){ stop("Data is not angular") }
    
    if(dim==2){ result <- dens_skewt_2d(data=x, rho=rho, alpha=alpha, mu=mu, c=c)}
    if(dim==3){ result <- dens_skewt_3d(data=x, rho=rho, alpha=alpha, mu=mu, c=c)}
    
    if(!vectorial){
      result <- prod(result)
    }
    
  }
  if(any(is.na(result))){
    result[is.na(result)] <- 0
  }
  if(log){
    if(any(result==0)){
      ind0 <- which(result==0)
      ind <- which(result!=0)
      result[ind0] <- -1e+300
      result[ind] <- log(result[ind])
      return(result)
    }else{
      return(log(result))
    }
  }else{
    return(result)
  }
}

####################################################################
####################################################################
####################################################################
### END Internal functions for PARAMETRIC and ANGULAR
####################################################################
####################################################################
####################################################################

####################################################################
####################################################################
####################################################################
### Internal functions for PARAMETRIC and NOT ANGULAR
####################################################################
####################################################################
####################################################################

dHuslerReiss <- function(x, lambda=1){
  if(any(is.na(x)))
    return(NA)
  .C("dHuslerReiss", as.double(x), as.double(lambda), out=double(1), NAOK=TRUE)$out
}

dmextst <- function(x, scale=1, shape=rep(0,length(x)), df=1){
  if(any(is.na(x)))
    return(NA)
  .C("dmextst", as.double(x), as.double(scale), as.double(df), as.double(shape), out=double(1), NAOK=TRUE)$out
}

####################################################################
####################################################################
####################################################################
### END Internal functions for PARAMETRIC and NOT ANGULAR
####################################################################
####################################################################
####################################################################

####################################################################
####################################################################
####################################################################
### Internal functions for NON-PARAMETRIC and ANGULAR
####################################################################
####################################################################
####################################################################

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

####################################################################
####################################################################
####################################################################
### END Internal functions for NON-PARAMETRIC and ANGULAR
####################################################################
####################################################################
####################################################################
