pExtDep <- function(q, type, method="Parametric", model, par, plot=TRUE, 
                    main, xlab, cex.lab, cex.axis, lwd,...){
  
  # Checking quantile vector or matrix
  if(is.vector(q)){
    dim <- length(q)
    q0 <- q
    if(is.matrix(par)){
      npar <- nrow(par)
      p <- vector(length=npar)
      if(any(q==0)){return(rep(NA,npar))}
    }else if(is.vector(par)){
      if(any(q==0)){return(NA)}
    }
  }else if(is.matrix(q)){
    dim <- ncol(q)
    qNA <- apply(q,1,function(x) any(x==0))
    
    q0 <- q
    nq0 <- nrow(q0)
    q <- q[!qNA,]
    if(is.matrix(q)){
      nq <- nrow(q)      
    }else if(is.vector(q)){
      nq <- 1
    }

    if(is.vector(par)){
      p <- vector(length=nq)      
    }else if(is.matrix(par)){
      p <- matrix(nrow=nrow(par), ncol=nq)
    }

  }else{
    stop("q must be a vector or a matrix")
  }
  
  # Checking method
  methods <- c("Parametric", "NonParametric")
  if(!any(method ==  methods)){ stop("Wrong method specified")}
  
  # Checking the type
  types <- c("lower", "inv.lower", "upper")
  if(method=="Parametric" && !any(type ==  types)){ stop("Wrong probability type specified")}
  
  if(method=="Parametric"){
    
    if(all(dim != c(2,3)) ){ stop("Dimension 2 or 3 only")}
    
    # Checking the model
    models <- c("HR", "ET", "EST")
    if(!any(model ==  models)){ stop("model wrongly specified")}
    
    # Checking model, par & dim are correctly specified
    if(is.vector(par)){
      if(!dim_ExtDep(model=model, par=par, dim=dim)){stop("Length of 'par' is incorrect")}      
    }else if(is.matrix(par)){
      if(!dim_ExtDep(model=model, par=par[1,], dim=dim)){stop("Length of 'par' is incorrect")} 
    }else{
      stop("'par' should be a vector or a matrix")
    }
    
    if(model == "HR"){
      if(is.vector(q)){
        if(is.vector(par)){
          p <- p.hr(q=q, par=par, type=type)          
        }else if(is.matrix(par)){
          p <- apply(par,1, p.hr, q=q, type=type)
        }
      }else if(is.matrix(q)){
        if(is.vector(par)){
          p <- apply(q,1, p.hr, par=par, type=type)
        }else if(is.matrix(par)){
          for(i in 1:nq){
            p[,i] <- apply(par,1, p.hr, q=q[i,], type=type)
          }
        }
      }
    }else if(model == "ET"){
      if(is.vector(q)){
        if(is.vector(par)){
          p <- p.et(q=q, par=par, type=type)          
        }else if(is.matrix(par)){
          p <- apply(par,1, p.et, q=q, type=type)
        }
        
      }else if(is.matrix(q)){
        if(is.vector(par)){
          p <- apply(q,1, p.et, par=par, type=type)
        }else if(is.matrix(par)){
          for(i in 1:nq){
            p[,i] <- apply(par,1, p.et, q=q[i,], type=type)
          }
        }
      }
    }else if(model == "EST"){
      if(is.vector(q)){
        if(is.vector(par)){
          p <- p.est(q=q, par=par, type=type)          
        }else if(is.matrix(par)){
          p <- apply(par,1, p.est, q=q, type=type)
        }
    
      }else if(is.matrix(q)){
        if(is.vector(par)){
          p <- apply(q,1, p.est, par=par, type=type)          
        }else if(is.matrix(par)){
          for(i in 1:nq){
            p[,i] <- apply(par,1, p.est, q=q[i,], type=type)
          }
        }
      }
    }
  }else if(method == "NonParametric"){
    
    if(is.vector(q)){
      if(is.vector(par)){
        p <- ph(w=q, beta=par)        
      }else if(is.matrix(par)){
        for(i in 1:npar){
          p[i] <- ph(w=q, beta=par[i,])
        }
      }
       
    }else if(is.matrix(q)){
      if(is.vector(par)){
        for(i in 1:nq){
          p[i] <- ph(w=q[i,], beta=par)   
        }        
      }else if(is.matrix(par)){
        for(i in 1:nq){
          for(j in 1:npar){
            p[j,i] <- ph(w=q[i,], beta=par[j,])
          }
        }
      }

    }
  }
  
  if(is.matrix(par)){
    if(plot){
      
      if(missing(main) ){main <- ""}
      if(missing(cex.lab)){cex.lab <- 1.4}
      if(missing(cex.axis)){cex.axis <- 1.4}
      if(missing(lwd)){lwd <- 2}
      
      if(is.vector(p)){
        
        if(missing(xlab)){
          if(dim==2){
            if(type=="lower"){
              xlab <- paste("P(X<",q[1],",Y<",q[2], ")", sep="")
            }else if(type=="inv.lower"){
              xlab <- paste("1-P(X<",q[1],",Y<",q[2], ")", sep="")
            }else if(type=="upper"){
              xlab <- paste("P(X>",q[1],",Y>",q[2], ")", sep="")              
            }
            
          }else if(dim==3){
            if(type=="lower"){
              xlab <- paste("P(X<",q[1],",Y<",q[2],",Z<",q[3], ")", sep="") 
            }else if(type=="inv.lower"){
              xlab <- paste("1-P(X<",q[1],",Y<",q[2],",Z<",q[3], ")", sep="") 
            }else if(type=="upper"){
              xlab <- paste("P(X>",q[1],",Y>",q[2],",Z>",q[3], ")", sep="")              
            }

          }
        }
        
        Ke <- density(p) # KDE
        Hi <- hist(p, prob=TRUE, col="lightgrey",
                   ylim=range(Ke$y), main=main, xlab=xlab ,
                   cex.lab=cex.lab, cex.axis=cex.axis, lwd=lwd, ...)
        p_ic <- quantile(p, probs=c(0.025, 0.5, 0.975))
        
        points(x=p_ic, y=c(0,0,0), pch=4, lwd=4)
        points(x=mean(p), y=0, pch=16, lwd=4)
        lines(Ke, lwd = 2, col = "dimgrey")
        
      }else if(is.matrix(p)){
        for(i in 1:nq){
          
          if(missing(xlab)){
            if(dim==2){
              if(type=="lower"){
                xlab <- paste("P(X<",q[i,1],",Y<",q[i,2], ")", sep="")
              }else if(type=="inv.lower"){
                xlab <- paste("1-P(X<",q[i,1],",Y<",q[i,2], ")", sep="")
              }else if(type=="upper"){
                xlab <- paste("P(X>",q[i,1],",Y>",q[i,2], ")", sep="")              
              }
              
            }else if(dim==3){
              if(type=="lower"){
                xlab <- paste("P(X<",q[i,1],",Y<",q[i,2],",Z<",q[i,3], ")", sep="") 
              }else if(type=="inv.lower"){
                xlab <- paste("1-P(X<",q[i,1],",Y<",q[i,2],",Z<",q[i,3], ")", sep="") 
              }else if(type=="upper"){
                xlab <- paste("P(X>",q[i,1],",Y>",q[i,2],",Z>",q[i,3], ")", sep="")              
              }
              
            }
          }
          
          Ke <- density(p[,i]) # KDE
          Hi <- hist(p[,i], prob=TRUE, col="lightgrey",
                     ylim=range(Ke$y), main=main, xlab=xlab ,
                     cex.lab=cex.lab, cex.axis=cex.axis, lwd=lwd, ...)
          p_ic <- quantile(p[,i], probs=c(0.025, 0.5, 0.975))
          
          points(x=p_ic, y=c(0,0,0), pch=4, lwd=4)
          points(x=mean(p[,i]), y=0, pch=16, lwd=4)
          lines(Ke, lwd = 2, col = "dimgrey")
        }
      }
    }
  }
  
  if(is.matrix(q0)){
    
    if(sum(qNA) !=0){
      p.temp <- p
      
      if(is.vector(par)){
        p <- vector(length=nq0)
        p[!qNA] <- p.temp
        p[qNA] <- rep(NA, sum(qNA))
      }else if(is.matrix(par)){
        p <- matrix(nrow=nrow(par), ncol=nq0)
        p[,!qNA] <- p.temp
        #p[,qNA] <- matrix(NA, nrow=nrow(par), ncol=sum(qNA))
      }

    }
    
  }
  
  return(p)
}

pFailure <- function(n, beta, u1, u2, mar1, mar2, type, plot, xlab, ylab, nlevels=10){
  
  if(!(type %in% c("and", "or", "both")) ){stop(" 'type' must be 'and', 'or' or 'both'.")}
  
  ### Simulate data
  
  w <- rh(n, beta)
  r <- runif(n)^(-1)
  y <- cbind(2*r*w, 2*r*(1-w))
  
  ### Compute ustar       
  
  u.star <- function(threshold, mar){
    if(mar[3] == 0){
      return( exp((threshold - mar[1])/mar[2]) )
    }else{
      return( (1 + mar[3]*(threshold-mar[1])/mar[2])^(1/mar[3]) )
    }
  }
  
  ustar1 <- u.star(u1, mar1)
  ustar2 <- u.star(u2, mar2)
  
  ustar <- expand.grid(ustar1, ustar2)
  
  if(type == "and"){
    AND <- apply(ustar, 1, function(u) sum(y[,1] > u[1] & y[,2] > u[2]) )/n
    AND <- matrix(AND, nrow=length(u1), ncol=length(u2), byrow=FALSE)
  }
  if(type == "or"){
    OR <- apply(ustar, 1, function(u) sum(y[,1] > u[1] | y[,2] > u[2]) )/n
    OR <- matrix(OR, nrow=length(u1), ncol=length(u2), byrow=FALSE)
  }
  if(type == "both"){
    and_or <- function(u, y){
      c(sum(y[,1] > u[1] & y[,2] > u[2]), sum(y[,1] > u[1] | y[,2] > u[2]))
    }
    
    AND_OR <- apply(ustar,1, and_or, y=y) / n
    AND <- matrix(AND_OR[1,], nrow=length(u1), ncol=length(u2), byrow=FALSE)
    OR <- matrix(AND_OR[2,], nrow=length(u1), ncol=length(u2), byrow=FALSE)
  }
  
  
  if(plot){
    
    op <- par(mai = c(0.8, 0.8, 0.35, 0.1), mgp = c(2.5, 1, 0), cex.axis=2, cex.lab=2, cex.main=2)
    on.exit(par(op))
    
    if(type == "and"){
      contour(u1, u2, AND, levels=round(seq(min(AND), max(AND), length=nlevels), 3), 
              xlab=xlab, ylab=ylab, main="P(X>x and Y>y)", labcex=1)
    }
    if(type == "or"){
      contour(u1, u2, OR, levels=round(seq(min(OR), max(OR), length=nlevels), 3), 
              xlab=xlab, ylab=ylab, main="P(X>x or Y>y)", labcex=1)
    }
    if(type == "both"){
      contour(u1, u2, AND, levels=round(seq(min(AND), max(AND), length=nlevels), 3), 
              xlab=xlab, ylab=ylab, main="P(X>x and Y>y)", labcex=1)
      
      contour(u1, u2, OR, levels=round(seq(min(OR), max(OR), length=nlevels), 3), 
              xlab=xlab, ylab=ylab, main="P(X>x or Y>y)", labcex=1)  
    }
    
  }
  
  if(type == "and"){
    return(list(AND=AND))
  }
  if(type == "or"){
    return(list(OR=OR))  }
  if(type == "both"){
    return(list(AND=AND, OR=OR))
  }
  
}

#########################################
#########################################
### Internal functions for pExtDep
#########################################
#########################################

## Husler-Reiss model

p.hr <- function(q, par, type){
  
  p.hr.2d <- function(q, par, type){
    
    I1 <- pnorm(par+log(q[2]/q[1])/(2*par))
    I2 <- pnorm(par+log(q[1]/q[2])/(2*par))
    
    if(type == "lower"){
      return( exp(-1/q[1] * I1 - 1/q[2] * I2 ))
    }else if(type == "inv.lower"){
      return( 1 - exp(-1/q[1] * I1 - 1/q[2] * I2 ))
    }else if(type == "upper"){
      return( 1 - exp(-1/q[1]) - exp(-1/q[2]) + exp(-1/q[1] * I1 - 1/q[2] * I2 ) )
    }
  }
  
  p.hr.3d <- function(q, par, type){
    
    lambda12 <- par[1]; lambda13 <- par[2]; lambda23 <- par[3];
    S1 <- matrix(c(4*lambda12^2, rep(2*(lambda12^2 + lambda13^2 - lambda23^2),2), 4*lambda13^2),nrow=2)
    S2 <- matrix(c(4*lambda12^2, rep(2*(lambda12^2 + lambda23^2 - lambda13^2),2), 4*lambda23^2),nrow=2)
    S3 <- matrix(c(4*lambda13^2, rep(2*(lambda13^2 + lambda23^2 - lambda12^2),2), 4*lambda23^2),nrow=2)
    
    I1 <- pmesn(x=c(log(q[2]/q[1])+2*lambda12^2, log(q[3]/q[1])+2*lambda13^2),scale = S1 )
    I2 <- pmesn(x=c(log(q[1]/q[2])+2*lambda12^2, log(q[3]/q[2])+2*lambda23^2),scale = S2 )
    I3 <- pmesn(x=c(log(q[1]/q[3])+2*lambda13^2, log(q[2]/q[3])+2*lambda23^2),scale = S3 )
    
    if(type == "lower"){
      return( exp(-1/q[1] * I1 - 1/q[2] * I2 - 1/q[3] * I3 ))
    }else if(type == "inv.lower"){
      return( 1 - exp(-1/q[1] * I1 - 1/q[2] * I2 - 1/q[3] * I3 ))
    }else if(type == "upper"){
      return( 1 - exp(-1/q[1]) - exp(-1/q[2]) - exp(-1/q[3])
              + p.hr.2d(q[1:2], lambda12, type="lower") + p.hr.2d(q[c(1,3)], lambda13, type="lower")
              + p.hr.2d(q[2:3], lambda23, type="lower") - exp(-1/q[1] * I1 - 1/q[2] * I2 - 1/q[3] * I3 ))
    }
  }
  
  if(any(par<=0)){stop('wrong value of parameters')}
  dim <- length(q)
  if(!dim_ExtDep(model="HR", par=par, dim=dim )){stop("Length of par is incorrect")}
  if(dim == 2){
    return( p.hr.2d(q=q, par=par, type=type))
  }else if(dim == 3){
    return( p.hr.3d(q=q, par=par, type=type))
  }
}

## Extremal-t model

p.et <- function(q, par, type){
  
  p.et.2d <- function(q, par, type){
    
    if( ( (par[1] <= -1) || (par[1] >=1) || (par[2] <= 0))==1){stop("Extremal-t parameters ill defined")}
    rho <- par[1]; nu <- par[2];
    
    I1 <- pt(sqrt((nu+1)/(1-rho^2))*((q[2]/q[1])^(1/nu)-rho), df=nu+1)
    I2 <- pt(sqrt((nu+1)/(1-rho^2))*((q[1]/q[2])^(1/nu)-rho), df=nu+1)
    
    if(type == "lower"){
      return( exp(-1/q[1] * I1 - 1/q[2] * I2 ))
    }else if(type == "inv.lower"){
      return( 1 - exp(-1/q[1] * I1 - 1/q[2] * I2 ))
    }else if(type == "upper"){
      return( 1/q[1] * (1-I1) + 1/q[2] * (1-I2) )
    }
  }
  
  p.et.3d <- function(q, par, type){
    
    if( (any(par[1:3] <= -1) || any(par[1:3] >=1) || (par[4] <= 0))==1){stop("Extremal-t parameters ill defined")}
    
    rho12 <- par[1]; rho13 <- par[2]; rho23 <- par[3]; nu <- par[4];
    R1 <- (rho23-rho12*rho13)/sqrt((1-rho12^2)*(1-rho13^2))
    R2 <- (rho13-rho12*rho23)/sqrt((1-rho12^2)*(1-rho23^2))
    R3 <- (rho12-rho13*rho23)/sqrt((1-rho13^2)*(1-rho23^2))
    if(sum(eigen(matrix(c(1,R1,R1,1),ncol=2))$values<0)>=1){stop("negative definite partial correlation matrix")}
    if(sum(eigen(matrix(c(1,R2,R2,1),ncol=2))$values<0)>=1){stop("negative definite partial correlation matrix")}
    if(sum(eigen(matrix(c(1,R3,R3,1),ncol=2))$values<0)>=1){stop("negative definite partial correlation matrix")}
    
    x1 <- sqrt((nu+1)/(1-rho12^2))*((q[2]/q[1])^(1/nu)-rho12)
    y1 <- sqrt((nu+1)/(1-rho13^2))*((q[3]/q[1])^(1/nu)-rho13)
    x2 <- sqrt((nu+1)/(1-rho12^2))*((q[1]/q[2])^(1/nu)-rho12)
    y2 <- sqrt((nu+1)/(1-rho23^2))*((q[3]/q[2])^(1/nu)-rho23)
    x3 <- sqrt((nu+1)/(1-rho13^2))*((q[1]/q[3])^(1/nu)-rho13)
    y3 <- sqrt((nu+1)/(1-rho23^2))*((q[2]/q[3])^(1/nu)-rho23)
    
    I1 <- pmest(x=c(x1,y1), scale=matrix(c(1,R1,R1,1),ncol=2), df=nu+1)
    I2 <- pmest(x=c(x2,y2), scale=matrix(c(1,R2,R2,1),ncol=2), df=nu+1)
    I3 <- pmest(x=c(x3,y3), scale=matrix(c(1,R3,R3,1),ncol=2), df=nu+1)
    
    if(type == "lower"){
      return( exp(-1/q[1] * I1 - 1/q[2] * I2 - 1/q[3] * I3 ))
    }else if(type == "inv.lower"){
      return( 1 - exp(-1/q[1] * I1 - 1/q[2] * I2 - 1/q[3] * I3 ))
    }else if(type == "upper"){
      return( 1 - exp(-1/q[1]) - exp(-1/q[2]) - exp(-1/q[3])
              + p.et.2d(q[1:2], c(rho12,nu), type="lower") + p.et.2d(q[c(1,3)], c(rho13,nu), type="lower")
              + p.et.2d(q[2:3], c(rho23,nu), type="lower") - exp(-1/q[1] * I1 - 1/q[2] * I2 - 1/q[3] * I3 ))
    }
  }
  
  dim <- length(q)
  if(!dim_ExtDep(model="ET", par=par, dim=dim )){stop("Length of par is incorrect")}
  if(dim == 2){
    return( p.et.2d(q=q, par=par, type=type))
  }else if(dim == 3){
    return( p.et.3d(q=q, par=par, type=type))
  }
}

## Extremal Skew-t

pmextst <- function(x, scale=1, shape=rep(0,length(x)), df=1){
  if(any(is.na(x)))
    return(NA)
  .C("pmextst", as.double(x), as.double(scale), as.double(df),
     as.double(shape), out=double(1), NAOK=TRUE)$out
}

p.est <- function(q, par, type){
  
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
  
  alpha_star_j <- function(rho,alpha,j){
    Alpha_tilde <- alpha_tilde(alpha,j)
    sj <- s_j(rho,j)
    return( Alpha_tilde %*% sj )
  }
  
  tau_star_j <- function(rho,alpha,nu,j){
    Sig <- Sigma(rho)
    Alpha_tilde <- alpha_tilde(alpha,j)
    return( sqrt(nu+1) * (alpha[j] + Sig[-j,j] %*% t(Alpha_tilde) )     )
  }
  
  nu_p <- function(rho,alpha,nu,j){
    Alpha_bar <- alpha_bar_j(rho,alpha,j)
    return( pt(sqrt(nu+1)*Alpha_bar,df=nu+1) * 2^(nu/2-1) * gamma((nu+1)/2) / sqrt(pi) )
  }
  
  x_bar <- function(x,rho,alpha,nu,j){
    return( x * nu_p(rho,alpha,nu,j) )
  }
  
  
  p.est.2d <- function(q, par, type){
    
    dim <- length(q)
    if(dim != 2){stop("Dimension should be 2")}
    
    rho <- par[1]
    alpha <- par[2:3]
    nu <- par[4]
    
    if( ( (rho <= -1) || (rho >=1) || (nu <= 0))==1){stop("Extremal-t parameters ill defined")}
    
    nu_1 <- nu_p(rho=rho, alpha=alpha, nu=nu, j=1)
    nu_2 <- nu_p(rho=rho, alpha=alpha, nu=nu, j=2)
    
    x1_bar <- q[1] * nu_1
    x2_bar <- q[2] * nu_2
    
    comp1_1 <- sqrt((nu+1)/(1-rho^2))*((x2_bar/x1_bar)^(1/nu)-rho)
    comp1_2 <- sqrt((nu+1)/(1-rho^2))*((x1_bar/x2_bar)^(1/nu)-rho)
    
    alpha_star_1 <- alpha_star_j(rho=rho, alpha=alpha, j=1)
    alpha_star_2 <- alpha_star_j(rho=rho, alpha=alpha, j=2)
    
    tau_star_1 <- tau_star_j(rho=rho, alpha=alpha, nu=nu, j=1)
    tau_star_2 <- tau_star_j(rho=rho, alpha=alpha, nu=nu, j=2)
    
    I1 <- pest(x=comp1_1,scale=1,shape=alpha_star_1,extended=tau_star_1,df=nu+1)
    I2 <- pest(x=comp1_2,scale=1,shape=alpha_star_2,extended=tau_star_2,df=nu+1)
    
    if(type == "lower"){
      return( exp(-1/q[1] * I1 - 1/q[2] * I2 ))
    }else if(type == "inv.lower"){
      return( 1 - exp(-1/q[1] * I1 - 1/q[2] * I2 ))
    }else if(type == "upper"){
      return( 1/q[1] * (1-I1) + 1/q[2] * (1-I2) )
    }
  }
  
  p.est.3d <- function(q, par, type){
    
    dim <- length(q)
    if(dim != 3){stop("Dimension should be 3")}
    
    rho <- par[1:3]
    alpha <- par[4:6]
    nu <- par[7]
    
    if( (any(rho <= -1) || any(rho >=1) || (nu <= 0))==1){stop("Extremal-t parameters ill defined")}
    
    Sigma_bar_1 <- Sigma_bar_j(rho=rho, j=1)
    Sigma_bar_2 <- Sigma_bar_j(rho=rho, j=2)
    Sigma_bar_3 <- Sigma_bar_j(rho=rho, j=3)
    
    if(sum(eigen(Sigma_bar_1)$values<0)>=1){stop("negative definite partial correlation matrix")}
    if(sum(eigen(Sigma_bar_2)$values<0)>=1){stop("negative definite partial correlation matrix")}
    if(sum(eigen(Sigma_bar_3)$values<0)>=1){stop("negative definite partial correlation matrix")}
    
    nu_1 <- nu_p(rho=rho, alpha=alpha, nu=nu, j=1)
    nu_2 <- nu_p(rho=rho, alpha=alpha, nu=nu, j=2)
    nu_3 <- nu_p(rho=rho, alpha=alpha, nu=nu, j=3)
    
    x1_bar <- q[1] * nu_1
    x2_bar <- q[2] * nu_2
    x3_bar <- q[3] * nu_3
    
    comp1_1 <- sqrt((nu+1)/(1-rho[1]^2))*((x2_bar/x1_bar)^(1/nu)-rho[1])
    comp2_1 <- sqrt((nu+1)/(1-rho[2]^2))*((x3_bar/x1_bar)^(1/nu)-rho[2])
    comp1_2 <- sqrt((nu+1)/(1-rho[1]^2))*((x1_bar/x2_bar)^(1/nu)-rho[1])
    comp2_2 <- sqrt((nu+1)/(1-rho[3]^2))*((x3_bar/x2_bar)^(1/nu)-rho[3])
    comp1_3 <- sqrt((nu+1)/(1-rho[2]^2))*((x1_bar/x3_bar)^(1/nu)-rho[2])
    comp2_3 <- sqrt((nu+1)/(1-rho[3]^2))*((x2_bar/x3_bar)^(1/nu)-rho[3])
    
    alpha_star_1 <- alpha_star_j(rho=rho, alpha=alpha, j=1)
    alpha_star_2 <- alpha_star_j(rho=rho, alpha=alpha, j=2)
    alpha_star_3 <- alpha_star_j(rho=rho, alpha=alpha, j=3)
    
    tau_star_1 <- tau_star_j(rho=rho, alpha=alpha, nu=nu, j=1)
    tau_star_2 <- tau_star_j(rho=rho, alpha=alpha, nu=nu, j=2)
    tau_star_3 <- tau_star_j(rho=rho, alpha=alpha, nu=nu, j=3)
    
    I1 <- pmest(x=c(comp1_1,comp2_1),scale=Sigma_bar_1,shape=alpha_star_1,extended=tau_star_1,df=nu+1)
    I2 <- pmest(x=c(comp1_2,comp2_2),scale=Sigma_bar_2,shape=alpha_star_2,extended=tau_star_2,df=nu+1)
    I3 <- pmest(x=c(comp1_3,comp2_3),scale=Sigma_bar_3,shape=alpha_star_3,extended=tau_star_3,df=nu+1)
    
    if(type == "lower"){
      return( exp(-1/q[1] * I1 - 1/q[2] * I2 - 1/q[3] * I3 ))
    }else if(type == "inv.lower"){
      return( 1 - exp(-1/q[1] * I1 - 1/q[2] * I2 - 1/q[3] * I3 ))
    }else if(type == "upper"){
      return( 1 - exp(-1/q[1]) - exp(-1/q[2]) - exp(-1/q[3])
              + p.est.2d(q[1:2], c(rho[1], alpha[1:2],nu), type="lower") + p.est.2d(q[c(1,3)], c(rho[2], alpha[c(1,3)],nu), type="lower")
              + p.est.2d(q[2:3], c(rho[3], alpha[2:3],nu), type="lower") - exp(-1/q[1] * I1 - 1/q[2] * I2 - 1/q[3] * I3 ))
    }
  }
  
  dim <- length(q)
  if(!dim_ExtDep(model="EST", par=par, dim=dim )){stop("Length of par is incorrect")}
  if(dim == 2){
    return( p.est.2d(q=q, par=par, type=type))
  }else if(dim == 3){
    return( p.est.3d(q=q, par=par, type=type))
  }
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


#########################################
#########################################
### END Internal functions for pExtDep
#########################################
#########################################
