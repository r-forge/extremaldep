plot_ExtDep <- function(object="angular", model, par, log=TRUE, data=NULL, 
                        contour=TRUE, style, labels, cex.dat=1, cex.lab=1, 
                        cex.cont=1, Q.fix, Q.range, Q.range0, cond=FALSE,...){
  
  Objects <- c("angular", "pickands", "returns")
  if(!any(object ==  Objects)){ stop("object is wrongly specified")}
  
  if(is.vector(par)){
    if(model=="AL"){
      dim <- length(par)-1
    }else{
      dim <- dim_ExtDep(model=model, par=par)      
    }
  }else if(is.matrix(par)){
    if(model=="AL"){
      dim <- ncol(par)-1
    }else{
      dim <- dim_ExtDep(model=model, par=par[1,])         
    }
  }
  
  if(missing(labels)){
    if(dim==2){
      labels="w"
    }
    if(dim==3){
      labels=c(expression(w[1]),expression(w[3]),expression(w[2]))
    }
  }
  
  if(object=="angular"){
    
    if(!is.vector(par)){stop("'par' must be a vector")}
    
    models <- c("PB", "HR", "ET", "EST", "TD", "AL")
    if(!any(model ==  models)){ stop("model wrongly specified")}
    
    if(dim == 2){
      AngDens2dPlot(model=model, para=par, log=log, data=data, style=style,
                    labels=labels, cex.lab=cex.lab, ...)
    }else if(dim == 3){
      AngDens3dPlot(model=model, para=par, log=log, data=data, contour=contour, 
                    labels=labels, cex.dat=cex.dat, cex.lab=cex.lab, cex.cont=cex.cont)
    }else{
      stop("Angular density plotted only for 2 or 3 dimensional models.")
    }
    
  }else if (object=="pickands"){
    
    if(!is.vector(par)){stop("'par' must be a vector")}
    
    models <- c("HR", "ET", "EST")
    if(!any(model ==  models)){ stop("model wrongly specified")}
    
    if(dim == 2){
      Pick2dPlot(model=model, para=par, label=labels, cex.lab=cex.lab,...)
    }else if(dim == 3){
      Pick3dPlot(model=model, para=par, contour=contour, 
                 labels=labels, cex.lab=cex.lab, cex.cont=cex.cont)
    }else{
      stop("Pickands dependence function plotted only for 2 or 3 dimensional models.")
    }
    
  }else if (object=="returns"){
    
    returns.plot(model=model, par=par, Q.fix=Q.fix, Q.range=Q.range, 
                 Q.range0=Q.range0, cond=cond, labels=labels, cex.lab=cex.lab, ...)
    
  }
  
}

#########################################
#########################################
### Internal functions
#########################################
#########################################

AngDens2dPlot <- function(model='PB', para=c(2,1), log=TRUE, data=NULL, style=NULL,
                          labels="w", cex.lab=1, ...){
  
  if(!is.null(data) && is.null(style)){stop("style must be provided ('hist' or 'ticks')" )}
  if(!is.null(data) && style=="hist" && log==TRUE){stop("style 'hist' only valid when log=FALSE")}
  
  s <- seq(from=0.0001, to=0.9999, length=200)
  h <- sapply(s, function(y) dExtDep(x=c(y,1-y), method="Parametric", model=model, par=para, angular=TRUE, c=0, log=log))/2 # Division by dim=2 to standardise
  
  if(model == "AL" || model == "ET" || model == "EST"){ # Model has mass in the corners
    pm0 <- dExtDep(x=c(0,1), method="Parametric", model=model, par=para, angular=TRUE, c=0, log=log)/2 
    pm1 <- dExtDep(x=c(1,0), method="Parametric", model=model, par=para, angular=TRUE, c=0, log=log)/2
    
    range.h <- range(c(h,pm0,pm1))
  }else{
    range.h <- range(h)
  }
  
  if(log){
    y.lab <- expression(log(h(w)))
  }else{
    y.lab <- expression(h(w))
  }
  
  if(!is.null(data)){
    if(style == "hist"){
      hist(data[,1], freq=FALSE, main="",xlab=labels, ylab=y.lab,...)
      lines(s,h, col="blue", lwd=2)
    }else if(style == "ticks"){
      plot(s,h, type="l", col="blue", lwd=2, xlab=labels, ylab=y.lab, ylim = range.h )
      points(data[,1], rep(range.h[1], length(data[,1])), pch="|")
    }
  }else{
    plot(s,h, type="l", col="blue", lwd=2, xlab=labels, ylab=y.lab, ylim = range.h )
  }
  
  if(model == "AL" || model == "ET" || model == "EST"){ # Model has mass in the corners
    points(0, dExtDep(x=c(0,1), method="Parametric", model=model, par=para, angular=TRUE, c=0, log=log),  pch=16, col="blue")
    points(1, dExtDep(x=c(1,0), method="Parametric", model=model, par=para, angular=TRUE, c=0, log=log),  pch=16, col="blue")
  }

}  
  
AngDens3dPlot <- function(model='PB', para=c(2,4,15,1), log=TRUE, data=NULL, contour=TRUE, 
                          labels=c(expression(w[1]),expression(w[3]),expression(w[2])),cex.dat=1, cex.lab=1, cex.cont=1){
  
  if(model=="AL"){para <- c(rep(1,3), para[1], rep(0,6), para[2:4])}
  
  RectTri<-function(x,model,para,log){
    if(is.vector(x)){x<-matrix(x,nrow=1)}
    n<-nrow(x)
    d<-ncol(x)
    ind<-(rowSums(x>0)==d)*(rowSums(x)<1) # Give which coordinates have both components greater than zero and the sum of the components is less than 1
    ind [ind == 0] <- NA
    x<-cbind(x,1-rowSums(x))
    lf<-vector("numeric")
    
    for(i in 1:nrow(x)){
      if(prod(x[i,]>=0)==1){
        lf[i] = dExtDep(x=x[i,], method="Parametric", model=model, par=para, angular=TRUE, c=0, log=log, vectorial=TRUE)
      } else {lf[i]=0}
    }
    
    f<-lf*ind    
    return(f)
  }
  
  EquiTri <-function(x,model,para,log){
    y <- x
    y[,1] <- x[,1]- x[,2]/sqrt(3)
    y[,2] <- x[,2]*2/sqrt(3)
    fx <- RectTri(x=y, model=model, para=para, log=log)/(sqrt(3)/2) ## adjust with Jacobian
    equil.ind <- !((sqrt(3)*x[,1] - x[,2] <= 0) & x[,1]<=1/2 | (sqrt(3)*x[,1] + x[,2] >= sqrt(3)) & x[,1]>=1/2)  ## indicator for equilateral triangle
    equil.ind [equil.ind == 0] <- NA
    return(fx*equil.ind)
    
  }
  
  x1 <- seq(0,1, length=301)
  x2 <- seq(0,1, length=301)
  xx <- as.matrix(expand.grid(x1,x2))
  
  equi <- EquiTri(x=xx, model= model, para=para, log=log)
  dec <- seq (from=0.1, to=0.9, by=0.1)
  
  if(!is.null(data)){
    quant=0
    for(i in 1:nrow(data)){
      quant[i] = dExtDep(x=data[i,], method="Parametric", model=model, par=para, angular=TRUE, c=0, log=log, vectorial=TRUE)
    }
    deciles <- c(min(equi,na.rm=TRUE),quantile(quant,dec,na.rm=TRUE)/(sqrt(3)/2),max(equi,na.rm=TRUE))
  }else{
    deciles <- c(min(equi,na.rm=TRUE),quantile(equi,dec,na.rm=TRUE),max(equi,na.rm=TRUE))
  }	
  
  image(x1,x2, matrix(equi, nrow=length(x1), ncol=length(x2)), asp=1, breaks=deciles, col=rev(heat.colors(length(deciles)-1)),
        axes=FALSE,xlab="",ylab="",xlim=c(-0.05,1.06),ylim=c(-0.08,0.92))
  if(!is.null(data)){ points(data[,1] + 0.5*data[,2],sqrt(3)/2*data[,2],pch=16,cex=cex.dat)}
  text(c(0.94,0.5,0.06),c(-0.05,0.905,-0.05),labels=labels,cex=cex.lab)
  segments(c(0,0,0.5),c(0,0,sqrt(3)/2),c(1,0.5,1),c(0,sqrt(3)/2,0))
  if(contour==TRUE){ contour(x1,x2, matrix(equi, nrow=length(x1), ncol=length(x2)),levels=deciles, 
                             labels=round(deciles,3), labcex=cex.cont, nlevels=15,add=TRUE)}
}


Pick2dPlot <- function(model='HR', para=2, label="w", cex.lab=1, ...){
 
  w <- seq(0.00001, .99999, length=100)
  pick <- vector(length=100)
  for(i in 1:100){pick[i] <- index.ExtDep(object="pickands", model=model, par=para, x=c(w[i],1-w[i]))}
  
  plot(w, pick, type="l", ylim=c(0.5, 1), xlab=label, ylab="A(w)", cex.lab=cex.lab,...)
  polygon(c(0, 0.5, 1), c(1, 0.5, 1), lwd=2, border = 'grey')
  
}  

Pick3dPlot <- function(model='HR', para=c(2,2,2), contour=TRUE, 
                       labels=c(expression(w[1]),expression(w[3]),expression(w[2])),cex.lab=1, cex.cont=1){
  
  RectTri<-function(x,model,para){
    if(is.vector(x)){x<-matrix(x,nrow=1)}
    n<-nrow(x)
    d<-ncol(x)
    ind<-(rowSums(x>0)==d)*(rowSums(x)<1) # Give which coordinates have both components greater than zero and the sum of the components is less than 1
    ind [ind == 0] <- NA
    x<-cbind(x,1-rowSums(x))
    lf<-vector("numeric")
    
    for(i in 1:nrow(x)){
      if(prod(x[i,]>=0)==1){
        lf[i] = index.ExtDep(object="pickands", model=model, par=para, x=x[i,])
      } else {lf[i]=0}
    }
    
    f<-lf*ind    
    return(f)
  }
  
  EquiTri <-function(x,model,para){
    y <- x
    y[,1] <- x[,1]- x[,2]/sqrt(3)
    y[,2] <- x[,2]*2/sqrt(3)
    fx <- RectTri(x=y, model=model, para=para)
    equil.ind <- !((sqrt(3)*x[,1] - x[,2] <= 0) & x[,1]<=1/2 | (sqrt(3)*x[,1] + x[,2] >= sqrt(3)) & x[,1]>=1/2)  ## indicator for equilateral triangle
    equil.ind [equil.ind == 0] <- NA
    return(fx*equil.ind)
    
  }
  
  x1 <- seq(0.0001,0.9999, length=301)
  x2 <- seq(0.0001,0.9999, length=301)
  xx <- as.matrix(expand.grid(x1,x2))
  
  equi <- EquiTri(x=xx, model= model, para=para)
  dec <- seq (from=0.1, to=0.9, by=0.1)
  deciles <- c(min(equi,na.rm=TRUE),quantile(equi,dec,na.rm=TRUE),max(equi,na.rm=TRUE))
  
  image(x1,x2, matrix(equi, nrow=length(x1), ncol=length(x2)), asp=1, breaks=deciles, col=rev(heat.colors(length(deciles)-1)),
        axes=FALSE,xlab="",ylab="",xlim=c(-0.05,1.06),ylim=c(-0.08,0.99))
  text(c(0.94,0.5,0.06),c(-0.05,0.905,-0.05),labels=labels,cex=cex.lab)
  segments(c(0,0,0.5),c(0,0,sqrt(3)/2),c(1,0.5,1),c(0,sqrt(3)/2,0))
  if(contour==TRUE){ contour(x1,x2, matrix(equi, nrow=length(x1), ncol=length(x2)),levels=deciles, 
                             labels=round(deciles,3), labcex=cex.cont, nlevels=15,add=TRUE)}
}

returns.plot <- function(model, par, Q.fix, Q.range, Q.range0, cond, labels, cex.lab, ...){
  
  # Checking the model
  models <- c("HR", "ET", "EST")
  if(!any(model ==  models)){ stop("model wrongly specified")}
  
  # Looking into par
  if(is.vector(par)){
    dim <- dim_ExtDep(model=model, par=par)
  }else if(is.matrix(par)){
    dim <- dim_ExtDep(model=model, par=par[1,])
    n <- nrow(par)
  }else{
    stop("'par' must be a vector or a matrix")
  }
  if(all(dim != c(2,3)) ){ stop("Dimension 2 or 3 only")}
  
  # Checking Q.fix
  if(is.vector(Q.fix)){
    nQ.fix <- length(Q.fix)
    if(nQ.fix != dim){stop("'Q.fix' must be of the same length as the dimension of the model")}
    if(!any(is.na(Q.fix))){
      stop("'Q.fix' must contain some NAs. See help!")
    }else if(all(is.na(Q.fix)) && dim==2){
      nnfix <- 2
      cond <- FALSE # 2 dimensions and none are fixed so cannot be conditional.
    }else if(all(is.na(Q.fix)) && dim==3){
      stop("'Q.fix' must not contain only NAs for 3 dimensional models. See help!")
    }else{  
      rl.var <- which(is.na(Q.fix))
      rl.cond <- which(!is.na(Q.fix))
      nnfix <- length(rl.var) # number of non-fixed components
    }
  }else{
    stop("'Q.fix' must be a vector")
  }

  # Checking Q.range / Q.range0
  if(is.vector(Q.range)){
    
    if(!is.vector(Q.range0)){stop("'Q.range' and 'Q.range0' must be of the same format")}
    if(length(Q.range) != length(Q.range0)){stop("'Q.range' and 'Q.range0' must be of the same dimension")}
    
    Q.range <- matrix(Q.range, ncol=1)
    nQ.range <- nrow(Q.range)
    if(nnfix != 1){stop("There should be one NAs in 'Q.fix'")}
  }else if(is.matrix(Q.range)){
    if(!is.matrix(Q.range0)){stop("'Q.range' and 'Q.range0' must be of the same format")}
    if(any(dim(Q.range) != dim(Q.range0))){stop("'Q.range' and 'Q.range0' must be of the same dimension")}
    nQ.range0 <- nrow(Q.range0)
    if(ncol(Q.range) > dim){stop("'Q.range' cannot have more columns than the model dimension")}
    if(ncol(Q.range) == dim && dim==3){stop("'Q.range' cannot have as many columns as the model dimension when the dimension is 3.")}
    if(nnfix != ncol(Q.range)){stop("There should be as many NAs in 'Q.fix' as columns in 'Q.range'")}
    Q.range <- as.matrix(expand.grid(Q.range[,1], Q.range[,2])) # Since limited to 3-dim model, when Q.range is a matrix then it can only have 2 columnns
    nQ.range <- nrow(Q.range)
  }else{
    stop("'Q.range' must be a vector or a matrix")
  }
  
  Qmat <- matrix(rep(Q.fix, nQ.range), nrow=nQ.range, ncol=dim, byrow=TRUE)
  Qmat[,rl.var] <- Q.range
  
  # Joint probability
  p1 <- pExtDep(q=Qmat, type="upper", method="Parametric", model=model, par=par, plot=FALSE) # (n x nQ.range) matrix
  
  # conditional probability (if cond=TRUE)
  if(cond){
    if(length(rl.cond)==1){
      p <- p1 / (1-exp(-Q.fix[rl.cond]^(-1)))  
    }else if(length(rl.cond)==2){
      # Probability of conditioning event
      sub.par <- subset.pars(model=model, par=par, sub=rl.cond)
      p2 <- pExtDep(q=Q.fix[rl.cond], type="upper", method="Parametric", model=model, par=sub.par, plot=FALSE) # vector of length(n)
      if(is.vector(par)){
        p <- p1/p2
      }else if(is.matrix(par)){
        p <- apply(p1,2, function(x) x/p2) # (n x nQ.range) matrix      
      }
      
    }
    
  }else{
    p <- p1
  }
  
  if(is.vector(par)){
    p.inv <- 1/p
  }else if(is.matrix(par)){
    p.inv <- apply(1/p, 2, func, conf=0.95) # ( 3 x nQ.range) matrix    
  }

  if(nnfix==1){ # Univariate return levels
    
    if(is.vector(par)){
      plot(p.inv, Q.range0, xlab="Return period 1/p", type="l", ylab=labels, lwd=2, cex.lab=cex.lab,...)
    }else if(is.matrix(par)){
      plot(p.inv[2,], Q.range0, xlab="Return period 1/p", type="l", ylab=labels, lwd=2, cex.lab=cex.lab,...)
      points(p.inv[1,], Q.range0, type="l", lwd=2, lty=2)
      points(p.inv[3,], Q.range0, type="l", lwd=2, lty=2)
    }
    
  }else if(nnfix==2){ # Bivariate return levels
    
    levels <- round(quantile(p.inv, probs=c(1:9)/10))
    
    if(is.vector(par)){
      rl <- matrix(p.inv, nrow=nQ.range0, ncol=nQ.range0)
      contour(x=Q.range0[,1], y=Q.range0[,2], z=rl, levels=levels, xlab=labels[1], ylab=labels[2], cex.lab=cex.lab,...)
    }else if(is.matrix(par)){
      rl.low <- matrix(p.inv[1,], nrow=nQ.range0, ncol=nQ.range0)
      rl.mean <- matrix(p.inv[2,], nrow=nQ.range0, ncol=nQ.range0)
      rl.upp <- matrix(p.inv[3,], nrow=nQ.range0, ncol=nQ.range0)
      
      contour(x=Q.range0[,1], y=Q.range0[,2], z=rl.mean, levels=levels, xlab=labels[1], ylab=labels[2], lwd=2, cex.lab=cex.lab,...)
      contour(x=Q.range0[,1], y=Q.range0[,2], z=rl.low, levels=levels, drawlabels=FALSE, add=TRUE, lty=2)
      contour(x=Q.range0[,1], y=Q.range0[,2], z=rl.upp, levels=levels, drawlabels=FALSE, add=TRUE, lty=3)
    }
    
  }
  
  return(p.inv)
  
}

subset.pars <- function(model, par, sub){
  
  # Find subset of parameters from 3-dimensions to 2-dimensions
  
  # Checking the model
  models <- c("HR", "ET", "EST")
  if(!any(model ==  models)){ stop("model wrongly specified")}
  
  if(is.vector(par)){
    dim <- dim_ExtDep(model=model, par=par)    
  }else if(is.matrix(par)){
    dim <- dim_ExtDep(model=model, par=par[1,])
  }

  if(length(unique(sub)) != length(sub)){stop("'sub' can only take unique values")}
  if( all(!is.finite(sub))){stop("Check 'sub' vector")}
  if(length(sub) > dim){stop("'sub' cannot be of greater length than the dimension of the model")}
  if( any((sub-floor(sub))!=0) ){stop("'sub' must be a vector of integers")}
  
  index <- matrix(c(1,2,1,3,2,3), nrow=2, byrow=FALSE)
  pair <- apply(index, 2, function(x) all(x==sub))
  pair.keep <- which(pair==TRUE)
  
  if(model == "HR"){
    if(is.vector(par)){
      par.new <- par[pair.keep]      
    }else if(is.matrix(par)){
      par.new <- par[,pair.keep]
    }
  }else if(model == "ET"){
    if(is.vector(par)){
      par.new <- par[c(pair.keep,4)]
    }else if(is.matrix(par)){
      par.new <- par[,c(pair.keep,4)]
    }  
  }else if(model== "EST"){
    if(is.vector(par)){
      par.new <- par[c(pair.keep, index[,pair.keep]+3, 7)]
    }else if(is.matrix(par)){
      par.new <- par[,c(pair.keep, index[,pair.keep]+3, 7)]
    }  
  }
  
  if(is.matrix(par)){
    par.new <- matrix(par.new, nrow=nrow(par))
  }
  
  return(par.new)

}


