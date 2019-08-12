#########################################################
### Authors: Giulia Marcon and Simone Padoan          ###
### Emails: giulia.marcon@phd.unibocconi.it,          ###
### simone.padoan@unibocconi.it                       ###
### Institution: Department of Decision Sciences,     ###
### University Bocconi of Milan                       ###
### File name: ConvexityConstraints.r                 ###
### Description:                                      ###
### This file contains a set of procedures            ###
### for computing the matrix of convexity constraints ###
### used in the quadprog procedure of optimization    ###
### Last change: 15/08/2016                           ###
#########################################################

## Permutation function taken from gtools package (hidden function)

# Main routine that provides a matrix of CONSTRAINTS that garanties
# the CONVEXITY. 

convexity <-function(d, k){
  
  if(d==2){
    R1 <- matrix(0,k-1,k+1)
    for (i in 1:(k-1))
      R1[i,i:(i+2)] <- c(1,-2,1)
  }
  
  else{
    ###########################
    #Subroutines for convexity:
    
    # 1)
    index <- function(k, d){
      beta.index <- expand.grid(rep(list(0:k),d))
      beta.index.sum <- rowSums(beta.index)
      restr <- which(beta.index.sum<=k)
      v <- as.matrix(beta.index[restr,ncol(beta.index):1])
      return(v)
    }
    
    # 2)
    lookUpTable <- function(kTilde, dTilde){
      nrows= (kTilde+1)^(dTilde)
      lookup <- matrix(nrow = nrows,ncol=2)
      feasibles =0;
      for(i in 0:nrows-1){
        if(isFeasible(i,kTilde+1)){
          lookup[i+1,1] = i
          lookup[i+1,2] = feasibles
          feasibles = feasibles+1;
        }
      }
      return(lookup)
    }
    
    # 3)
    isFeasible <- function (number, base){
      dividend = number
      sum =0
      while(dividend>=base){
        sum = sum + dividend%%base
        dividend = floor(dividend/base)
      }
      sum =sum+dividend
      return(sum<=(base-1)) 
    }
    
    # 4)
    frombaseKToBase10 <-function(row, kTilde){
      dTilde <- length(row)
      baseNum <- kTilde+1
      index <- 0;
      for(i in dTilde:1){
        index = index + row[i]*baseNum^(dTilde-i)
      }
      return(index)
    }
    
    # 5)
    numberToModifyInCompleteTable <- function(numberBase10, power_direction, kTilde, order){
      return(numberBase10 + ((kTilde+1)^(power_direction-1))*order)
    }
    
    # 6)
    diff1inDirections <-function (lookUp,numberBase10, k, vLength, power_direction1, power_direction2) {
      resultVector = vector(length = vLength,mode = "numeric")
      
      numberAfterFirstDirection = numberToModifyInCompleteTable(numberBase10, power_direction1, k,1)
      indexToBeModified11Complete = numberToModifyInCompleteTable(numberAfterFirstDirection, power_direction2, k,1) + 1
      index11 = lookUp[indexToBeModified11Complete,2]
      
      indexToBeModified01Complete = numberToModifyInCompleteTable(numberBase10, power_direction2, k,1) + 1
      index01 = lookUp[indexToBeModified01Complete,2]
      
      indexToBeModified10Complete = numberToModifyInCompleteTable(numberBase10, power_direction1, k,1) + 1
      index10 = lookUp[indexToBeModified10Complete,2]
      
      index0 = lookUp[numberBase10 +1,2]
      
      resultVector[index11] = 1;
      resultVector[index10] =  resultVector[index10]-1;
      resultVector[index01] =  resultVector[index01] -1;
      resultVector[index0] = 1;
      return(resultVector)
    }
    ##########################################
    
    dTilde = d-1
    kTilde = k-2
    v = index(k, dTilde)
    vSub = index(kTilde, dTilde)
    vLengthSub =  dim(vSub)[1]
    vLength = dim(v)[1]
    lookUp = lookUpTable(k, dTilde)
    
    directions = 1: (dTilde)
    perm <- permutations(2,dTilde-1,c(1,-1),repeats.allowed =TRUE);
    dimPerm <-  dim(perm)[1]
    
    R1 = matrix(data= 0, nrow = dTilde*dimPerm*vLengthSub, ncol=vLength)
    for(vIndex in 1:vLengthSub){
      numberBase10 = frombaseKToBase10(vSub[vIndex,],k)
      for(i in directions){
        diff1is = matrix(nrow = dTilde, ncol=vLength)
        diff2s = matrix(nrow = dTilde, ncol=vLength)
        diff2s = diff1inDirections(lookUp,numberBase10, k, vLength, i, i)
        spuriousDirections = directions[-i]
        for(j in spuriousDirections){
          diff1is[j,] =diff1inDirections(lookUp,numberBase10, k, vLength, i, j)
        }
        for(t in 1:dimPerm){
          begin =(vIndex-1)*dimPerm*dTilde
          rowR1 = diff2s +  perm[t,] %*%diff1is[-i,]
          R1[begin +(dTilde-i)*dimPerm+t,] = rowR1
        }
      }
    }
  }
  return(R1)
}



# Main routine that provides a matrix of CONSTRAINTS that garanties
# the CONVEXITY. OLD VERSION
convexity_old <- function(v, d){
  
  combn_eqk <- v
  q <- nrow(v)
  k <- max(v)
  
  if(d==2){
    Aconvx <- matrix(0,k-1,k+1)
    for (i in 1:(k-1))
      Aconvx[i,i:(i+2)] <- c(1,-2,1)
  }
  
  if(d==3){
    B <- matrix(c(0,  1,  rep(0,k-1),  -1,  -1,  rep(0,k-1),  1,
                  2,  -1,	rep(0,k-1),	-3,	 1,	rep(0,k-1),	1,
                  0,	-1,	1,	rep(0,k-2),	 1,	-1,	rep(0,k),
                  2,	-3,	1,	rep(0,k-2),	-1,	 1,	rep(0,k)),
                4,ncol=3+2*k,byrow=T)
    ri <- 0
    co <- 0
    for (i in 1:(k-1)){
      ri <- ri + k-i
      co <- c(co, (i+(i-1)*k):(i*k-1))
    }
    righe <- seq(1,ri*4,by=4)
    colonne <- co[-1]
    
    Aconvx <- matrix(0,length(righe)*nrow(B),(k+1)^2)
    for (i in 1:length(righe)){
      Aconvx[righe[i]:(righe[i]+nrow(B)-1),colonne[i]:(colonne[i]+ncol(B)-1)] <- B
    }
    restr <- labels(v)[[1]]
    restr <- as.numeric(restr)
    
    Aconvx <- Aconvx[,restr]
  }
  
  if(d>3){  
    
    sumvb <- apply(v,1,sum)
    leqk_2 <- which(sumvb<=k-2)
    combn_leqk_2 <- as.matrix(v[leqk_2,])
    
    ##########################
    # Subroutines of Convexity:
    
    # 1) Difference operators for convexity constraints
    DELTAl0 <- function(d,combn_leqk_2,diffl,diff0,l){
      diffll <- diffl0 <- diffl
      diff0l <- diff00 <- diff0
      
      for(m in 1:(d-1)){
        for (i in 1:nrow(combn_leqk_2)){
          diffll[i,l+1,m] <- diffl[i,l+1,m] + 1
          diffl0[i,1,m] <- diffl[i,1,m] + 1 
          diff0l[i,l+1,m] <- diff0[i,l+1,m] + 1
          diff00[i,1,m] <- diff0[i,1,m] + 1
        }
      }
      return(list(diffll=diffll,diffl0=diffl0,
                  diff0l=diff0l,diff00=diff00))
    }
    
    # 2) 
    DELTAl0d <- function(d,combn_leqk_2,diffl,diff0,l){
      diffll <- diffl0 <- diffl
      diff0l <- diff00 <- diff0
      
      for(m in 1:(d-1)){
        for (i in 1:nrow(combn_leqk_2)){
          diffll[i,l,m] <- diffl[i,l,m] + 1
          diffl0[i,1,m] <- diffl[i,1,m]  
          diff0l[i,l,m] <- diff0[i,l,m] + 1
          diff00[i,1,m] <- diff0[i,1,m]
        }
      }
      return(list(diffll=diffll,diffl0=diffl0,
                  diff0l=diff0l,diff00=diff00))
    }
    
    # 3)
    DELTAl0DELTAl0 <- function(d,combn_leqk_2,diffl,diff0,l){
      deltal0 <- NULL
      for (l in 1:(d-1)){
        deltal0[[l]] <- DELTAl0(d=d,combn_leqk_2,diffl,diff0,l=l)
      }
      return(deltal0)
    }
    
    # 4)
    DELTAl0DELTAl0d <- function(d,combn_leqk_2,diffl,diff0,l){
      deltal0 <- NULL
      for (l in 1:(d-1)){
        deltal0[[l]] <- DELTAl0d(d=d,combn_leqk_2,diffl,diff0,l=l)
      }
      return(deltal0)
    }
    
    
    # DELTA_i0
    diffl <- diff0 <- array(rep(combn_leqk_2,length(leqk_2)),dim=c(length(leqk_2),d-1,d-1))
    
    for (l in 1:(d-1)){
      for (i in 1:length(leqk_2)){
        diffl[i,l,l] <- combn_leqk_2[i,l] + 1
        diff0[i,1,l] <- combn_leqk_2[i,1]   
      }
    }
    
    # DELTAl0 DELTAl0
    DIFF <- DELTAl0DELTAl0d(d=d,combn_leqk_2,diffl,diff0,l=l)
    
    # + and - combinations
    perm <- permutations(2,(d-2),c('+','-'),repeats.allowed=TRUE)
    piuomeno <- NULL
    piuomeno[[1]] <- cbind(rep("+",nrow(perm)),perm)
    piuomeno[[d-1]] <- cbind(perm,rep("+",nrow(perm)))
    
    for(l in 1:(d-3)){
      piuomeno[[l+1]] <- cbind(perm[,1:l],rep("+",nrow(perm)),perm[,(l+1):(d-2)])
    }
    
    rowA <- nrow(perm)*nrow(combn_leqk_2)*(d-1)
    Aconvx <- matrix(0,rowA,q)
    
    ind <-0
    f <- 0
    o <- 0
    g <- 0
    
    for (l in 1:(d-1)){
      g <- (l-1) * nrow(perm)  
      
      for(m in 1:(d-1)){
        for(i in 1:nrow(combn_leqk_2)){
          f <- (i-1+o) * nrow(perm) 
          
          for(h in 1:(nrow(perm))){
            for(j in 1:nrow(combn_eqk)){
              ind <- h+f+g
              
              if(l==m){
                
                if(all(DIFF[[l]]$diffll[i,,m]==combn_eqk[j,]))
                  Aconvx[ind,j] <- Aconvx[ind,j] + 1
                if(all(DIFF[[l]]$diffl0[i,,m]==combn_eqk[j,]))
                  Aconvx[ind,j] <- Aconvx[ind,j] - 1
                if(all(DIFF[[l]]$diff0l[i,,m]==combn_eqk[j,]))
                  Aconvx[ind,j] <- Aconvx[ind,j] - 1
                if(all(DIFF[[l]]$diff00[i,,m]==combn_eqk[j,]))
                  Aconvx[ind,j] <- Aconvx[ind,j] + 1
                
              }
              else{
                
                pm <- switch(piuomeno[[m]][h,l],"+"= 1, "-"= -1)
                if(all(DIFF[[l]]$diffll[i,,m]==combn_eqk[j,]))
                  Aconvx[ind,j] <- Aconvx[ind,j] + 1 * pm
                if(all(DIFF[[l]]$diffl0[i,,m]==combn_eqk[j,]))
                  Aconvx[ind,j] <- Aconvx[ind,j] - 1 * pm
                if(all(DIFF[[l]]$diff0l[i,,m]==combn_eqk[j,]))
                  Aconvx[ind,j] <- Aconvx[ind,j] - 1 * pm
                if(all(DIFF[[l]]$diff00[i,,m]==combn_eqk[j,]))
                  Aconvx[ind,j] <- Aconvx[ind,j] + 1 * pm
              }
            }
          }
        }
      }
      o <- o + 2
    }
  }
  return(Aconvx)
}

###################################################################

########################### START CONSTRAINTS #####################
########################### ONLY BIVARIATE CASE ###################

constraints <- function(k){
  kp <- k+1
  km <- k-1
  # Create the A, b0 matrix for constraints
  # R1) Convexity
  R1 <- matrix(0,km,kp)
  for (i in 1:km)
    R1[i,i:(i+2)] <- c(1,-2,1)
  r1 <- rep(0,nrow(R1))
  
  # R2) upper bound
  # A(ei) = 1
  R2 <- matrix(0,4,kp)
  R2[1,1] <- R2[2,kp] <- 1
  R2[3,1]<- R2[4,kp] <- -1
  r2 <- c(1,1,-1,-1)
  
  # R3) lower bound
  # A(w) >= max(w)
  R3 <- matrix(0,2,kp)
  R3[1,2] <- R3[2,k] <- 1
  r3 <- rep(1-1/k,2)
  
  return(list(R=rbind(R1, R2, R3), r = c(r1, r2, r3) ))
  
}
########################### ONLY BIVARIATE CASE ###################
########################## END CONSTRAINTS ########################
