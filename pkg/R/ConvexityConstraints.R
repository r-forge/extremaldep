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
### Last change: 05/10/2014                           ###
#########################################################


# Main routine that provides a matrix of CONSTRAINTS that garanties
# the CONVEXITY. 
convexity <- function(v, d){
  
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


