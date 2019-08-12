###
### Generate the random variable from clover distribution
###

simulationsClover = function(n){

  R1 = numeric(n);
  Theta = numeric(n);
  ind = numeric(n);
  
  for(j in 1:n){
    U1=runif(1,0,1);
    U2=runif(1,0,1);
    
    R1[j] = ((1-U1)^(-2)-1)^(1/2);
    f = function(t){2*t/pi+3*sin(4*t)/(10*pi)-U2}
    Theta[j] = uniroot(f,lower=0,upper=pi/2)$root;
    ind[j] = 1;
  }
  
  Z = numeric(2*n);
  dim(Z) = c(n,2);
  Z[,1] = R1*cos(Theta);
  Z[,2] = R1*sin(Theta); #    Z is simulated from clover distribution
  Z }

###
### Generate the random variable from asymmetric distribution
###

# auxiliary functions 

fr = function(r){0.30746/(r^(5/12)*(r+1))}
fth = function(th){0.15739/(th^(3/4)*((th+1)^(7/12)))}

Mr = function(r){if(r <= 1){ret = r^(-5/12)}
  if(r > 1){ret = r^(-17/12)}
  0.30746*ret}

Mth = function(th){if(th <= 1){ret = th^(-3/4)}
  if(th > 1){ret = th^(-4/3)}
  0.15739*ret}

fr = Vectorize(fr, vectorize.args = "r")
fth = Vectorize(fth, vectorize.args = "th")

Mr = Vectorize(Mr, vectorize.args = "r")
Mth = Vectorize(Mth, vectorize.args = "th")


simData=function(n){
  
  # Step 1: simulate data from m (its df)
  
  Ur = runif(2*n, 0, 1)
  Uth = runif(2*n, 0, 1)
  
  small_r = which(Ur <= 5/12)
  big_r = which(Ur > 5/12)
  
  small_th = which(Uth <= 4/7)
  big_th = which(Uth > 4/7)
  
  Tr = c(1:(2*n))*0
  Tth = c(1:(2*n))*0
  
  Tr[small_r] = ((12/5)*Ur[small_r])^(12/7)
  Tr[big_r] = ((12/7)*(1-Ur[big_r]))^(-12/5)
  
  Tth[small_th] = ((7/4)*Uth[small_th])^4
  Tth[big_th] = ((7/3)*(1-Uth[big_th]))^(-3)
  
  # Step 2: accept if M(T)xU <= f(T)
  
  Umr = runif(2*n, 0, 1)
  Umth = runif(2*n, 0, 1)
  
  accept = which(Mr(Tr)*Umr <= fr(Tr) & Mth(Tth)*Umth <= fth(Tth))
  
  L = length(accept)
  
  if(L>=n){
    R = c(1:n)*0
    TH = c(1:n)*0
    
    R = Tr[accept[1:n]]
    TH = Tth[accept[1:n]]
    
    # i=1
    # while(i<=n)
    # {prod[i] = R[i]*TH[i]; i+1}
    
    X1 = (R/(TH+1))^(1/3)
    X2 = ((R*TH)/(TH+1))^(1/4)
    ret = list()
    ret$rez = cbind(X1, X2)
    ret$L = L}
  
  if(L<n){print('try again :-)')}
  ret}


##  Simulating the data 

# in R Console type:
# X = simData(n) 
#		     X gives two outputs: 
#		     1. n data points (they are matched
#		        here in the way that the point is rejected
#                   if either of 'r' or 'theta' is a bloop)
#		     2. the number of the accepted data

# X = simData(n)$rez 
#		     this gives only the n data points



