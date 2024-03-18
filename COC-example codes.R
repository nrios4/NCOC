###############################
### Part I: Main functions for NCOC 
###############################
KLD_pen <- function(xmat, ymat, B, lambda, theta, gam_vec,eps=10^-10){
## This function calculates the penalized KLD with augmented lagrangian, which corresponds to Eq. 6 of the main text
  
  n <- nrow(ymat)
  q <- ncol(ymat)
  y2 <- xmat%*%t(B)
  
  if(any(y2 < eps)){
    neginds <- which(y2 < eps)
    y2[neginds] <- eps  ## otherwise log(y2) contains NaNs
  }
  
  part1 <- -(1/n)*sum(ymat*log(y2)) 
  Btilde2 <- (B - 1/q)^2
  part2 <- lambda*sum(sqrt(colSums(Btilde2)))  
  part3 <- theta/2*sum((colSums(B)-1)^2)
  part4 <- sum(gam_vec*(colSums(B)-1))
  res <- part1 + part2 + part3 + part4
  
  return(res)
}



prox_gradient <- function(xmat, ymat, B0, gam_vec, lambda = 0.01,
                          theta = 0.5,step_size = 1,tol = 1e-5, max_iter = 1000,eps=10^-10){
## This function performs a proximal gradient update for B matrix, which corresponds to Eq. 9 of the main text
  
  q <- ncol(ymat)
  n <- nrow(ymat)
  p <- ncol(xmat)
  s <- 1
  delta <- 100000
  KLD0 <- KLD_pen(xmat, ymat, B0, lambda, theta, gam_vec)
  
  while(delta > tol & s < max_iter){
    
    xB <- xmat%*%t(B0)  # n by q matrix
    if(any(xB < eps)){
      neginds <- which(xB < eps)
      xB[neginds] <- eps  ## xB is used in denominator, otherwise it contains NaNs
    }
    # Bmin1sums <- colSums(B0-1) 
    Bmin1sums <- colSums(B0) - 1 # check this..
    # B2sums <- colSums((B0 - 1/q)^2)
    
    # copy previous value of B
    B1 <- B0
    # ts <- step_size/sqrt(s)
    ts <- step_size/log(s+1)
    
    for(j in 1:p){
      
      gvec <- numeric(q)
      
      for(k in 1:q){
        gvec[k] <- -(1/n)*sum(xmat[,j]*ymat[,k]/xB[,k]) + theta*Bmin1sums[j] + gam_vec[j]
      }
      
      # update column j of B1
      # print(gvec)
      norm1 <- sqrt(sum((B0[,j] - ts*gvec - 1/q)^2))
      # newcol <- max(0,1 - ts/(lambda*norm1))*(B0[,j] - ts*gvec - 1/q) + 1/q
      newcol <- pmax(0,1 - lambda*ts/norm1)*(B0[,j] - ts*gvec - 1/q) + 1/q
      B1[,j] <- pmax(0, newcol)
    }
    
    # Update the KLD to check for convergence
    KLD1 <- KLD_pen(xmat, ymat, B1, lambda, theta, gam_vec)
    delta <- abs(KLD1 - KLD0)
    # print(c(KLD1,KLD0,delta,s))
    
    # update for next iteration
    KLD0 <- KLD1
    B0 <- B1
    s <- s + 1
  }
  
  return(B1)  
}


######################### Regression analysis function
reg_coc <- function(xmat, ymat, B0, lambda = 0.01,
                              theta = 1, step_size = 0.01,tol = 1e-5, max_iter = 1000){
  
  q <- ncol(ymat)
  n <- nrow(ymat)
  p <- ncol(xmat)
  if(missing(B0)){
    B0 <- matrix(1/q, nrow = q, ncol = p) # initial B0 matrix
  }
  
  s <- 1
  gam_vec0 <- numeric(p)  # gammas
  delta <- 100000
  KLD0 <- KLD_pen(xmat, ymat,B0, lambda, theta, gam_vec0)
  
  
  while(delta > tol & s < max_iter){
    
    # ts <- step_size/sqrt(s)
    # print(paste("Iteration", s, sep = " "))
    # Update B
    B1 <- prox_gradient(xmat, ymat, B0, gam_vec0, lambda, theta, step_size)
    #print("updated B matrix")
    #print(round(B1,3))
    
    # Update gamma
    gam_vec1 <- gam_vec0 + theta*(colSums(B1) - 1)
    #print("updated gamma vector")
    
    # Update the KLD to check for convergence
    KLD1 <- KLD_pen(xmat, ymat, B1, lambda, theta, gam_vec1)
    delta <- abs(KLD1 - KLD0)
    #print(c(KLD1,KLD0,delta,s))
    
    
    # update for next iteration
    KLD0 <- KLD1
    B0 <- B1
    gam_vec0 <- gam_vec1
    s <- s + 1
  }
  
  ## due to machine epsilon, sometimes we need normalizaitons£»
  ## yet this inexact solution is well accomendated in proposition 1 
  for(j in 1:ncol(B0)){ B0[,j]=B0[,j]/colSums(B0)[j]}

  return(B0)
  
}

############## CV function to select tuning parameter lambda
cv_coc = function(xmat,ymat,B0,n_fold = 5, lambda_vec){
  n = nrow(xmat)
  q = ncol(ymat)
  p = ncol(xmat)
  if(missing(B0)){
    B0 = matrix(1/q, nrow = q, ncol = p) # initial B0 matrix
  }

  KLD_vec = numeric(length(lambda_vec))
  # randomly shuffle data
  rand_seq = sample(1:n) 
  # create indices for n_fold folds
  folds = cut(rand_seq, breaks = n_fold, labels = FALSE) 
 
for(i in 1:length(lambda_vec)){
    lambda = lambda_vec[i]
    KLD_CV = 0
    for(j in 1:n_fold){
      test_inds = which(folds==j)
      X_test = xmat[test_inds,]
      Y_test = ymat[test_inds,]
      X_train = xmat[-test_inds,]
      Y_train = ymat[-test_inds,]
      B_CV <- reg_coc(X_train, Y_train, B0, step_size = 0.01, theta =1, lambda = lambda)

      ypred1=X_test%*%t(B_CV)
      div2=rep(NA,nrow(X_test))
      for(a in 1:nrow(X_test)){
        P=Y_test[a,]
        Q=ypred1[a,]
        xpq=rbind(P,Q)
        div2[a]=KL(xpq,unit="log2",est.prob = "empirical")
      }
      KLD_CV = KLD_CV + sum(div2) 
    }
    KLD_vec[i] = KLD_CV   
  }
  
  plot(lambda_vec, KLD_vec, type = "l")
  return(lambda_vec[which.min(KLD_vec)])
}


##########################################
### Part II: Example codes to implement NCOC 
##########################################
## These example codes may take 3 minutes
set.seed(2024)
n=50  
p=50
q=10
lambda_vec = c(10^-4,10^-3,10^-2,0.1,1)
library(codalm)
library(dirmult)
library(philentropy)

X=matrix(NA,n,p)
for(i in 1:n){X[i,]=rdirichlet(1,alpha=rep(1/p,p))}
Bnull=B=matrix(1/q,nrow=q,ncol=p)
pstar=sample(c(1:p),5,replace=F)
for(j in 1:p) {
if(j %in% pstar){B[,j]=rdirichlet(1,alpha=rep(1,q))}
}
Ey=X%*%t(B)
y=Ey
for(i in 1:nrow(y)){y[i,]=rdirichlet(1,alpha=(10*Ey[i,]))}

oldB <- t(codalm(y, X))
lambda0 = cv_coc(X,y,Bnull,n_fold = 5,lambda_vec)
lambda0
newB <-  reg_coc(X, y, step_size =  0.01, Bnull , theta =1, lambda =lambda0) 
FofDR=norm((oldB-B),type="F");FofDR
FofCOC=norm((newB-B),type="F");FofCOC


