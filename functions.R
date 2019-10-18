library(MASS)
library(MCMCpack)
library(mvtnorm)
library(Matrix)
library(Matching)

logit<-function(x){log(x/(1-x))}
expit<-function(x){exp(x)/(1+exp(x))}
#Draw samples from a truncated normal:
rtnorm<-function(n,mu,sigma,lower,upper){ 
  lp<-pnorm(lower,mu,sigma) 
  up<-pnorm(upper,mu,sigma)  
  qnorm(runif(n,lp,up),mu,sigma) 
}


MI_MNAR <- function(y, X, A, sd.beta, sd.alpha, sd.gamma,a, b, mu0.X, cov0.X, eps=1,iters=1000){
  
  n <- length(y)
  low<-ifelse(A==1,0,-Inf)
  high<-ifelse(A==1,Inf,0)
  p0 <- ncol(X)
  p.out <- p0 + 1
  p.A <- p0 + 1
  p.miss <- p0 + 2
  
  # Initial values
  miss <- is.na(X)
  mis.ind <- which(apply(miss,2,sum)!=0)
  for (ind in mis.ind){
    X[miss[,ind],ind] <- mean(X[,ind],na.rm=TRUE)
  }
  nmis.X <- length(mis.ind)
  R <- as.matrix(1- miss[,mis.ind],n,nmis.X)
  lowR<-ifelse(R==1,0,-Inf)
  highR<-ifelse(R==1,Inf,0)
  
  miss.y <- is.na(y)
  Ry <- 1-miss.y
  lowRy<-ifelse(Ry==1,0,-Inf)
  highRy<-ifelse(Ry==1,Inf,0)
  nm.y       <- sum(miss.y)
  y[miss.y] <- mean(y,na.rm=TRUE)
  
  
  beta0 <- rep(0,p.out)
  beta1 <- rep(0,p.out)
  alpha <- rep(0,p.A)
  gamma <- matrix(0,p.miss,nmis.X)
  gammay <- rep(0,p.miss)
  tau0 <- 1
  tau1 <- 1
  mu.X <- rep(0,p0)
  Q0 <- diag(p0)
  
  # Stores samples
  keep.beta0 <- matrix(0,iters,p.out)
  keep.beta1 <- matrix(0,iters,p.out)
  keep.alpha <- matrix(0,iters,p.A)
  keep.gamma <- array(0,c(iters,p.miss,nmis.X))
  keep.gammay <- matrix(0,iters,p.miss)
  keep.sigma20 <- rep(0,iters)
  keep.sigma21 <- rep(0,iters)
  keep.X <- array(0,c(iters,n,p0))
  keep.y <- array(0,c(iters,n))
  keep.muX <- matrix(0,iters,p0)
  keep.covX <- array(0,c(iters,p0,p0))
  keep.Q <- array(0,c(iters,p0,p0))
  
  for(iter in 1:iters){
    if(nm.y>0){
      y[miss.y] <- as.vector(rnorm(nm.y,cbind(1,X[miss.y,])%*%beta0,sqrt(1/tau0)))*
        (1-A[miss.y]) + 
        as.vector(rnorm(nm.y,cbind(1,X[miss.y,])%*%beta1,sqrt(1/tau1)))*A[miss.y]
    }
    
    X.out0 <- cbind(1,X)[A==0,]
    X.out1 <- cbind(1,X)[A==1,]
    y0 <- y[A==0]
    y1 <- y[A==1]
    X.A <- cbind(1,X)
    X.miss <- cbind(1,X,A)
    
    # Update beta0
    MMM <- tau0*t(X.out0)%*%y0
    VVV <- solve(tau0*t(X.out0)%*%X.out0+diag(p.out)/sd.beta^2)
    beta0 <- VVV%*%MMM + t(chol(VVV))%*%rnorm(p.out)
    
    # Update beta1
    MMM <- tau1*t(X.out1)%*%y1
    VVV <- solve(tau1*t(X.out1)%*%X.out1+diag(p.out)/sd.beta^2)
    beta1 <- VVV%*%MMM + t(chol(VVV))%*%rnorm(p.out)
    
    # Update tau0=1/sigma20
    tau0 <- rgamma(1,length(y0)/2+a,sum((y0-X.out0%*%beta0)^2)/2+b)
    # Update tau1=1/sigma21
    tau1 <- rgamma(1,length(y1)/2+a,sum((y1-X.out1%*%beta1)^2)/2+b)
    
    # Update latent probit variable z
    z <- rtnorm(n,X.A%*%alpha,1,low,high)
    
    # Update alpha
    MMM <- t(X.A)%*%z
    VVV <- solve(t(X.A)%*%X.A+diag(p.A)/sd.alpha^2)
    alpha <- VVV%*%MMM + t(chol(VVV))%*%rnorm(p.A)
    
    # Update latent probit variable zR
    zR <- matrix(0,n,nmis.X)
    for (i in 1:nmis.X){
      zR[,i] <- rtnorm(n,X.miss%*%gamma[,i],1,lowR[,i],highR[,i])
      # Update gamma
      MMM <- t(X.miss)%*%zR[,i]
      VVV <- solve(t(X.miss)%*%X.miss+diag(p.miss)/sd.gamma^2)
      gamma[,i] <- VVV%*%MMM + t(chol(VVV))%*%rnorm(p.miss)
    }
    
    
    # Update latent probit variable zRy
    zRy <- rtnorm(n,X.miss%*%gammay,1,lowRy,highRy)
    # Update gammay
    MMM <- t(X.miss)%*%zRy
    VVV <- solve(t(X.miss)%*%X.miss+diag(p.miss)/sd.gamma^2)
    gammay <- VVV%*%MMM + t(chol(VVV))%*%rnorm(p.miss)
    
    # Update mean and covariance of X ~mvtnorm(mu.X,solve(Q0))
    MMM <- Q0%*%colSums(X) + solve(cov0.X)%*%mu0.X
    VVV <- solve(n*Q0 + solve(cov0.X))
    mu.X <- VVV%*%MMM + t(chol(VVV))%*%rnorm(p0)
    
    SSS <- sweep(X,2,mu.X,"-")
    SSS <- t(SSS)%*%SSS
    Q0  <- rwish(n+p0+eps,solve(SSS+diag(p0)*(p0+eps))) 
    cov.X <- solve(Q0)
    
    #Imput X
    for (ind in mis.ind){
      X.obs <- X[,-ind]
      X.mis <- X[,ind]
      beta0.obs <- beta0[-(ind+1)]
      beta1.obs <- beta1[-(ind+1)]
      beta0.mis <- beta0[ind+1]
      beta1.mis <- beta1[ind+1]
      alpha.obs <- alpha[-(ind+1)]
      alpha.mis <- alpha[ind+1]
      gamma.obs <- gamma[-(ind+1),]
      gamma.mis <- gamma[ind+1,]
      gammay.obs <- gammay[-(ind+1)]
      gammay.mis <- gammay[ind+1]
      mean.cond <- mu.X[ind] + as.vector(cov.X[ind,-ind]%*%solve(cov.X[-ind,-ind])%*%(X.obs-mu.X[-ind]))
      var.cond <- as.vector(cov.X[ind,ind]-cov.X[ind,-ind]%*%solve(cov.X[-ind,-ind])%*%cov.X[-ind,ind])
      
      MMM <- as.vector(alpha.mis*(z-cbind(1,X.obs)%*%alpha.obs) + 
                         beta0.mis*tau0*(y-cbind(1,X.obs)%*%beta0.obs)*(1-A) +
                         beta1.mis*tau1*(y-cbind(1,X.obs)%*%beta1.obs)*A +
                         mean.cond/var.cond + 
                         t(gamma.mis%*%t(zR-cbind(1,X.obs,A)%*%gamma.obs)) + 
                         gammay.mis*(zRy-cbind(1,X.obs,A)%*%gammay.obs))
      VVV <- 1/(alpha.mis^2 + sum(gamma.mis^2) + gammay.mis^2 +
                  beta0.mis^2*tau0*(1-A)+beta1.mis^2*tau1*A + as.vector(1/var.cond))
      impute <- VVV*MMM + sqrt(VVV)*rnorm(n)
      #impute <- rnorm(n,mean.cond,sqrt(var.cond))
      X[miss[,ind],ind] <- impute[miss[,ind]]
    }
    keep.beta0[iter,] <- beta0
    keep.beta1[iter,] <- beta1
    keep.alpha[iter,] <- alpha
    keep.gamma[iter,,] <- gamma
    keep.gammay[iter,] <- gammay
    keep.sigma20[iter] <- 1/tau0
    keep.sigma21[iter] <- 1/tau1
    keep.X[iter,,] <- X
    keep.y[iter,] <- y
    keep.muX[iter,] <- mu.X
    keep.covX[iter,,]<-cov.X
    keep.Q[iter,,]<-Q0
  }
  out <- list(beta0=keep.beta0,beta1=keep.beta1,alpha=keep.alpha,
              gamma=keep.gamma,gammay=keep.gammay,
              sigma20=keep.sigma20,sigma21=keep.sigma21,
              X=keep.X,y=keep.y,muX=keep.muX,covX=keep.covX,Q=keep.Q)
  return(out)
}
