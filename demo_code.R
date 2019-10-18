source("functions.R")
source("mibs.R")

nimpute <- 10 # Multiple imputation times

nbs <- 1000  # Number of Bootstrap samples

# Simulate one dataset
n <- 2000
mu.X <- c(0,0)
cov.X <- matrix(c(1,0.3,0.3,1),2,2)
px <- 2
x <- mvrnorm(n,mu.X,cov.X)
x1 <- x[,1]
x2 <- x[,2]
alpha <- c(-0.2,0.4,0.3)
A <- as.vector(ifelse(cbind(1,x)%*%alpha+rnorm(n)>0,1,0))
beta0 <- c(1,3,1.5)
beta1 <- c(0,2,1)
y0 <- as.vector(cbind(1,x)%*%beta0 + rnorm(n))
y1 <- as.vector(cbind(1,x)%*%beta1 + rnorm(n))
y <- y0*(1-A) + y1*A

gamma <- c(0.8,0,1,0,0)
gammay <- c(1,0.5,0.5,0,0.2)
obs.loc.x2 <- as.vector(ifelse(cbind(1,x,y,A)%*%gamma+rnorm(n)>0,1,0))
obs.loc.y <- as.vector(ifelse(cbind(1,x,y,A)%*%gammay+rnorm(n)>0,1,0))
x2[obs.loc.x2==0] <- NA
y[obs.loc.y==0] <- NA
X <- cbind(x1,x2)

#Point estimate, variance, confidence interval of MI estimators
result.reg <- mibs(y,X,A,method='regression',nimpute=10,nbs=1000)
result.IPW <- mibs(y,X,A,method='IPW',nimpute=10,nbs=1000)
result.AIPW <- mibs(y,X,A,method='AIPW',nimpute=10,nbs=1000)
result.matching <- mibs(y,X,A,method='matching',nimpute=10,nbs=1000)

