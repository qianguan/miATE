mibs <- function(y,X,A,method,nimpute=10,nbs=1000,nsamps=100,M=1,iters=5000,burn=2000){
  
  n <- nrow(X)
  px <- ncol(X)
  
  miss <- is.na(X)
  mis.ind <- which(apply(miss,2,sum)!=0)
  nmis.X <- length(mis.ind)
  obs.loc <- as.matrix(1-miss[,mis.ind],n,nmis.X)
  fit <- MI_MNAR(y,X,A,100,100,100,0.01,0.01,rep(0,px),diag(px),iters=iters)
  imputes <- sample((burn+1):iters,nimpute,replace=FALSE)
  samps <- sample((burn+1):iters,nsamps,replace=FALSE)
  alpha.mle <- apply(fit$alpha[samps,],2,mean)
  gamma.mle <- apply(array(fit$gamma[samps,,],c(nsamps,dim(fit$gamma)[2:3])),c(2,3),mean)
  gammay.mle <- apply(fit$gammay[samps,],2,mean)
  beta0.mle <- apply(fit$beta0[samps,],2,mean)
  beta1.mle <- apply(fit$beta1[samps,],2,mean)
  sigma20.mle <- mean(fit$sigma20[samps])
  sigma21.mle <- mean(fit$sigma21[samps])
  muX.mle <- apply(fit$muX[samps,],2,mean)
  covX.mle <- apply(fit$covX[samps,,],c(2,3),mean)
  Q.mle <- apply(fit$Q[samps,,],c(2,3),mean)
  
  trt_effect <- rep(0,nimpute)
  var_effect <- rep(0,nimpute)
  iid <- matrix(0,n,nimpute)
  
  if(method=="regression"){
    mu.deriv <- rep(0,px+1)
    score.deriv0 <- matrix(0,px+1,px+1)
    score.deriv1 <- matrix(0,px+1,px+1)
    for(i in 1:nsamps){
      X.impute <- fit$X[samps[i],,]
      mu.deriv <- mu.deriv + 1/n*rep(1,n)%*%cbind(1,X.impute)/nsamps
      score.deriv0 <- score.deriv0 + 
        1/n*t(cbind(1,X.impute[A==0,]))%*%cbind(1,X.impute[A==0,])/nsamps
      score.deriv1 <- score.deriv1 + 
        1/n*t(cbind(1,X.impute[A==1,]))%*%cbind(1,X.impute[A==1,])/nsamps
    }
    for(i in 1:nimpute){
      X.impute <- fit$X[imputes[i],,]
      y.impute <- as.vector(fit$y[imputes[i],])
      designX <- cbind(1,X.impute)
      
      beta0H <- coef(lm(y.impute[A==0]~X.impute[A==0,]))
      beta1H <- coef(lm(y.impute[A==1]~X.impute[A==1,]))
      score0 <- cbind(1,X.impute)*as.vector((y.impute-cbind(1,X.impute)%*%beta0.mle)*(1-A))
      score1 <- cbind(1,X.impute)*as.vector((y.impute-cbind(1,X.impute)%*%beta1.mle)*A)
      
      iid[,i] <- cbind(1,X.impute)%*%beta1.mle - cbind(1,X.impute)%*%beta0.mle -
        t(mu.deriv%*%solve(score.deriv1)%*%t(score1) -
            mu.deriv%*%solve(score.deriv0)%*%t(score0))
      
      trt_effect[i] <- mean(cbind(1,X.impute)%*%(beta1H-beta0H))
      # cov.beta0 <- summary(lm(y.impute[A==0]~X.impute[A==0,]))$cov.unscaled
      # cov.beta1 <- summary(lm(y.impute[A==1]~X.impute[A==1,]))$cov.unscaled
      # var_effect[i] <- sum(cbind(1,X.impute)%*%(cov.beta0+cov.beta1)%*%t(cbind(1,X.impute)))/n^2 +
      # mean((cbind(X.impute)%*%(beta1H[2:3]-beta0H[2:3]))^2)/n
      var_effect[i] <- var(iid[,i])/n
    }
  }else if(method=="IPW"){
    ipw.deriv <- rep(0,px+1)
    info_mat <- matrix(0,px+1,px+1)
    for(i in 1:nsamps){
      X.impute <- fit$X[samps[i],,]
      designX <- cbind(1,X.impute)
      y.impute <- as.vector(fit$y[samps[i],])
      
      ps <- as.vector(pnorm(designX%*%alpha.mle))
      ps.deriv <- as.vector(dnorm(designX%*%alpha.mle))*designX
      ipw.deriv <- ipw.deriv + apply((A*y.impute/ps^2+(1-A)*y.impute/(1-ps)^2)*
                                       ps.deriv,2,mean)/nsamps
      q <- 2*A-1
      lambda <- as.vector(dnorm(q*designX%*%alpha.mle)*q/pnorm(q*designX%*%alpha.mle))
      # info_mat <- info_mat + 
      #   t(lambda*(as.vector(designX%*%alpha.mle)+lambda)*designX)%*%designX/n/nsamps
      info_mat <- info_mat + t(1/(ps*(1-ps))*ps.deriv)%*%ps.deriv/n/nsamps
    }
    
    for(i in 1:nimpute){
      X.impute <- fit$X[imputes[i],,]
      y.impute <- as.vector(fit$y[imputes[i],])
      designX <- cbind(1,X.impute)
      
      ps <- as.vector(pnorm(designX%*%alpha.mle))
      ps.deriv <- as.vector(dnorm(designX%*%alpha.mle))*designX
      score.alpha <- (A-ps)/(ps*(1-ps))*ps.deriv
      
      iid[,i] <- A*y.impute/ps-(1-A)*y.impute/(1-ps) - 
        t(ipw.deriv%*%solve(info_mat)%*%t(score.alpha))
      
      alphaH <- coef(glm(A~X.impute,family=binomial(link='probit')))
      psH <- pnorm(cbind(1,X.impute)%*%alphaH)
      trt_effect[i] <- mean(A*y.impute/psH-(1-A)*y.impute/(1-psH))
      var_effect[i] <- var(iid[,i])/n
    }
  }else if(method=='AIPW'){
    mu.deriv0 <- rep(0,px+1)
    mu.deriv1 <- rep(0,px+1)
    score.deriv0 <- matrix(0,px+1,px+1)
    score.deriv1 <- matrix(0,px+1,px+1)
    aipw.deriv <- rep(0,px+1)
    info_mat <- matrix(0,px+1,px+1)
    for(i in 1:nsamps){
      X.impute <- fit$X[samps[i],,]
      designX <- cbind(1,X.impute)
      y.impute <- as.vector(fit$y[samps[i],])
      
      mu0 <- as.vector(designX%*%beta0.mle)
      mu1 <- as.vector(designX%*%beta1.mle)
      ps <- as.vector(pnorm(designX%*%alpha.mle))
      ps.deriv <- as.vector(dnorm(designX%*%alpha.mle))*designX
      mu.deriv0 <- mu.deriv0 + apply((1-(1-A)/(1-ps))*designX,2,mean)/nsamps
      mu.deriv1 <- mu.deriv1 + apply((1-A/ps)*designX,2,mean)/nsamps
      score.deriv0 <- score.deriv0 + 
        1/n*t(designX[A==0,])%*%designX[A==0,]/nsamps
      score.deriv1 <- score.deriv1 + 
        1/n*t(designX[A==1,])%*%designX[A==1,]/nsamps
      aipw.deriv <- aipw.deriv + apply((A*(y.impute-mu1)/ps^2+
                      (1-A)*(y.impute-mu0)/(1-ps)^2)*ps.deriv,2,mean)/nsamps
      q <- 2*A-1
      lambda <- as.vector(dnorm(q*designX%*%alpha.mle)*
                            q/pnorm(q*designX%*%alpha.mle))
      # info_mat <- info_mat + 
      #   t(lambda*(as.vector(designX%*%alpha.mle)+lambda)*designX)%*%designX/n/nsamps
      info_mat <- info_mat + t(1/(ps*(1-ps))*ps.deriv)%*%ps.deriv/n/nsamps
    }
    
    for(i in 1:nimpute){
      X.impute <- fit$X[imputes[i],,]
      y.impute <- as.vector(fit$y[imputes[i],])
      designX <- cbind(1,X.impute)
      
      mu0 <- as.vector(designX%*%beta0.mle)
      mu1 <- as.vector(designX%*%beta1.mle)
      score0 <- designX*as.vector((y.impute-designX%*%beta0.mle)*(1-A))
      score1 <- designX*as.vector((y.impute-designX%*%beta1.mle)*A)
      ps <- as.vector(pnorm(designX%*%alpha.mle))
      ps.deriv <- as.vector(dnorm(designX%*%alpha.mle))*designX
      score.alpha <- (A-ps)/(ps*(1-ps))*ps.deriv
      
      iid[,i] <- (A*y.impute/ps+(1-A/ps)*mu1)-((1-A)*y.impute/(1-ps)+(1-(1-A)/(1-ps))*mu0)-
        t(aipw.deriv%*%solve(info_mat)%*%t(score.alpha)) +
        t(mu.deriv0%*%solve(score.deriv0)%*%t(score0))-
        t(mu.deriv1%*%solve(score.deriv1)%*%t(score1))
      
      beta0H <- coef(lm(y.impute[A==0]~X.impute[A==0,]))
      beta1H <- coef(lm(y.impute[A==1]~X.impute[A==1,]))
      alphaH <- coef(glm(A~X.impute,family=binomial(link='probit')))
      psH <- pnorm(cbind(1,X.impute)%*%alphaH)
      trt_effect[i] <- mean((A*y.impute/psH+(1-A/psH)*designX%*%beta1H)-
                              ((1-A)*y.impute/(1-psH)+(1-(1-A)/(1-psH))*designX%*%beta0H))
      var_effect[i] <- var(iid[,i])/n
    }
  }else if(method=='matching'){
    coef.mu0 <- rep(0,px+1)
    coef.mu1 <- rep(0,px+1)
    for(i in 1:nsamps){
      X.impute <- fit$X[samps[i],,]
      y.impute <- as.vector(fit$y[samps[i],])
      designX <- cbind(1,X.impute)
      coef.mu0 <- coef.mu0 + coef(lm(y.impute[A==0]~X.impute[A==0,]))/nsamps
      coef.mu1 <- coef.mu1 + coef(lm(y.impute[A==1]~X.impute[A==1,]))/nsamps
    }

    for(i in 1:nimpute){
      X.impute <- fit$X[imputes[i],,]
      y.impute <- as.vector(fit$y[imputes[i],])
      designX <- cbind(1,X.impute)
      
      out <- Match(y.impute,A,X.impute,estimand="ATE",M=M,distance.tolerance=0,ties=FALSE,Weight=2)
      mdata <- out$mdata
      KMi <- table(c(out$index.treated,out$index.control))-M
      mu.pred.comp <- as.vector(designX%*%coef.mu1)*(1-A)+as.vector(designX%*%coef.mu0)*A
      mu.pred <- as.vector(designX%*%coef.mu0)*(1-A)+as.vector(designX%*%coef.mu1)*A
      iid[,i] <- (2*A-1)*(y.impute-mu.pred.comp+KMi/M*(y.impute-mu.pred))
      
      coef.mu0.impute <- coef(lm(y.impute[A==0]~X.impute[A==0,]))
      coef.mu1.impute <- coef(lm(y.impute[A==1]~X.impute[A==1,]))
      mu.pred.comp.impute <- as.vector(designX%*%coef.mu1.impute)*(1-A)+
        as.vector(designX%*%coef.mu0.impute)*A
      mu.pred.impute <- as.vector(designX%*%coef.mu0.impute)*(1-A)+
        as.vector(designX%*%coef.mu1.impute)*A
      
      trt_effect[i] <- mean((2*A-1)*(y.impute-mu.pred.comp.impute+
                                       KMi/M*(y.impute-mu.pred.impute)))
      var_effect[i] <- var(iid[,i])/n
    }
    # else{
    #   stop("The method is not supported")
    #   }
  }
  effect_est_rubin <- mean(trt_effect)
  var_est_rubin <- mean(var_effect)+(1+1/nimpute)*var(trt_effect)
  df <- (nimpute-1)/((1+1/nimpute)*var(trt_effect)/var_est_rubin)^2
  rubin_lower <- effect_est_rubin - qt(0.975,df)*sqrt(var_est_rubin)
  rubin_upper <- effect_est_rubin + qt(0.975,df)*sqrt(var_est_rubin)
  
  score.obs <- matrix(0,n,(px+1)*3+2+px+px*(px+1)/2+(px+2)*(nmis.X+1))
  Gamma1 <- matrix(0,n,(px+1)*3+2+px+px*(px+1)/2+(px+2)*(nmis.X+1))
  iid.samp <- matrix(0,n,nsamps)
  for(i in 1:nsamps){
    X.impute <- fit$X[samps[i],,]
    y.impute <- as.vector(fit$y[samps[i],])
    designX <- cbind(1,X.impute)
    
    if(method=="regression"){
      score0 <- cbind(1,X.impute)*as.vector((y.impute-cbind(1,X.impute)%*%beta0.mle)*(1-A))
      score1 <- cbind(1,X.impute)*as.vector((y.impute-cbind(1,X.impute)%*%beta1.mle)*A)
      
      iid.samp[,i] <- cbind(1,X.impute)%*%beta1.mle - cbind(1,X.impute)%*%beta0.mle -
        t(mu.deriv%*%solve(score.deriv1)%*%t(score1) -
            mu.deriv%*%solve(score.deriv0)%*%t(score0))
      
    }else if(method=='IPW'){
      ps <- as.vector(pnorm(designX%*%alpha.mle))
      ps.deriv <- as.vector(dnorm(designX%*%alpha.mle))*designX
      score.alpha <- (A-ps)/(ps*(1-ps))*ps.deriv
      
      iid.samp[,i] <-  A*y.impute/ps-(1-A)*y.impute/(1-ps) - 
        t(ipw.deriv%*%solve(info_mat)%*%t(score.alpha))
      
    }else if(method=="AIPW"){
      mu0 <- as.vector(designX%*%beta0.mle)
      mu1 <- as.vector(designX%*%beta1.mle)
      score0 <- designX*as.vector((y.impute-designX%*%beta0.mle)*(1-A))
      score1 <- designX*as.vector((y.impute-designX%*%beta1.mle)*A)
      ps <- as.vector(pnorm(designX%*%alpha.mle))
      ps.deriv <- as.vector(dnorm(designX%*%alpha.mle))*designX
      score.alpha <- (A-ps)/(ps*(1-ps))*ps.deriv
      
      iid.samp[,i] <- (A*y.impute/ps+(1-A/ps)*mu1)-((1-A)*y.impute/(1-ps)+(1-(1-A)/(1-ps))*mu0)-
        t(aipw.deriv%*%solve(info_mat)%*%t(score.alpha)) +
        t(mu.deriv0%*%solve(score.deriv0)%*%t(score0))-
        t(mu.deriv1%*%solve(score.deriv1)%*%t(score1))
      
    }else if(method=="matching"){
      out <- Match(y.impute,A,X.impute,estimand="ATE",M=M,distance.tolerance=0,ties=FALSE,Weight=2)
      mdata <- out$mdata
      KMi <- table(c(out$index.treated,out$index.control))-M
      mu.pred.comp <- as.vector(designX%*%coef.mu1)*(1-A)+as.vector(designX%*%coef.mu0)*A
      mu.pred <- as.vector(designX%*%coef.mu0)*(1-A)+as.vector(designX%*%coef.mu1)*A
      iid.samp[,i] <- (2*A-1)*(y.impute-mu.pred.comp+KMi/M*(y.impute-mu.pred))
    }
    # else{
    #   stop("The method is not supported")
    # }
    
    score.beta0 <- designX*as.vector((y.impute-designX%*%beta0.mle)/sigma20.mle*(1-A))
    score.sigma20 <- (-1/(2*sigma20.mle)+(y.impute-designX%*%beta0.mle)^2/(2*sigma20.mle^2))*(1-A)
    score.beta1 <- designX*as.vector((y.impute-designX%*%beta1.mle)/sigma21.mle*A)
    score.sigma21 <- (-1/(2*sigma21.mle)+(y.impute-designX%*%beta1.mle)^2/(2*sigma21.mle^2))*A
    # Q.mle <- solve(covX.mle)
    score.muX <- t(Q.mle%*%(t(X.impute)-muX.mle))
    score.covX <-- array(0,c(n,px,px))
    for(j in 1:n){
      temp <- -1/2*(Q.mle-Q.mle%*%(X.impute[j,]-muX.mle)%*%t(X.impute[j,]-muX.mle)%*%Q.mle)
      score.covX[j,,] <- 2*temp-diag(diag(temp))
    }
    Falpha <- pnorm(designX%*%alpha.mle)
    score.alpha <- as.vector(dnorm(designX%*%alpha.mle)/(Falpha*(1-Falpha))*(A-Falpha))*designX
    Fgammay <- pnorm(cbind(designX,A)%*%gammay.mle)
    Fgammay[which(Fgammay==1)] <- 0.999
    Fgammay[which(Fgammay==0)] <- 0.001
    score.gammay <- as.vector(dnorm(cbind(designX,A)%*%gammay.mle)/
                                (Fgammay*(1-Fgammay))*(obs.loc.y-Fgammay))*
      cbind(designX,A)
    score <- cbind(score.beta0,score.sigma20,score.beta1,score.sigma21,score.alpha,
                   score.gammay,score.muX,matrix(as.vector(score.covX),n,px*px)[,c(1,2,4)])
    
    for(j in 1:nmis.X){
      Fgamma <- pnorm(cbind(designX,A)%*%gamma.mle[,j])
      Fgamma[which(Fgamma==1)] <- 0.999
      Fgamma[which(Fgamma==0)] <- 0.001
      score.gamma <- as.vector(dnorm(cbind(designX,A)%*%gamma.mle[,j])/
                                 (Fgamma*(1-Fgamma))*(obs.loc[,j]-Fgamma))*
        cbind(designX,A)
      score <- cbind(score,score.gamma)
    }
    Gamma1 <- Gamma1 + as.vector(iid.samp[,i])*score/nsamps
    score.obs <- score.obs + score/nsamps
  }
  
  iid.expect <- apply(iid.samp,1,mean)
  Gamma <- apply(Gamma1 - score.obs*as.vector(iid.expect),2,mean)
  info <- t(score.obs)%*%score.obs/n
  part3 <- as.vector(Gamma%*%solve(info)%*%t(score.obs))
  
  effect_bs1 <- rep(0,nbs)
  effect_bs2 <- rep(0,nbs)
  for(bs in 1:nbs){
    wt1 <- rbinom(n*nimpute, 1, (sqrt(5)+1)/(2*sqrt(5))) * (-(sqrt(5)-1)/2)
    wt1[wt1==0] <- (sqrt(5)+1)/2
    wt1 <- matrix(wt1,n,nimpute)
    wt2 <- rbinom(n, 1, (sqrt(5)+1)/(2*sqrt(5))) * (-(sqrt(5)-1)/2)
    wt2[wt2==0] <- (sqrt(5)+1)/2
    
    
    effect_bs1[bs] <- mean((iid-iid.expect)*wt1)
    effect_bs2[bs] <- mean((iid.expect+part3-mean(trt_effect))*wt2)
  }
  
  effect_est_bs <- mean(trt_effect)-mean(effect_bs1+effect_bs2)
  var_est_bs <- var(effect_bs1+effect_bs2)
  
  bs_lower <- mean(trt_effect)-quantile(effect_bs1+effect_bs2,0.975)
  bs_upper <- mean(trt_effect)-quantile(effect_bs1+effect_bs2,0.025)
  bs_lower_wald <- mean(trt_effect)-1.96*sqrt(var_est_bs)
  bs_upper_wald <- mean(trt_effect)+1.96*sqrt(var_est_bs)
  result <- list(ACE.hat=effect_est_rubin,ACE.SE=sqrt(var_est_bs),
                 ci.quantile=c(bs_lower,bs_upper),
                 ci.wald=c(bs_lower_wald,bs_upper_wald))
}