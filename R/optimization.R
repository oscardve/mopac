# MoPAC barrier method implementation functions
# version: 1.0
# date: May 2018
# description: performs minimization of control gene variance through the barrier method
# author:  Oscar Villarreal
# affiliation: University of Texas MD Anderson Cancer Center. Laboratory of Dr. Han Xu
# contact: oscardvillarreal AT gmail.com

do.differentiate <- function(mup,mun,w0,Di,Dijp,Dijn,dDik,dDpijk,dDnijk,t,tau){
  # Calculate gradient for sgRNA rank weights:
  gradt <- tensorA::to.tensor(-c(0, (1/w0[-1] - 1/(1-sum(w0[-1]))) / t), c(k=length(w0)))
  gradp <- tensorA::mul.tensor(dDpijk,"i",Dijp/Di^2,"i",by=c("j","k")) -
    tensorA::mul.tensor(Dijp^2/Di^3,"i",dDik,"i",by=c("j","k"))
  gradn <- tensorA::mul.tensor(dDnijk,"i",Dijn/Di^2,"i",by=c("j","k")) -
    tensorA::mul.tensor(Dijn^2/Di^3,"i",dDik,"i",by=c("j","k"))
  gradient <- gradt + 2*(mup-mun)^2 * (tensorA::mean.tensor(gradp,"j") + tensorA::mean.tensor(gradn,"j"))
  # Calculate hessian for sgRNA weights:
  if(length(w0)==2) {
    hesst <- tensorA::to.tensor(cbind(0,rbind(0,(1/t)*(diag(as.matrix(1/w0[-1]^2))+(1/(1-sum(w0[-1])))^2))))
  } else hesst <- tensorA::to.tensor(cbind(0,rbind(0,(1/t)*(diag(1/w0[-1]^2)+(1/(1-sum(w0[-1])))^2))))
  names(hesst) <- c("k'","k*")
  hess1p <- tensorA::einstein.tensor( tensorA::mark(dDpijk,"'","k") , tensorA::mark(dDpijk/(ncol(Dijp)*Di^2),"*","k") )
  hess1n <- tensorA::einstein.tensor( tensorA::mark(dDnijk,"'","k") , tensorA::mark(dDnijk/(ncol(Dijn)*Di^2),"*","k") )
  hessp <- -2* tensorA::einstein.tensor( tensorA::mark(tensorA::mean.tensor(dDpijk*Dijp,"j"),"'","k") , tensorA::mark(dDik/Di^3,"*","k") )
  hessn <- -2* tensorA::einstein.tensor( tensorA::mark(tensorA::mean.tensor(dDnijk*Dijn,"j"),"'","k") , tensorA::mark(dDik/Di^3,"*","k") )
  hess2p <- hessp + tensorA::mark(t(hessp),c("k'","k*"))
  hess2n <- hessn + tensorA::mark(t(hessn),c("k'","k*"))
  hess3p <- 3* tensorA::einstein.tensor( tensorA::mark(tensorA::mean.tensor(Dijp^2,"j")*dDik,"'","k") , tensorA::mark(dDik/Di^4,"*","k") )
  hess3n <- 3* tensorA::einstein.tensor( tensorA::mark(tensorA::mean.tensor(Dijn^2,"j")*dDik,"'","k") , tensorA::mark(dDik/Di^4,"*","k") )
  hessian <- 2*(mup-mun)^2 * (hess1p+hess1n+hess2p+hess2n+hess3p+hess3n) + hesst
  # Calculate Newton step and Newton decrement:
  step <- -solve(hessian[-1,-1,drop=F],gradient[-1,drop=F],i="k'",j="k"); names(step) <- "k"
  lambda2 <- -gradient[-1]%*%(step*tau)
  return(list(Step=step,Lambda2=lambda2))
}

do.optimize <- function(xp,xn,mup=-1,mun=0,t=1000,mu=100,eps=1E-9,alpha=0.1,beta=0.5){
  # Dataset constants:
  mode <- dim(xp)[3] #number of guides per gene
  xpt <- tensorA::to.tensor(xp); names(xpt) <- c("i","j","k")
  xnt <- tensorA::to.tensor(xn); names(xnt) <- c("i","j","k")
  Dik <- tensorA::mean.tensor(xpt,"j") - tensorA::mean.tensor(xnt,"j")
  Dpijk <- xpt - tensorA::mean.tensor(xpt,"j")
  Dnijk <- xnt - tensorA::mean.tensor(xnt,"j")
  dDik <- Dik - Dik[[k=1]]
  dDpijk <- Dpijk - Dpijk[[k=1]]
  dDnijk <- Dnijk - Dnijk[[k=1]]
  # Initial guess of the parameters:
  variance <- variancep <- variancen <- numeric()
  w0 <- tensorA::to.tensor(rep(1,mode)/mode, c(k=mode)); w1 <- w0; w <- t(w0)
  A0 <- (mup-mun)/(tensorA::einstein.tensor(Dik,w0)); A <- t(A0)
  B0 <- mup - A0 * tensorA::mean.tensor(tensorA::einstein.tensor(xpt,w0),"j"); B <- t(B0)
  # Barrier method with stopping criterion -> (# of log-barriers)/t > eps :
  while(length(w0)/t > eps){
    # Get Newton step and Newton decrement:
    Di <- Dik[[k=1]] + tensorA::einstein.tensor(dDik,w0)
    Dijp <- Dpijk[[k=1]] + tensorA::einstein.tensor(dDpijk,w0)
    Dijn <- Dnijk[[k=1]] + tensorA::einstein.tensor(dDnijk,w0)
    diff <- do.differentiate(mup,mun,w0,Di,Dijp,Dijn,dDik,dDpijk,dDnijk,t,tau=1)
    step <- diff$Step; lambda2 <- diff$Lambda2
    while(lambda2/2 > eps){
      tau <- 1
      w1[-1] <- w0[-1] + tau*step
      while(min(w1[-1])<0){
        tau <- tau * beta
        w1[-1] <- w0[-1] + tau*step
      }
      while(1-sum(w1[2:mode])<0){
        tau <- tau * beta
        w1[-1] <- w0[-1] + tau*step
      }
      current <- (1/t) * (log(1-sum(w0[-1])) + log(prod(w0[-1]))) + (mup-mun)^2 *
        (tensorA::mean.tensor(tensorA::margin.tensor((Dijp/Di)^2,"i"),"j") +
           tensorA::mean.tensor(tensorA::margin.tensor((Dijn/Di)^2,"i"),"j"))
      w1[1] <- 1-sum(w1[-1])
      w1 <- tensorA::to.tensor(w1, c(k=mode))
      Di <- Dik[[k=1]] + tensorA::einstein.tensor(dDik,w1)
      Dijp <- Dpijk[[k=1]] + tensorA::einstein.tensor(dDpijk,w1)
      Dijn <- Dnijk[[k=1]] + tensorA::einstein.tensor(dDnijk,w1)
      future <- (1/t) * (log(1-sum(w1[-1])) + log(prod(w1[-1]))) + (mup-mun)^2 *
        (tensorA::mean.tensor(tensorA::margin.tensor((Dijp/Di)^2,"i"),"j") +
           tensorA::mean.tensor(tensorA::margin.tensor((Dijn/Di)^2,"i"),"j"))
      while(future > (current - alpha*tau*lambda2)){
        tau <- tau * beta
        w1[-1] <- w0[-1] + tau*step
        w1[1] <- 1-sum(w1[-1])
        w1 <- tensorA::to.tensor(w1, c(k=mode))
        Di <- Dik[[k=1]] + tensorA::einstein.tensor(dDik,w1)
        Dijp <- Dpijk[[k=1]] + tensorA::einstein.tensor(dDpijk,w1)
        Dijn <- Dnijk[[k=1]] + tensorA::einstein.tensor(dDnijk,w1)
        future <- (1/t) * (log(1-sum(w1[-1])) + log(prod(w1[-1]))) + (mup-mun)^2 *
          (tensorA::mean.tensor(tensorA::margin.tensor((Dijp/Di)^2,"i"),"j") +
             tensorA::mean.tensor(tensorA::margin.tensor((Dijn/Di)^2,"i"),"j"))
      }
      # Update parameters and variance:
      w0 <- w1
      w <- rbind(w,t(w0))
      A0 <- (mup-mun)/(tensorA::einstein.tensor(Dik,w0))
      B0 <- mup - A0 * tensorA::mean.tensor(tensorA::einstein.tensor(xpt,w0),"j")
      A <- rbind(A,t(A0))
      B <- rbind(B,t(B0))
      variancepi <- (mup-mun)^2 * tensorA::mean.tensor(tensorA::margin.tensor((Dijp/Di)^2,"i"),"j")
      varianceni <- (mup-mun)^2 * tensorA::mean.tensor(tensorA::margin.tensor((Dijn/Di)^2,"i"),"j")
      variancep <- append(variancep, variancepi)
      variancen <- append(variancen, varianceni)
      variance <- append(variance, variancepi + varianceni)
      # Get Newton step and Newton decrement:
      Di <- Dik[[k=1]] + tensorA::einstein.tensor(dDik,w0)
      Dijp <- Dpijk[[k=1]] + tensorA::einstein.tensor(dDpijk,w0)
      Dijn <- Dnijk[[k=1]] + tensorA::einstein.tensor(dDnijk,w0)
      diff <- do.differentiate(mup,mun,w0,Di,Dijp,Dijn,dDik,dDpijk,dDnijk,t,tau)
      step <- diff$Step; lambda2 <- diff$Lambda2
    }
    t <- t*mu
  }
  return(list(Weights=w,A=A,B=B,Variance=variance,Variancep=variancep,Variancen=variancen))
}
