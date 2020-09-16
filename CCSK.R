# CCSK.R: Contains CCSK survival estimator function.
# Author: Yuri Ahuja
# Last Updated: 9/15/2020
#
# The Censor-Time Current Status Kernel Estimator (CCSK) is a consistent, non-parametric
# survival curve estimator that boosts efficiency over existing non-parametric estimators
# by (1) utilizing unlabeled patients in a semi-supervised fashion, and (2) leveraging
# information-dense engineered EHR features to maximize unlabeled set imputation precision
# See Ahuja et al. (2020) Submitted to ... for details


library(survival)
library(parallel)
library(doParallel)
library(foreach)
library(Matrix)
library(pracma)

expit <- function(x){
  ifelse(x>=100,1,exp(x)/(1+exp(x)))
}

logit <- function(x){
  log(x/(1-x))
}

# Default kernel function = standard normal PDF
Knorm <- function(t0,t,b=1){
  dnorm(abs(t-t0),sd=b)
}

# Kernel-Smoothed Cox/Breslow estimator of C|Z0
g <- function(C,Z=NULL,K=Knorm,b=1){
  N <- length(C)
  if (is.null(Z)){Z <- matrix(1,nrow=length(C),ncol=1)}
  beta_C.Z <- coxph(Surv(C)~Z)$coefficients; beta_C.Z[is.na(beta_C.Z)] <- 0
  denom <- sapply(1:N,function(j){sum(exp(as.matrix(Z[C >= C[j],]) %*% beta_C.Z))})
  hC0 <- 1 / denom
  hC0k <- foreach(c1=C, .combine=c, .export='K') %dopar% {
    hC0 %*% sapply(C,function(c2){K(c1,c2,b)})
  }
  HC0 <- foreach(c=C, .combine=c, .export='trapz') %dopar% {
    lt <- which(C <= c)
    lt <- lt[order(C[lt])]
    trapz(C[lt],hC0k[lt])
  }
  SC0 <- exp(-HC0)
  h <- exp(Z %*% beta_C.Z)
  
  as.vector(h * hC0k * SC0^h)
}

## CCSK Estimator
# Inputs:
# Delta: Labeled set current status labels (I(T<C)); C: Labeled set censoring times;
# t0.all: Times at which to estimate survival; C.UL: Unlabeled set censoring times;
# filter: Labeled set binary filter indicators; filter.UL: Unlabeled set filter indicators;
# Z0: Labeled set baseline feature matrix; Z0.UL: Unlabeled set baseline feature matrix;
# Zehr: Labeled set EHR-derived feature matrix; Zehr.UL: Unlabeled set EHR-derived feature matrix;
# K: Kernel function (defaults to standard normal) with bandwidth b (default set heuristically)
# ghat: N^1/4-consistent pdf estimator of C|Z0 (defaults to Kernel-Smoothed Cox/Breslow estimator)
# Outputs:
# S_hat: Survival function estimates at times t0.all
# StdErrs: Asymptotically consistent standard error estimates corresponding to S_hat
ccsk <- function(Delta,C,t0.all,C.UL=NULL,filter=NULL,filter.UL=NULL,Z0=NULL,Z0.UL=NULL,
                     Zehr=NULL,Zehr.UL=NULL,K=Knorm,b=NULL,ghat=NULL){
  Ctot <- c(C,C.UL)
  N <- length(C)
  Ntot <- length(Ctot)
  if (is.null(filter)){filter <- rep(TRUE,N)}
  if (is.null(filter.UL)){filter.UL <- rep(TRUE,length(C.UL))}
  filtertot <- c(filter,filter.UL)
  if (!is.null(Z0)){Z0 <- as.matrix(Z0)} else{Z0 <- matrix(1,N,1)}
  if (is.null(Z0.UL)){Z0.UL <- rep(1,length(C.UL))}
  Z0tot <- rbind(Z0,as.matrix(Z0.UL))
  if (!is.null(Zehr)){Zehr <- as.matrix(Zehr)}
  Zehrtot <- Zehr; if (!is.null(Zehr.UL)){Zehrtot <- rbind(Zehrtot,as.matrix(Zehr.UL))}
  Cfp <- Ctot[filtertot]
  Z0fp <- Z0tot[filtertot,]
  Zehrfp <- Zehrtot[filtertot,]
  Nfp <- sum(filtertot)
  
  if (is.null(b)){
    b <- N^(-1/3) * min(sd(Cfp), (quantile(Cfp,0.75)-quantile(Cfp,0.25))/1.34)
  }
  nu <- Ntot^(-1/4) * min(sd(Ctot), (quantile(Ctot,0.75)-quantile(Ctot,0.25))/1.34)
  
  if (is.null(ghat)){
    ghat <- g(Ctot,Z0tot,K,nu)
  }
  
  # Estimate S_T
  Kmat <- sapply(t0.all,function(t){
    sapply(Ctot,function(c){
      K(c,t,b)
    })
  })
  beta_T.Z <- sapply(1:length(t0.all),function(i){
    glm(Delta~cbind(Z0,Zehr), family=binomial, weights=filter*Kmat[1:N,i])$coef
  })
  beta_T.Z[is.na(beta_T.Z)] <- 0
  
  Fi_hat <- sapply(1:length(t0.all),function(i){
    t <- t0.all[i]
    expit(cbind(1,Z0tot,Zehrtot) %*% beta_T.Z[,i])
  })
  Fi_hat[!filtertot,] <- 0
  S_hat <- 1 - (colSums(Fi_hat*Kmat/ghat) / colSums(Kmat/ghat))
  
  # Estimate Std Errs
  tau_sq <- 1 / (2*sqrt(pi))
  StdErrs <- sapply(1:length(t0.all),function(i){
    tryCatch({
      W <- cbind(1,Z0fp,Zehrfp)
      gTW <- expit(W %*% beta_T.Z[,i])
      gTWprime <- gTW * (1-gTW)
      if (sum(gTWprime)==0){NA}
      else{
        At <- matrix(0,ncol(W),ncol(W))
        fps <- which(filtertot)
        for (j in 1:Nfp){
          At <- At + gTWprime[j] * Kmat[fps[j],i] * (W[j,] %*% t(W[j,]))
        }
        At <- At / sum(Kmat[fps,i])
        ABA <- solve(At) / (b * sum(Kmat[1:N,i]) / tau_sq)
        
        gTWprime <- Fi_hat[filtertot,i] * (1-Fi_hat[filtertot,i])
        P <- t(W) %*% (Kmat[filtertot,i]*gTWprime/ghat[filtertot]) / sum(Kmat[,i]/ghat)
        ABAP <- c(t(P) %*% ABA %*% P)
        
        sqrt(ABAP)
      }      
    }, error=function(e){NA})
  })
  
  return(list('S_hat'=S_hat, 'StdErrs'=StdErrs))
}
