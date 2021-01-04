`%do%` <- foreach::`%do%`
`%dopar%` <- foreach::`%dopar%`
utils::globalVariables("i")

#' @importFrom stats binomial dnorm glm quantile sd
NULL

# SCORNET.R: Contains SCORNET survival estimator function.
# Author: Yuri Ahuja
# Last Updated: 12/15/2020
#
# Semi-Supervised Calibration of Risk with Noisy Event Times (SCORNET) is a consistent, non-parametric
# survival curve estimator that boosts efficiency over existing non-parametric estimators
# by (1) utilizing unlabeled patients in a semi-supervised fashion, and (2) leveraging
# information-dense engineered EHR features to maximize unlabeled set imputation precision
# See Ahuja et al. (2020) BioArxiv for details

expit <- function(x){
  1/(1+exp(-x))
}

logit <- function(x){
  log(x/(1-x))
}

# Default kernel function: standard normal PDF
Knorm <- function(t0,t,b=1){
  dnorm(abs(t-t0),sd=b)
}

# Kernel-Smoothed Cox/Breslow estimator of C|Z0
estimate.h <- function(C,Z=NULL,b=NULL,nCores=1){
  N <- length(C)
  if (is.null(Z)){Z <- matrix(1,nrow=length(C),ncol=1)}
  if (is.null(b)){b <- N^(-1/4) * min(sd(C), (quantile(C,0.75)-quantile(C,0.25))/1.34)}
  
  beta_C.Z <- survival::coxph(survival::Surv(C)~Z)$coefficients; beta_C.Z[is.na(beta_C.Z)] <- 0
  
  denom <- sapply(1:N,function(i){sum(exp(as.matrix(Z[C>=C[i],]) %*% beta_C.Z))})
  hC0 <- 1 / denom
  hC0k <- c(kernelSmoothen(hC0,C,b))
  
  if (nCores == 1){
    HC0 <- foreach::foreach(i=1:N, .combine=c) %do% {sum(hC0[C<=C[i]])}
  }
  else{
    HC0 <- foreach::foreach(i=1:N, .combine=c) %dopar% {sum(hC0[C<=C[i]])}
  }
  SC0 <- exp(-HC0)
  h <- exp(Z %*% beta_C.Z)
  
  as.vector(h * hC0k * SC0^h)
}

#' SCORNET Estimator
#' @param Delta Labeled set current status labels (I(T<C))
#' @param C Labeled set censoring times
#' @param t0.all Times at which to estimate survival 
#' @param C.UL Unlabeled set censoring times
#' @param filter Labeled set binary filter indicators
#' @param filter.UL Unlabeled set filter indicators
#' @param Z0 Labeled set baseline feature matrix
#' @param Z0.UL Unlabeled set baseline feature matrix
#' @param Zehr Labeled set EHR-derived feature matrix
#' @param Zehr.UL Unlabeled set EHR-derived feature matrix
#' @param K Kernel function (defaults to standard normal) 
#' @param b bandwidth (optional)
#' @param bexp bandwidth exponent (must be between -1/5 and -1/3, defaults to -1/4)
#' @param fc N^1/4-consistent pdf estimator of C|Z0 (defaults to Kernel-Smoothed Cox/Breslow estimator)
#' @param nCores Number of cores to use for parallelization (defaults to 1)
#' @return S_hat: Survival function estimates at times t0.all; StdErrs: Asymptotically consistent standard error estimates corresponding to S_hat
#' @export
scornet <- function(Delta, C, t0.all, C.UL = NULL, filter = NULL, filter.UL = NULL, Z0 = NULL, Z0.UL = NULL,
                    Zehr = NULL, Zehr.UL = NULL, K = Knorm, b = NULL, bexp = -1/4, fc = NULL, nCores = 1) {
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
  
  if (nCores > 1){
    logfile <- "SCORNET.log"
    writeLines(c(""), file(logfile,'w'))
    clust <- parallel::makeCluster(nCores, outfile=logfile)
    doParallel::registerDoParallel(clust)
  }
  
  if (is.null(b)){
    b <- N^bexp * min(sd(Cfp), (quantile(Cfp,0.75)-quantile(Cfp,0.25))/1.34)
  }
  nu <- Ntot^(-1/4) * min(sd(Ctot), (quantile(Ctot,0.75)-quantile(Ctot,0.25))/1.34)
  
  
  # STEP 1: Estimate conditional censoring density f(C|Z0)
  
  if (is.null(fc)){
    fc <- estimate.h(Ctot,Z0tot,nu,nCores)
  }
  

  # STEP 2: Train imputation model P(T<=t|Z,Z0)
  
  Kmat1 <- outer(C,t0.all,function(x,y){K(x,y,b)})
  weights1 <- filter * Kmat1 / fc[1:N]
  beta_T.Z <- sapply(1:length(t0.all),function(i){
    glm(Delta~cbind(Z0,Zehr), family=quasibinomial, weights=weights1[,i])$coef
  })
  beta_T.Z[is.na(beta_T.Z)] <- 0
  

  # STEP 3: Estimate marginal survival function S_T(t)
  
  Kmat2 <- outer(Ctot,t0.all,function(x,y){K(x,y,nu)})
  weights2 <- Kmat2 / fc
  Fi_hat <- expit(cbind(1,Z0tot,Zehrtot) %*% beta_T.Z)
  Fi_hat[!filtertot,] <- 0
  S_hat <- 1 - (colSums(Fi_hat*weights2) / colSums(weights2))
  
  
  # Estimate standard errors for S(t)
  
  filterprop <- apply(weights2,2,function(wt){(wt %*% filtertot) / sum(wt)})
  StdErrs <- sapply(1:length(t0.all),function(i){
    tryCatch({
      V_t <- sum((Fi_hat[1:N,i] - Delta)^2 * weights1[,i]^2) * filterprop[i] / N / (b*N)
      sqrt(V_t)
    }, error=function(e){NA})
  })
  
  
  return(list('S_hat'=S_hat, 'StdErrs'=StdErrs))
}
