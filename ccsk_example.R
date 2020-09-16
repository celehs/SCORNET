source("CCSK.R")

sim <- function(N){
  dat <- data.frame('ID'=1:N)
  dat$Z0 <- runif(N,-1,1)
  dat$C <- 10*(rexp(N)*exp(-0.75*dat$Z))^(2/3)
  dat$T <- 15*(rexp(N)*exp(-0.75*dat$Z))^(2/5)
  dat$X <- pmin(dat$T,dat$C)
  dat$Delta <- dat$T <= dat$C
  dat$filter <- as.logical(rbinom(N,1,0.98)*dat$Delta + rbinom(N,1,0.12)*(1-dat$Delta))
  dat$Zehr <- pmin(dat$T+rnorm(N,0,2),dat$C)
  return(dat)
}

N <- 5000
n <- 200
dat <- sim(N)
t0.all <- seq(quantile(dat$C,.1),quantile(dat$C,.9),length.out=100)

ccsk_est <- ccsk(dat$Delta[1:n],dat$C[1:n],t0.all,dat$C[-c(1:n)],dat$filter[1:n],dat$filter[-c(1:n)],
                 dat$Z0[1:n],dat$Z0[-c(1:n)],dat$Zehr[1:n],dat$Zehr[-c(1:n)])

plot(t0.all,ccsk_est$S_hat,type='l',xlab='Time',ylab='Survival')
plot(t0.all,ccsk_est$StdErrs,type='l',xlab='Time',ylab='Std. Err.')

