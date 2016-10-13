summaryMEclustnet <-
function(fit, Y)
{
  # Compute AICM
  N = nrow(Y)
  meanll = mean(fit$LLstore,na.rm=TRUE)
  varll = var(fit$LLstore,na.rm=TRUE)
  AICM = 2*(meanll - varll)
  BICM = meanll - varll*(log(N) - 1)
  nu = ncol(fit$betastore) + (fit$G*fit$d) + fit$G + (fit$G*dim(fit$taustore)[2])
  BICMCMC = 2*max(fit$LLstore) - nu*log(N)

  betamean = apply(fit$betastore,2,mean,na.rm=TRUE)
  betasd = apply(fit$betastore,2,sd,na.rm=TRUE)


  taumean = apply(fit$taustore, c(1,2), mean, na.rm=T)
  tausd = apply(fit$taustore, c(1,2), sd, na.rm=T)


  mumean = apply(fit$mustore, c(1,2), mean, na.rm=T)
  musd = apply(fit$mustore, c(1,2), sd, na.rm=T)

  sigma2mean = apply(fit$sigma2store, 2, mean, na.rm=T)
  sigma2sd = apply(fit$sigma2store, 2, sd, na.rm=T)


  Kmode = as.numeric(apply(fit$Kstore,2,function(v){names(sort(-table(v)))[1]}))
  zmean = apply(fit$zstore,c(1:fit$d),mean,na.rm=TRUE)

  list(AICM = AICM, BICM = BICM, BICMCMC = BICMCMC, betamean = betamean, betasd = betasd, taumean = taumean, tausd = tausd, mumean = mumean, musd = musd, sigma2mean = sigma2mean, sigma2sd = sigma2sd, Kmode = Kmode, zmean = zmean)
}
