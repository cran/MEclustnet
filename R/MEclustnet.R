MEclustnet <-
function(Y, covars, link.vars = c(1:ncol(covars)), mix.vars = c(1:ncol(covars)), G=2, d=2, itermax = 10000, uphill = 100, burnin = 1000, thin = 10, rho.input = 1, verbose=TRUE, ...)
{
  n = nrow(Y)
  y = c(Y)
  res = formatting.covars(covars, link.vars, mix.vars, n)
  x.link = res[[1]]; x.mix = res[[2]]
  L = ncol(x.link)-1     # no of covariates in link probs
  r = ncol(x.mix)-1      # no of covariates in mixing proportions
  p = L+1; s = r+1       # dimension of coefficient vectors
  delete = seq(from=1, to=n*n, by=(n+1))
  y = y[-delete]
  n.tilde = length(y)
  x.link  = x.link[(-delete),]

  # Start with latent locations from LPCM
  if(verbose) cat("Obtaining starting values via the ergmm function in the latentnet package....\n")
  latentnet.fit = suppressWarnings(ergmm(as.network(Y) ~ euclidean(d=d, G=G)))
  z =  zMAP = latentnet.fit$mcmc.pmode$Z
  K = latentnet.fit$mcmc.pmode$Z.K
  mu = latentnet.fit$mcmc.pmode$Z.mean
  sigma2 = latentnet.fit$mcmc.pmode$Z.var


  delta = c(-as.matrix(dist(z)))[-delete]
  # Obtain starting values for beta from fitting a glm with binary response y
  fitglm = glm(y~x.link+delta-1,family="binomial")
  beta = matrix(coef(fitglm)[1:p],1,p)

  tau = matrix(0, G, s)
  if(G>1)
  {
    if(r>0)
    {
      tau[2:G,]  = coef(multinom(K~cbind(1,x.mix[,2:s])-1, trace=FALSE))/n
    }else{
      tau[2:G,] = coef(multinom(as.factor(K)~matrix(1, n, 1)-1, trace=FALSE))
    }
  }
  lambda = calclambda(tau, x.mix)
  perms = permutations(G)
  m = rep(0,G)
  m = calcm(m, G, K)                    # Vector of counts in each cluster

  Id = diag(1, d)                # d dimensional identity matrix
  pis = calcpis(beta, x.link, delta, n.tilde)  # Vector of link probabilities

  #### Set up required prior parameters ####
  epsilon<-matrix(latentnet.fit$prior$beta.mean, 1, p)            # Mean of betas MVN prior
  psi<-matrix(latentnet.fit$prior$beta.var*diag(p), p, p)        # Covariance matrix of betas MVN prior
  psi.inv = solve(psi)
  sigma02<-latentnet.fit$prior$Z.var      # Scale parameter of Inv-Chi squared prior for sigma2_g
  alpha<-latentnet.fit$prior$Z.var.df/2    # dof parameter of Inv-Chi squared prior for sigma2_g
  omega2<-latentnet.fit$prior$Z.mean.var     # Variance of MVN prior for mu_g
  Sigmag<-matrix(latentnet.fit$prior$beta.var*diag(s), s, s)       # Covariance matrix of MVN prior for tau_g
  #Sigmag = matrix(2*diag(s), s, s)    # Covariance matrix of MVN prior for tau_g
  Sigmag.inv = solve(Sigmag)
  #gammag = matrix(latentnet.fit$prior$beta.mean, 1, s)             # Mean vector of MVN prior for tau_g
  gammag = matrix(1, 1, s)

  store.dim = (itermax-(burnin+uphill))/thin
  zstore = array(NA, c(n, d, store.dim))             # Each sheet is a sampled value
  betastore = matrix(NA, store.dim, p)                # Each row is a sampled value
  Kstore = matrix(NA, store.dim, n)                   # Each row is a sampled value
  mustore = array(NA, c(G, d, store.dim))             # Each sheet is a sampled value
  sigma2store = matrix(NA, store.dim, G)              # Each row is a sampled value
  lambdastore = array(NA, c(n, G, store.dim))         # Each sheet is a sampled value
  taustore = array(NA, c(G, s, store.dim))            # Each sheet is a smapled value
  LLstore = rep(NA,store.dim)                         # Each value is a computed log-likelihood for current parameters.
  countbeta = countz = counttau = 0                 # Record acceptance rates

  if((sum(link.vars)==1 && sum(mix.vars)==1) == F){
    #### Run Metropolis/Gibbs algorithm ####
    if(verbose){
      cat("Fitting the mixture of expert latent position cluster model....\n")
      flush.console()
      pbar = txtProgressBar(min = 1, max = itermax, style = 3)
    }
    for(iter in 1:itermax)
    {
      if(verbose) setTxtProgressBar(pbar, iter)

      if (iter<=uphill){rho = 0.01}else{rho = rho.input}

      res = updatez(n, z, x.link, delta, beta, y, mu, K, sigma2, Id, pis, iter, uphill, countz, delete, d, n.tilde)
      z = res[[1]]
      delta = res[[2]]
      pis = res[[3]]
      countz = res[[4]]
      if(iter == uphill)     # At end of uphill steps store zMAP estimates
      {
        zMAP = z
        zMAP = sweep(zMAP,2,apply(zMAP,2,mean), "-")     # Reference configuration centred at the origin
        z = zMAP
      }
      # Once have zMAP estimates, perform rotation.
      if(iter > uphill){ z = invariant(z, zMAP)}

      res = updatebeta(beta, p, x.link, delta, y, epsilon, psi, psi.inv, pis, countbeta, rho, n.tilde)
      beta = res[[1]]
      countbeta = res[[2]]
      pis = calcpis(beta, x.link, delta, n.tilde)

      res = updatetau(G, x.mix, lambda, Sigmag, Sigmag.inv, K, gammag, tau, counttau, rho)
      tau = res[[1]]
      lambda = res[[2]]
      counttau = res[[3]]

      mu = updatemu(G, z, K, m, sigma2, omega2, Id, mu, d)
      # Store centred, rotated mu's as muMAP for label switching reference
      if(iter == uphill){ muMAP = mu }

      sigma2 = updatesigma2(G, alpha, m, d, sigma02, z, K, mu, sigma2)
      K = updateK(G, K, z, mu, sigma2, Id, lambda)
      m = calcm(m, G, K)

      if((iter/thin) == round(iter/thin) & iter > (uphill+burnin))               # If it's a thinned iteration
      {
        res = labelswitch(mu, sigma2,lambda, tau, K, G, d, perms, muMAP, iter, uphill, burnin, thin, s, x.mix)
        mu = res[[1]]
        sigma2 = res[[2]]
        lambda = res[[3]]
        tau = res[[4]]
        K = res[[5]]
        m = calcm(m , G, K)

        store.ind = (iter-burnin-uphill)/thin
        zstore[,,store.ind] = z                   # Store thinned,rotated samples
        betastore[store.ind,] = beta
        Kstore[store.ind,] = K
        mustore[,,store.ind] = mu
        sigma2store[store.ind,] = sigma2
        lambdastore[,,store.ind] = lambda
        taustore[,,store.ind] = tau
        LLstore[store.ind] = calcloglikelihood(pis,y)
      } #if iter/thin
      #plot(z, col=K)
    } # iter

    if(verbose)  close(pbar)
    res = list(zstore = zstore, betastore = betastore, Kstore = Kstore, mustore = mustore, sigma2store = sigma2store,lambdastore = lambdastore, taustore = taustore, LLstore = LLstore, G = G, d = d, countbeta = countbeta, counttau = counttau)
    res
  }else{
    if(verbose){
      cat("You have requested to fit the latent position cluster model with no covariates.\n")
      cat("Thus an ergmm object from the latentnet package has been be returned.")
    }
    latentnet.fit
  }
} # End MEclustnet
