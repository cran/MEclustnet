plotMEclustnet <-
function(fit, Y, link.vars, mix.vars)
{
  n = nrow(Y)

  zmean = apply(fit$zstore,c(1:2),mean,na.rm=TRUE)
  Kmode = as.numeric(apply(fit$Kstore,2,function(v){names(sort(-table(v)))[1]}))

  plot(zmean,type="n",xlim=c(-7,7), ylim=c(-6,6),col=Kmode,xlab="Dimension 1",ylab="Dimension 2", main = paste("Link variables:", c(link.vars), "Mix variables:", mix.vars))

  mat = matrix(1:n,n,n)
  edges = cbind(mat[Y==1],t(mat)[Y==1])
  for(m in 1:nrow(edges))
  {
    segments(zmean[edges[m,1],1],zmean[edges[m,1],2],zmean[edges[m,2],1],zmean[edges[m,2],2],col="lightgray",lwd=0.75)
  }

  zvar = array(NA,c(2,2,n))
  for(i in 1:n)
  {
    zvar[,,i] = var(t(fit$zstore[i,,]),na.rm=TRUE)
    points(ellipse(zvar[,,i],centre=zmean[i,],level=0.5),type="l",col=Kmode[i])
  }
  points(zmean, col=Kmode, pch=Kmode+1, bg=Kmode, cex=1)
}
