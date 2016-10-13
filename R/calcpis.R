calcpis <-
function(beta, x.link, delta, n.tilde)
{
  beta.tilde = matrix(c(beta,1), nrow=1)
  x.link.tilde = matrix(c(x.link, delta),ncol=n.tilde, byrow=T)
  res = exp(beta.tilde%*%x.link.tilde)
  c(res/(1 + res))
}
