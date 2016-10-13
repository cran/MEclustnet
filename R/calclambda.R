calclambda <-
function(tau, x.mix)
{
  res = tau%*%t(x.mix)
  res = t(exp(res-apply(res,2,max)))
  sweep(res, 1, rowSums(res), "/")
}
