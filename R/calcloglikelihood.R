calcloglikelihood <-
function(pis, y)
{
  sum((y*log(pis))+((1-y)*log(1-pis)))
}
