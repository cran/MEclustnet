formatting.covars <-
function(covars, link.vars, mix.vars, n)
{
  # covariate data for link regression. First column is 1.
  covars.link = as.data.frame(covars[,link.vars], nrow=n)
  x.link = as.data.frame(matrix(1, n^2, length(link.vars)))
  colnames(x.link) = colnames(covars.link)
  if(length(link.vars) > 1){
    for(j in 2:ncol(covars.link))
    {
      if(is.factor(covars.link[,j]))
      {
        x.link[,j] = as.numeric(matrix(as.numeric(covars.link[,j]),n,n)-t(matrix(as.numeric(covars.link[,j]),n,n)) == 0)   # Gives a 1 if the covariates are the same
      }else{
        x.link[,j] = c(matrix(covars.link[,j],n,n)-t(matrix(covars.link[,j],n,n)))
      }
    }}

  # Covariate data for mixing proportions. First column is 1.
  covars.mix = as.data.frame(covars[,mix.vars], nrow=n)
  nl = sapply(covars.mix, nlevels)-1
  nl.total = sum(nl[nl!=-1]) + sum(nl==-1)
  x.mix = as.data.frame(matrix(1,n,nl.total))
  colnames(x.mix)[1] = colnames(covars.mix)[1]
  counter = 2
  if(length(mix.vars) > 1){
    for(j in 2:ncol(covars.mix))
    {
      if(is.factor(covars.mix[,j]))
      {
        for(i in 1:nl[j])
        {
          x.mix[,counter] = as.numeric(covars.mix[,j])==i
          colnames(x.mix)[counter] = colnames(covars.mix)[j]
          counter = counter+1
        }
      }else{
        x.mix[,counter] = covars.mix[,j]
        colnames(x.mix)[counter] = colnames(covars.mix)[j]
        counter = counter+1
      }
    }}

  list(x.link = as.matrix(x.link), x.mix = as.matrix(x.mix))
}
