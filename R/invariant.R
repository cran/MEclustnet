invariant <-
function(z, zMAP)
{
  procrustes(zMAP, z, scale = FALSE)$Yrot
}
