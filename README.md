cd_BAR =  function(y, x, lambda, max.iter = 100, eps = 1e-8)
{
  n =  length(y)
  p =  ncol(x)
  
  # standardize
  u =  y - mean(y) #scale(y,scale=F)
  z =  t(t(x) - apply(x, 2, mean))
  norm.z =  apply(z^2, 2, mean)
  z =  t(t(z)/sqrt(norm.z)) #scale(x)
  
  # initialize beta
  init = rep(1,p)
  beta = init
  
  # residual 
  resid =  (u - z %*% beta)
  
  # start update
  for (t in 1:max.iter)
  {
    new.beta =  beta
    for (j in 1:p)                                  
    {
      zj = crossprod(z[,j], resid)/n + beta[j]     
      new.beta[j] = ST.BAR(zj, lambda, n)         
      resid =  resid - z[,j] * (new.beta[j] - beta[j]) 
    }     
    if (max(abs(beta - new.beta)) < eps) break
    beta =  new.beta
  }
  
  # transform back
  return(beta / sqrt(norm.z))
}
