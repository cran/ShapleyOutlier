shapley_base <- function(x, mu, Sigma, inverted = FALSE){
  if(inverted){
    Sigma_inv <- Sigma
  } else{
    Sigma_inv <- solve(Sigma)
  }
  return((x-mu)*crossprod((x-mu),Sigma_inv))
}
