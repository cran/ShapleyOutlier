min_MD <- function(x, mu, Sigma_inv = NULL, S = numeric()){
  mu_tilde <- numeric()
  p <- length(x)
  if(p != length(S)){
    for(j in 1:p){
      sub <- unique(c(S, j))
      mu_tilde[j] <- (x[j]  - solve(Sigma_inv[sub,sub])[,which(sub == j)]%*%Sigma_inv[sub,]%*%(x - mu))
    }
  } else {
    mu_tilde <- mu
  }
  return(mu_tilde)
}
