initial_estimates <- function(X, method = "cellMCD", alpha = 0.5){
  if(method == "cellMCD"){
    call_estimator <- cellWise::cellMCD(X, alpha = alpha, checkPars = list(coreOnly = TRUE, silent = TRUE))
    mu <- call_estimator$mu
    Sigma <- call_estimator$S
  } else if(method == "MCD"){
    call_estimator <- robustbase::covMcd(X, alpha = alpha)
    mu <- call_estimator$center
    Sigma <- call_estimator$cov
  } else{
    stop("Input -- method -- must be either cellMCD or MCD")
  }
  return(list("mu" = mu, "Sigma" = Sigma))
}
