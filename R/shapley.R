#' Decomposition of squared Mahalanobis distance using Shapley values.
#'
#' @description The \code{shapley} function computes a \eqn{p}-dimensional vector containing the decomposition of the
#' squared Mahalanobis distance of \code{x} (with respect to \code{mu} and \code{Sigma})
#' into outlyingness contributions of the individual variables \insertCite{Mayrhofer2022}{ShapleyOutlier}.
#' The value of the \eqn{j}-th coordinate of this vector represents the
#' average marginal contribution of the \eqn{j}-th variable to the squared Mahalanobis distance of
#' the individual observation \code{x}.\cr
#' If \code{cells} is provided, Shapley values of \code{x} are computed with respect to a local reference point,
#' that is based on a cellwise prediction of each coordinate, using the information of the regular cells of \code{x}, see \insertCite{Mayrhofer2022}{ShapleyOutlier}.\cr
#' If \code{x} is a \eqn{n \times p} matrix, a \eqn{n \times p} matrix is returned, containing the decomposition for each row.
#'
#' @param x Data vector with \eqn{p} entries or data matrix with \eqn{n \times p} entries containing only numeric entries.
#' @param mu Either \code{NULL} (default) or mean vector of \code{x}. If NULL, \code{method} is used for parameter estimation.
#' @param Sigma Either \code{NULL} (default) or covariance matrix \eqn{p \times p} of \code{x}. If NULL, \code{method} is used for parameter estimation.
#' @param inverted Logical. If \code{TRUE}, \code{Sigma} is supposed to contain the inverse of the covariance matrix.
#' @param method Either "cellMCD" (default) or "MCD". Specifies the method used for parameter estimation if \code{mu} and/or \code{Sigma} are not provided.
#' @param check Logical. If \code{TRUE} (default), inputs are checked before running the function
#' and an error message is returned if one of the inputs is not as expected.
#' @param cells Either \code{NULL} (default) or a vector/matrix of the same dimension as \code{x},
#' indicating the outlying cells. The matrix must contain only zeros and ones, or \code{TRUE}/\code{FALSE}.
#'
#' @references
#' \insertRef{Mayrhofer2022}{ShapleyOutlier}
#'
#' @return
#' \item{\code{phi}}{A \eqn{p}-dimensional vector (or a \eqn{n \times p} matrix) containing the Shapley values (outlyingness-scores) of \code{x}.}
#' \item{\code{mu_tilde}}{A \eqn{p}-dimensional vector (or a \eqn{n \times p} matrix) containing the
#' alternative reference points based on the regular cells of the original observations.}
#' \item{\code{non_centrality}}{The non-centrality parameters for the Chi-Squared distribution, given by \code{mahlanobis(mu_tilde, mu, Sigma)}}
#' @export
#'
#' @examples
#' ## Without outlying cells as input in the 'cells' argument#'
#' # Single observation
#' p <- 5
#' mu <- rep(0,p)
#' Sigma <- matrix(0.9, p, p); diag(Sigma) = 1
#' Sigma_inv <- solve(Sigma)
#' x <- c(0,1,2,2.3,2.5)
#' shapley(x, mu, Sigma)
#' phi <- shapley(x, mu, Sigma_inv, inverted = TRUE)
#' plot(phi)
#'
#' # Multiple observations
#' library(MASS)
#' set.seed(1)
#' n <- 100; p <- 10
#' mu <- rep(0,p)
#' Sigma <- matrix(0.9, p, p); diag(Sigma) = 1
#' X <- mvrnorm(n, mu, Sigma)
#' X_clean <- X
#' X[sample(1:(n*p), 100, FALSE)] <- rep(c(-5,5),50)
#' call_shapley <- shapley(X, mu, Sigma)
#' plot(call_shapley, subset = 20)
#'
#'
#' ## Giving outlying cells as input in the 'cells' argument
#' # Single observation
#' p <- 5
#' mu <- rep(0,p)
#' Sigma <- matrix(0.9, p, p); diag(Sigma) = 1
#' Sigma_inv <- solve(Sigma)
#' x <- c(0,1,2,2.3,2.5)
#' call_shapley <- shapley(x, mu, Sigma_inv, inverted = TRUE,
#' method = "cellMCD", check = TRUE, cells = c(1,1,0,0,0))
#' plot(call_shapley)
#'
#' # Multiple observations
#' library(MASS)
#' set.seed(1)
#' n <- 100; p <- 10
#' mu <- rep(0,p)
#' Sigma <- matrix(0.9, p, p); diag(Sigma) = 1
#' X <- mvrnorm(n, mu, Sigma)
#' X_clean <- X
#' X[sample(1:(n*p), 100, FALSE)] <- rep(c(-5,5),50)
#' call_shapley <- shapley(X, mu, Sigma, cells = (X_clean - X)!=0)
#' plot(call_shapley, subset = 20)
shapley <- function(x, mu = NULL, Sigma = NULL, inverted = FALSE, method = "cellMCD", check = TRUE, cells = NULL){
  if(check){# Only check inputs if check == TRUE, use check = FALSE to avoid overhead of checks
    if(is.null(mu) | is.null(Sigma)){# If mu and Sigma are not provided, try to estimate them
      if(is.null(dim(x))){# For vectors, estimation is not possible
        stop("If mu and/or Sigma are not provided, x must be a matrix n times p matrix (n must be larger than p).")
      } else {
        if(nrow(x) <= ncol(x)){# Check if there are enough observations for estimation
          stop("If mu and/or Sigma are not provided, n must be larger than p")
        }
        message(paste("Parameters mu and/or Sigma are not provided --", method, "-- is used for estimation."))
        estimates <- initial_estimates(x, method = method)
        inverted = FALSE
        if(is.null(mu)){
          mu <- estimates$mu
        }
        Sigma <- estimates$Sigma
      }
    }
  }
  if(is.null(dim(x))){# VECTOR
    if(check){
      p <- length(x)
      check_vector(x)
      check_vector(mu, p)
      check_matrix(Sigma, p)
    }
    if(is.null(cells)){
      S = NULL
    } else{
      if(check){
        check_indicator(cells, p)
      }
      S = which(cells == 1)
    }
    if(inverted){# Check if matrix must be inverted
      Sigma_inv <- Sigma
    } else{
      Sigma_inv <- solve(Sigma)
    }
    # ACTUAL COMPUTATION OF SHAPLEY VALUE for a single instance
    if(is.null(S)|length(S) == 0){
      mu_tilde <- mu
      phi <- shapley_base(x, mu, Sigma_inv, inverted = TRUE)
      non_centrality <- 0
    } else{
      mu_tilde <- min_MD(x = x, mu = mu, Sigma_inv = Sigma_inv, S = S)
      phi <- shapley_base(x, mu_tilde, Sigma_inv, inverted = TRUE)
      non_centrality <-   mahalanobis(mu_tilde, mu, Sigma_inv, inverted = TRUE)
    }
  } else{ # MATRIX
    x <- as.matrix(x)
    p <- ncol(x)
    n <- nrow(x)
    check_data(x)
    check_vector(mu, p)
    check_matrix(Sigma, p)
    if(!is.null(cells)){
      check_indicator_matrix(cells, p = p, n = n)
    }
    if(inverted){# Check if matrix must be inverted
      Sigma_inv <- Sigma
    } else{
      Sigma_inv <- solve(Sigma)
    }
    phi <- array(dim = dim(x), dimnames = dimnames(x))
    mu_tilde <- array(dim = dim(x))
    non_centrality <- numeric()
    for(i in 1:nrow(x)){
      call_shapley <- shapley(x[i,], mu, Sigma_inv, inverted = TRUE, method = method, check = FALSE, cells = cells[i,])
      phi[i,] <- call_shapley$phi
      mu_tilde[i,] <- call_shapley$mu_tilde
      non_centrality[i] <- call_shapley$non_centrality
    }
  }
  res <- new_shapley(phi, mu_tilde, non_centrality)
  return(res)
}

