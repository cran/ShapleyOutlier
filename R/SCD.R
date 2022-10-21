#' Detecting cellwise outliers using Shapley values.
#'
#' @description The \code{SCD} function indicates outlying cells for
#' a data vector with \eqn{p} entries or data matrix with \eqn{n \times p} entries containing only numeric entries \code{x}
#' for a given center \code{mu} and covariance matrix \code{Sigma} using the Shapley value \insertCite{Mayrhofer2022}{ShapleyOutlier}.
#'
#' @param Sigma_inv Either \code{NULL} (default) or Sigma's inverse \eqn{p \times p} matrix.
#' If \code{NULL}, the inverse of \code{Sigma} is computed using \code{solve(Sigma)}.
#' @param step_size Numeric. Step size for the imputation of outlying cells, with \code{step_size} \eqn{\in [0,1]}. Defaults to \eqn{0.1}.
#' @param min_deviation Numeric. Detection threshold, with \code{min_deviation} \eqn{\in [0,1]}. Defaults to \eqn{0.2}
#' @param max_step Either \code{NULL} (default) or an integer. The maximum number of steps in each iteration.
#' If \code{NULL}, \code{max_step} \eqn{= p}.
#' @param max_iter Integer. The maximum number of iterations.
#' @param q Numeric. The quantile of the Chi-squared distribution for detection and imputation of outliers. Defaults to \eqn{0.99}.
#' @inheritParams shapley
#'
#' @return A list of class \code{shapley_algorithm} (\code{\link{new_shapley_algorithm}}) containing the following:
#' \item{\code{x}}{A \eqn{p}-dimensional vector (or a \eqn{n \times p} matrix) containing the imputed data.}
#' \item{\code{phi}}{A \eqn{p}-dimensional vector (or a \eqn{n \times p} matrix) containing the Shapley values (outlyingness-scores) of \code{x}; see \code{\link{shapley}}.}
#' \item{\code{x_original}}{A \eqn{p}-dimensional vector (or a \eqn{n \times p} matrix) containing the original data.}
#' \item{\code{x_history}}{The path of how the original data vector was modified.}
#' \item{\code{phi_history}}{The Shapley values corresponding to \code{x_history}.}
#' \item{\code{S_history}}{The indices of the outlying cells in each iteration.}
#'
#' @references
#' \insertAllCited{}
#'
#' @export
#'
#' @examples
#' p <- 5
#' mu <- rep(0,p)
#' Sigma <- matrix(0.9, p, p); diag(Sigma) = 1
#' Sigma_inv <- solve(Sigma)
#' x <- c(0,1,2,2.3,2.5)
#' SCD_x <- SCD(x = x, mu = mu, Sigma = Sigma)
#' plot(SCD_x)
#'
#' library(MASS)
#' set.seed(1)
#' n <- 100; p <- 10
#' mu <- rep(0,p)
#' Sigma <- matrix(0.9, p, p); diag(Sigma) = 1
#' X <- mvrnorm(n, mu, Sigma)
#' X[sample(1:(n*p), 100, FALSE)] <- rep(c(-5,5),50)
#' SCD_X <- SCD(X, mu, Sigma)
#' plot(SCD_X, subset = 20)
SCD <- function(x, mu, Sigma, Sigma_inv = NULL, step_size = 0.1, min_deviation = 0, max_step = NULL, max_iter = 1000, q = 0.99,
                method = "cellMCD", check = TRUE, cells = NULL){
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
    p <- length(x)
    if(check){
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
    if(is.null(Sigma_inv)){
      Sigma_inv <- solve(Sigma)
    }
    check_matrix(Sigma_inv, p)
    # ACTUAL COMPUTATION SCD
      x_tilde <- x
      x_final <- x
      x_history <- x_tilde
      S_history <- list()

      #md of original observation
      phi <- shapley_base(x, mu, Sigma_inv, inverted = TRUE)
      phi_start <- phi
      phi_history <- numeric()
      md <- sum(phi)

      if(md > qchisq(q, p)){
        S_history <- append(S_history, list(S))
        iter <- 0
        while((md > qchisq(q, p))&(iter < max_iter)){
          iter <- iter + 1
          S_old <- numeric()
          count <- 0
          while(((length(S_old) != length(S))|(length(S_old) == 0))&(count <= p)){
            count <- count + 1
            S_old <- S
            phi <- shapley_base(x_tilde, mu, Sigma_inv, inverted = TRUE)
            S <- unique(c(S, which(abs(phi - max(phi))<10^-5)))
            S_history <- append(S_history, list(S))
          }
          correction <- coordinate_correction(x_tilde = x_tilde, mu = mu, mu_tilde = mu, Sigma_inv = Sigma_inv,
                                              phi = phi, S = S, step_size = step_size, max_step = max_step)
          x_tilde <- correction$x
          x_history <- rbind(x_history, correction$x_history)
          phi_history <- rbind(phi_history, correction$phi_history)

          phi <- shapley_base(x_tilde, mu, Sigma_inv, inverted = TRUE)
          md <- sum(phi)
        }
        if(min_deviation>0){
          difference <- abs(x_tilde - x)
          weak_deviation <- which(difference < max(difference)*min_deviation)
          x_tilde[weak_deviation] <- x[weak_deviation]
        }
        x_final <- x_tilde
      }
  } else{ # MATRIX
    p <- ncol(x)
    n <- nrow(x)
    check_data(x)
    check_vector(mu, p)
    check_matrix(Sigma, p)
    if(!is.null(cells)){
      check_indicator_matrix(cells, p = p, n = n)
    }
    if(is.null(Sigma_inv)){
      Sigma_inv <- solve(Sigma)
    }
    check_matrix(Sigma_inv, p)
    phi <- array(dim = dim(x))
    mu_tilde <- array(dim = dim(x))
    non_centrality <- numeric()
    for(i in 1:nrow(x)){
      call_shapley <- shapley(x[i,], mu, Sigma_inv, inverted = TRUE, method = method, check = FALSE, cells = cells[i,])
      phi[i,] <- call_shapley$phi
      mu_tilde[i,] <- call_shapley$mu_tilde
      non_centrality[i] <- call_shapley$non_centrality
    }

    x <- as.matrix(x)
    x_final <- array(NA, dim = dim(x), dimnames = dimnames(x))
    phi_start <- array(NA, dim = dim(x), dimnames = dimnames(x))
    x_history <- list()
    phi_history <- list()
    S_history <- list()
    for(i in 1:nrow(x)){
      call_SCD <- SCD(x = x[i,], mu = mu, Sigma = Sigma, Sigma_inv = Sigma_inv, step_size = step_size,
                      min_deviation = min_deviation, max_step = max_step, max_iter = max_iter, q = q, check = FALSE, cells = cells[i,])
      x_final[i,] <- call_SCD$x
      phi_start[i,] <- call_SCD$phi
      x_history[[i]] <- call_SCD$x_history
      phi_history[[i]] <- call_SCD$phi_history
      S_history[[i]] <- call_SCD$S_history
    }
  }
  res <- new_shapley_algorithm(x = x_final, phi = phi_start, x_original = x, x_history = x_history, phi_history = phi_history, S_history = S_history)
  return(res)
}
