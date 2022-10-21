#' Detecting cellwise outliers using Shapley values based on local outlyingness.
#'
#' @description The \code{MOE} function indicates outlying cells for
#' a data vector with \eqn{p} entries or data matrix with \eqn{n \times p} entries containing only numeric entries \code{x}
#' for a given center \code{mu} and covariance matrix \code{Sigma} using the Shapley value.
#' It is a more sophisticated alternative to the \code{\link{SCD}} algorithm,
#' which uses the information of the regular cells to derive an alternative reference point  \insertCite{Mayrhofer2022}{ShapleyOutlier}.
#'
#' @param local Logical. If TRUE (default), the non-central Chi-Squared distribution is used to determine the cutoff value based on \code{mu_tilde}.
#' @param check_outlyingness Logical. If TRUE (default), the outlyingness is rechecked after applying \code{min_deviation}.
#' @inheritParams SCD
#'
#' @return A list of class \code{shapley_algorithm} (\code{\link{new_shapley_algorithm}}) containing the following:
#' \item{\code{x}}{A \eqn{p}-dimensional vector (or a \eqn{n \times p} matrix) containing the imputed data.}
#' \item{\code{phi}}{A \eqn{p}-dimensional vector (or a \eqn{n \times p} matrix) containing the Shapley values (outlyingness-scores) of \code{x}; see \code{\link{shapley}}.}
#' \item{\code{mu_tilde}}{A \eqn{p}-dimensional vector (or a \eqn{n \times p} matrix) containing the alternative reference points based on the regular cells of the original observations.}
#' \item{\code{x_original}}{A \eqn{p}-dimensional vector (or a \eqn{n \times p} matrix) containing the original data.}
#' \item{\code{x_original}}{The non-centrality parameters for the Chi-Squared distribution}
#' \item{\code{x_history}}{A list with \eqn{n} elements, each containing the path of how the original data vector was modified.}
#' \item{\code{phi_history}}{A list with \eqn{n} elements, each containing the Shapley values corresponding to \code{x_history}.}
#' \item{\code{mu_tilde_history}}{A list with \eqn{n} elements, each containing the alternative reference points corresponding to \code{x_history}.}
#' \item{\code{S_history}}{A list with \eqn{n} elements, each containing the indices of the outlying cells in each iteration.}
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
#' MOE_x <- MOE(x = x, mu = mu, Sigma = Sigma)
#' plot(MOE_x)
#'
#' library(MASS)
#' set.seed(1)
#' n <- 100; p <- 10
#' mu <- rep(0,p)
#' Sigma <- matrix(0.9, p, p); diag(Sigma) = 1
#' X <- mvrnorm(n, mu, Sigma)
#' X[sample(1:(n*p), 100, FALSE)] <- rep(c(-5,5),50)
#' MOE_X <- MOE(X, mu, Sigma)
#' plot(MOE_X, subset = 20)
MOE <- function(x, mu, Sigma, Sigma_inv = NULL, step_size = 0.1, min_deviation = 0, max_step = NULL,
                local = TRUE, max_iter = 1000, q = 0.99, check_outlyingness = FALSE, check = TRUE, cells = NULL, method = "cellMCD"){
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

    # ACTUAL COMPUTATION MOE
    x_tilde <- x
    x_final <- x
    mu_tilde <- mu
    d <- numeric(length = p)
    non_centrality <- 0

    #md of original observation
    phi <- shapley_base(x, mu, Sigma_inv, inverted = TRUE)
    md <- sum(phi)

    #historys
    x_history <- numeric()
    phi_history <- numeric()
    mu_tilde_history <- numeric()
    S_history <- list()

    if(md > qchisq(q, p)){
      iter <- 0
      conv_crit <- TRUE
      while(conv_crit&(iter < max_iter)){
        iter <- iter + 1
        S_old <- numeric()
        count <- 0
        while(((length(S_old) != length(S))|(length(S_old) == 0))&(count <= p)){
          count <- count + 1
          S_old <- S
          mu_tilde <- min_MD(x = x, mu = mu, Sigma_inv = Sigma_inv , S = S)
          phi <- shapley_base(x_tilde, mu_tilde, Sigma_inv, inverted = TRUE)
          S <- unique(c(S, which(abs(phi - max(phi))<10^-5)))
        }

        correction <- coordinate_correction(x_tilde = x_tilde, mu = mu, mu_tilde = mu_tilde, Sigma_inv = Sigma_inv,
                                            phi = phi, S = S, step_size = step_size, max_step = max_step)
        x_tilde <- correction$x
        phi <- correction$phi
        d <- d + correction$d

        x_history <- rbind(x_history, correction$x_history)
        phi_history <- rbind(phi_history, correction$phi_history)
        mu_tilde_history <- rbind(mu_tilde_history, correction$mu_tilde_history)
        S_history <- append(S_history, list(S))
        if(local){
          non_centrality <- mahalanobis(mu_tilde, mu, Sigma_inv, inverted = TRUE)
        }
        conv_crit <- sum(phi) > qchisq(q, p, ncp = non_centrality)
      }
      x_final <- x_tilde
      x_final[S] <- mu_tilde[S]
      # Detection threshold
      if(min_deviation>0){
        difference <- abs(d)/sqrt(diag(Sigma))
        weak_deviation <- which(difference < max(difference)*min_deviation)
        x_tilde_weak <- x_tilde
        x_tilde_weak[weak_deviation] <- x[weak_deviation]

        S_weak <- which(abs(x_tilde_weak - x) > 0)
        mu_tilde_weak <- min_MD(x = x, mu = mu, Sigma_inv = Sigma_inv , S = S_weak)
        x_tilde_weak[S_weak] <- mu_tilde_weak[S_weak]
        phi_weak <- shapley_base(x_tilde_weak, mu_tilde_weak, Sigma_inv, inverted = TRUE)
        if(local){
          non_centrality_weak <- mahalanobis(mu_tilde_weak, mu, Sigma_inv, inverted = TRUE)
        }
        #check if x_tilde_weak is outlying or not (FALSE --> outlier, TRUE --> regular)
        if(check_outlyingness & (sum(phi_weak) < qchisq(q, p, ncp = non_centrality_weak))){
          phi <- phi_weak
          x_final <- x_tilde_weak
          mu_tilde <- mu_tilde_weak
          non_centrality <- non_centrality_weak

          x_history <- rbind(x_history, x_tilde_weak)
          phi_history <- rbind(phi_history, phi_weak)
          mu_tilde_history <- rbind(mu_tilde_history, mu_tilde_weak)
          S_history <- append(S_history, list(S_weak))
        } else{ # if outlyingness of changed observation is not checked, overwrite
          phi <- phi_weak
          x_final <- x_tilde_weak
          mu_tilde <- mu_tilde_weak
          non_centrality <- non_centrality_weak

          x_history <- rbind(x_history, x_tilde_weak)
          phi_history <- rbind(phi_history, phi_weak)
          mu_tilde_history <- rbind(mu_tilde_history, mu_tilde_weak)
          S_history <- append(S_history, list(S_weak))
        }
      }
      phi <- shapley_base(x, mu_tilde, Sigma_inv, inverted = TRUE)
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

    x <- as.matrix(x)
    x_final <- array(NA, dim = dim(x), dimnames = dimnames(x))
    phi <- array(NA, dim = dim(x), dimnames = dimnames(x))
    mu_tilde <- array(NA, dim = dim(x), dimnames = dimnames(x))
    x_history <- list()
    phi_history <- list()
    S_history <- list()
    mu_tilde_history <- list()
    non_centrality <- numeric(length = nrow(x))
    names(non_centrality) <- rownames(x)
    for(i in 1:nrow(x)){
      call_MOE <- MOE(x = x[i,], mu = mu, Sigma = Sigma, Sigma_inv = Sigma_inv, step_size = step_size,
                 min_deviation = min_deviation, max_step = max_step, local = local, max_iter = max_iter, q = q,
                 check_outlyingness = check_outlyingness, check = FALSE, cells = cells[i,], method = method)
      x_final[i,] <- call_MOE$x
      phi[i,] <- call_MOE$phi
      mu_tilde[i,] <- call_MOE$mu_tilde
      non_centrality[i] <- call_MOE$non_centrality

      x_history[[i]] <- call_MOE$x_history
      phi_history[[i]] <- call_MOE$phi_history
      mu_tilde_history[[i]] <- call_MOE$mu_tilde_history
      S_history[[i]] <- call_MOE$S_history
    }
  }
  res <- new_shapley_algorithm(x = x_final, phi = phi, x_original = x, mu_tilde = mu_tilde, non_centrality = non_centrality,
                               x_history = x_history, phi_history = phi_history, S_history = S_history, mu_tilde_history = mu_tilde_history)
  return(res)
}
