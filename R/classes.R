#' Class constructor for class \code{shapley}.
#'
#' This function creates an object of class \code{shapley} that is returned by the \code{\link{shapley}} function.
#'
#' @param phi A \eqn{p}-dimensional vector (or a \eqn{n \times p} matrix) containing the
#' Shapley values (outlyingness-scores) of a \eqn{p}-dimensional data vector (or a \eqn{n \times p} data matrix).
#' @param mu_tilde Optional. A \eqn{p}-dimensional vector (or a \eqn{n \times p} matrix) containing the
#' alternative reference points based on the regular cells of the original observations.
#' @param non_centrality Optional. The non-centrality parameters for the Chi-Squared distribution,
#' which are given by \code{mahlanobis(mu_tilde, mu, Sigma)}.
#'
#' @return Named list of class \code{shapley}, containing the input parameters.
#' @export
new_shapley <- function(phi = numeric(), mu_tilde = NULL, non_centrality = NULL){
  stopifnot(is.numeric(phi)|is.null(phi))
  stopifnot(is.numeric(mu_tilde)|is.null(mu_tilde))
  stopifnot(is.numeric(non_centrality)|is.null(non_centrality))

  shapley <- list("phi" = phi, "mu_tilde" = mu_tilde, "non_centrality" = non_centrality)
  class(shapley) <- c("shapley", class(shapley))
  return(shapley)
}

#' Print function for class \code{shapley}.
#'
#' @param x List of class \code{shapley}.
#' @param ... Optional arguments passed to methods.
#'
#' @return Prints the list entries of \code{x} that are not \code{NULL}.
#' @export
print.shapley <- function(x,...){
  class(x) <- NULL
  if(is.null(x$mu_tilde)&is.null(x$non_centrality)){
    print(x$phi)
  } else{
    print(x)
  }
}

#' Class constructor for class \code{shapley_algorithm}.
#'
#' This function creates an object of class \code{shapley_algorithm} that is returned
#' by the \code{\link{SCD}} and \code{\link{MOE}} functions.
#'
#'
#' @param x A \eqn{p}-dimensional vector (or a \eqn{n \times p} matrix) containing the imputed data.
#' @param x_original A \eqn{p}-dimensional vector (or a \eqn{n \times p} matrix) containing the original data.
#' @param x_history Optional. A list with \eqn{n} elements, each containing the path of how the original data vector was modified.
#' @param phi_history Optional. A list with \eqn{n} elements, each containing the Shapley values corresponding to \code{x_history}.
#' @param mu_tilde_history Optional. A list with \eqn{n} elements, each containing the alternative reference points corresponding to \code{x_history}.
#' @param S_history Optional. A list with \eqn{n} elements, each containing the indices of the outlying cells in each iteration.
#' @inheritParams new_shapley
#'
#' @return  Named list of class \code{shapley_algorithm}, containing the input parameters.
#' @export
new_shapley_algorithm <- function(x = numeric(), phi = numeric(), x_original = numeric(),
                                  mu_tilde = NULL, non_centrality = NULL,
                                  x_history = NULL, phi_history = NULL, mu_tilde_history = NULL, S_history = NULL){
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(phi))
  stopifnot(is.numeric(x_original))

  shapley_algorithm <- list("x" = x, "phi" = phi, "x_original" = x_original,
                            "mu_tilde" = mu_tilde, "non_centrality" = non_centrality,
                            "x_history" = x_history, "phi_history" = phi_history, "mu_tilde_history" = mu_tilde_history, "S_history" = S_history)
  class(shapley_algorithm) <- c("shapley_algorithm", class(shapley_algorithm))
  return(shapley_algorithm)
}

#' Print function for class \code{shapley_algorithm}.
#'
#' @param x List of class \code{shapley_algorithm}.
#' @param ... Optional arguments passed to methods.
#'
#' @return Prints the imputed data and the Shapley values.
#' @export
print.shapley_algorithm <- function(x,...){
  print(list("x" = x$x, "phi" = x$phi))
}


#'  Class constructor for class \code{shapley_interaction}.
#'
#' This function creates an object of class \code{shapley_interaction} that is returned
#' by the \code{\link{shapley_interaction}} function.
#'
#' @param PHI A \eqn{p \times p} matrix containing the decomposition of the squared Mahalanobis distance
#' of a \eqn{p}-dimensional numeric vector into outlyingness scores for pairs of variables.
#'
#' @return Matrix of class \code{shapley_interaction}, containing input matrix \code{PHI}.
#' @export
new_shapley_interaction <- function(PHI = numeric()){
  stopifnot(is.numeric(PHI))

  shapley_interaction <- PHI
  class(shapley_interaction) <- c("shapley_interaction", class(shapley_interaction))
  return(shapley_interaction)
}

#' Print function for class \code{shapley_interaction}.
#'
#' @param x Matrix of class \code{shapley_interaction}.
#' @param ... Optional arguments passed to methods.
#'
#' @return Prints the Shapley interaction indices.
#' @export
print.shapley_interaction <- function(x,...){
  class(x) <- NULL
  print(x)
}
