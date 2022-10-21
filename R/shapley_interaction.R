#' Decomposition of squared Mahalanobis distance using Shapley interaction indices.
#'
#' @description The \code{shapley_interaction} function computes a \eqn{p \times p} matrix
#' containing pairwise outlyingness scores based on Shapley interaction indices.
#' It decomposes the squared Mahalanobis distance of \code{x} (with respect to \code{mu} and \code{Sigma})
#' into outlyingness contributions of pairs of variables \insertCite{Mayrhofer2022}{ShapleyOutlier}.
#'
#' @param x Data vector with \eqn{p} entries containing only numeric entries.
#' @inheritParams shapley
#'
#' @return A \eqn{p \times p} matrix containing the decomposition of the squared Mahalanobis distance of \code{x}
#' into outlyingness scores for pairs of variables with respect to \code{mu} and \code{Sigma}.
#' @export
#'
#' @references
#' \insertRef{Mayrhofer2022}{ShapleyOutlier}
#'
#' @examples
#' p <- 5
#' mu <- rep(0,p)
#' Sigma <- matrix(0.9, p, p); diag(Sigma) = 1
#' Sigma_inv <- solve(Sigma)
#' x <- c(0,1,2,2.3,2.5)
#' shapley_interaction(x, mu, Sigma)
#' PHI <- shapley_interaction(x, mu, Sigma_inv, inverted = TRUE)
#' plot(PHI)
shapley_interaction <- function(x, mu, Sigma, inverted = FALSE){
  p <- length(x)
  check_vector(x)
  check_vector(mu, p)
  check_matrix(Sigma, p)

  if(inverted){
    Sigma_inv <- Sigma
  } else{
    Sigma_inv <- solve(Sigma)
  }
  interaction_matrix <- 2* tcrossprod((x-mu)) * Sigma_inv
  interaction_values <- interaction_matrix - diag(rowSums(interaction_matrix/2))
  return(new_shapley_interaction(PHI = interaction_values))
}
