% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shapley.R
\name{shapley}
\alias{shapley}
\title{Decomposition of squared Mahalanobis distance using Shapley values.}
\usage{
shapley(
  x,
  mu = NULL,
  Sigma = NULL,
  inverted = FALSE,
  method = "cellMCD",
  check = TRUE,
  cells = NULL
)
}
\arguments{
\item{x}{Data vector with \eqn{p} entries or data matrix with \eqn{n \times p} entries containing only numeric entries.}

\item{mu}{Either \code{NULL} (default) or mean vector of \code{x}. If NULL, \code{method} is used for parameter estimation.}

\item{Sigma}{Either \code{NULL} (default) or covariance matrix \eqn{p \times p} of \code{x}. If NULL, \code{method} is used for parameter estimation.}

\item{inverted}{Logical. If \code{TRUE}, \code{Sigma} is supposed to contain the inverse of the covariance matrix.}

\item{method}{Either "cellMCD" (default) or "MCD". Specifies the method used for parameter estimation if \code{mu} and/or \code{Sigma} are not provided.}

\item{check}{Logical. If \code{TRUE} (default), inputs are checked before running the function
and an error message is returned if one of the inputs is not as expected.}

\item{cells}{Either \code{NULL} (default) or a vector/matrix of the same dimension as \code{x},
indicating the outlying cells. The matrix must contain only zeros and ones, or \code{TRUE}/\code{FALSE}.}
}
\value{
\item{\code{phi}}{A \eqn{p}-dimensional vector (or a \eqn{n \times p} matrix) containing the Shapley values (outlyingness-scores) of \code{x}.}
\item{\code{mu_tilde}}{A \eqn{p}-dimensional vector (or a \eqn{n \times p} matrix) containing the
alternative reference points based on the regular cells of the original observations.}
\item{\code{non_centrality}}{The non-centrality parameters for the Chi-Squared distribution, given by \code{mahlanobis(mu_tilde, mu, Sigma)}}
}
\description{
The \code{shapley} function computes a \eqn{p}-dimensional vector containing the decomposition of the
squared Mahalanobis distance of \code{x} (with respect to \code{mu} and \code{Sigma})
into outlyingness contributions of the individual variables \insertCite{Mayrhofer2022}{ShapleyOutlier}.
The value of the \eqn{j}-th coordinate of this vector represents the
average marginal contribution of the \eqn{j}-th variable to the squared Mahalanobis distance of
the individual observation \code{x}.\cr
If \code{cells} is provided, Shapley values of \code{x} are computed with respect to a local reference point,
that is based on a cellwise prediction of each coordinate, using the information of the regular cells of \code{x}, see \insertCite{Mayrhofer2022}{ShapleyOutlier}.\cr
If \code{x} is a \eqn{n \times p} matrix, a \eqn{n \times p} matrix is returned, containing the decomposition for each row.
}
\examples{
## Without outlying cells as input in the 'cells' argument#'
# Single observation
p <- 5
mu <- rep(0,p)
Sigma <- matrix(0.9, p, p); diag(Sigma) = 1
Sigma_inv <- solve(Sigma)
x <- c(0,1,2,2.3,2.5)
shapley(x, mu, Sigma)
phi <- shapley(x, mu, Sigma_inv, inverted = TRUE)
plot(phi)

# Multiple observations
library(MASS)
set.seed(1)
n <- 100; p <- 10
mu <- rep(0,p)
Sigma <- matrix(0.9, p, p); diag(Sigma) = 1
X <- mvrnorm(n, mu, Sigma)
X_clean <- X
X[sample(1:(n*p), 100, FALSE)] <- rep(c(-5,5),50)
call_shapley <- shapley(X, mu, Sigma)
plot(call_shapley, subset = 20)


## Giving outlying cells as input in the 'cells' argument
# Single observation
p <- 5
mu <- rep(0,p)
Sigma <- matrix(0.9, p, p); diag(Sigma) = 1
Sigma_inv <- solve(Sigma)
x <- c(0,1,2,2.3,2.5)
call_shapley <- shapley(x, mu, Sigma_inv, inverted = TRUE,
method = "cellMCD", check = TRUE, cells = c(1,1,0,0,0))
plot(call_shapley)

# Multiple observations
library(MASS)
set.seed(1)
n <- 100; p <- 10
mu <- rep(0,p)
Sigma <- matrix(0.9, p, p); diag(Sigma) = 1
X <- mvrnorm(n, mu, Sigma)
X_clean <- X
X[sample(1:(n*p), 100, FALSE)] <- rep(c(-5,5),50)
call_shapley <- shapley(X, mu, Sigma, cells = (X_clean - X)!=0)
plot(call_shapley, subset = 20)
}
\references{
\insertRef{Mayrhofer2022}{ShapleyOutlier}
}
