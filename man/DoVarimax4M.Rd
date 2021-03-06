% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/DoVarimax4M.R
\name{DoVarimax4M}
\alias{DoVarimax4M}
\title{Performs Varimax rotation of principal components or factors.}
\usage{
DoVarimax4M(X, maxit = 100L, tol = 1e-04, normalize = TRUE,
  verbose = FALSE)
}
\arguments{
\item{X}{A matrix. Should be the loadings matrix from eigenvalue or
  singular value decomposition.}

\item{maxit}{Integer. Maximum number of iterations allowed before accepting
 a solution to the rotation. Defaults to 100 iterations.}

\item{tol}{Convergence criterion. Defaults to 0.0001.}

\item{normalize}{Logical scalar. Should Kaiser normalization be used?
  Defaults to TRUE.}

\item{verbose}{Logical scalar. Should text output on the attempted rotation
  be displayed? Defaults to FALSE.}
}
\value{
List containing
\item{Y}{A matrix. Rotated loadings from \code{X}.}
\item{G}{A matrix. History of attempted solutions for the rotation problem.}
}
\description{
This function performs a Varimax rotation on factor loadings
from an eigenvalue or singular value decomposition. According to the original
MATLAB code by Jurgen Kayser, it ``emulates algorithms described by Harman
(1967, pp. 304-308) as implemented in BMDP-4M (Dixon, 1992, pp. 602-603).''
}
\examples{
## from the help of stats::varimax, using factor analysis
fa <- factanal( ~., 2, data = swiss)
stats::varimax(loadings(fa), normalize = FALSE)
DoVarimax4M(loadings(fa), normalize = FALSE)$Y

## using the iris dataset and principal components analysis
iris.pca <- prcomp(iris[, 1:4])
stats::varimax(iris.pca$rotation)
DoVarimax4M(iris.pca$rotation)$Y
}
\seealso{
\code{\link[stats]{varimax}} from package \code{stats} and
  \code{\link[GPArotation]{Varimax}} from package \code{GPArotation} for
  other implementations of Varimax rotations.
}

