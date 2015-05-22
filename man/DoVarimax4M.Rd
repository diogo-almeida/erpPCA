\name{DoVarimax4M}
\alias{DoVarimax4M}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
DoVarimax4M(X, maxit = 100, tol = 1e-04, normalize = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{maxit}{
%%     ~~Describe \code{maxit} here~~
}
  \item{tol}{
%%     ~~Describe \code{tol} here~~
}
  \item{normalize}{
%%     ~~Describe \code{normalize} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, maxit = 100, tol = 1e-04, normalize = TRUE) 
{
    repmat <- function(a, n, m) {
        kronecker(matrix(1, n, m), a)
    }
    cat("--------- Varimax Rotation (4M) -------------\n")
    p <- nrow(X)
    m <- ncol(X)
    cat("Matrix rows:                     ", p, "\n")
    cat("Matrix columns:                  ", m, "\n")
    cat("Max. # of iterations:            ", maxit, "\n")
    cat("Convergence criterion:           ", tol, "\n")
    if (normalize) {
        cat("Kaiser's Normalization:          ", "Yes", "\n")
        h <- sqrt(colSums(t(X)^2))
        H <- repmat(h, 1, m)
        Y <- X/H
        Y[is.nan(Y)] <- X[is.nan(Y)]
    }
    else {
        cat("Kaiser's Normalization:          ", "No", "\n")
        Y <- X
    }
    g <- ComputeSimplicityCriterion(Y)
    it <- 0
    G <- matrix(0, nrow = maxit, ncol = 3)
    G[1, ] <- c(it, g, tol)
    Gold <- g
    YY <- t(Y)
    cat("     #    SimplicityCriterion     Convergence\n")
    cat(sprintf("\%6d     \%18.8f \%14.8f", it, g, tol))
    cat("\n")
    for (it in 1:maxit) {
        for (i in 1:(m - 1)) {
            for (j in (i + 1):m) {
                t.park <- Arg(sum(complex(real = Y[, i], imaginary = Y[, 
                  j])^4)/p - (sum(complex(real = Y[, i], imaginary = Y[, 
                  j])^2)/p)^2)/4
                XY <- matrix(c(Y[, i], Y[, j]), ncol = 2, byrow = F) \%*\% 
                  matrix(c(cos(t.park), -sin(t.park), sin(t.park), 
                    cos(t.park)), nrow = 2, byrow = T)
                Y[, i] <- XY[, 1]
                Y[, j] <- XY[, 2]
            }
        }
        g <- ComputeSimplicityCriterion(Y)
        cat(sprintf("\%6d     \%18.8f \%14.8f", it, g, Gold - g))
        cat("\n")
        if ((Gold - g) < tol) {
            if (Gold < g) {
                Y <- YY
            }
            else {
                G[it + 1, ] <- c(it, g, (Gold - g))
            }
            break
        }
        YY <- Y
        G[it + 1, ] = c(it, g, (Gold - g))
        Gold <- g
    }
    if (normalize) {
        Y <- Y * H
    }
    cat("---------------------------------------------")
    output <- list(Y = Y, G = G[1:it, ])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
