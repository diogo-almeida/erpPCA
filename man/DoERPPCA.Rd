\name{DoERPPCA}
\alias{DoERPPCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
DoERPPCA(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
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
function (X) 
{
    require(GPArotation)
    repmat <- function(a, n, m) {
        kronecker(matrix(1, n, m), a)
    }
    if (is.list(X)) {
        cat("Your data needs to be input as a matrix. You input a list.\n")
        cat("I'm converting your data to the appropriate format now.\n")
        X <- as.matrix(X)
        if (is.double(X)) {
            cat("Your data is now a matrix!\n")
        }
    }
    D <- cov(X)
    eigenD <- eigen(D)
    EM <- eigenD$vectors
    EV <- eigenD$values
    UL <- EM \%*\% sqrt(diag(EV))
    rk <- sum(svd(cor(X))$d > 1e-04)
    u <- sort(EV, decreasing = T)[1:rk]
    LU <- UL[, 1:rk]
    s <- matrix(1, ncol = rk)
    s[abs(apply(LU, 2, max)) < abs(apply(LU, 2, min))] <- -1
    LU <- LU * repmat(s, dim(LU)[1], 1)
    RL <- DoVarimax4M(LU)$Y
    EVr <- colSums(RL * RL)
    r <- sort(EVr, decreasing = T)
    rx <- order(EVr, decreasing = T)
    LR <- RL[, rx]
    s <- matrix(1, ncol = ncol(LR))
    s[abs(apply(LR, 2, max)) < abs(apply(LR, 2, min))] <- -1
    LR <- LR * repmat(s, ncol(X), 1)
    tv <- sum(EV)
    VT <- matrix(c(u, (100 * u)/tv, r, (100 * r)/tv), byrow = F, 
        nrow = length(u))
    FSCFr <- LR \%*\% solve(t(LR) \%*\% LR)
    FSCFr <- FSCFr * repmat(sqrt(diag(D)), 1, rk)
    mu <- apply(X, 2, mean)
    sigma <- apply(X, 2, sd)
    Xc <- sweep(X, 2, mu)
    FSr = matrix(0, nrow = nrow(X), ncol = ncol(X))
    for (n in 1:nrow(X)) {
        for (m in 1:rk) {
            FSr[n, m] <- sum((Xc[n, ]/sigma) * FSCFr[, m])
        }
    }
    factor.labels <- paste("Factor", 1:nrow(VT), sep = ".")
    header.labels <- c("Eigenvalue.Before", "Var.Explained.Before", 
        "Eigenvalue.After", "Var.Explained.After")
    rownames(VT) <- factor.labels
    colnames(VT) <- header.labels
    output <- list(Unrotated = LU, Rotated = LR, Factor.Scores = FSr, 
        Variance = VT)
    return(output)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
