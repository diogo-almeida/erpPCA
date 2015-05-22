#' @title Performs Varimax rotation of principal components or factors.
#' 
#' @description This function performs a Varimax rotation on factor loadings from an
#' eigenvalue or singular value decomposition. According to the original
#' MATLAB code by J\"{u}rgen Kayser, it ```emulates algorithms described by Harman
#' (1967, pp. 304-308) as implemented in BMDP-4M (Dixon, 1992, pp. 602-603).'''
#' 
#' @param X A matrix. Should be the loadings matrix from eigenvalue or 
#'   singular value decomposition.
#' @param maxit Integer. Maximum number of iterations allowed before accepting
#'  a solution to the rotation.
#' @param tol Convergence criterion.
#' @param normalize Logical scalar. Should Kaiser normalization be used? 
#'   Defaults to TRUE.
#' @param verbose Logical scalar. Should text output on the attempted rotation 
#'   be displayed? Defaults to FALSE.
#'   
#' @return List containing the rotated loadings (Y) and the history of 
#'   attempted solutions for the rotation problem (G).
#'   
#' @seealso \code{\link[stats]{varimax}} and \code{\link[GPArotation]{Varimax}}
#'   for other implementations of Varimax rotations.
#'   
#' @export
#'   
#' @examples
#' ## from the help of stats::varimax, using factor analysis
#' fa <- factanal( ~., 2, data = swiss)
#' stats::varimax(loadings(fa), normalize = FALSE)
#' DoVarimax4M(loadings(fa), normalize = FALSE)$Y
#' 
#' ## using the iris dataset and principal components analysis
#' iris.pca <- prcomp(iris[, 1:4])
#' stats::varimax(iris.pca$rotation)
#' DoVarimax4M(iris.pca$rotation)$Y
DoVarimax4M <- function(X, maxit = 100, tol = 1e-4, normalize = TRUE, 
                        verbose = FALSE) {
  
  # Auxiliary function definitions
  # ------------------------------
  # repmat (emulate repmat from matlab)
  # ComputeSimplicityCriterion
  # ------------------------------
  repmat <- function(a, n, m) {
    kronecker(matrix(1, n, m), a)
  }
  #
  ComputeSimplicityCriterion <- function(Y) {
    g <- 0
    ncol.Y <- ncol(Y)
    nrow.Y <- nrow(Y)
    for (i in 1:ncol.Y) {
      for (j in 1:ncol.Y) {
        if (i != j) {
          g <- g + sum( (Y[, i]^2) * (Y[, j]^2) ) -
            ((1/nrow.Y) * (sum(Y[, i]^2) * sum(Y[, j]^2)))
        }
      }
    }
    G <- g
    G
  } 
  
  # Varimax4M starts here.
  if (verbose) {
    cat("--------- Varimax Rotation (4M) -------------\n")
  }
  p <- nrow(X)
  m <- ncol(X)
  
  if (verbose) {
    cat("Matrix rows:                     ", p,     "\n")
    cat("Matrix columns:                  ", m,     "\n")
    cat("Max. # of iterations:            ", maxit, "\n")
    cat("Convergence criterion:           ", tol,   "\n")
  }  
  if (normalize) {
    if (verbose) {
      # Use Kaiser's normalization across the rows
      cat("Kaiser's Normalization:          ", "Yes",   "\n")
    }
    h <- sqrt(colSums(t(X)^2))   # communality column vector
    H <- repmat(h, 1, m)  # communality normalization matrix
    Y <- X / H               # normalize X by rows
    Y[is.nan(Y)] <- X[is.nan(Y)]
  } else {
    if (verbose) {
      cat("Kaiser's Normalization:          ", "No",   "\n")
    }
    Y <- X
  }
  
  # Compute rotation criterion for input matrix
  # SimplicityG function
  g <- ComputeSimplicityCriterion(Y)
  
  it <- 0
  G <- matrix(0, nrow = maxit, ncol = 3)
  G[1, ] <- c(it, g, tol)
  Gold <- g    # previous simplicity criterion
  YY <- t(Y)   # rotated matrix at begining of current iteration
  if (verbose) {
    cat("     #    SimplicityCriterion     Convergence\n")
    cat(sprintf("%6d     %18.8f %14.8f", it, g, tol))
    cat("\n")
  }
  for (it in 1:maxit) {
    for (i in 1:(m - 1)) {
      for (j in (i + 1):m) {
        # computes the rotation angle phi as the angle in the complex plane 
        # (Park, 2003):
        t.park <- Arg(
          sum(complex(real = Y[, i], imaginary = Y[, j])^4) / p -
            (sum(complex(real = Y[, i], imaginary = Y[, j])^2) / p)^2
        ) / 4
        
        # rotate the two vectors
        XY <- matrix(c(Y[,i], Y[, j]), ncol = 2, byrow = F) %*% matrix(
          c(cos(t.park), -sin(t.park), sin(t.park), cos(t.park)), nrow = 2, 
          byrow = T
        )
        
        # replace the two columns in the matrix
        Y[, i] <- XY[, 1]
        Y[, j] <- XY[, 2]      
      }
    }
    # Compute rotation criterion for the iteration
    g <- ComputeSimplicityCriterion(Y)
    if (verbose) {
      cat(sprintf("%6d     %18.8f %14.8f", it, g, Gold - g))
      cat("\n")
    }
    if ((Gold - g) < tol) {
      if (Gold < g) {     # if the previous solution was better
        Y <- YY           # report the previous one
      } else {                  
        G[it + 1, ] <- c(it, g, (Gold - g))
      }
      break
    }    
    
    YY <- Y
    G[it + 1, ] = c(it, g, (Gold - g))
    Gold <- g
  }
  
  if (normalize) {
    Y <- Y * H             # reverse Kaiser's normalization
  }
  if (verbose) {
    cat("---------------------------------------------\n")
  }
  output <- list(Y = Y, G = G[1:it, ])
  output
}
